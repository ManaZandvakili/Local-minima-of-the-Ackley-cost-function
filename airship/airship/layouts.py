"""Propeller layouts: the JSON exchange format between the optimizer and
everything downstream (simulator, allocation, PID, RL, evaluation).

JSON fields (all lengths in metres, thrust in newtons):
    name          unique layout name, e.g. "sphere_pareto_03"
    shape         "sphere" or "ellipsoid"
    semi_axes_m   [a, b, c]
    positions_m   N x 3 mounting points, body frame (origin at hull centre, z up)
    directions    N x 3 unit thrust directions
    u_max_N       maximum thrust per propeller (scalar, or a list of N values)
    source        "optimizer" | "random" | "hand" | "perturbed" | "failure"
    objectives    dict of objective values (may be empty)
    seed          int or null
    notes         dict of free-form provenance (optional)
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from .geometry import Hull, unit
from .wrench import build_M

SOURCES = ("optimizer", "random", "hand", "perturbed", "failure")
_ON_HULL_TOL_M = 1e-6
_UNIT_TOL = 1e-6


@dataclass(frozen=True)
class Layout:
    name: str
    shape: str
    semi_axes_m: tuple
    positions_m: np.ndarray
    directions: np.ndarray
    u_max_N: object = 2.5
    source: str = "optimizer"
    objectives: dict = field(default_factory=dict)
    seed: int | None = None
    notes: dict = field(default_factory=dict)

    def __post_init__(self):
        p = np.atleast_2d(np.asarray(self.positions_m, dtype=float))
        d = np.atleast_2d(np.asarray(self.directions, dtype=float))
        if p.shape != d.shape or p.ndim != 2 or p.shape[1] != 3:
            raise ValueError(f"{self.name}: positions {p.shape} and directions {d.shape} must both be N x 3")
        norms = np.linalg.norm(d, axis=1)
        if np.any(np.abs(norms - 1.0) > _UNIT_TOL):
            raise ValueError(f"{self.name}: directions must be unit vectors (norms {norms.round(6)})")
        if self.source not in SOURCES:
            raise ValueError(f"{self.name}: source must be one of {SOURCES}, got {self.source!r}")
        u = np.asarray(self.u_max_N, dtype=float)
        if u.ndim not in (0, 1) or (u.ndim == 1 and u.shape[0] != p.shape[0]) or np.any(u < 0):
            raise ValueError(f"{self.name}: u_max_N must be >= 0, a scalar or one value per propeller")
        hull = Hull(self.shape, tuple(float(x) for x in self.semi_axes_m))
        off = hull.surface_distance(p)
        if np.any(off > _ON_HULL_TOL_M):
            raise ValueError(f"{self.name}: propellers {np.flatnonzero(off > _ON_HULL_TOL_M).tolist()} "
                             f"are off the hull surface (max {off.max():.3g} m)")
        object.__setattr__(self, "positions_m", p)
        object.__setattr__(self, "directions", d)
        object.__setattr__(self, "semi_axes_m", hull.semi_axes_m)

    @property
    def n_props(self):
        return self.positions_m.shape[0]

    @property
    def hull(self):
        return Hull(self.shape, self.semi_axes_m)

    @property
    def u_max(self):
        """Per-propeller maximum thrust as an N-vector [N]."""
        return np.broadcast_to(np.asarray(self.u_max_N, dtype=float), (self.n_props,)).copy()

    def M(self):
        """6 x N wrench matrix in physical units."""
        return build_M(self.positions_m, self.directions)

    def to_dict(self):
        u = np.asarray(self.u_max_N, dtype=float)
        return {
            "name": self.name,
            "shape": self.shape,
            "semi_axes_m": [float(x) for x in self.semi_axes_m],
            "positions_m": self.positions_m.tolist(),
            "directions": self.directions.tolist(),
            "u_max_N": float(u) if u.ndim == 0 else u.tolist(),
            "source": self.source,
            "objectives": {k: float(v) for k, v in self.objectives.items()},
            "seed": self.seed,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            name=d["name"], shape=d["shape"], semi_axes_m=tuple(d["semi_axes_m"]),
            positions_m=np.asarray(d["positions_m"], dtype=float),
            directions=np.asarray(d["directions"], dtype=float),
            u_max_N=d.get("u_max_N", 2.5), source=d.get("source", "optimizer"),
            objectives=dict(d.get("objectives") or {}), seed=d.get("seed"),
            notes=dict(d.get("notes") or {}),
        )


def save(layout, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(layout.to_dict(), indent=2) + "\n")
    return path


def load(path):
    return Layout.from_dict(json.loads(Path(path).read_text()))


def remove_prop(layout, i):
    """Layout with propeller i removed (a failed motor whose column is dropped from M)."""
    if not 0 <= i < layout.n_props:
        raise IndexError(f"propeller {i} out of range for {layout.n_props} propellers")
    keep = np.arange(layout.n_props) != i
    u = np.asarray(layout.u_max_N, dtype=float)
    return replace(
        layout,
        name=f"{layout.name}_minus{i}",
        positions_m=layout.positions_m[keep],
        directions=layout.directions[keep],
        u_max_N=float(u) if u.ndim == 0 else u[keep],
        source="failure",
        objectives={},
        notes={"parent": layout.name, "removed_index": int(i)},
    )


def _rotate_small(d, angle_rad, rng):
    """Rotate unit vector d by angle_rad about a random axis perpendicular to d."""
    a = rng.normal(size=3)
    a -= a.dot(d) * d
    a = unit(a)
    return unit(d * np.cos(angle_rad) + np.cross(a, d) * np.sin(angle_rad))


def perturb(layout, pos_sigma_m, ang_sigma_deg, rng, name=None):
    """Manufacturing-tolerance copy of a layout.

    Each mounting point gets Gaussian noise (std pos_sigma_m per axis) and is
    projected back onto the hull. Each thrust direction is tilted by a
    Gaussian angle (std ang_sigma_deg) about a random perpendicular axis.
    """
    hull = layout.hull
    p = hull.project(layout.positions_m + rng.normal(scale=pos_sigma_m, size=layout.positions_m.shape))
    tilts = np.deg2rad(rng.normal(scale=ang_sigma_deg, size=layout.n_props))
    d = np.array([_rotate_small(di, t, rng) for di, t in zip(layout.directions, tilts)])
    return replace(
        layout,
        name=name or f"{layout.name}_pert",
        positions_m=p, directions=d, source="perturbed", objectives={},
        notes={"parent": layout.name, "pos_sigma_m": float(pos_sigma_m),
               "ang_sigma_deg": float(ang_sigma_deg)},
    )
