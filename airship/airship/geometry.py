"""Hull shapes, candidate mounting points, surface normals and gondola gates.

Ported from reference/Final.py (fibonacci_directions, fibonacci_ellipsoid,
surface_normals, project_onto_surface, gondola_*_hard_check), changed only
to work in metres. Final.py stores positions normalised by the sphere radius
(1.575 m), so  positions_m = positions_normalised * 1.575.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

GONDOLA_DOWN = np.array([0.0, 0.0, -1.0])
GONDOLA_KEEPOUT_DEG = 35.0


def unit(v, eps=1e-12):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, eps)


def fibonacci_directions(K):
    """K unit vectors on a spherical Fibonacci lattice (same order as Final.py)."""
    i = np.arange(K)
    ga = 2 - (1 + 5**0.5) / 2
    theta = 2 * np.pi * ((i * ga) % 1.0)
    z = 1 - (2 * i + 1) / K
    r = np.sqrt(np.maximum(1 - z * z, 0))
    return unit(np.stack([r * np.cos(theta), r * np.sin(theta), z], axis=-1))


@dataclass(frozen=True)
class Hull:
    """An ellipsoidal hull (a sphere when all semi-axes are equal)."""
    shape: str                      # "sphere" or "ellipsoid"
    semi_axes_m: tuple              # (a, b, c) along body x, y, z [m]

    def __post_init__(self):
        a = np.asarray(self.semi_axes_m, dtype=float)
        if a.shape != (3,) or np.any(a <= 0):
            raise ValueError(f"semi_axes_m must be three positive numbers, got {self.semi_axes_m}")
        if self.shape == "sphere" and not np.allclose(a, a[0]):
            raise ValueError("a sphere needs equal semi-axes")
        if self.shape not in ("sphere", "ellipsoid"):
            raise ValueError(f"unknown shape {self.shape!r}")

    @classmethod
    def sphere(cls, radius_m):
        return cls("sphere", (radius_m,) * 3)

    @property
    def abc(self):
        return np.asarray(self.semi_axes_m, dtype=float)

    def candidate_points(self, K):
        """K Fibonacci mounting points on the hull surface [m]."""
        return fibonacci_directions(K) * self.abc

    def to_unit_sphere(self, positions_m):
        """Map hull points to unit-sphere directions (divide by the semi-axes)."""
        return unit(np.atleast_2d(positions_m) / self.abc)

    def normals(self, positions_m):
        """Outward unit surface normals at hull points."""
        p = np.atleast_2d(positions_m)
        return unit(p / self.abc**2)

    def project(self, positions_m):
        """Radial projection onto the hull surface (in unit-sphere space)."""
        return self.to_unit_sphere(positions_m) * self.abc

    def surface_distance(self, positions_m):
        """|p - project(p)| for each point [m]; 0 for points on the hull."""
        p = np.atleast_2d(positions_m)
        return np.linalg.norm(p - self.project(p), axis=1)

    # ── gondola keep-out gates (same rule as Final.py F10/F12/F13) ──────
    def gondola_position_violation(self, positions_m, keepout_deg=GONDOLA_KEEPOUT_DEG):
        """True for each point within keepout_deg of the bottom pole."""
        cos_ang = self.to_unit_sphere(positions_m) @ GONDOLA_DOWN
        return cos_ang > np.cos(np.deg2rad(keepout_deg))


def gondola_direction_violation(directions, keepout_deg=GONDOLA_KEEPOUT_DEG):
    """True for each thrust direction within keepout_deg of straight down."""
    return unit(np.atleast_2d(directions)) @ GONDOLA_DOWN > np.cos(np.deg2rad(keepout_deg))
