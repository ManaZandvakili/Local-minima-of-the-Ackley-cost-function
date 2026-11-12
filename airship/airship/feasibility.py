"""Bounded-thrust feasibility and capacity (replaces Final.py's rank-6 test).

The rank-6 test only guarantees full control if thrust can be negative. With
real propellers each thrust is bounded,  u_min <= u <= u_max  (u_min = 0 for
push-only), and the question becomes: how large a wrench can the layout
produce along a given direction e?

    capacity(e) = max lambda   s.t.  M u = lambda e,   u_min <= u <= u_max

(one small LP). A layout can produce every wrench direction exactly when
capacity > 0 along all 12 signed axes (+-Fx, ..., +-Mz): the reachable cone
then contains +-e_i for every i, which spans the whole wrench space. So the
12 axis capacities are both the feasibility test and the "capability vector"
handed to the RL agent.

Mixed (non-axis) directions are compared in a unit-consistent space: torque
rows divided by a characteristic length L (the hull radius), so capacities
there are in newtons.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog

from .wrench import DOF_NAMES, scale_torque_rows

AXIS_LABELS = tuple(f"{s}{n}" for n in DOF_NAMES for s in "+-")   # +Fx, -Fx, +Fy, ...


def _bounds(n, u_min, u_max):
    lo = np.broadcast_to(np.asarray(u_min, dtype=float), (n,))
    hi = np.broadcast_to(np.asarray(u_max, dtype=float), (n,))
    return lo, hi


def capacity(M, direction, u_min=0.0, u_max=1.0):
    """Largest lambda >= 0 with M u = lambda * direction and bounded u.

    `direction` need not be unit length; lambda is measured in multiples of it.
    Returns 0.0 when only lambda = 0 is achievable.
    """
    M = np.asarray(M, dtype=float)
    e = np.asarray(direction, dtype=float)
    n = M.shape[1]
    lo, hi = _bounds(n, u_min, u_max)
    if np.any(lo > 0) or np.any(hi < lo):
        raise ValueError("bounds must satisfy u_min <= 0 <= u_max (u = 0 must be allowed)")
    c = np.zeros(n + 1); c[-1] = -1.0                       # maximise lambda
    A_eq = np.hstack([M, -e[:, None]])
    bounds = list(zip(lo, hi)) + [(0.0, None)]
    res = linprog(c, A_eq=A_eq, b_eq=np.zeros(M.shape[0]), bounds=bounds, method="highs")
    if res.status != 0:
        raise RuntimeError(f"capacity LP failed: {res.message}")
    return max(float(res.x[-1]), 0.0)


def axis_capacities(M, u_min=0.0, u_max=1.0):
    """Capacities along the 12 signed axes, in physical units (N or N*m).

    Order: +Fx, -Fx, +Fy, -Fy, +Fz, -Fz, +Mx, -Mx, +My, -My, +Mz, -Mz.
    """
    out = np.empty(12)
    for k in range(12):
        e = np.zeros(6); e[k // 2] = 1.0 if k % 2 == 0 else -1.0
        out[k] = capacity(M, e, u_min, u_max)
    return out


def sample_directions(n, seed=0):
    """n unit directions in 6-D wrench space (seeded, uniform on the sphere)."""
    v = np.random.default_rng(seed).normal(size=(n, 6))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


@dataclass
class Requirements:
    """Minimum capacity per signed axis; a value of None means not required."""
    min_force_N: float = 1.0                 # +-Fx, +-Fy, +-Fz
    min_yaw_Nm: float = 1.0                  # +-Mz
    min_roll_pitch_Nm: float | None = 0.25   # +-Mx, +-My (lower: the gondola pendulum helps)
    required_dofs: tuple = DOF_NAMES

    def minima(self):
        """12-vector of minimum capacities (np.nan where not required)."""
        per_dof = {"Fx": self.min_force_N, "Fy": self.min_force_N, "Fz": self.min_force_N,
                   "Mx": self.min_roll_pitch_Nm, "My": self.min_roll_pitch_Nm, "Mz": self.min_yaw_Nm}
        out = []
        for name in DOF_NAMES:
            v = per_dof[name] if name in self.required_dofs else None
            out += [np.nan if v is None else float(v)] * 2
        return np.array(out)

    @classmethod
    def from_config(cls, d):
        d = dict(d or {})
        if "required_dofs" in d:
            d["required_dofs"] = tuple(d["required_dofs"])
        return cls(**d)


def shortfall(caps, req):
    """Sum over required axes of max(0, 1 - cap/min); 0 means feasible.

    Smooth and graded (each missing axis adds up to 1), so it can be used as an
    optimizer penalty without a cliff.
    """
    m = req.minima()
    use = ~np.isnan(m)
    return float(np.sum(np.clip(1.0 - caps[use] / m[use], 0.0, 1.0)))


@dataclass
class FeasibilityReport:
    feasible: bool
    axis_capacity: np.ndarray        # 12, physical units
    capacity_min_N: float            # worst over required axes and sampled directions (scaled space, N)
    shortfall: float

    def as_dict(self):
        return {
            "feasible": bool(self.feasible),
            "axis_capacity": {k: float(v) for k, v in zip(AXIS_LABELS, self.axis_capacity)},
            "capacity_min_N": float(self.capacity_min_N),
            "shortfall": float(self.shortfall),
        }


def assess(M, u_min, u_max, req, length_m, n_sampled=50, seed=0):
    """Full check: 12 axis capacities plus n_sampled mixed directions.

    Mixed directions are drawn in the torque-scaled space (rows / length_m), and
    only within the required DOFs. capacity_min_N is the smallest capacity over
    the required axes (torques divided by length_m) and the sampled directions.
    """
    caps = axis_capacities(M, u_min, u_max)
    sf = shortfall(caps, req)
    Ms = scale_torque_rows(M, length_m)
    scaled_axis = caps / np.repeat([1, 1, 1, length_m, length_m, length_m], 2)
    m = req.minima()
    worst = [scaled_axis[~np.isnan(m)].min()] if np.any(~np.isnan(m)) else []
    if n_sampled:
        req_idx = [DOF_NAMES.index(d) for d in req.required_dofs]
        dirs = np.zeros((n_sampled, 6))
        dirs[:, req_idx] = sample_directions(n_sampled, seed)[:, :len(req_idx)]
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        worst.append(min(capacity(Ms, e, u_min, u_max) for e in dirs))
    cap_min = float(min(worst)) if worst else float("nan")
    return FeasibilityReport(sf == 0.0 and cap_min > 0.0, caps, cap_min, sf)
