"""Wrench (allocation) matrix in physical units.

    w = M u,    w = [Fx, Fy, Fz, Mx, My, Mz]   [N, N, N, N*m, N*m, N*m]
                u = per-propeller thrust        [N]

    M = [ d_1        ...  d_N        ]   force rows   [N per N]
        [ r_1 x d_1  ...  r_N x d_N  ]   torque rows  [N*m per N]

Ported from Final.py's build_M_matrix. There, positions are normalised by the
hull radius R, so its torque rows are in units of (N * R); the simulator then
multiplies them by R. Building M from positions in metres gives the same
matrix directly (checked in tests/test_wrench.py).

The one-way-thrust capability metrics come in Phase 1 (feasibility.py).
"""
from __future__ import annotations

import numpy as np

DOF_NAMES = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")


def build_M(positions_m, directions):
    """6 x N wrench matrix from positions [m] and unit thrust directions."""
    r = np.atleast_2d(np.asarray(positions_m, dtype=float))
    d = np.atleast_2d(np.asarray(directions, dtype=float))
    if r.shape != d.shape or r.shape[1] != 3:
        raise ValueError(f"positions {r.shape} and directions {d.shape} must both be N x 3")
    return np.vstack([d.T, np.cross(r, d).T])


def scale_torque_rows(M, length_m):
    """Divide the torque rows by a characteristic length so all rows are in N.

    Use this before comparing singular values or condition numbers, which
    otherwise mix N and N*m. The paper should state the length used
    (default: the hull radius).
    """
    Ms = np.array(M, dtype=float, copy=True)
    Ms[3:] /= float(length_m)
    return Ms
