"""The paper's five objectives (and the authority / conditioning metrics used
for selection and regularisation), ported from reference/Final.py.

To keep the objective values identical to the paper, they are evaluated on
positions divided by the characteristic length L (the hull radius, 1.575 m),
i.e. on the torque-scaled wrench matrix M / L-scaled rows. This is the
stated choice of characteristic length for the paper text.

    paper name                code name in Final.py   key here
    J1 propeller count        J3                      count     N / (K/2)
    J2 control effort         J4                      effort    sum_j w_j ||M+ e_j||_inf
    J3 position symmetry      J5                      pos_sym   ||mean(r_i / L)||
    J4 direction symmetry     J6                      dir_sym   ||mean(d_i)||
    J5 geometric suitability  J8                      geometric spacing + angle diversity
                                                                + gondola soft penalty + downwash
Differences from Final.py:
  - no random noise is added before the SVD (Final.py adds 1e-9 * randn);
  - the excluded objectives (authority-as-cost, residual, fault tolerance,
    aerodynamic tangential) are not computed during optimisation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..geometry import GONDOLA_DOWN
from ..wrench import build_M, scale_torque_rows

OBJECTIVES = ("count", "effort", "pos_sym", "dir_sym", "geometric")
PAPER_KEYS = {"count": "J1_count", "effort": "J2_effort", "pos_sym": "J3_position_symmetry",
              "dir_sym": "J4_direction_symmetry", "geometric": "J5_geometric"}


@dataclass(frozen=True)
class GeometricRules:
    """Constants of the geometric term (values from Final.py)."""
    prop_diameter_m: float = 7 * 0.0254      # Final.py: 7-inch; spacing rule = 2.5 diameters
    spacing_diameters: float = 2.5
    ang_min_deg: float = 40.0
    gondola_keepout_deg: float = 35.0
    downwash_cone_deg: float = 15.0
    length_m: float = 1.575                  # characteristic length (hull radius)

    def d_min(self, hull):
        """Minimum spacing in normalised units, scaled by the geometric-mean semi-axis (Final.py I6)."""
        d = self.spacing_diameters * self.prop_diameter_m / self.length_m
        if hull.shape == "sphere":
            return d
        return d * float(np.prod(hull.abc / self.length_m)) ** (1.0 / 3.0)


def _gondola_soft_penalty(p, d, keepout_deg):
    """Final.py gondola_keepout_penalty (positions normalised by their own length)."""
    k = np.deg2rad(keepout_deg)
    ang = np.arccos(np.clip((p / np.linalg.norm(p, axis=1, keepdims=True)) @ GONDOLA_DOWN, -1, 1))
    pen = float(np.sum(np.maximum(k - ang, 0)))
    pen += float(np.sum(np.maximum(d @ GONDOLA_DOWN - np.cos(k), 0)))
    return pen


def _downwash(p, d, cone_deg):
    """Final.py downwash_interference_metric."""
    n = p.shape[0]
    if n < 2:
        return 0.0
    ct = np.cos(np.deg2rad(cone_deg))
    r = p[None, :, :] - p[:, None, :]                 # r[i, j] = p_j - p_i
    dist = np.linalg.norm(r, axis=2)
    np.fill_diagonal(dist, np.inf)
    cos_a = np.einsum("ik,ijk->ij", d, r) / dist
    return float(np.sum(np.where(cos_a > ct, (cos_a - ct) / (1.0 - ct), 0.0)))


def geometric(p_norm, d, hull, rules):
    """Final.py geometric_constraints_metric on normalised positions."""
    n = p_norm.shape[0]
    if n < 2:
        return 0.0
    d_min = rules.d_min(hull)
    i, j = np.triu_indices(n, 1)
    dij = np.linalg.norm(p_norm[i] - p_norm[j], axis=1)
    ang = np.arccos(np.clip(np.sum(d[i] * d[j], axis=1), -1, 1))
    pen = float(np.sum(np.maximum(d_min - dij, 0)) + np.sum(np.maximum(np.deg2rad(rules.ang_min_deg) - ang, 0)))
    pen += _gondola_soft_penalty(p_norm, d, rules.gondola_keepout_deg)
    pen += _downwash(p_norm, d, rules.downwash_cone_deg)
    return pen


def evaluate_objectives(positions_m, directions, hull, K, dof_weights, rules):
    """The five objectives plus authority and condition number of the scaled M."""
    p = np.asarray(positions_m, dtype=float)
    d = np.asarray(directions, dtype=float)
    Ms = scale_torque_rows(build_M(p, d), rules.length_m)
    S = np.linalg.svd(Ms, compute_uv=False)
    A = np.linalg.pinv(Ms, rcond=1e-9)
    pn = p / rules.length_m
    return {
        "count": p.shape[0] / (K / 2.0),
        "effort": float(np.sum(np.asarray(dof_weights) * np.abs(A).max(axis=0))),
        "pos_sym": float(np.linalg.norm(pn.mean(axis=0))),
        "dir_sym": float(np.linalg.norm(d.mean(axis=0))),
        "geometric": geometric(pn, d, hull, rules),
        "authority": float(S.sum()),
        "planar_authority": float(np.linalg.svd(Ms[[0, 1, 5]], compute_uv=False).sum()),
        "cond": float(S.max() / max(S.min(), 1e-9)),
        "rank": int(np.sum(S > 1e-6)),
    }


def conditioning_penalty(cond, threshold=100.0):
    """Final.py's conditioning term."""
    return 0.01 * cond / threshold if cond <= threshold else 0.5 * np.log10(cond / threshold)
