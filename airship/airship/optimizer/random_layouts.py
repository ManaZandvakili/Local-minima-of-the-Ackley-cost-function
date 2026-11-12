"""Matched-N random feasible layouts (Final.py's sample_random_at_n_feasible),
with the rank-6 test replaced by the bounded-thrust capacity test."""
from __future__ import annotations

import numpy as np

from ..feasibility import axis_capacities, shortfall
from ..geometry import gondola_direction_violation
from ..wrench import build_M


def sample_matched(target_n, hull, candidates_m, req, u_max, rules, n_target=20,
                   max_attempts=20000, seed=42):
    """Up to n_target random layouts with exactly target_n propellers.

    Each passes the same hard constraints as optimizer layouts: gondola position
    and direction gates, minimum spacing, and the capacity requirements with
    0 <= u <= u_max. Returns a list of (positions_m, directions) and the number
    of attempts used.
    """
    rng = np.random.default_rng(seed)
    P = np.asarray(candidates_m, dtype=float)
    safe = np.flatnonzero(~hull.gondola_position_violation(P))
    if len(safe) < target_n:
        return [], 0
    cos_t = np.cos(np.deg2rad(rules.gondola_keepout_deg))
    phi_max = np.arccos(-cos_t)                       # never point into the gondola cone
    d_min = rules.d_min(hull)
    out, attempts = [], 0
    while len(out) < n_target and attempts < max_attempts:
        attempts += 1
        idx = rng.choice(safe, size=target_n, replace=False)
        p = P[idx]
        pn = p / rules.length_m
        dist = np.linalg.norm(pn[:, None] - pn[None], axis=2)
        np.fill_diagonal(dist, np.inf)
        if dist.min() < d_min:
            continue
        theta = rng.uniform(0.0, 2 * np.pi, target_n)
        phi = rng.uniform(0.0, phi_max, target_n)
        d = np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=-1)
        d[np.sum(d * hull.normals(p), axis=1) < 0] *= -1
        if gondola_direction_violation(d).any():
            continue
        if shortfall(axis_capacities(build_M(p, d), 0.0, u_max), req) > 0:
            continue
        out.append((p.copy(), d.copy()))
    return out, attempts
