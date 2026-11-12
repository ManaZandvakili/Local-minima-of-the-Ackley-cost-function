"""Pareto filtering and layout selection (Final.py's pareto_filter and pick rule)."""
from __future__ import annotations

import numpy as np

# Final.py's primary Pareto front: authority (maximised), effort, count
PARETO_KEYS = ("authority", "effort", "count")


def pareto_filter(points):
    """Indices of non-dominated rows (all columns minimised). Same rule as Final.py."""
    pts = np.asarray(points, dtype=float)
    n = pts.shape[0]
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        others = np.delete(np.arange(n), i)
        dom = np.all(pts[others] <= pts[i], axis=1) & np.any(pts[others] < pts[i], axis=1)
        keep[i] = not dom.any()
    return np.flatnonzero(keep)


def pareto_points(objs):
    """Rows of (-authority, effort, count) for a list of objective dicts."""
    return np.array([[-o["authority"], o["effort"], o["count"]] for o in objs])


def knee_index(points):
    """Pareto point closest to the utopia point after min-max normalisation."""
    pts = np.asarray(points, dtype=float)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    return int(np.argmin(np.linalg.norm((pts - lo) / span, axis=1)))


def paper_rule_order(objs, idx):
    """Pareto indices sorted by planar authority, highest first (Final.py's pick order
    for the path-following comparison; the first one with enough matched random
    layouts is used)."""
    return sorted(idx, key=lambda k: -objs[k]["planar_authority"])
