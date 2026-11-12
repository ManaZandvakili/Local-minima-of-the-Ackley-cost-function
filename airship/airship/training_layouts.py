"""The set of layouts the RL agent trains on, and the held-out set (brief, section 8).

Training: Pareto layouts from the optimizer runs, matched random layouts,
reference and hand-designed layouts, small manufacturing perturbations of
each, and single-propeller-failure variants (kept only if the layout can still
push both ways in x and y and yaw both ways). Held out: about 25 % of the
random layouts and one Pareto layout, to test generalisation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .feasibility import axis_capacities
from .layouts import load, perturb, remove_prop

ROOT = Path(__file__).resolve().parent.parent
LAYOUT_DIR = ROOT / "results" / "layouts"
PLANAR_AXES = [0, 1, 2, 3, 10, 11]          # +-Fx, +-Fy, +-Mz in AXIS_LABELS order


@dataclass
class LayoutPool:
    train: list = field(default_factory=list)     # Layout objects
    holdout: list = field(default_factory=list)
    knee: object = None
    caps: dict = field(default_factory=dict)      # name -> 12 axis capacities (physical units)

    def capability(self, layout):
        if layout.name not in self.caps:
            self.caps[layout.name] = axis_capacities(layout.M(), layout.u_min, layout.u_max)
        return self.caps[layout.name]

    def summary(self):
        kinds = {}
        for L in self.train:
            kinds[L.source] = kinds.get(L.source, 0) + 1
        return {"train": len(self.train), "holdout": len(self.holdout), "train_by_source": kinds,
                "holdout_names": [L.name for L in self.holdout]}


def _unique(L, prefix):
    """Copy of a layout with a run-prefixed name (names repeat across runs)."""
    from dataclasses import replace
    return replace(L, name=f"{prefix}:{L.name}")


def build_pool(cfg):
    """Build the training / held-out pool from configs/train_rl.yaml['layouts']."""
    lc = cfg["layouts"]
    rng = np.random.default_rng(lc["split_seed"])
    pool = LayoutPool()
    base_train = []

    pareto = []
    for run in lc["pareto_runs"]:
        pareto += [_unique(load(p), run) for p in sorted((LAYOUT_DIR / run).glob("sphere_s*_t*.json"))]
    knee_src = load(LAYOUT_DIR / lc["pareto_runs"][0] / "sphere_knee.json")
    pool.knee = _unique(knee_src, lc["pareto_runs"][0])
    knee_name = f"{lc['pareto_runs'][0]}:{knee_src.notes.get('same_as', '')}"
    candidates = [i for i, L in enumerate(pareto) if L.name != knee_name]
    held = set(rng.choice(candidates, size=min(lc["holdout_pareto"], len(candidates)), replace=False).tolist())
    for i, L in enumerate(pareto):
        (pool.holdout if i in held else base_train).append(L)

    for run in lc["random_runs"]:
        rand = [_unique(load(p), run) for p in sorted((LAYOUT_DIR / run).glob("sphere_random_N*_*.json"))]
        idx = rng.permutation(len(rand))
        n_hold = int(round(lc["holdout_random_frac"] * len(rand)))
        for k, i in enumerate(idx):
            (pool.holdout if k < n_hold else base_train).append(rand[i])

    for name in lc.get("include_references", []):
        base_train.append(load(LAYOUT_DIR / f"{name}.json"))
    base_train += [load(p) for p in sorted(LAYOUT_DIR.glob(lc.get("hand_designed_glob", "_none_")))]

    pool.train = list(base_train)
    failures = []
    for L in base_train:
        for k in range(lc.get("perturbations_per_layout", 0)):
            P = perturb(L, lc["perturb_pos_sigma_m"], lc["perturb_ang_sigma_deg"], rng, name=f"{L.name}_pert{k}")
            if np.all(pool.capability(P)[PLANAR_AXES] > 1e-6):
                pool.train.append(P)
        if lc.get("failure_variants", False):
            for i in range(L.n_props):
                failures.append((L, i))
    # single-failure variants: a random subset (otherwise they would dominate the pool)
    n_fail = min(lc.get("max_failure_variants", len(failures)), len(failures))
    for j in rng.permutation(len(failures)):
        if n_fail <= 0:
            break
        F = remove_prop(*failures[j])
        if np.all(pool.capability(F)[PLANAR_AXES] > 1e-6):
            pool.train.append(F)
            n_fail -= 1
    for L in pool.train + pool.holdout + [pool.knee]:
        pool.capability(L)
    return pool
