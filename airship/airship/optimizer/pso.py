"""Inner PSO for propeller layouts, ported from Final.py's InnerPSO.

Kept from Final.py: the particle encoding (per candidate point: selection
logit, thrust azimuth, thrust polar angle), the sigmoid >= 0.5 selection,
flipping thrust to point outward, the gondola hard gates, the inertia /
cognitive / social update, the paper's five objectives with per-trial
Dirichlet weights, and the conditioning term.

Changed (Phase 1 of the brief):
  - the rank-6 gate and penalty are replaced by the bounded-thrust capacity
    test (feasibility.shortfall), a graded penalty with no cliff;
  - no "force the top-6 propellers" fallback: the capacity test decides
    (with push-only thrust at least 7 propellers are needed);
  - positions in metres; objectives on the torque-scaled matrix (L = hull radius);
  - no random noise before the SVD; every random number comes from the seed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..feasibility import axis_capacities, shortfall
from ..geometry import gondola_direction_violation
from ..wrench import build_M
from .objectives import OBJECTIVES, GeometricRules, conditioning_penalty, evaluate_objectives


@dataclass
class PSOConfig:
    swarm_size: int = 30
    iters: int = 80
    w_start: float = 0.90
    w_end: float = 0.40
    c_cog: float = 1.6
    c_soc: float = 1.6
    v_max_frac: float = 0.2
    sel_min: float = -5.0
    sel_max: float = 5.0
    threshold: float = 0.5
    penalty_inactive: float = 1e6
    penalty_feasibility: float = 1e4
    cond_threshold: float = 100.0
    dof_weights: tuple = (2.0, 2.0, 0.3, 0.3, 0.3, 2.0)

    @classmethod
    def from_config(cls, d):
        d = dict(d or {})
        if "dof_weights" in d:
            d["dof_weights"] = tuple(d["dof_weights"])
        return cls(**d)


@dataclass
class Evaluation:
    score: float
    mask: np.ndarray
    positions_m: np.ndarray
    directions: np.ndarray
    objectives: dict = field(default_factory=dict)
    axis_capacity: np.ndarray | None = None
    shortfall: float = np.inf
    gate: str = ""                       # "" or the reason the layout was rejected outright


def _sph(theta, phi):
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=-1)


class LayoutPSO:
    def __init__(self, hull, candidates_m, weights, cfg, req, u_max, rules, seed):
        self.hull, self.P, self.cfg, self.req, self.rules = hull, np.asarray(candidates_m), cfg, req, rules
        self.K = self.P.shape[0]
        self.normals = hull.normals(self.P)
        self.u_max = float(u_max)
        self.w = np.asarray(weights, dtype=float)
        if self.w.shape != (len(OBJECTIVES),):
            raise ValueError(f"need one weight per objective {OBJECTIVES}")
        self.rng = np.random.default_rng(seed)
        dim = 3 * self.K
        self.lb = np.empty(dim); self.ub = np.empty(dim)
        self.lb[0::3], self.ub[0::3] = cfg.sel_min, cfg.sel_max
        self.lb[1::3], self.ub[1::3] = 0.0, 2 * np.pi
        self.lb[2::3], self.ub[2::3] = 0.0, np.pi
        self.v_max = cfg.v_max_frac * (self.ub - self.lb)
        self.X = self.rng.uniform(self.lb, self.ub, size=(cfg.swarm_size, dim))
        self.V = np.zeros_like(self.X)
        self.pbest = self.X.copy()
        self.pbest_val = np.full(cfg.swarm_size, np.inf)
        self.gbest, self.gbest_val, self.gbest_eval = None, np.inf, None
        self.n_evals = 0

    def decode(self, x):
        mask = 1.0 / (1.0 + np.exp(-x[0::3])) >= self.cfg.threshold
        d = _sph(x[1::3], x[2::3])[mask]
        d[np.sum(d * self.normals[mask], axis=1) < 0] *= -1        # thrust points outward
        return mask, self.P[mask], d

    def evaluate(self, x):
        self.n_evals += 1
        mask, p, d = self.decode(x)
        if p.shape[0] == 0:
            return Evaluation(self.cfg.penalty_inactive, mask, p, d, gate="no propellers")
        if self.hull.gondola_position_violation(p).any():
            return Evaluation(self.cfg.penalty_inactive, mask, p, d, gate="gondola position")
        if gondola_direction_violation(d).any():
            return Evaluation(self.cfg.penalty_inactive, mask, p, d, gate="gondola direction")
        obj = evaluate_objectives(p, d, self.hull, self.K, self.cfg.dof_weights, self.rules)
        score = float(np.dot(self.w, [obj[k] for k in OBJECTIVES]))
        caps = axis_capacities(build_M(p, d), 0.0, self.u_max)
        sf = shortfall(caps, self.req)
        score += self.cfg.penalty_feasibility * sf + conditioning_penalty(obj["cond"], self.cfg.cond_threshold)
        return Evaluation(score, mask, p, d, obj, caps, sf)

    def _update(self, i):
        ev = self.evaluate(self.X[i])
        if ev.score < self.pbest_val[i]:
            self.pbest_val[i], self.pbest[i] = ev.score, self.X[i].copy()
            if ev.score < self.gbest_val:
                self.gbest_val, self.gbest, self.gbest_eval = ev.score, self.X[i].copy(), ev

    def run(self):
        cfg, conv = self.cfg, []
        for it in range(cfg.iters):
            w = cfg.w_start - (cfg.w_start - cfg.w_end) * it / max(cfg.iters - 1, 1)
            for i in range(cfg.swarm_size):
                self._update(i)
            rp, rg = self.rng.random(self.X.shape), self.rng.random(self.X.shape)
            self.V = w * self.V + cfg.c_cog * rp * (self.pbest - self.X) + cfg.c_soc * rg * (self.gbest - self.X)
            self.V = np.clip(self.V, -self.v_max, self.v_max)
            self.X = np.clip(self.X + self.V, self.lb, self.ub)
            conv.append(self.gbest_val)
        for i in range(cfg.swarm_size):
            self._update(i)
        return self.gbest_eval, np.array(conv)


def trial_weights(num_trials, rng):
    """Per-trial objective weights, Dirichlet(0.5) over the five objectives.

    Same distribution as Final.py's _scalarisation_weights: Dirichlet(0.5) over
    nine objectives, masked to the active five and renormalised, is Dirichlet(0.5)
    over the five. (The count term's 1/(K/2) scale is already inside `count`.)
    """
    return rng.dirichlet(np.full(len(OBJECTIVES), 0.5), size=num_trials)


def run_trials(hull, candidates_m, cfg, req, u_max, rules, num_trials, seed, log=print):
    """num_trials PSO runs with seeded weights; returns one dict per trial."""
    ss = np.random.SeedSequence(seed)
    w_rng, *trial_ss = [np.random.default_rng(s) for s in ss.spawn(num_trials + 1)]
    W = trial_weights(num_trials, w_rng)
    out = []
    for t in range(num_trials):
        pso = LayoutPSO(hull, candidates_m, W[t], cfg, req, u_max, rules, seed=trial_ss[t])
        ev, conv = pso.run()
        out.append({"trial": t, "weights": W[t], "eval": ev, "convergence": conv, "n_evals": pso.n_evals})
        log(f"  seed {seed} trial {t + 1}/{num_trials}: N={ev.positions_m.shape[0]} "
            f"score={ev.score:.4f} shortfall={ev.shortfall:.3f}")
    return out
