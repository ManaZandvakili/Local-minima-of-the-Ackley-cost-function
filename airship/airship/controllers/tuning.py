"""Per-layout PID gain tuning with PSO, following Final.py's tune_pd_gains
(swarm, iterations, inertia, coefficients, velocity clamp, early stop).

Every layout is tuned with identical bounds, cost, scenario and budget, so
tuning effort is fair. Cost (the paper's J_tune plus a yaw term, because the
6-DOF controller also has yaw gains to tune):

    J = w_pos * RMS position error [m] + w_yaw * RMS yaw error [rad]
        + w_effort * integral of sum(u^2) dt [N^2 s]   (+ a large penalty if it diverges)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..closed_loop import run_episode
from .pid import GAIN_NAMES, PIDController, PIDGains


@dataclass
class TuningConfig:
    lower: tuple
    upper: tuple
    swarm: int = 12
    iters: int = 30
    w: float = 0.7
    c1: float = 1.5
    c2: float = 1.5
    vmax_frac: float = 0.2
    conv_tol: float = 1e-3
    conv_window: int = 20
    w_pos: float = 5.0
    w_yaw: float = 1.0
    w_effort: float = 1e-3
    divergence_penalty: float = 1e3
    seed: int = 17

    @classmethod
    def from_config(cls, d):
        d = dict(d)
        lo = tuple(float(d["bounds"][k][0]) for k in GAIN_NAMES)
        hi = tuple(float(d["bounds"][k][1]) for k in GAIN_NAMES)
        d.pop("bounds")
        return cls(lower=lo, upper=hi, **d)


def make_pid(gains, params, control_hz=10):
    return PIDController(gains, params.mass_kg + params.added_mass_kg,
                         np.diag(params.I_o) + np.diag(params.J_a), dt=1.0 / control_hz)


def episode_cost(res, tc, n_total):
    """Tuning cost of one episode; a diverged run costs more the earlier it failed."""
    c = tc.w_pos * res["rms_pos_m"] + tc.w_yaw * res["rms_yaw_rad"] + tc.w_effort * res["effort_N2s"]
    if not res["survived"]:
        c += tc.divergence_penalty * (2.0 - res["steps"] / n_total)
    return float(c)


def tune(layout, params, scenario, tc, control_hz=10):
    """Returns (best PIDGains, best cost, cost history)."""
    rng = np.random.default_rng(tc.seed)
    n_total = int(round(scenario.duration_s * control_hz))
    lo, hi = np.array(tc.lower), np.array(tc.upper)
    d = lo.size

    def cost(v):
        return episode_cost(run_episode(layout, params, make_pid(PIDGains.from_vector(v), params), scenario),
                            tc, n_total)

    X = rng.uniform(lo, hi, size=(tc.swarm, d))
    V = np.zeros_like(X)
    vmax = tc.vmax_frac * (hi - lo)
    pbest, pval = X.copy(), np.array([cost(x) for x in X])
    g = int(np.argmin(pval))
    gbest, gval = pbest[g].copy(), float(pval[g])
    hist = [gval]
    for _ in range(tc.iters):
        r1, r2 = rng.random((tc.swarm, d)), rng.random((tc.swarm, d))
        V = np.clip(tc.w * V + tc.c1 * r1 * (pbest - X) + tc.c2 * r2 * (gbest - X), -vmax, vmax)
        X = np.clip(X + V, lo, hi)
        for i in range(tc.swarm):
            c = cost(X[i])
            if c < pval[i]:
                pval[i], pbest[i] = c, X[i].copy()
                if c < gval:
                    gval, gbest = c, X[i].copy()
        hist.append(gval)
        if len(hist) > tc.conv_window and max(hist[-tc.conv_window:]) - min(hist[-tc.conv_window:]) < tc.conv_tol:
            break
    return PIDGains.from_vector(gbest), gval, hist
