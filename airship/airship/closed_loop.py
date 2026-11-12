"""Closed-loop episodes: controller -> allocation -> 6-DOF simulator.

The same function runs PID now and the RL policy later, so every controller
is scored identically.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .allocation import Allocator
from .dynamics import Simulator, euler_from_quat
from .wind import WindConfig, WindModel


@dataclass
class Scenario:
    name: str = "station"
    duration_s: float = 60.0
    wind: WindConfig = field(default_factory=WindConfig)
    wind_seed: int = 0
    init_pos: tuple = (0.0, 0.0, 0.0)
    init_rpy: tuple = (0.0, 0.0, 0.0)
    yaw_ref: float = 0.0
    yaw_step_time_s: float | None = None      # heading command changes once (tests yaw authority)
    yaw_step_to: float = 0.0
    disturbance_torque: tuple = (0.0, 0.0, 0.0)
    fail_prop: int | None = None
    fail_time_s: float = 20.0

    @classmethod
    def from_config(cls, d):
        d = dict(d)
        d["wind"] = WindConfig.from_config(d.get("wind"))
        for k in ("init_pos", "init_rpy", "disturbance_torque"):
            if k in d:
                d[k] = tuple(d[k])
        if "init_rpy_deg" in d:
            d["init_rpy"] = tuple(np.deg2rad(d.pop("init_rpy_deg")))
        for k in ("yaw_ref", "yaw_step_to"):
            if k + "_deg" in d:
                d[k] = float(np.deg2rad(d.pop(k + "_deg")))
        return cls(**d)

    def yaw_at(self, t):
        if self.yaw_step_time_s is not None and t >= self.yaw_step_time_s:
            return self.yaw_step_to
        return self.yaw_ref


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def run_episode(layout, params, controller, scenario, physics_hz=100, control_hz=10, alloc_w=None,
                log=False, divergence_m=10.0, max_tilt_deg=60.0):
    """Fly one episode; returns a metrics dict (and per-step logs if log=True)."""
    sim = Simulator(params, layout, WindModel(scenario.wind, scenario.wind_seed, scenario.duration_s + 5),
                    physics_hz, control_hz, tau_dist=scenario.disturbance_torque)
    sim.reset(position=scenario.init_pos, rpy=scenario.init_rpy)
    alloc = Allocator(layout.M(), layout.u_min, layout.u_max, **({"W": alloc_w} if alloc_w is not None else {}))
    controller.reset()
    dt = 1.0 / control_hz
    n = int(round(scenario.duration_s * control_hz))
    pos_err, yaw_err, tilt, effort, sat, gap_frac = (np.zeros(n) for _ in range(6))
    rec = {"t": [], "x": [], "u": [], "wrench_cmd": [], "wrench_ach": [], "wind": []} if log else None
    failed, diverged = False, False
    for k in range(n):
        t = k * dt
        if scenario.fail_prop is not None and not failed and t >= scenario.fail_time_s:
            sim.fail_prop(scenario.fail_prop)
            alloc.set_bounds(sim.u_lo, sim.u_hi)
            failed = True
        yaw_ref = scenario.yaw_at(t)
        tau_cmd = controller.act(sim.x, np.zeros(3), yaw_ref)
        res = alloc.solve(tau_cmd)
        if hasattr(controller, "notify"):
            controller.notify(res)
        x = sim.step(res.u)
        roll, pitch, yaw = euler_from_quat(x[3:7])
        pos_err[k] = np.linalg.norm(x[:3])
        yaw_err[k] = wrap(yaw - yaw_ref)
        tilt[k] = np.degrees(np.arccos(np.clip(np.cos(roll) * np.cos(pitch), -1, 1)))
        effort[k] = np.sum(x[13:] ** 2) * dt
        live = sim.u_hi - sim.u_lo > 0
        sat[k] = np.sum(res.saturated) / max(np.sum(live), 1)
        gap_frac[k] = np.linalg.norm(res.gap) / max(np.linalg.norm(tau_cmd), 1e-6)
        if log:
            rec["t"].append(sim.t); rec["x"].append(x); rec["u"].append(x[13:].copy())
            rec["wrench_cmd"].append(tau_cmd); rec["wrench_ach"].append(res.achieved); rec["wind"].append(sim.last_wind)
        if not np.all(np.isfinite(x)) or pos_err[k] > divergence_m or tilt[k] > max_tilt_deg:
            diverged = True
            n = k + 1
            break
    out = {
        "rms_pos_m": float(np.sqrt(np.mean(pos_err[:n] ** 2))),
        "max_pos_m": float(np.max(pos_err[:n])),
        "rms_yaw_rad": float(np.sqrt(np.mean(yaw_err[:n] ** 2))),
        "max_tilt_deg": float(np.max(tilt[:n])),
        "effort_N2s": float(np.sum(effort[:n])),
        "saturation_pct": float(100 * np.mean(sat[:n])),
        "gap_pct": float(100 * np.mean(gap_frac[:n])),
        "survived": not diverged,
        "steps": n,
    }
    if scenario.fail_prop is not None and not diverged:
        k0 = int(round(scenario.fail_time_s * control_hz))
        within = np.flatnonzero(pos_err[k0:n] > 0.5)
        out["recovery_s"] = float((within[-1] + 1) * dt) if within.size else 0.0
    if log:
        out["log"] = {k: np.array(v) for k, v in rec.items()}
    return out
