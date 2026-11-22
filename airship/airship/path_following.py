"""Circular path following (the conference paper's validation, secondary task in the brief).

The reference is the paper's circle (Final.py, _path_reference): radius 50 m,
0.5 m/s, starting at the origin heading along +x, the nose pointing along the
path. Same simulator, allocation and controllers as station keeping; the
controller is given the moving reference (position, velocity, acceleration,
heading and heading rate).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .allocation import Allocator
from .dynamics import Simulator, euler_from_quat
from .wind import WindModel


@dataclass
class CirclePath:
    radius_m: float = 50.0
    speed_m_s: float = 0.5
    duration_s: float = 650.0          # one lap is 2*pi*50/0.5 = 628 s

    @property
    def omega(self):
        return self.speed_m_s / self.radius_m

    def ref(self, t):
        """Reference position, velocity, acceleration (world), heading and heading rate at time t."""
        R, V, w = self.radius_m, self.speed_m_s, self.omega
        s, c = np.sin(w * t), np.cos(w * t)
        p = np.array([R * s, R * (1.0 - c), 0.0])
        v = np.array([V * c, V * s, 0.0])
        a = np.array([-V * w * s, V * w * c, 0.0])
        return p, v, a, w * t, w


def _wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def run_path_episode(layout, params, controller, path, wind_cfg, wind_seed, physics_hz=100, control_hz=10,
                     log=False, divergence_m=10.0, max_tilt_deg=60.0):
    """Fly the circle once; returns metrics (and per-second logs if log=True).

    Error = distance from the reference point at the same time (as in the paper)."""
    sim = Simulator(params, layout, WindModel(wind_cfg, wind_seed, path.duration_s + 5), physics_hz, control_hz)
    sim.reset()
    alloc = Allocator(layout.M(), layout.u_min, layout.u_max)
    controller.reset()
    dt = 1.0 / control_hz
    n = int(round(path.duration_s * control_hz))
    err, yaw_err, tilt, effort, sat = (np.zeros(n) for _ in range(5))
    rec = {"t": [], "x": [], "ref": [], "u": []} if log else None
    diverged = False
    for k in range(n):
        t = k * dt
        p_ref, v_ref, a_ref, yaw_ref, yaw_rate = path.ref(t)
        tau = controller.act(sim.x, p_ref, yaw_ref, v_ref=v_ref, a_ref=a_ref, yaw_rate_ref=yaw_rate)
        res = alloc.solve(tau)
        if hasattr(controller, "notify"):
            controller.notify(res)
        x = sim.step(res.u)
        p_next = path.ref(t + dt)[0]
        roll, pitch, yaw = euler_from_quat(x[3:7])
        err[k] = np.linalg.norm(x[:3] - p_next)
        yaw_err[k] = _wrap(yaw - path.ref(t + dt)[3])
        tilt[k] = np.degrees(np.arccos(np.clip(np.cos(roll) * np.cos(pitch), -1, 1)))
        effort[k] = np.sum(x[13:] ** 2) * dt
        live = sim.u_hi - sim.u_lo > 0
        sat[k] = np.sum(res.saturated) / max(np.sum(live), 1)
        if log and k % control_hz == control_hz - 1:
            rec["t"].append(t + dt); rec["x"].append(x.copy()); rec["ref"].append(p_next); rec["u"].append(x[13:].copy())
        if not np.all(np.isfinite(x)) or err[k] > divergence_m or tilt[k] > max_tilt_deg:
            diverged = True
            n = k + 1
            break
    out = {"rms_err_m": float(np.sqrt(np.mean(err[:n] ** 2))), "max_err_m": float(np.max(err[:n])),
           "rms_yaw_rad": float(np.sqrt(np.mean(yaw_err[:n] ** 2))), "effort_N2s": float(np.sum(effort[:n])),
           "saturation_pct": float(100 * np.mean(sat[:n])), "survived": not diverged,
           "completed": (not diverged) and float(np.max(err[:n])) < 3.0, "steps": n}
    if log:
        out["log"] = {k: np.array(v) for k, v in rec.items()}
        out["log"]["err_1hz"] = err[control_hz - 1:n:control_hz].copy()
    return out
