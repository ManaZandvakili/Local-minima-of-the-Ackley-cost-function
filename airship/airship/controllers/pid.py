"""Cascaded 6-DOF PID baseline (brief, section 7).

Position loop (world frame): PID with integral action and anti-windup ->
desired force, rotated into the body frame. Attitude loop: quaternion-error PD
on roll/pitch (the gondola pendulum already helps) and PID on yaw -> desired
torque. Gains are accelerations per unit error, scaled by the effective mass
(m + m_a) and inertia, so one gain set means the same thing on any airship.
"""
from __future__ import annotations

from dataclasses import astuple, dataclass, fields

import numpy as np

from ..dynamics import quat_to_rot

GAIN_NAMES = ("kp_xy", "ki_xy", "kd_xy", "kp_z", "ki_z", "kd_z",
              "kp_rp", "kd_rp", "kp_yaw", "ki_yaw", "kd_yaw")


@dataclass
class PIDGains:
    kp_xy: float = 0.4
    ki_xy: float = 0.05
    kd_xy: float = 1.2
    kp_z: float = 0.4
    ki_z: float = 0.05
    kd_z: float = 1.2
    kp_rp: float = 1.0
    kd_rp: float = 1.0
    kp_yaw: float = 0.5
    ki_yaw: float = 0.05
    kd_yaw: float = 1.0

    @classmethod
    def from_vector(cls, v):
        return cls(*[float(x) for x in v])

    def vector(self):
        return np.array(astuple(self))

    def as_dict(self):
        return {f.name: float(getattr(self, f.name)) for f in fields(self)}


def _qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def attitude_error(q, yaw_ref):
    """Small-angle rotation vector (body frame) from the current attitude to level-at-yaw_ref."""
    q_ref = np.array([np.cos(yaw_ref / 2), 0.0, 0.0, np.sin(yaw_ref / 2)])
    qe = _qmul(np.array([q[0], -q[1], -q[2], -q[3]]), q_ref)
    if qe[0] < 0:
        qe = -qe
    s = np.linalg.norm(qe[1:])
    if s < 1e-12:
        return np.zeros(3)
    return 2.0 * np.arctan2(s, qe[0]) * qe[1:] / s


class PIDController:
    def __init__(self, gains, mass_eff_kg, inertia_diag, dt=0.1, i_max_pos=5.0, i_max_yaw=2.0):
        self.g = gains
        self.m = float(mass_eff_kg)
        self.I = np.asarray(inertia_diag, dtype=float)
        self.dt = dt
        self.i_max_pos, self.i_max_yaw = i_max_pos, i_max_yaw
        self.kp = np.array([gains.kp_xy, gains.kp_xy, gains.kp_z])
        self.ki = np.array([gains.ki_xy, gains.ki_xy, gains.ki_z])
        self.kd = np.array([gains.kd_xy, gains.kd_xy, gains.kd_z])
        self.reset()

    def reset(self):
        self.i_pos = np.zeros(3)
        self.i_yaw = 0.0
        self.saturated = False

    def notify(self, alloc):
        """Anti-windup: stop integrating while the allocation cannot deliver the force."""
        f_cmd = alloc.achieved[:3] + alloc.gap[:3]
        self.saturated = np.linalg.norm(alloc.gap[:3]) > 0.05 * np.linalg.norm(f_cmd) + 0.05

    def act(self, x, p_ref, yaw_ref, v_ref=np.zeros(3), a_ref=np.zeros(3), yaw_rate_ref=0.0):
        """Wrench command (body frame) for state x and the reference."""
        R = quat_to_rot(x[3:7])
        e = np.asarray(p_ref) - x[:3]
        v_world = R @ x[7:10]
        if not self.saturated:
            self.i_pos = np.clip(self.i_pos + e * self.dt, -self.i_max_pos, self.i_max_pos)
        f_world = self.m * (self.kp * e + self.ki * self.i_pos + self.kd * (np.asarray(v_ref) - v_world)
                            + np.asarray(a_ref))
        w = x[10:13]
        ea = attitude_error(x[3:7], yaw_ref)
        if not self.saturated:
            self.i_yaw = float(np.clip(self.i_yaw + ea[2] * self.dt, -self.i_max_yaw, self.i_max_yaw))
        g = self.g
        torque = self.I * np.array([g.kp_rp * ea[0] - g.kd_rp * w[0],
                                    g.kp_rp * ea[1] - g.kd_rp * w[1],
                                    g.kp_yaw * ea[2] + g.ki_yaw * self.i_yaw - g.kd_yaw * (w[2] - yaw_rate_ref)])
        return np.concatenate([R.T @ f_world, torque])
