"""6-DOF rigid-body airship model (Fossen-style), integrated with RK4.

Frames: world z up; body frame origin at the hull centre (centre of buoyancy),
z up, gondola at -z.

State: [p (3, world), q (4, w-x-y-z, body->world), v (3, body), w (3, body),
        u (N, propeller thrusts)]

Equations (body frame, generalised velocity nu = [v; w], nu_r = nu - [R^T wind; 0]):

    (M_rb + M_a) nu_dot = tau_thrust + tau_dist + g(eta)
                          - C_rb(nu) nu - C_a(nu_r) nu_r - D(nu_r)

    M_rb = [[m I, -m S(r_g)], [m S(r_g), I_o]]     I_o: inertia about the hull centre
    M_a  = diag(m_a, m_a, m_a, J_a)                 sphere: m_a = 0.5 rho V, J_a ~ 0
    C(nu) nu for a symmetric M: [w x n1; v x n1 + w x n2], n = M nu
    D: quadratic drag 0.5 rho C_D A |v_r| v_r at the hull centre, plus linear and
       quadratic rotational damping
    g: buoyancy rho V g up at the hull centre, weight m g down at the CG (r_g)
    tau_thrust = M_phys u,   motor lag  u_dot = (u_cmd - u) / tau_m

Simplification: the wind is held constant over each physics step and its
acceleration is not fed back as a force (no Froude-Krylov term); the added-mass
Coriolis term uses the air-relative velocity, as for a constant current.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

G = 9.81


def _skew(a):
    return np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])


@njit(cache=True)
def quat_to_rot(q):
    w, x, y, z = q[0], q[1], q[2], q[3]
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


@njit(cache=True)
def _cross(a, b):
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


@njit(cache=True)
def _coriolis(M, nu):
    n = M @ nu
    out = np.empty(6)
    out[:3] = _cross(nu[3:], n[:3])
    out[3:] = _cross(nu[:3], n[:3]) + _cross(nu[3:], n[3:])
    return out


@njit(cache=True)
def derivative(x, u_cmd, wind_w, Minv, Mrb, Ma, weight, buoy, r_g, A, tau_m, k_drag,
               d_lin, d_quad, tau_dist):
    q = x[3:7]
    v = x[7:10]
    w = x[10:13]
    u = x[13:]
    R = quat_to_rot(q)
    nu = x[7:13].copy()
    v_air = R.T @ wind_w
    nu_r = nu.copy()
    nu_r[:3] = v - v_air

    tau = A @ u
    tau[3:] += tau_dist
    # buoyancy at the origin, weight at the CG
    up = R.T @ np.array([0.0, 0.0, 1.0])
    tau[:3] += (buoy - weight) * up
    tau[3:] += _cross(r_g, -weight * up)
    # aerodynamic drag on the air-relative velocity, rotational damping
    sp = np.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2)
    tau[:3] -= k_drag * sp * nu_r[:3]
    tau[3:] -= d_lin * w + d_quad * np.abs(w) * w
    tau -= _coriolis(Mrb, nu) + _coriolis(Ma, nu_r)

    dx = np.empty_like(x)
    dx[:3] = R @ v
    dx[3] = -0.5 * (q[1] * w[0] + q[2] * w[1] + q[3] * w[2])
    dx[4] = 0.5 * (q[0] * w[0] + q[2] * w[2] - q[3] * w[1])
    dx[5] = 0.5 * (q[0] * w[1] + q[3] * w[0] - q[1] * w[2])
    dx[6] = 0.5 * (q[0] * w[2] + q[1] * w[1] - q[2] * w[0])
    dx[7:13] = Minv @ tau
    dx[13:] = (u_cmd - u) / tau_m
    return dx


@njit(cache=True)
def integrate(x, u_cmd, winds, dt, u_lo, u_hi, Minv, Mrb, Ma, weight, buoy, r_g, A, tau_m,
              k_drag, d_lin, d_quad, tau_dist):
    """len(winds) RK4 steps with a fixed thrust command; returns the new state."""
    x = x.copy()
    for k in range(winds.shape[0]):
        wv = winds[k]
        k1 = derivative(x, u_cmd, wv, Minv, Mrb, Ma, weight, buoy, r_g, A, tau_m, k_drag, d_lin, d_quad, tau_dist)
        k2 = derivative(x + 0.5 * dt * k1, u_cmd, wv, Minv, Mrb, Ma, weight, buoy, r_g, A, tau_m, k_drag, d_lin, d_quad, tau_dist)
        k3 = derivative(x + 0.5 * dt * k2, u_cmd, wv, Minv, Mrb, Ma, weight, buoy, r_g, A, tau_m, k_drag, d_lin, d_quad, tau_dist)
        k4 = derivative(x + dt * k3, u_cmd, wv, Minv, Mrb, Ma, weight, buoy, r_g, A, tau_m, k_drag, d_lin, d_quad, tau_dist)
        x = x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        x[3:7] /= np.sqrt(np.sum(x[3:7] ** 2))
        for i in range(u_lo.shape[0]):
            x[13 + i] = min(max(x[13 + i], u_lo[i]), u_hi[i])
    return x


@dataclass
class DynamicsParams:
    mass_kg: float
    buoyancy_kg: float
    r_g: np.ndarray            # CG position [m], body frame
    I_o: np.ndarray            # 3x3 inertia about the hull centre
    added_mass_kg: float
    J_a: np.ndarray            # 3x3 added rotational inertia
    k_drag: float              # 0.5 * rho * C_D * A
    d_lin: np.ndarray          # rotational damping, linear (3)
    d_quad: np.ndarray         # rotational damping, quadratic (3)
    tau_m: float               # motor time constant [s]

    def matrices(self):
        m, S = self.mass_kg, _skew(self.r_g)
        Mrb = np.zeros((6, 6))
        Mrb[:3, :3] = m * np.eye(3); Mrb[:3, 3:] = -m * S; Mrb[3:, :3] = m * S; Mrb[3:, 3:] = self.I_o
        Ma = np.zeros((6, 6))
        Ma[:3, :3] = self.added_mass_kg * np.eye(3); Ma[3:, 3:] = self.J_a
        return Mrb, Ma

    @classmethod
    def from_config(cls, cfg, mass_props=None):
        from .massprops import from_config as mass_from_config
        mp = mass_props or mass_from_config(cfg)
        R = cfg["hull"]["diameter_m"] / 2.0
        rho = cfg["air"]["density_kg_m3"]
        dyn = cfg["dynamics"]
        return cls(mass_kg=mp.mass_kg, buoyancy_kg=mp.buoyancy_kg, r_g=mp.cg_m, I_o=mp.inertia_kgm2,
                   added_mass_kg=mp.added_mass_kg, J_a=mp.added_inertia_kgm2,
                   k_drag=0.5 * rho * dyn["drag_coeff"] * np.pi * R**2,
                   d_lin=np.asarray(dyn["rot_damping_lin"], dtype=float),
                   d_quad=np.asarray(dyn["rot_damping_quad"], dtype=float),
                   tau_m=float(cfg["propulsion"]["motor_time_constant_s"]))


def quat_from_euler(roll, pitch, yaw):
    cr, sr = np.cos(roll / 2), np.sin(roll / 2)
    cp, sp = np.cos(pitch / 2), np.sin(pitch / 2)
    cy, sy = np.cos(yaw / 2), np.sin(yaw / 2)
    return np.array([cr * cp * cy + sr * sp * sy, sr * cp * cy - cr * sp * sy,
                     cr * sp * cy + sr * cp * sy, cr * cp * sy - sr * sp * cy])


def euler_from_quat(q):
    w, x, y, z = q
    roll = np.arctan2(2 * (w * x + y * z), 1 - 2 * (x * x + y * y))
    pitch = np.arcsin(np.clip(2 * (w * y - z * x), -1, 1))
    yaw = np.arctan2(2 * (w * z + x * y), 1 - 2 * (y * y + z * z))
    return roll, pitch, yaw


class Simulator:
    """One airship with one layout in one wind field.

    step(u_cmd) holds the thrust command for one control period (physics_hz /
    control_hz RK4 steps) and returns the new state.
    """

    def __init__(self, params, layout, wind, physics_hz=100, control_hz=10, tau_dist=None):
        self.p, self.layout, self.wind = params, layout, wind
        self.dt = 1.0 / physics_hz
        self.n_sub = int(round(physics_hz / control_hz))
        self.Mrb, self.Ma = params.matrices()
        self.Minv = np.linalg.inv(self.Mrb + self.Ma)
        self.A = np.ascontiguousarray(layout.M())
        self.u_lo = layout.u_min.astype(float)
        self.u_hi = layout.u_max.astype(float)
        self.tau_dist = np.zeros(3) if tau_dist is None else np.asarray(tau_dist, dtype=float)
        self.weight = params.mass_kg * G
        self.buoy = params.buoyancy_kg * G
        self.t = 0.0
        self.x = None
        self.reset()

    @property
    def n_props(self):
        return self.A.shape[1]

    def reset(self, position=(0, 0, 0), rpy=(0, 0, 0), v_body=(0, 0, 0), w_body=(0, 0, 0)):
        self.x = np.concatenate([np.asarray(position, float), quat_from_euler(*rpy),
                                 np.asarray(v_body, float), np.asarray(w_body, float), np.zeros(self.n_props)])
        self.t = 0.0
        return self.x.copy()

    def fail_prop(self, i):
        """Motor i stops: its thrust drops to zero and stays there."""
        self.u_lo[i] = self.u_hi[i] = 0.0
        self.x[13 + i] = 0.0

    def step(self, u_cmd, winds=None):
        u = np.clip(np.asarray(u_cmd, dtype=float), self.u_lo, self.u_hi)
        if winds is None:
            winds = self.wind.sample(self.n_sub, self.dt)
        self.last_wind = winds[-1]
        self.x = integrate(self.x, u, winds, self.dt, self.u_lo, self.u_hi, self.Minv, self.Mrb, self.Ma,
                           self.weight, self.buoy, self.p.r_g, self.A, self.p.tau_m, self.p.k_drag,
                           self.p.d_lin, self.p.d_quad, self.tau_dist)
        self.t += self.n_sub * self.dt
        return self.x.copy()

    # convenient views of the state
    @property
    def position(self):
        return self.x[:3]

    @property
    def quaternion(self):
        return self.x[3:7]

    @property
    def rot(self):
        return quat_to_rot(self.x[3:7])

    @property
    def v_body(self):
        return self.x[7:10]

    @property
    def w_body(self):
        return self.x[10:13]

    @property
    def thrust(self):
        return self.x[13:]

    def energy(self):
        """Kinetic (incl. added mass) + potential energy [J], for still air."""
        nu = self.x[7:13]
        ke = 0.5 * nu @ (self.Mrb + self.Ma) @ nu
        z_cg = self.x[2] + (self.rot @ self.p.r_g)[2]
        return ke + self.weight * z_cg - self.buoy * self.x[2]
