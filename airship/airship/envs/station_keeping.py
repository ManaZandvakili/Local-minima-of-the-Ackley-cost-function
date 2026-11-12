"""Station-keeping environment (brief, section 8).

Task: hold position (0, 0, 0) and a commanded heading for 60 s; the heading
command changes once per episode at a random time.

Action: 6 numbers in [-1, 1] -> body-frame wrench (scaled per axis) -> the
same bounded allocation as the PID -> propeller thrust commands.

Observation (body frame, normalised; 41 numbers):
    position error (3), linear velocity (3), gravity direction (3),
    yaw error as sin/cos (2), angular velocity (3), clipped integral of the
    position error (3), previous action (6), previous wrench gap achieved -
    requested (6), layout capability vector +-capacity per axis (12).
    The wind is not observed.
"""
from __future__ import annotations

from dataclasses import replace

import gymnasium as gym
import numpy as np

from ..allocation import Allocator
from ..dynamics import Simulator, euler_from_quat, quat_to_rot
from ..wind import WindConfig, WindModel

OBS_DIM = 41


def _wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def make_obs(p, v, w, R, yaw_ref, integral, a_prev, gap, caps, env_cfg, scale):
    """The 41-number observation (shared by the env and the RL controller wrapper)."""
    yaw = np.arctan2(R[1, 0], R[0, 0])
    ey = _wrap(yaw_ref - yaw)
    return np.concatenate([
        R.T @ (-p) / env_cfg["pos_scale_m"], v / env_cfg["vel_scale_m_s"], R.T @ np.array([0.0, 0.0, -1.0]),
        [np.sin(ey), np.cos(ey)], w / env_cfg["rate_scale_rad_s"],
        R.T @ integral / env_cfg["integral_clip_m_s"], a_prev, gap / scale,
        caps / np.repeat(scale, 2)]).astype(np.float32)


class StationKeepingEnv(gym.Env):
    metadata = {"render_modes": []}

    def __init__(self, pool, cfg, params, stage=None, seed=0, eval_layouts=None, fixed=None):
        """pool: LayoutPool; cfg: configs/train_rl.yaml; params: nominal DynamicsParams.

        stage: curriculum stage dict (layouts / wind / randomize / failures).
        eval_layouts / fixed: for evaluation, a fixed list of layouts and a dict of
        fixed episode settings (see _sample_episode).
        """
        super().__init__()
        self.pool, self.cfg, self.p0 = pool, cfg, params
        e = cfg["env"]
        self.scale = np.asarray(e["action_scale"], dtype=float)
        self.dt = 1.0 / e["control_hz"]
        self.n_steps = int(round(e["episode_s"] * e["control_hz"]))
        self.stage = stage or {"layouts": "all", "wind": "full", "randomize": True, "failures": True}
        self.eval_layouts, self.fixed = eval_layouts, fixed or {}
        self.rng = np.random.default_rng(seed)
        self.action_space = gym.spaces.Box(-1.0, 1.0, shape=(6,), dtype=np.float32)
        self.observation_space = gym.spaces.Box(-np.inf, np.inf, shape=(OBS_DIM,), dtype=np.float32)
        self.sim = None

    # ── curriculum hook (called through VecEnv.env_method) ─────────────
    def set_stage(self, stage):
        self.stage = dict(stage)

    # ── episode setup ──────────────────────────────────────────────────
    def _pick_layout(self):
        if "layout" in self.fixed:
            return self.fixed["layout"]
        if self.eval_layouts:
            return self.eval_layouts[self.rng.integers(len(self.eval_layouts))]
        if self.stage["layouts"] == "knee":
            return self.pool.knee
        return self.pool.train[self.rng.integers(len(self.pool.train))]

    def _wind(self):
        if "wind" in self.fixed:
            return self.fixed["wind"], self.fixed.get("wind_seed", 0)
        w, r = self.cfg["wind"], self.rng
        kind = self.stage["wind"]
        if kind == "calm":
            return WindConfig(turbulence_sigma_m_s=(0, 0, 0)), 0
        speed = r.uniform(*w["mean_speed_m_s"])
        cfg = WindConfig(mean_speed_m_s=speed, mean_dir_deg=r.uniform(0, 360), turbulence_sigma_m_s=(0, 0, 0))
        if kind == "full":
            s = r.uniform(*w["turbulence_sigma_m_s"])
            cfg = replace(cfg, turbulence_sigma_m_s=(s, s, s / 2), turbulence_tau_s=w["turbulence_tau_s"],
                          gust_rate_per_min=r.uniform(*w["gust_rate_per_min"]),
                          gust_amplitude_m_s=tuple(w["gust_amplitude_m_s"]),
                          gust_duration_s=tuple(w["gust_duration_s"]))
        return cfg, int(r.integers(2**31))

    def _params(self):
        p = self.p0
        if not self.stage.get("randomize", False) or self.fixed.get("nominal", False):
            return p, np.zeros(3), 0.0
        d, r = self.cfg["randomization"], self.rng
        s = 1 + r.uniform(-d["mass_frac"], d["mass_frac"])
        net = r.uniform(-d["net_lift_kg"], d["net_lift_kg"])
        q = replace(p, mass_kg=p.mass_kg * s, buoyancy_kg=p.mass_kg * s + net, I_o=p.I_o * s,
                    added_mass_kg=p.added_mass_kg * (1 + r.uniform(-d["added_mass_frac"], d["added_mass_frac"])),
                    k_drag=p.k_drag * (1 + r.uniform(-d["drag_frac"], d["drag_frac"])),
                    r_g=p.r_g * (1 + r.uniform(-d["h_cg_frac"], d["h_cg_frac"])),
                    tau_m=r.uniform(*d["motor_tau_s"]))
        dist = r.normal(scale=d["disturbance_torque_Nm"], size=3)
        delay = float(r.random() < d["action_delay_prob"])
        return q, dist, delay

    def reset(self, seed=None, options=None):
        super().reset(seed=seed)
        if seed is not None:
            self.rng = np.random.default_rng(seed)
        e, r = self.cfg["env"], self.rng
        self.layout = self._pick_layout()
        wcfg, wseed = self._wind()
        params, dist, self.delay = self._params()
        self.sim = Simulator(params, self.layout, WindModel(wcfg, wseed, e["episode_s"] + 5),
                             e["physics_hz"], e["control_hz"], tau_dist=dist)
        off = self.fixed.get("init_pos")
        if off is None:
            v = r.normal(size=3); v[2] *= 0.25
            off = v / np.linalg.norm(v) * r.uniform(0, e["init_offset_m"])
        yaw0 = self.fixed.get("init_yaw", r.uniform(-1, 1) * np.deg2rad(e["init_yaw_deg"]))
        self.sim.reset(position=off, rpy=(0.0, 0.0, yaw0))
        self.alloc = Allocator(self.layout.M(), self.layout.u_min, self.layout.u_max)
        self.yaw_ref = 0.0
        self.yaw_change_step = self.fixed.get("yaw_change_step", int(r.integers(100, self.n_steps - 100)))
        self.yaw_new = self.fixed.get("yaw_new", r.uniform(-1, 1) * np.deg2rad(e["heading_change_deg"]))
        fail_p = self.cfg["randomization"]["midflight_failure_prob"] if self.stage.get("failures") else 0.0
        self.fail_step, self.fail_prop = None, None
        if "fail_prop" in self.fixed:
            self.fail_prop, self.fail_step = self.fixed["fail_prop"], self.fixed.get("fail_step", 200)
        elif r.random() < fail_p:
            self.fail_prop, self.fail_step = int(r.integers(self.layout.n_props)), int(r.integers(50, self.n_steps - 50))
        self.caps = self.pool.capability(self.layout)
        self.k = 0
        self.integral = np.zeros(3)
        self.a_prev = np.zeros(6)
        self.a_filt = np.zeros(6)
        self.wrench_prev = np.zeros(6)
        self.gap = np.zeros(6)
        self.stats = {"pos2": 0.0, "yaw2": 0.0, "effort": 0.0, "sat": 0.0, "gap": 0.0}
        self.terms = {}
        return self._obs(), {}

    # ── observation ────────────────────────────────────────────────────
    def _obs(self):
        e, x = self.cfg["env"], self.sim.x
        n = self.cfg["randomization"]["sensor_noise"] if self.stage.get("randomize") and not self.fixed.get("nominal") else None
        p, v, w, q = x[:3].copy(), x[7:10].copy(), x[10:13].copy(), x[3:7]
        R = quat_to_rot(q)
        if n:
            p += self.rng.normal(scale=n["pos_m"], size=3)
            v += self.rng.normal(scale=n["vel_m_s"], size=3)
            w += self.rng.normal(scale=n["rate_rad_s"], size=3)
            dth = self.rng.normal(scale=n["att_rad"], size=3)
            R = R @ np.array([[1, -dth[2], dth[1]], [dth[2], 1, -dth[0]], [-dth[1], dth[0], 1]])
        return make_obs(p, v, w, R, self.yaw_ref, self.integral, self.a_prev, self.gap, self.caps, e, self.scale)

    # ── step ───────────────────────────────────────────────────────────
    def step(self, action):
        e, rw = self.cfg["env"], self.cfg["reward"]
        a = np.clip(np.asarray(action, dtype=float), -1.0, 1.0)
        tau_f = e.get("action_filter_tau_s", 0.0)
        if tau_f > 0:
            alpha = self.dt / (tau_f + self.dt)
            self.a_filt = self.a_filt + alpha * (a - self.a_filt)
            a_used = self.a_filt
        else:
            a_used = a
        wrench = a_used * self.scale
        applied = self.wrench_prev if self.delay else wrench
        self.wrench_prev = wrench
        if self.fail_step is not None and self.k == self.fail_step:
            self.sim.fail_prop(self.fail_prop)
            self.alloc.set_bounds(self.sim.u_lo, self.sim.u_hi)
        if self.k == self.yaw_change_step:
            self.yaw_ref = self.yaw_new
        res = self.alloc.solve(applied)
        x = self.sim.step(res.u)
        self.k += 1
        self.gap = res.achieved - applied
        self.integral = np.clip(self.integral + (-x[:3]) * self.dt, -e["integral_clip_m_s"], e["integral_clip_m_s"])

        roll, pitch, yaw = euler_from_quat(x[3:7])
        ep = np.linalg.norm(x[:3])
        ey = _wrap(self.yaw_ref - yaw)
        u = x[13:]
        n_live = max(np.sum(self.sim.u_hi - self.sim.u_lo > 0), 1)
        umax = np.max(np.abs(self.layout.u_max))
        # every penalty is bounded: q(z) = z^2 / (1 + z^2) is quadratic near 0 and < 1, so with
        # r_alive > sum of the weights, staying alive always pays (no incentive to crash early)
        def q(z2):
            return z2 / (1.0 + z2)
        terms = {
            "pos": -rw["w_p"] * q((ep / rw["pos_norm_m"]) ** 2),
            "yaw": -rw["w_yaw"] * q((ey / rw["yaw_norm_rad"]) ** 2),
            "att": -rw["w_att"] * q((roll ** 2 + pitch ** 2) / rw["att_norm_rad"] ** 2),
            "rate": -rw["w_w"] * q(float(x[10:13] @ x[10:13]) / rw["rate_norm_rad_s"] ** 2),
            "effort": -rw["w_u"] * float(u @ u) / (n_live * umax ** 2),
            "gap": -rw["w_gap"] * q(float(np.sum((self.gap / self.scale) ** 2))),
            "smooth": -rw["w_da"] * q(float(np.sum((a - self.a_prev) ** 2))),
            "alive": rw["r_alive"],
        }
        reward = float(sum(terms.values()))
        for k_, v_ in terms.items():
            self.terms[k_] = self.terms.get(k_, 0.0) + v_
        self.a_prev = a
        s = self.stats
        s["pos2"] += ep ** 2; s["yaw2"] += ey ** 2; s["effort"] += float(u @ u) * self.dt
        s["sat"] += np.sum(res.saturated) / n_live
        s["gap"] += np.linalg.norm(res.gap) / max(np.linalg.norm(applied), 1e-6)

        tilt = np.degrees(np.arccos(np.clip(np.cos(roll) * np.cos(pitch), -1, 1)))
        terminated = bool(ep > e["terminate_pos_m"] or tilt > e["terminate_tilt_deg"] or not np.all(np.isfinite(x)))
        truncated = self.k >= self.n_steps
        info = {}
        if terminated or truncated:
            k = self.k
            info["episode_metrics"] = {
                "rms_pos_m": float(np.sqrt(s["pos2"] / k)), "rms_yaw_rad": float(np.sqrt(s["yaw2"] / k)),
                "effort_N2s": float(s["effort"]), "saturation_pct": float(100 * s["sat"] / k),
                "gap_pct": float(100 * s["gap"] / k),
                "survived": not terminated, "steps": k, "layout": self.layout.name,
            }
            info["reward_terms"] = {k_: v_ / k for k_, v_ in self.terms.items()}
        return self._obs(), reward, terminated, truncated, info
