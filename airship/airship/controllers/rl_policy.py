"""Run a trained SB3 policy as a controller in closed_loop.run_episode, so RL and
PID are flown and scored by exactly the same code.

The wrapper rebuilds the environment's observation (same function), applies
the saved VecNormalize statistics, and returns the deterministic action scaled
to a wrench. No sensor noise, delay or randomisation (evaluation).
"""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np

from ..dynamics import quat_to_rot
from ..envs.station_keeping import make_obs


class RLController:
    def __init__(self, model, obs_rms, env_cfg, layout_caps, clip_obs=10.0, epsilon=1e-8, dt=0.1):
        self.model, self.obs_rms, self.env_cfg = model, obs_rms, env_cfg
        self.scale = np.asarray(env_cfg["action_scale"], dtype=float)
        self.caps = np.asarray(layout_caps, dtype=float)
        self.clip, self.eps, self.dt = clip_obs, epsilon, dt
        self.reset()

    @classmethod
    def load(cls, run_dir, layout_caps, checkpoint_steps=None):
        """Final model, or the checkpoint saved after `checkpoint_steps` training steps."""
        import yaml
        from stable_baselines3 import PPO
        run_dir = Path(run_dir)
        cfg = yaml.safe_load(open(run_dir / "config.yaml"))["train_rl"]
        if checkpoint_steps is None:
            model_path, vn_path = run_dir / "model", run_dir / "vecnormalize.pkl"
        else:
            c = run_dir / "checkpoints"
            model_path = c / f"ppo_{checkpoint_steps}_steps"
            vn_path = c / f"ppo_vecnormalize_{checkpoint_steps}_steps.pkl"
        model = PPO.load(model_path, device="cpu")
        vn = pickle.load(open(vn_path, "rb"))
        return cls(model, vn.obs_rms, cfg["env"], layout_caps, clip_obs=vn.clip_obs, epsilon=vn.epsilon,
                   dt=1.0 / cfg["env"]["control_hz"])

    def reset(self):
        self.integral = np.zeros(3)
        self.a_prev = np.zeros(6)
        self.gap = np.zeros(6)
        self.first = True

    def notify(self, alloc):
        self.gap = -alloc.gap                         # achieved - requested, as in the environment

    def act(self, x, p_ref, yaw_ref, **_):
        if not self.first:
            self.integral = np.clip(self.integral + (np.asarray(p_ref) - x[:3]) * self.dt,
                                    -self.env_cfg["integral_clip_m_s"], self.env_cfg["integral_clip_m_s"])
        self.first = False
        obs = make_obs(x[:3] - np.asarray(p_ref), x[7:10], x[10:13], quat_to_rot(x[3:7]), yaw_ref,
                       self.integral, self.a_prev, self.gap, self.caps, self.env_cfg, self.scale)
        o = np.clip((obs - self.obs_rms.mean) / np.sqrt(self.obs_rms.var + self.eps), -self.clip, self.clip)
        a, _ = self.model.predict(o.astype(np.float32), deterministic=True)
        a = np.clip(a, -1.0, 1.0)
        self.a_prev = a
        return a * self.scale
