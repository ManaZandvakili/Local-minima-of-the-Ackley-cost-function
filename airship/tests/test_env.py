"""Phase 4: environment checks, and that the RL controller wrapper flies exactly
like the environment (so RL is scored by the same code as PID)."""
from types import SimpleNamespace

import numpy as np
import pytest

from airship.closed_loop import Scenario, run_episode
from airship.config import load_config
from airship.dynamics import DynamicsParams
from airship.envs.station_keeping import OBS_DIM, StationKeepingEnv
from airship.layouts import load
from airship.training_layouts import LAYOUT_DIR, LayoutPool
from airship.wind import WindConfig

CFG = load_config("train_rl")
PARAMS = DynamicsParams.from_config(load_config("sphere"))
LAY = load(LAYOUT_DIR / "sphere_paper_best.json")


def _pool():
    pool = LayoutPool(train=[LAY], holdout=[], knee=LAY)
    pool.capability(LAY)
    return pool


class _Linear:
    """A fixed random linear policy standing in for a trained network."""
    def __init__(self, seed=0):
        self.W = np.random.default_rng(seed).normal(scale=0.3, size=(6, OBS_DIM))

    def predict(self, obs, deterministic=True):
        return np.tanh(self.W @ obs), None


def test_env_api_and_termination():
    env = StationKeepingEnv(_pool(), CFG, PARAMS, seed=0)
    obs, _ = env.reset(seed=3)
    assert obs.shape == (OBS_DIM,) and np.all(np.isfinite(obs))
    total, done, steps = 0.0, False, 0
    while not done:
        obs, r, term, trunc, info = env.step(np.zeros(6))
        total += r; steps += 1; done = term or trunc
    assert "episode_metrics" in info and "reward_terms" in info
    assert steps <= 600


def test_rl_wrapper_matches_env():
    from airship.controllers.rl_policy import RLController
    wind = WindConfig(mean_speed_m_s=2.0, mean_dir_deg=40.0, turbulence_sigma_m_s=(0.3, 0.3, 0.15))
    fixed = {"layout": LAY, "wind": wind, "wind_seed": 5, "init_pos": np.array([1.0, -0.5, 0.2]),
             "init_yaw": 0.1, "yaw_change_step": 300, "yaw_new": 0.5, "nominal": True}
    stage = {"layouts": "all", "wind": "full", "randomize": False, "failures": False}
    env = StationKeepingEnv(_pool(), CFG, PARAMS, stage=stage, seed=0, fixed=fixed)
    pol = _Linear()
    obs, _ = env.reset()
    xs = []
    for _ in range(100):
        obs, *_ = env.step(pol.predict(obs)[0])
        xs.append(env.sim.x[:3].copy())
    sc = Scenario(duration_s=10.0, wind=wind, wind_seed=5, init_pos=(1.0, -0.5, 0.2), init_rpy=(0, 0, 0.1),
                  yaw_ref=0.0, yaw_step_time_s=30.0, yaw_step_to=0.5)
    rms = SimpleNamespace(mean=np.zeros(OBS_DIM), var=np.ones(OBS_DIM) - 1e-8)
    ctrl = RLController(pol, rms, CFG["env"], _pool().capability(LAY))
    res = run_episode(LAY, PARAMS, ctrl, sc, log=True)
    n = len(res["log"]["x"])                  # the random policy may crash the airship early
    assert n >= 40
    np.testing.assert_allclose(res["log"]["x"][:, :3], np.array(xs)[:n], atol=1e-9)
