"""Phase 4: train the station-keeping policy with PPO (Stable-Baselines3).

One policy flies every layout: it outputs a wrench, and the layout only lives
inside the allocation. Curriculum (configs/train_rl.yaml): calm air with the
knee layout -> calm air, all layouts -> steady wind -> turbulence, gusts,
failures and domain randomisation. Each reward term is logged separately to
TensorBoard so reward hacking is visible.

Usage (from the airship/ folder):
    python scripts/train_rl.py --seed 0                    # full budget from the config
    python scripts/train_rl.py --seed 0 --steps 200000     # quick check

Outputs: results/rl/<timestamp>_seed<k>/ with checkpoints, the final model,
VecNormalize statistics, the config, curriculum log and TensorBoard logs.
"""
import argparse
import json
import pickle
import sys
import time
from collections import deque
from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import torch                                                        # noqa: E402
from stable_baselines3 import PPO                                   # noqa: E402
from stable_baselines3.common.callbacks import BaseCallback, CheckpointCallback  # noqa: E402
from stable_baselines3.common.vec_env import SubprocVecEnv, VecMonitor, VecNormalize  # noqa: E402

from airship.config import load_config                              # noqa: E402
from airship.dynamics import DynamicsParams                         # noqa: E402
from airship.envs.station_keeping import StationKeepingEnv          # noqa: E402
from airship.training_layouts import build_pool                     # noqa: E402


def make_env(pool, cfg, params, stage, seed):
    def _f():
        return StationKeepingEnv(pool, cfg, params, stage=stage, seed=seed)
    return _f


class CurriculumAndTerms(BaseCallback):
    """Logs episode metrics and reward terms; advances the curriculum."""

    def __init__(self, stages, out_dir, window=50, min_steps_per_stage=100_000):
        super().__init__()
        self.stages, self.out_dir = stages, Path(out_dir)
        self.i, self.window, self.min_steps = 0, window, min_steps_per_stage
        self.stage_start = 0
        self.metrics, self.terms = deque(maxlen=window), deque(maxlen=window)
        self.log = []

    def _on_step(self):
        for info in self.locals["infos"]:
            if "episode_metrics" in info:
                self.metrics.append(info["episode_metrics"]); self.terms.append(info["reward_terms"])
        return True

    def _on_rollout_end(self):
        if not self.metrics:
            return
        m = list(self.metrics)
        rms = float(np.mean([x["rms_pos_m"] for x in m]))
        surv = float(np.mean([x["survived"] for x in m]))
        self.logger.record("episode/rms_pos_m", rms)
        self.logger.record("episode/survival", surv)
        self.logger.record("episode/rms_yaw_rad", float(np.mean([x["rms_yaw_rad"] for x in m])))
        self.logger.record("episode/saturation_pct", float(np.mean([x["saturation_pct"] for x in m])))
        for k in self.terms[0]:
            self.logger.record(f"reward_terms/{k}", float(np.mean([t[k] for t in self.terms])))
        self.logger.record("curriculum/stage", self.i)
        st = self.stages[self.i]
        thr = st.get("advance_rms_m")
        self.log.append({"steps": self.num_timesteps, "stage": st["name"], "rms_pos_m": rms, "survival": surv})
        if (thr is not None and self.i < len(self.stages) - 1 and len(m) >= self.window
                and self.num_timesteps - self.stage_start >= self.min_steps and rms < thr and surv > 0.95):
            self.i += 1
            self.stage_start = self.num_timesteps
            self.metrics.clear(); self.terms.clear()
            self.training_env.env_method("set_stage", self.stages[self.i])
            print(f"[curriculum] step {self.num_timesteps}: -> stage {self.i} {self.stages[self.i]['name']}", flush=True)
        json.dump({"stage": self.i, "log": self.log}, open(self.out_dir / "curriculum.json", "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--steps", type=int, help="override total steps")
    ap.add_argument("--n-envs", type=int)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    cfg = load_config("train_rl")
    pc = cfg["ppo"]
    n_envs = a.n_envs or pc["n_envs"]
    total = a.steps or pc["total_steps"]
    out = ROOT / "results" / "rl" / f"{datetime.now():%Y%m%d_%H%M%S}{a.tag}_seed{a.seed}"
    out.mkdir(parents=True)
    torch.set_num_threads(1)

    params = DynamicsParams.from_config(load_config("sphere"))
    t0 = time.time()
    pool = build_pool(cfg)
    print(f"layout pool: {pool.summary()['train']} train, {pool.summary()['holdout']} held out "
          f"({time.time() - t0:.0f} s)", flush=True)
    pickle.dump(pool, open(out / "pool.pkl", "wb"))
    yaml.safe_dump({"train_rl": cfg, "sphere": load_config("sphere"), "seed": a.seed, "n_envs": n_envs,
                    "total_steps": total, "pool": pool.summary()}, open(out / "config.yaml", "w"))

    stages = cfg["curriculum"]
    venv = SubprocVecEnv([make_env(pool, cfg, params, stages[0], a.seed * 1000 + i) for i in range(n_envs)])
    venv = VecNormalize(VecMonitor(venv), norm_obs=True, norm_reward=True, gamma=pc["gamma"])
    act = {"tanh": torch.nn.Tanh, "relu": torch.nn.ReLU}[pc["activation"]]
    model = PPO("MlpPolicy", venv, n_steps=pc["n_steps"], batch_size=pc["batch_size"], n_epochs=pc["n_epochs"],
                learning_rate=pc["learning_rate"], gamma=pc["gamma"], gae_lambda=pc["gae_lambda"],
                clip_range=pc["clip_range"], ent_coef=pc["ent_coef"], seed=a.seed, verbose=0,
                policy_kwargs={"net_arch": list(pc["net_arch"]), "activation_fn": act},
                tensorboard_log=str(out / "tb"), device="cpu")
    cb = [CurriculumAndTerms(stages, out),
          CheckpointCallback(save_freq=max(500_000 // n_envs, 1), save_path=str(out / "checkpoints"),
                             name_prefix="ppo", save_vecnormalize=True)]
    t0 = time.time()
    model.learn(total_timesteps=total, callback=cb, progress_bar=False)
    model.save(out / "model")
    venv.save(str(out / "vecnormalize.pkl"))
    print(f"done: {total} steps in {(time.time() - t0) / 60:.1f} min -> {out}", flush=True)


if __name__ == "__main__":
    main()
