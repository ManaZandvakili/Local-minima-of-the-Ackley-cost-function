"""Videos: the three controllers side by side, same layout, same wind.

For each propeller type (reversible, one-way) and its optimized knee layout:
  <case>_knee_gusty.mp4          60 s in 3 m/s gusty wind (wind toward 120 deg), heading change at 30 s
  <case>_knee_motor_failure.mp4  45 s in 2 m/s wind, motor 0 fails at 15 s
Panels: classic controller (PID) | residual learning | plain learning.
Each panel's title gives its RMS position error. Videos play at 2x speed.

Usage (from the airship/ folder):
    python scripts/make_controller_videos.py
Writes results/videos/.
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.closed_loop import Scenario, run_episode            # noqa: E402
from airship.config import load_config                           # noqa: E402
from airship.controllers.pid import PIDGains                     # noqa: E402
from airship.controllers.rl_policy import RLController           # noqa: E402
from airship.controllers.tuning import make_pid                  # noqa: E402
from airship.dynamics import DynamicsParams                      # noqa: E402
from airship.feasibility import axis_capacities                  # noqa: E402
from airship.layouts import load                                 # noqa: E402

RL = ROOT / "results" / "rl"
CASES = {
    "reversible": {"pid": ROOT / "results" / "pid" / "20261117_153423",
                   "residual": RL / "20261117_155720_bidir_residual_seed0",
                   "plain": RL / "20261117_155720_bidir_plain_seed0",
                   "label": "reversible propellers"},
    "oneway": {"pid": ROOT / "results" / "pid" / "20261117_203431",
               "residual": RL / "20261117_205532_oneway_residual_seed0",
               "plain": RL / "20261117_205532_oneway_plain_seed0",
               "label": "one-way propellers"},
}
GUSTY = {"name": "gusty", "duration_s": 60.0,
         "wind": {"mean_speed_m_s": 3.0, "mean_dir_deg": 120.0, "turbulence_sigma_m_s": [0.5, 0.5, 0.25],
                  "gust_rate_per_min": 2.0},
         "wind_seed": 501, "init_pos": [1.0, -0.6, 0.3], "yaw_step_time_s": 30.0, "yaw_step_to_deg": 45.0}
FAILURE = {"name": "motor_failure", "duration_s": 45.0,
           "wind": {"mean_speed_m_s": 2.0, "mean_dir_deg": 60.0, "turbulence_sigma_m_s": [0.3, 0.3, 0.15]},
           "wind_seed": 777, "fail_prop": 0, "fail_time_s": 15.0}
NAMES = {"pid": "Classic controller (PID)", "residual": "Residual learning", "plain": "Plain learning"}


def make(job):
    from airship.viz.animate import animate
    case, scen_name = job
    c = CASES[case]
    params = DynamicsParams.from_config(load_config("sphere"))
    g = json.load(open(c["pid"] / "gains" / "sphere_knee.json"))
    lay = load(ROOT / g["path"])
    gains = PIDGains(**g["gains"])
    caps = axis_capacities(lay.M(), lay.u_min, lay.u_max)
    sc = Scenario.from_config(GUSTY if scen_name == "gusty" else FAILURE)
    runs, lines = [], []
    for key in ("pid", "residual", "plain"):
        ctrl = (make_pid(gains, params) if key == "pid" else
                RLController.load(c[key], caps, None, base_pid=make_pid(gains, params)))
        r = run_episode(lay, params, ctrl, sc, log=True)
        status = f"RMS {r['rms_pos_m']:.2f} m" if r["survived"] else f"lost control after {r['steps'] / 10:.0f} s"
        runs.append((r["log"], lay, f"{NAMES[key]}\n{status}"))
        lines.append(f"  {NAMES[key]}: {status}")
    out = ROOT / "results" / "videos" / f"{case}_knee_{scen_name}.mp4"
    animate(runs, out)
    return f"{out.relative_to(ROOT)} ({lay.n_props} propellers, {c['label']})\n" + "\n".join(lines)


if __name__ == "__main__":
    (ROOT / "results" / "videos").mkdir(parents=True, exist_ok=True)
    jobs = [(c, s) for c in CASES for s in ("gusty", "motor_failure")]
    with Pool(4) as p:
        for line in p.imap(make, jobs):
            print(line, flush=True)
