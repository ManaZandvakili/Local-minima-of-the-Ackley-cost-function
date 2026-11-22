"""Path-following videos on the paper's 50 m circle (3 m/s gusty wind toward 0 deg).

Per propeller type:
  path_<case>_knee_3controllers.mp4   classic | residual learning | plain learning, optimized knee layout
  path_<case>_knee_vs_random.mp4      classic controller: optimized knee vs a typical random layout
                                      (the median of the 20 random layouts with the same number of propellers)
One-way also:
  path_oneway_capped9_vs_random.mp4   classic controller: 9-propeller capped knee vs a typical random 9-propeller layout
One frame per simulated second, 20 frames per second (20x speed).

Usage (from the airship/ folder):  python scripts/make_path_videos.py
"""
import json
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.config import load_config                           # noqa: E402
from airship.controllers.pid import PIDGains                     # noqa: E402
from airship.controllers.rl_policy import RLController           # noqa: E402
from airship.controllers.tuning import make_pid                  # noqa: E402
from airship.dynamics import DynamicsParams                      # noqa: E402
from airship.feasibility import axis_capacities                  # noqa: E402
from airship.layouts import load                                 # noqa: E402
from airship.path_following import CirclePath, run_path_episode  # noqa: E402
from airship.wind import WindConfig                              # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from path_following import CASES                                  # noqa: E402

WIND = WindConfig(mean_speed_m_s=3.0, mean_dir_deg=0.0, turbulence_sigma_m_s=(0.5, 0.5, 0.25), gust_rate_per_min=2.0)
SEED = 7
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
NAMES = {"PID": "Classic controller", "residual": "Residual learning", "plain": "Plain learning"}


def setup(case, label):
    c = CASES[case]
    g = json.load(open(c["pid"] / "gains" / f"{label}.json"))
    lay = load(ROOT / g["path"])
    return lay, PIDGains(**g["gains"])


def controller(case, name, lay, gains, params):
    if name == "PID":
        return make_pid(gains, params)
    caps = axis_capacities(lay.M(), lay.u_min, lay.u_max)
    return RLController.load(CASES[case][name], caps, None, base_pid=make_pid(gains, params))


def fly(case, label, name, params, log=True):
    lay, gains = setup(case, label)
    r = run_path_episode(lay, params, controller(case, name, lay, gains, params), CirclePath(), WIND, SEED, log=log)
    return r, lay


def title(name, label, lay, r):
    what = "optimized" if "random" not in label else "random"
    status = f"RMS error {r['rms_err_m']:.2f} m" if r["survived"] else f"lost the path at t = {r['steps'] / 10:.0f} s"
    return f"{NAMES[name]}, {what} layout ({lay.n_props} propellers)\n{status}"


def typical_random(case, prefix, n, params):
    """The random layout whose classic-controller error is the median of the 20."""
    labels = sorted(f.stem for f in (CASES[case]["pid"] / "gains").glob(f"{prefix}sphere_random_N{n}_*.json"))
    errs = [(fly(case, lab, "PID", params, log=False)[0]["rms_err_m"], lab) for lab in labels]
    errs.sort()
    return errs[len(errs) // 2][1]


def make(job):
    from airship.viz.animate import animate_path
    kind, case = job
    params = DynamicsParams.from_config(load_config("sphere"))
    out = ROOT / "results" / "videos"
    if kind == "3controllers":
        runs = []
        for name, col in (("PID", BLUE), ("residual", ORANGE), ("plain", AQUA)):
            r, lay = fly(case, "sphere_knee", name, params)
            runs.append((r["log"], lay, title(name, "sphere_knee", lay, r), col))
        path = out / f"path_{case}_knee_3controllers.mp4"
    else:
        knee = "maxN11_sphere_knee" if kind == "capped" else "sphere_knee"
        prefix = "maxN11_" if kind == "capped" else ""
        r1, lay1 = fly(case, knee, "PID", params)
        rand = typical_random(case, prefix, lay1.n_props, params)
        r2, lay2 = fly(case, rand, "PID", params)
        runs = [(r1["log"], lay1, title("PID", knee, lay1, r1), BLUE),
                (r2["log"], lay2, title("PID", rand, lay2, r2), ORANGE)]
        path = out / (f"path_{case}_capped{lay1.n_props}_vs_random.mp4" if kind == "capped"
                      else f"path_{case}_knee_vs_random.mp4")
    animate_path(runs, CirclePath(), path)
    return f"{path.relative_to(ROOT)}: " + " | ".join(t.replace(chr(10), ", ") for _, _, t, _ in runs)


if __name__ == "__main__":
    (ROOT / "results" / "videos").mkdir(parents=True, exist_ok=True)
    jobs = [("3controllers", "reversible"), ("3controllers", "oneway"), ("vs_random", "reversible"),
            ("vs_random", "oneway"), ("capped", "oneway")]
    with Pool(len(jobs)) as p:
        for line in p.imap_unordered(make, jobs):
            print(line, flush=True)
