"""Motor-failure comparison: PID vs one or more RL runs (plain or residual).

For each chosen layout, every propeller in turn fails at 15 s in 2 m/s wind
with light turbulence (45 s episode). Reports survival and RMS position error
after the failure, per controller.

Usage (from the airship/ folder):
    python scripts/compare_failures.py --pid results/pid/<run> \
        --rl results/rl/<run1> results/rl/<run2> --layouts sphere_knee maxN11_sphere_knee
Writes <first rl run>/eval/failures.md and failures.csv.
"""
import argparse
import csv
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

FAIL_T = 15.0


def scenario(i, seed):
    return Scenario.from_config({"name": f"fail{i}", "duration_s": 45.0,
                                 "wind": {"mean_speed_m_s": 2.0, "mean_dir_deg": 60.0,
                                          "turbulence_sigma_m_s": [0.3, 0.3, 0.15]},
                                 "wind_seed": seed, "init_pos": [0.5, -0.3, 0.2],
                                 "fail_prop": i, "fail_time_s": FAIL_T})


def job(args):
    label, gfile, rl_dirs, i = args
    g = json.load(open(gfile))
    params = DynamicsParams.from_config(load_config("sphere"))
    lay = load(ROOT / g["path"])
    gains = PIDGains(**g["gains"])
    caps = axis_capacities(lay.M(), lay.u_min, lay.u_max)
    ctrls = [("PID", make_pid(gains, params))]
    for d in rl_dirs:
        ctrls.append((Path(d).name, RLController.load(d, caps, None, base_pid=make_pid(gains, params))))
    rows = []
    for name, c in ctrls:
        r = run_episode(lay, params, c, scenario(i, 900 + i), log=True)
        x, t = r["log"]["x"], r["log"]["t"]
        after = t >= FAIL_T
        rms_after = float(np.sqrt(np.mean(np.sum(x[after, :3] ** 2, axis=1)))) if after.any() else np.nan
        rows.append({"layout": label, "n_props": lay.n_props, "failed_prop": i, "controller": name,
                     "survived": bool(r["survived"]), "rms_after_m": rms_after,
                     "max_pos_m": float(r["max_pos_m"])})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pid", required=True)
    ap.add_argument("--rl", nargs="+", required=True)
    ap.add_argument("--layouts", nargs="+", default=["sphere_knee", "maxN11_sphere_knee",
                                                     "sphere_starlab_6thruster"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    pid_dir = Path(a.pid).resolve()
    rl_dirs = [str(Path(d).resolve()) for d in a.rl]
    jobs = []
    for lab in a.layouts:
        gfile = pid_dir / "gains" / f"{lab}.json"
        n = load(ROOT / json.load(open(gfile))["path"]).n_props
        jobs += [(lab, str(gfile), rl_dirs, i) for i in range(n)]
    with Pool(a.workers) as p:
        rows = [r for rs in p.map(job, jobs) for r in rs]
    out = Path(rl_dirs[0]) / "eval"
    out.mkdir(exist_ok=True)
    with open(out / "failures.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    names = ["PID"] + [Path(d).name for d in rl_dirs]
    L = ["# Motor-failure comparison\n",
         f"Each propeller in turn fails at {FAIL_T:.0f} s (2 m/s wind, light turbulence, 45 s). "
         f"PID gains `{pid_dir.name}`. RMS position error after the failure: median over the failed "
         "propeller, survivors only.\n",
         "| Layout | N | " + " | ".join(f"{n} survived" for n in names) + " | "
         + " | ".join(f"{n} RMS after [m]" for n in names) + " |",
         "|---|---|" + "---|" * (2 * len(names))]
    for lab in a.layouts:
        rs = [r for r in rows if r["layout"] == lab]
        surv, rms = [], []
        for n in names:
            c = [r for r in rs if r["controller"] == n]
            surv.append(f"{sum(r['survived'] for r in c)}/{len(c)}")
            ok = [r["rms_after_m"] for r in c if r["survived"]]
            rms.append(f"{np.median(ok):.3f}" if ok else "-")
        L.append(f"| {lab} | {rs[0]['n_props']} | " + " | ".join(surv) + " | " + " | ".join(rms) + " |")
    text = "\n".join(L) + "\n"
    (out / "failures.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
