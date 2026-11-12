"""Phase 3: tune the PID baseline for each layout (identical procedure for all),
then fly the acceptance scenario (2 m/s steady wind) and a turbulent check.

Usage (from the airship/ folder):
    python scripts/tune_pid.py --run results/layouts/phase1_20261112_015323 --workers 4
    python scripts/tune_pid.py --layouts path/a.json path/b.json

Outputs go to a new timestamped folder:
    results/pid/<timestamp>/gains/<layout>.json   tuned gains, cost, history
    results/pid/<timestamp>/summary.csv, report.md
"""
import argparse
import csv
import json
import sys
import time
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.closed_loop import Scenario, run_episode            # noqa: E402
from airship.config import load_config                           # noqa: E402
from airship.controllers.tuning import TuningConfig, make_pid, tune  # noqa: E402
from airship.dynamics import DynamicsParams                      # noqa: E402
from airship.layouts import load                                 # noqa: E402

CHECK_SCENARIO = {"name": "turbulent_3ms", "duration_s": 60.0,
                  "wind": {"mean_speed_m_s": 3.0, "mean_dir_deg": 120.0, "turbulence_sigma_m_s": [0.5, 0.5, 0.25],
                           "gust_rate_per_min": 2.0},
                  "wind_seed": 101, "init_pos": [0.0, 0.0, 0.0]}


def layout_set(run_dir):
    """The layouts flown in Phases 3-5: chosen optimized layouts, the minimum-N Pareto
    layout, matched random layouts, and the reference layouts."""
    run_dir = Path(run_dir)
    paths = [run_dir / "sphere_knee.json", run_dir / "sphere_paper_pick.json"]
    pareto = [load(p) for p in sorted(run_dir.glob("sphere_s*_t*.json"))]
    min_n = min(pareto, key=lambda L: (L.n_props, L.objectives.get("J2_effort", 0)))
    paths.append(run_dir / f"{min_n.name}.json")
    paths += sorted(run_dir.glob("sphere_random_N*_*.json"))
    paths += [ROOT / "results" / "layouts" / "sphere_paper_best.json",
              ROOT / "results" / "layouts" / "sphere_starlab_6thruster.json"]
    return [p for p in paths if p.exists()]


def work(args):
    path, out_dir = args
    cfg, pcfg = load_config("sphere"), load_config("pid")
    params = DynamicsParams.from_config(cfg)
    lay = load(path)
    tc = TuningConfig.from_config(pcfg["tuning"])
    t0 = time.time()
    gains, cost, hist = tune(lay, params, Scenario.from_config(pcfg["tuning_scenario"]), tc)
    out = {"layout": lay.name, "path": str(Path(path).relative_to(ROOT)), "n_props": lay.n_props,
           "source": lay.source, "gains": gains.as_dict(), "tuning_cost": cost, "history": hist,
           "minutes": (time.time() - t0) / 60}
    for key, sc in (("acceptance", pcfg["acceptance_scenario"]), ("turbulent", CHECK_SCENARIO)):
        r = run_episode(lay, params, make_pid(gains, params), Scenario.from_config(sc))
        out[key] = {k: v for k, v in r.items() if k != "log"}
    json.dump(out, open(Path(out_dir) / "gains" / f"{lay.name}.json", "w"), indent=2)
    print(f"{lay.name:28s} N={lay.n_props:2d} cost={cost:.3f} "
          f"accept RMS={out['acceptance']['rms_pos_m']:.3f} m  turb RMS={out['turbulent']['rms_pos_m']:.3f} m  "
          f"({out['minutes']:.1f} min)", flush=True)
    return out


def group(name):
    if name.startswith("sphere_random_N"):
        return name.rsplit("_", 1)[0]
    return name


def report(out_dir, rows):
    L = [f"# PID tuning `{out_dir.name}`\n",
         "Same bounds, cost, tuning scenario (1.5 m/s wind + light turbulence, 1 m offset, 45° heading step) "
         "and PSO budget for every layout.\n",
         "| Layout | N | Tuning cost | 2 m/s steady: RMS pos [m] | effort [N²s] | 3 m/s turbulent: RMS pos [m] | "
         "max pos [m] | effort [N²s] | saturation % | survived |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    groups = {}
    for r in rows:
        groups.setdefault(group(r["layout"]), []).append(r)
    for g, rs in groups.items():
        if len(rs) == 1:
            r = rs[0]; a, t = r["acceptance"], r["turbulent"]
            L.append(f"| {g} | {r['n_props']} | {r['tuning_cost']:.3f} | {a['rms_pos_m']:.3f} | {a['effort_N2s']:.0f} | "
                     f"{t['rms_pos_m']:.3f} | {t['max_pos_m']:.2f} | {t['effort_N2s']:.0f} | "
                     f"{t['saturation_pct']:.1f} | {t['survived']} |")
        else:
            def med(f):
                return np.median([f(r) for r in rs]), np.percentile([f(r) for r in rs], 25), np.percentile([f(r) for r in rs], 75)
            cols = [med(lambda r: r["tuning_cost"]), med(lambda r: r["acceptance"]["rms_pos_m"]),
                    med(lambda r: r["acceptance"]["effort_N2s"]), med(lambda r: r["turbulent"]["rms_pos_m"]),
                    med(lambda r: r["turbulent"]["max_pos_m"]), med(lambda r: r["turbulent"]["effort_N2s"]),
                    med(lambda r: r["turbulent"]["saturation_pct"])]
            fmt = [".3f", ".3f", ".0f", ".3f", ".2f", ".0f", ".1f"]
            cells = [f"{m:{f}} [{a:{f}}-{b:{f}}]" for (m, a, b), f in zip(cols, fmt)]
            surv = sum(r["turbulent"]["survived"] for r in rs)
            L.append(f"| {g} (median [IQR] of {len(rs)}) | {rs[0]['n_props']} | " + " | ".join(cells) + f" | {surv}/{len(rs)} |")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", help="phase-1 layout folder (results/layouts/phase1_...)")
    ap.add_argument("--layouts", nargs="*")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    paths = [Path(p) for p in a.layouts] if a.layouts else layout_set(a.run)
    out_dir = ROOT / "results" / "pid" / f"{datetime.now():%Y%m%d_%H%M%S}"
    (out_dir / "gains").mkdir(parents=True)
    yaml.safe_dump({"pid": load_config("pid"), "check_scenario": CHECK_SCENARIO,
                    "layouts": [str(p) for p in paths]}, open(out_dir / "config.yaml", "w"))
    print(f"tuning {len(paths)} layouts -> {out_dir}", flush=True)
    with Pool(a.workers) as pool:
        rows = list(pool.imap(work, [(p, out_dir) for p in paths]))
    with open(out_dir / "summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["layout", "n_props", "tuning_cost", "accept_rms_pos_m", "accept_effort_N2s",
                    "turb_rms_pos_m", "turb_max_pos_m", "turb_effort_N2s", "turb_saturation_pct", "turb_survived"])
        for r in rows:
            a_, t = r["acceptance"], r["turbulent"]
            w.writerow([r["layout"], r["n_props"], r["tuning_cost"], a_["rms_pos_m"], a_["effort_N2s"],
                        t["rms_pos_m"], t["max_pos_m"], t["effort_N2s"], t["saturation_pct"], t["survived"]])
    text = report(out_dir, rows)
    (out_dir / "report.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
