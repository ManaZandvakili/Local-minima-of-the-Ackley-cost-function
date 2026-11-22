"""Circular path following (the conference paper's test) with every layout and controller.

The paper's circle: radius 50 m, 0.5 m/s, 650 s (about one lap), nose along the
path. Wind: the paper's 1.5 m/s and the harder 3 m/s used for station keeping,
each from 3 directions, with turbulence (0.5 m/s) and gusts (2 per minute).
Controllers: the tuned classic controller (PID), residual learning and plain
learning (trained for station keeping; they see their error relative to the
moving target, so they are used as they are).

Usage (from the airship/ folder):
    python scripts/path_following.py --case reversible
    python scripts/path_following.py --case oneway
    python scripts/path_following.py --case reversible --pid results/pid/<run> \
        --residual results/rl/<residual run> --plain results/rl/<plain run>     # your own runs
Writes results/path/<timestamp>_<case>/: eval.csv, report.md, path_<layout>.png
(optimized knee and 9/8-propeller capped knee vs their 20 random layouts).
To redo only the figures:  python scripts/path_following.py --case oneway --figures-only results/path/<folder>
"""
import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

# one thread per worker: the learning controllers are 7x faster than when the workers fight over cores
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")

from airship.config import load_config                           # noqa: E402
from airship.controllers.pid import PIDGains                     # noqa: E402
from airship.controllers.rl_policy import RLController           # noqa: E402
from airship.controllers.tuning import make_pid                  # noqa: E402
from airship.dynamics import DynamicsParams                      # noqa: E402
from airship.feasibility import axis_capacities                  # noqa: E402
from airship.layouts import load                                 # noqa: E402
from airship.path_following import CirclePath, run_path_episode  # noqa: E402
from airship.wind import WindConfig                              # noqa: E402

RL = ROOT / "results" / "rl"
CASES = {
    "reversible": {"pid": ROOT / "results" / "pid" / "20261117_153423",
                   "residual": RL / "20261117_155720_bidir_residual_seed0",
                   "plain": RL / "20261117_155720_bidir_plain_seed0", "label": "reversible propellers"},
    "oneway": {"pid": ROOT / "results" / "pid" / "20261117_203431",
               "residual": RL / "20261117_205532_oneway_residual_seed0",
               "plain": RL / "20261117_205532_oneway_plain_seed0", "label": "one-way propellers"},
}
WINDS = (1.5, 3.0)
DIRECTIONS = (0.0, 120.0, 240.0)
SEED = 7                                    # the paper's path seed
CONTROLLERS = ("PID", "residual", "plain")
OPTIMIZED = ("sphere_knee", "sphere_paper_pick", "maxN11_sphere_knee")


def wind(speed, direction):
    return WindConfig(mean_speed_m_s=speed, mean_dir_deg=direction, turbulence_sigma_m_s=(0.5, 0.5, 0.25),
                      gust_rate_per_min=2.0)


def group(label):
    m = re.match(r"(.*sphere_random_N\d+)_\d+$", label)
    return m.group(1) if m else label


def job(args):
    c, gfile = args
    g = json.load(open(gfile))
    params = DynamicsParams.from_config(load_config("sphere"))
    lay = load(ROOT / g["path"])
    gains = PIDGains(**g["gains"])
    caps = axis_capacities(lay.M(), lay.u_min, lay.u_max)
    rows = []
    for name in CONTROLLERS:
        for ws in WINDS:
            for d in DIRECTIONS:
                ctrl = (make_pid(gains, params) if name == "PID" else
                        RLController.load(c[name], caps, None, base_pid=make_pid(gains, params)))
                r = run_path_episode(lay, params, ctrl, CirclePath(), wind(ws, d), SEED)
                rows.append({"layout": g["layout"], "n_props": lay.n_props, "controller": name, "wind_m_s": ws,
                             "wind_dir_deg": d, **{k: r[k] for k in ("rms_err_m", "max_err_m", "rms_yaw_rad",
                                                                     "effort_N2s", "saturation_pct", "survived",
                                                                     "completed")}})
    return rows


def summarize(rows):
    """{(wind, controller, layout): dict(n, rms, effort, done)} medians over directions."""
    per = {}
    for r in rows:
        per.setdefault((r["wind_m_s"], r["controller"], r["layout"]), []).append(r)
    return {k: {"n": v[0]["n_props"], "rms": float(np.median([x["rms_err_m"] for x in v])),
                "effort": float(np.median([x["effort_N2s"] for x in v])),
                "done": float(np.mean([x["completed"] for x in v]))} for k, v in per.items()}


def pairs(summ, ws, ctrl):
    out = []
    for opt in OPTIMIZED:
        if (ws, ctrl, opt) not in summ:
            continue
        o = summ[(ws, ctrl, opt)]
        prefix = "maxN11_" if opt.startswith("maxN11_") else ""
        pat = re.compile(rf"^{prefix}sphere_random_N{o['n']}_\d+$")
        rand = [s for (w, c, lab), s in summ.items() if w == ws and c == ctrl and pat.match(lab)]
        if rand:
            out.append((o["n"], opt, o, rand))
    return sorted(out)


def report(case, summ, out):
    c = CASES[case]
    L = [f"# Path following: {c['label']}\n",
         "The paper's circle: radius 50 m, 0.5 m/s, 650 s (about one lap), nose along the path. Wind from 3 "
         "directions with turbulence (0.5 m/s) and gusts. Error = distance from the reference point; values = "
         "median over the 3 directions, then median over the 20 random layouts. \"Completed\" = error stayed "
         "under 3 m for the whole lap (the paper's rule), share of flights.\n",
         f"PID gains `{c['pid'].name}`; learning runs `{c['residual'].name}`, `{c['plain'].name}`.\n"]
    for ws in WINDS:
        L += [f"\n## Wind {ws} m/s\n", "| Layout | N | PID | Residual | Plain | Completed (PID / residual / plain) |",
              "|---|---|---|---|---|---|"]
        for n, opt, _, _ in pairs(summ, ws, "PID"):
            for kind in ("optimized", "random"):
                cells, done = [], []
                for ctrl in CONTROLLERS:
                    pr = {p[0]: p for p in pairs(summ, ws, ctrl)}.get(n)
                    if kind == "optimized":
                        cells.append(f"{pr[2]['rms']:.3f} m"); done.append(f"{pr[2]['done']:.0%}")
                    else:
                        v = [s["rms"] for s in pr[3]]
                        cells.append(f"{np.median(v):.3f} m [{np.percentile(v, 25):.2f}-{np.percentile(v, 75):.2f}]")
                        done.append(f"{np.mean([s['done'] for s in pr[3]]):.0%}")
                name = f"{opt} (optimized)" if kind == "optimized" else f"random, {n} propellers (20)"
                L.append(f"| {name} | {n} | " + " | ".join(cells) + " | " + " / ".join(done) + " |")
        for ref in ("sphere_paper_best", "sphere_starlab_6thruster"):
            if (ws, "PID", ref) in summ:
                v = [summ[(ws, ctrl, ref)] for ctrl in CONTROLLERS]
                L.append(f"| {ref} | {v[0]['n']} | " + " | ".join(f"{s['rms']:.3f} m" for s in v) + " | "
                         + " / ".join(f"{s['done']:.0%}" for s in v) + " |")
    (out / "report.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def fly_logged(args):
    """3 m/s toward 0 deg with logs, every controller, one layout (for the figures)."""
    c, gfile = args
    g = json.load(open(gfile))
    params = DynamicsParams.from_config(load_config("sphere"))
    lay = load(ROOT / g["path"])
    gains = PIDGains(**g["gains"])
    caps = axis_capacities(lay.M(), lay.u_min, lay.u_max)
    cv = {"n": lay.n_props}
    for name in CONTROLLERS:
        ctrl = (make_pid(gains, params) if name == "PID" else
                RLController.load(c[name], caps, None, base_pid=make_pid(gains, params)))
        r = run_path_episode(lay, params, ctrl, CirclePath(), wind(3.0, 0.0), SEED, log=True)
        cv[name] = {"err": r["log"]["err_1hz"], "xy": r["log"]["x"][:, :2], "survived": r["survived"]}
    return g["layout"], cv


def figure(case, opt, curves, out):
    """Error around the circle and top view: one optimized layout vs the 20 random ones with the same N."""
    import matplotlib.pyplot as plt
    knee = curves[opt]
    n = knee["n"]
    prefix = "maxN11_" if opt.startswith("maxN11_") else ""
    rand = [cv for lab, cv in curves.items() if re.match(rf"^{prefix}sphere_random_N{n}_\d+$", lab)]
    path = CirclePath()
    deg = lambda m: np.degrees(path.omega * np.arange(1, m + 1))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6))
    top = 0
    for j, ctrl in enumerate(CONTROLLERS):
        ax, ax2 = axes[0, j], axes[1, j]
        lost = 0
        for i, cv in enumerate(rand):
            e = cv[ctrl]["err"]
            ax.plot(deg(len(e)), e, color="0.55", lw=0.6, alpha=0.45,
                    label=f"random {n}-propeller layouts (20)" if i == 0 else None)
            if not cv[ctrl]["survived"]:
                lost += 1
                ax.plot(deg(len(e))[-1], min(e[-1], 10), "x", color="#eb6834", ms=6,
                        label="random layout lost the path" if lost == 1 else None)
            top = max(top, float(np.max(e)))
        e = knee[ctrl]["err"]
        ax.plot(deg(len(e)), e, color="#2a78d6", lw=1.5, label=f"optimized ({opt})")
        if not knee[ctrl]["survived"]:
            ax.plot(deg(len(e))[-1], min(e[-1], 10), "X", color="#2a78d6", ms=9)
        name = {"PID": "Classic controller", "residual": "Residual learning", "plain": "Plain learning"}[ctrl]
        ax.set_title(f"{name}  (random layouts lost: {lost} of {len(rand)})", fontsize=11)
        ax.set_xlabel("Position around the circle [deg]"); ax.set_xlim(0, 360); ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel("Tracking error [m]"); ax.legend(fontsize=8, loc="upper right")
        th = np.linspace(0, 2 * np.pi, 400)
        ax2.plot(path.radius_m * np.sin(th), path.radius_m * (1 - np.cos(th)), "k--", lw=0.8, label="reference")
        worst = max(rand, key=lambda cv: (not cv[ctrl]["survived"], np.sqrt(np.mean(cv[ctrl]["err"] ** 2))))
        ax2.plot(*worst[ctrl]["xy"].T, color="#eb6834", lw=1.0, label="worst random layout")
        ax2.plot(*knee[ctrl]["xy"].T, color="#2a78d6", lw=1.2, label="optimized")
        ax2.set_aspect("equal"); ax2.grid(alpha=0.3); ax2.set_xlabel("x [m]")
        if j == 0:
            ax2.set_ylabel("y [m]"); ax2.legend(fontsize=8, loc="center")
    for ax in axes[0]:
        ax.set_ylim(0, min(top, 10) * 1.05)
    fig.suptitle(f"Path following, {CASES[case]['label']}, {n} propellers: 50 m circle at 0.5 m/s, "
                 f"3 m/s gusty wind (toward 0 deg). x = lost the path (more than 10 m off or tipped over)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out / f"path_{opt}.png", dpi=160)
    plt.close(fig)


def figures(case, out, workers):
    c = CASES[case]
    files = {f.stem: f for f in (c["pid"] / "gains").glob("*.json")}
    want = [lab for lab in ("sphere_knee", "maxN11_sphere_knee") if lab in files]
    ns = {lab: load(ROOT / json.load(open(files[lab]))["path"]).n_props for lab in want}
    for lab in list(ns):
        prefix = "maxN11_" if lab.startswith("maxN11_") else ""
        want += [r for r in files if re.match(rf"^{prefix}sphere_random_N{ns[lab]}_\d+$", r)]
    with Pool(workers) as p:
        curves = dict(p.map(fly_logged, [(c, str(files[lab])) for lab in want]))
    for opt in ("sphere_knee", "maxN11_sphere_knee"):
        if opt in curves:
            figure(case, opt, curves, out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=list(CASES), required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--figures-only", metavar="DIR", help="redo only the figures in an existing results folder")
    ap.add_argument("--pid", help="tuned PID folder (results/pid/<run>); default: the run listed in CASES")
    ap.add_argument("--residual", help="residual learning run (results/rl/<run>); default: CASES")
    ap.add_argument("--plain", help="plain learning run (results/rl/<run>); default: CASES")
    a = ap.parse_args()
    for k in ("pid", "residual", "plain"):
        if getattr(a, k):
            CASES[a.case][k] = Path(getattr(a, k)).resolve()
    c = CASES[a.case]
    if a.figures_only:
        figures(a.case, Path(a.figures_only), a.workers)
        return
    out = ROOT / "results" / "path" / f"{datetime.now():%Y%m%d_%H%M%S}_{a.case}"
    out.mkdir(parents=True)
    jobs = [(c, str(f)) for f in sorted((c["pid"] / "gains").glob("*.json"))]
    with Pool(a.workers) as p:
        res = p.map(job, jobs)
    rows = [r for rs in res for r in rs]
    with open(out / "eval.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    report(a.case, summarize(rows), out)
    figures(a.case, out, a.workers)
    print(f"-> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
