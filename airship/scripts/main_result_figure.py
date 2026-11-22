"""Main result of the journal paper: optimized vs random layouts with the same number of
propellers, flown by the classic controller, residual learning and plain learning.

Reversible propellers, the 8-propeller layout from the optimizer run limited to 11 motors
(the most the airship can lift) vs its 20 random layouts. Two tests in 3 m/s gusty wind:
station keeping (results/rl/<run>/eval/eval.csv) and the paper's 50 m circle
(results/path/<run>/eval.csv). Error = median over the 3 wind directions; random = median
of the 20 layouts with the interquartile range.

Usage (from the airship/ folder):
    python scripts/main_result_figure.py
    python scripts/main_result_figure.py --residual results/rl/<residual run> --plain results/rl/<plain run> \
        --circle results/path/<run>          # after steps 5 and 6 of the README
Writes results/figures/main_result.png / .svg and prints the numbers.
"""
import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
STATION = {"PID": ROOT / "results/rl/20261117_155720_bidir_residual_seed0/eval/eval.csv",
           "residual": ROOT / "results/rl/20261117_155720_bidir_residual_seed0/eval/eval.csv",
           "plain": ROOT / "results/rl/20261117_155720_bidir_plain_seed0/eval/eval.csv"}
CIRCLE = ROOT / "results/path/20261122_053908_reversible/eval.csv"
OPT, RAND = "maxN11_sphere_knee", re.compile(r"^maxN11_sphere_random_N8_\d+$")
CTRL = (("PID", "Classic (PID)"), ("residual", "Residual learning"), ("plain", "Plain learning"))
BLUE, GREY, INK, MUTED = "#2a78d6", "#8c8c8c", "#222222", "#666666"


def station(ctrl):
    rows = [r for r in csv.DictReader(open(STATION[ctrl])) if r["controller"] == ("PID" if ctrl == "PID" else "RL")]
    per = {}
    for r in rows:
        per.setdefault(r["layout"], []).append((float(r["rms_pos_m"]), r["survived"] == "True"))
    return per


def circle(ctrl):
    per = {}
    for r in csv.DictReader(open(CIRCLE)):
        if r["controller"] == ctrl and r["wind_m_s"] == "3.0":
            per.setdefault(r["layout"], []).append((float(r["rms_err_m"]), r["completed"] == "True"))
    return per


def summary(per):
    opt = per[OPT]
    rand = [v for k, v in per.items() if RAND.match(k)]
    med = [np.median([e for e, _ in v]) for v in rand]
    return {"opt": np.median([e for e, _ in opt]), "opt_ok": np.mean([ok for _, ok in opt]),
            "rand": np.median(med), "q1": np.percentile(med, 25), "q3": np.percentile(med, 75),
            "rand_ok": np.mean([ok for v in rand for _, ok in v])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--residual", help="residual learning run with eval/eval.csv (holds the PID rows too)")
    ap.add_argument("--plain", help="plain learning run with eval/eval.csv")
    ap.add_argument("--circle", help="path-following results folder with eval.csv")
    a = ap.parse_args()
    global CIRCLE
    if a.residual:
        STATION["PID"] = STATION["residual"] = Path(a.residual) / "eval" / "eval.csv"
    if a.plain:
        STATION["plain"] = Path(a.plain) / "eval" / "eval.csv"
    if a.circle:
        CIRCLE = Path(a.circle) / "eval.csv"
    data = {"Holding position": {c: summary(station(c)) for c, _ in CTRL},
            "50 m circle": {c: summary(circle(c)) for c, _ in CTRL}}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, (test, d) in zip(axes, data.items()):
        for i, (c, name) in enumerate(CTRL):
            s = d[c]
            ax.errorbar(i + 0.12, s["rand"], yerr=[[s["rand"] - s["q1"]], [s["q3"] - s["rand"]]], fmt="o",
                        color=GREY, ms=8, capsize=4, lw=1.5, label="Random layouts (median, middle half)" if i == 0 else None)
            ax.plot(i - 0.12, s["opt"], "o", color=BLUE, ms=9, label="Optimized layout" if i == 0 else None)
            ax.annotate(f"{s['opt']:.2f}", (i - 0.12, s["opt"]), xytext=(-8, 0), textcoords="offset points",
                        ha="right", va="center", fontsize=9, color=INK)
            ax.annotate(f"{s['rand']:.2f}", (i + 0.12, s["rand"]), xytext=(8, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=9, color=MUTED)
            if test == "50 m circle":
                ax.text(i, -0.12, f"laps finished\n{s['opt_ok']:.0%} vs {s['rand_ok']:.0%}", transform=ax.get_xaxis_transform(),
                        ha="center", va="top", fontsize=8.5, color=MUTED)
            print(f"{test:16s} {name:18s} optimized {s['opt']:.3f} m | random {s['rand']:.3f} m "
                  f"[{s['q1']:.2f}-{s['q3']:.2f}] | finished/survived {s['opt_ok']:.0%} vs {s['rand_ok']:.0%}")
        ax.set_xticks(range(3), [n for _, n in CTRL])
        ax.set_xlim(-0.6, 2.6)
        ax.set_ylim(0, None)
        ax.set_title(f"{test}, 3 m/s gusty wind", fontsize=11, color=INK, loc="left")
        ax.set_ylabel("RMS distance from target [m]", color=MUTED)
        ax.grid(axis="y", color="#e6e6e6"); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle("8 reversible propellers: the optimized layout beats random layouts with every controller,\n"
                 "and both learning controllers beat the classic controller", fontsize=11.5, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = ROOT / "results" / "figures" / "main_result"
    fig.savefig(out.with_suffix(".png"), dpi=200)
    fig.savefig(out.with_suffix(".svg"))
    print(f"-> {out.relative_to(ROOT)}.png/.svg")


if __name__ == "__main__":
    main()
