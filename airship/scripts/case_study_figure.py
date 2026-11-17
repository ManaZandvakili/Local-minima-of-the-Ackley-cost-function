"""Case-study figure: optimized vs matched random layouts, one-way vs reversible propellers.

Reads the eval.csv of one RL run per case (it holds the PID rows too, flown on the
same wind seeds) and plots, for each propeller count, the optimized layout next to
the median of its 20 matched random layouts (whiskers: interquartile range).
Rows: controllers (PID, RL). Columns: cases. Two figures: RMS position error and
control effort (never two scales on one axis).

Usage (from the airship/ folder):
    python scripts/case_study_figure.py \
        --case "One-way propellers=results/rl/<oneway residual run>" \
        --case "Reversible propellers=results/rl/<reversible residual run>" \
        --rl-name "Residual RL" --out results/figures/case_study
"""
import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402
import matplotlib.ticker                                           # noqa: E402,F401
import numpy as np                                                # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OPTIMIZED = ("maxN11_sphere_knee", "sphere_knee", "sphere_paper_pick")
C_OPT, C_RAND = "#2a78d6", "#eb6834"          # reference palette slots 1-2 (validated pair)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def per_layout(rl_dir):
    """{(controller, layout): dict(n, rms, effort, surv)}: median over wind directions."""
    rows = list(csv.DictReader(open(Path(rl_dir) / "eval" / "eval.csv")))
    out = {}
    for r in rows:
        out.setdefault((r["controller"], r["layout"]), []).append(r)
    return {k: {"n": int(v[0]["n_props"]),
                "rms": float(np.median([float(x["rms_pos_m"]) for x in v])),
                "effort": float(np.median([float(x["effort_N2s"]) for x in v])),
                "surv": float(np.mean([x["survived"] == "True" for x in v]))} for k, v in out.items()}


def pairs(data, ctrl):
    """[(N, optimized stats, [random stats])] sorted by N."""
    res = []
    for opt in OPTIMIZED:
        if (ctrl, opt) not in data:
            continue
        o = data[(ctrl, opt)]
        prefix = "maxN11_" if opt.startswith("maxN11_") else ""
        pat = re.compile(rf"^{prefix}sphere_random_N{o['n']}_\d+$")
        rand = [s for (c, lab), s in data.items() if c == ctrl and pat.match(lab)]
        if rand:
            res.append((o["n"], o, rand))
    return sorted(res, key=lambda t: t[0])


def panel(ax, prs, key, log):
    """Log scale: dots (bars on a log axis mislead). Linear scale: bars from zero."""
    w = 0.34 if not log else 0.18
    labels = []
    for i, (n, o, rand) in enumerate(prs):
        v = np.array([s[key] for s in rand])
        med, lo, hi = np.median(v), np.percentile(v, 25), np.percentile(v, 75)
        xo, xr = i - w / 2 - 0.01, i + w / 2 + 0.01
        if log:
            ax.plot([xo], [o[key]], "o", ms=9, color=C_OPT, mec="white", mew=1.5, zorder=4)
            ax.plot([xr], [med], "o", ms=9, color=C_RAND, mec="white", mew=1.5, zorder=4)
        else:
            ax.bar(xo, o[key], w, color=C_OPT, zorder=3)
            ax.bar(xr, med, w, color=C_RAND, zorder=3)
        ax.errorbar(xr, med, yerr=[[med - lo], [hi - med]], fmt="none", ecolor=INK2,
                    elinewidth=1.2, capsize=3, zorder=3)
        surv = np.mean([s["surv"] for s in rand])
        labels.append(f"N = {n}" + (f"\n{surv:.0%} of random\nkept control" if surv < 0.999 else ""))
    ax.set_xticks(np.arange(len(prs)), labels)
    ax.set_xlim(-0.6, len(prs) - 0.4)
    if log:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda y, _: f"{y:g}"))
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_yticks([0.2, 0.5, 1, 2, 5])
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)


def figure(cases, key, ylabel, log, rl_name, path):
    ctrls = [("PID", "PID"), ("RL", rl_name)]
    fig, axes = plt.subplots(len(ctrls), len(cases), figsize=(4.2 * len(cases), 3.2 * len(ctrls)),
                             sharey=True, squeeze=False)
    for j, (cname, data) in enumerate(cases):
        for i, (ctrl, cl) in enumerate(ctrls):
            ax = axes[i, j]
            panel(ax, pairs(data, ctrl), key, log)
            if i == 0:
                ax.set_title(cname, fontsize=11, color=INK, loc="left")
            if j == 0:
                ax.set_ylabel(f"{cl}\n{ylabel}", fontsize=10, color=INK)
    if log:
        handles = [plt.Line2D([], [], ls="none", marker="o", ms=9, color=c) for c in (C_OPT, C_RAND)]
    else:
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (C_OPT, C_RAND)]
    fig.legend(handles, ["Optimized layout", "Random layouts, same N (median, IQR)"], loc="upper center",
               ncol=2, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "svg"):
        fig.savefig(f"{path}.{ext}", dpi=200, facecolor="white")
    plt.close(fig)
    print(f"{path}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", action="append", required=True, help='"Title=results/rl/<run>"')
    ap.add_argument("--rl-name", default="RL")
    ap.add_argument("--out", default="results/figures/case_study")
    a = ap.parse_args()
    cases = []
    for c in a.case:
        title, d = c.split("=", 1)
        cases.append((title, per_layout(ROOT / d)))
    out = ROOT / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    figure(cases, "rms", "RMS position error [m]", True, a.rl_name, f"{out}_rms")
    figure(cases, "effort", "Control effort [N²s]", False, a.rl_name, f"{out}_effort")


if __name__ == "__main__":
    main()
