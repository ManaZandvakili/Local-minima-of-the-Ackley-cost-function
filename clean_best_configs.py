"""
"Optimization results" slide: the PSO best configuration for each hull,
drawn as shapes only (hull, gondola keep-out, propellers, thrust arrows),
with NO grid, axes, tick numbers, labels, titles or legend.  Put the
"Sphere - N = 16" style labels underneath in PowerPoint.

It runs the same optimisation as Final.py (same PSO, seeds, Pareto pick
of the best-authority solution), with the paper's cost function, i.e.
direction symmetry (J6) ACTIVE.  Results are cached in best_config_*.npz,
so re-plotting after the first run is instant.

Usage (Final.py must be in the same folder, or pass --final):
    python clean_best_configs.py                  # all three hulls
    python clean_best_configs.py --shape sphere   # just one (sphere|cigar|oblate)
    python clean_best_configs.py --plot-only      # re-plot from the cache

Outputs (transparent PNGs, 300 dpi):
    best_config_sphere.png, best_config_cigar.png, best_config_oblate.png
    best_configs_all.png   (the three side by side)
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401

K           = 30
NUM_TRIALS  = 20
USE_PAPER_COST = True     # True: J6 (direction symmetry) active, as in the paper
SHAPES = {                # name: (Final.py shape, shape_params, label)
    "sphere": ("sphere",    None,            "Sphere"),
    "cigar":  ("ellipsoid", (1.0, 1.0, 2.5), "Cigar (ellipsoid)"),
    "oblate": ("ellipsoid", (2.0, 2.0, 0.5), "Flat Disk (oblate)"),
}
# Same look as Final.py's plot_config_3d
HULL_COLOR, HULL_ALPHA = "skyblue", 0.10
KO_COLOR,   KO_ALPHA   = "salmon", 0.25
PROP_COLOR, ARROW_COLOR = "red", "darkorange"
VIEW = dict(elev=30, azim=-60)            # matplotlib's default 3-D view


def load_final(path):
    spec = importlib.util.spec_from_file_location("final_module", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["final_module"] = mod
    spec.loader.exec_module(mod)           # main() is not run on import
    if USE_PAPER_COST:
        mod.ACTIVE_OBJ_MASK[5] = True      # J6 back on (paper's 5 objectives)
    return mod


def optimise(mod, key):
    shape, params, label = SHAPES[key]
    if shape == "sphere":
        cand = mod.fibonacci_sphere(K, radius=1.0)
    else:
        cand = mod.fibonacci_ellipsoid(K, *params)
    weights = np.ones(9); weights[0] = 2.0; weights /= weights.sum()   # as in main()

    print(f"\n=== {label}: running {NUM_TRIALS} PSO trials ===")
    (_, all_best, _, _, _, best_pos_list, best_dirs_list, _, _) = \
        mod.run_multiple_trials(cand, weights, num_trials=NUM_TRIALS,
                                shape=shape, shape_params=params)
    pts1, idx1 = mod.pareto_filter(all_best[:, [0, 3, 2]])          # as in _run_geometry
    best_k = idx1[int(np.argmin(pts1[:, 0]))]
    pos, dirs = best_pos_list[best_k], best_dirs_list[best_k]
    print(f"    >>> {label}: best-authority solution has {len(pos)} props")
    np.savez(f"best_config_{key}.npz", positions=pos, directions=dirs)
    return pos, dirs


def draw(ax, key, pos, dirs):
    shape, params, _ = SHAPES[key]
    a, b, c = (1.0, 1.0, 1.0) if shape == "sphere" else params
    arrow_scale = 0.30 if shape == "sphere" else 0.55

    u = np.linspace(0, 2*np.pi, 60); v = np.linspace(0, np.pi, 30)
    ax.plot_surface(a*np.outer(np.cos(u), np.sin(v)), b*np.outer(np.sin(u), np.sin(v)),
                    c*np.outer(np.ones_like(u), np.cos(v)),
                    alpha=HULL_ALPHA, color=HULL_COLOR, linewidth=0)
    phi_ko = np.deg2rad(35.0)
    u2 = np.linspace(0, 2*np.pi, 40); p2 = np.linspace(np.pi - phi_ko, np.pi, 12)
    ax.plot_surface(a*np.outer(np.cos(u2), np.sin(p2)), b*np.outer(np.sin(u2), np.sin(p2)),
                    c*np.outer(np.ones_like(u2), np.cos(p2)),
                    alpha=KO_ALPHA, color=KO_COLOR, linewidth=0)
    ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=90, c=PROP_COLOR, zorder=5)
    for p, d in zip(pos, dirs):
        ax.quiver(p[0], p[1], p[2], *(d*arrow_scale), color=ARROW_COLOR,
                  linewidth=1.8, arrow_length_ratio=0.35)

    ax.set_box_aspect((a, b, c))
    ax.set_xlim(-a*1.15, a*1.15); ax.set_ylim(-b*1.15, b*1.15); ax.set_zlim(-c*1.15, c*1.15)
    ax.view_init(**VIEW)
    ax.set_axis_off()                       # no panes, grid, ticks or labels


def save_single(key, pos, dirs):
    fig = plt.figure(figsize=(7, 7))
    draw(fig.add_subplot(111, projection="3d"), key, pos, dirs)
    fig.savefig(f"best_config_{key}.png", dpi=300, bbox_inches="tight",
                pad_inches=0.0, transparent=True)
    plt.close(fig)
    print(f"Saved best_config_{key}.png  ({len(pos)} props)")


def save_all(results):
    keys = [k for k in SHAPES if k in results]
    fig = plt.figure(figsize=(6*len(keys), 6.5))
    for i, k in enumerate(keys):
        draw(fig.add_subplot(1, len(keys), i + 1, projection="3d"), k, *results[k])
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)
    fig.savefig("best_configs_all.png", dpi=300, bbox_inches="tight",
                pad_inches=0.0, transparent=True)
    print("Saved best_configs_all.png  ("
          + ", ".join(f"{k}: {len(results[k][0])} props" for k in keys) + ")")
    plt.show()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", default="Final.py", help="path to Final.py")
    ap.add_argument("--shape", choices=list(SHAPES), help="run only this hull")
    ap.add_argument("--plot-only", action="store_true", help="re-plot from the cache")
    args = ap.parse_args()

    keys = [args.shape] if args.shape else list(SHAPES)
    mod = None
    results = {}
    for k in keys:
        cache = f"best_config_{k}.npz"
        if os.path.exists(cache):
            d = np.load(cache, allow_pickle=True)
            results[k] = (d["positions"], d["directions"])
            print(f"Loaded {cache}  ({len(results[k][0])} props)")
        elif args.plot_only:
            print(f"[skip] {cache} not found - run without --plot-only first")
            continue
        else:
            if mod is None:
                mod = load_final(args.final)
            results[k] = optimise(mod, k)
        save_single(k, *results[k])

    # combined figure uses whatever caches exist for all three
    for k in SHAPES:
        if k not in results and os.path.exists(f"best_config_{k}.npz"):
            d = np.load(f"best_config_{k}.npz", allow_pickle=True)
            results[k] = (d["positions"], d["directions"])
    if results:
        save_all(results)


if __name__ == "__main__":
    main()
