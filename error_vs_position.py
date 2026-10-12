"""
Closed-loop path-following slide: tracking error as a function of the
position around the circle (deg), PSO-optimised layout vs. the matched
feasible random layouts.

It repeats Final.py's validation exactly (same PSO trials and seeds, same
Pareto choice of the PSO layout, same matched-N random population, same
PD-gain tuning, same 8 wind directions and gust seeds), with the paper's
cost function (J6 direction symmetry ACTIVE), and records the error at
every time step instead of only the RMS.

Position around the circle: theta = (V / R) * t, in degrees (0 deg = start).
For each layout the error is averaged over the 8 wind directions.

    PSO            : solid line
    Random median  : dashed line (the configuration Final.py calls the median)
    Random IQR     : shaded band, 25th-75th percentile of the 20 random layouts

Usage (Final.py in the same folder, or pass --final):
    python error_vs_position.py                  # all three hulls (slow, ~1 h each)
    python error_vs_position.py --shape sphere   # one hull
    python error_vs_position.py --plot-only      # re-plot from the cached results

Outputs:  error_vs_position_<shape>.png  and  error_vs_position_all.png
Cache:    error_vs_position_<shape>.npz
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

K, NUM_TRIALS, N_ANGLES = 30, 20, 8
N_RAND_TARGET, N_RAND_MIN = 20, 15
SHAPES = {                 # key: (Final.py shape, params, label)
    "sphere": ("sphere",    None,            "Sphere"),
    "cigar":  ("ellipsoid", (1.0, 1.0, 2.5), "Cigar"),
    "oblate": ("ellipsoid", (2.0, 2.0, 0.5), "Oblate"),
}
C_PSO, C_RAND, C_BAND = "#1f4e9c", "#d62728", "#b9c0c7"
ONE_LAP = True             # True: x-axis 0-360 deg (first lap only)


# ── Run Final.py's pipeline and keep the error time series ─────────────
def load_final(path):
    spec = importlib.util.spec_from_file_location("final_module", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["final_module"] = mod
    spec.loader.exec_module(mod)
    mod.ACTIVE_OBJ_MASK[5] = True          # paper cost: J6 active
    return mod


def run_shape(mod, key):
    shape, params, label = SHAPES[key]
    cand = (mod.fibonacci_sphere(K, radius=1.0) if shape == "sphere"
            else mod.fibonacci_ellipsoid(K, *params))
    weights = np.ones(9); weights[0] = 2.0; weights /= weights.sum()

    print(f"\n=== {label}: {NUM_TRIALS} PSO trials ===", flush=True)
    (_, all_best, _, _, _, best_pos_list, best_dirs_list, _, _) = \
        mod.run_multiple_trials(cand, weights, num_trials=NUM_TRIALS,
                                shape=shape, shape_params=params)
    _, idx1 = mod.pareto_filter(all_best[:, [0, 3, 2]])

    def planar_auth(p, d):
        M = mod.build_M_matrix(p, d)
        return float(np.linalg.svd(M[[0, 1, 5], :], full_matrices=False)[1].sum())

    ranked = sorted(((k, len(best_pos_list[k]),
                      planar_auth(best_pos_list[k], best_dirs_list[k])) for k in idx1),
                    key=lambda x: -x[2])
    chosen, rand_samples = None, []
    for k, n_k, _ in ranked:
        c = mod.sample_random_at_n_feasible(n_k, cand, shape, params,
                                            n_target=N_RAND_TARGET, seed=42 + n_k)
        if len(c) >= N_RAND_MIN:
            chosen, rand_samples = k, c
            break
    if chosen is None:
        raise RuntimeError("could not build a matched random population")
    pso_pos, pso_dirs = best_pos_list[chosen], best_dirs_list[chosen]
    n = len(pso_pos)
    print(f"    matched N = {n}; tuning PD gains for PSO + {len(rand_samples)} random",
          flush=True)

    seeds = [mod.SIM_PATH_SEED + i*1009 for i in range(N_ANGLES)]
    pso_res = mod.simulate_multi_wind(pso_pos, pso_dirs, n_angles=N_ANGLES,
                                      gains=mod.tune_pd_gains(pso_pos, pso_dirs),
                                      per_angle_seeds=seeds)
    rand_res = [mod.simulate_multi_wind(rp, rd, n_angles=N_ANGLES,
                                        gains=mod.tune_pd_gains(rp, rd),
                                        per_angle_seeds=seeds)
                for rp, rd in rand_samples]

    curves = lambda res: np.array([r["error_curve"] for r in res["runs"]])  # (8, T)
    r_rms = np.array([r["rms_mean"] for r in rand_res])
    med = int(np.argsort(r_rms)[len(r_rms)//2])           # same as Final.py
    t = np.arange(curves(pso_res).shape[1]) * mod.SIM_DT
    theta = np.degrees(mod.SIM_CRUISE_V / mod.SIM_PATH_R * t)

    out = dict(theta=theta, n=n,
               pso=curves(pso_res),                                   # (8, T)
               rand_all=np.array([curves(r).mean(0) for r in rand_res]),  # (20, T)
               rand_med=curves(rand_res[med]),                        # (8, T)
               pso_rms=pso_res["rms_mean"], rand_med_rms=rand_res[med]["rms_mean"])
    np.savez(f"error_vs_position_{key}.npz", **out)
    print(f"    {label}: N={n}  RMS  PSO={out['pso_rms']:.3f} m  "
          f"random median={out['rand_med_rms']:.3f} m", flush=True)
    return out


# ── Plot ──────────────────────────────────────────────────────────────
def draw(ax, key, d, show_ylabel=True, show_legend=True):
    th = d["theta"]
    keep = th <= 360.0 if ONE_LAP else np.ones_like(th, bool)
    th = th[keep]
    pso = d["pso"][:, keep].mean(0)
    med = d["rand_med"][:, keep].mean(0)
    q25, q75 = np.percentile(d["rand_all"][:, keep], [25, 75], axis=0)

    ax.fill_between(th, q25, q75, color=C_BAND, alpha=0.55, lw=0,
                    label="Random layouts (25–75%)")
    ax.plot(th, med, color=C_RAND, lw=2.0, ls="--",
            label=f"Median random  (RMS {float(d['rand_med_rms']):.3f} m)")
    ax.plot(th, pso, color=C_PSO, lw=2.4,
            label=f"PSO optimized  (RMS {float(d['pso_rms']):.3f} m)")

    ax.set_xlim(0, th.max())
    ax.set_xticks(np.arange(0, th.max() + 1, 45))
    # scale to the steady part; the start-up transient (first ~10 deg) may clip
    steady = th > 15.0
    ax.set_ylim(0, 1.25*max(pso[steady].max(), med[steady].max(), q75[steady].max()))
    ax.set_xlabel("Position around the circle  [deg]", fontsize=13)
    if show_ylabel:
        ax.set_ylabel("Tracking error  [m]", fontsize=13)
    ax.set_title(f"{SHAPES[key][2]}  (N = {int(d['n'])})", fontsize=15,
                 fontweight="bold")
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=11)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if show_legend:
        ax.legend(fontsize=10.5, frameon=False, loc="upper right")


def save_plots(results):
    for k, d in results.items():
        fig, ax = plt.subplots(figsize=(8, 4.8))
        draw(ax, k, d)
        fig.tight_layout()
        fig.savefig(f"error_vs_position_{k}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved error_vs_position_{k}.png")

    keys = [k for k in SHAPES if k in results]
    fig, axes = plt.subplots(1, len(keys), figsize=(6.2*len(keys), 4.8),
                             sharey=True, squeeze=False)
    for i, k in enumerate(keys):
        draw(axes[0, i], k, results[k], show_ylabel=(i == 0))
    fig.text(0.5, -0.02, "Error averaged over 8 wind directions;  "
             "random = 20 feasible layouts with the same N",
             ha="center", fontsize=11, color="#555555")
    fig.tight_layout()
    fig.savefig("error_vs_position_all.png", dpi=300, bbox_inches="tight")
    print("Saved error_vs_position_all.png")
    plt.show()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", default="Final.py")
    ap.add_argument("--shape", choices=list(SHAPES))
    ap.add_argument("--plot-only", action="store_true")
    args = ap.parse_args()

    keys = [args.shape] if args.shape else list(SHAPES)
    mod, results = None, {}
    for k in keys:
        cache = f"error_vs_position_{k}.npz"
        if os.path.exists(cache):
            results[k] = dict(np.load(cache))
            print(f"Loaded {cache}")
        elif not args.plot_only:
            mod = mod or load_final(args.final)
            results[k] = run_shape(mod, k)
    for k in SHAPES:                     # include other cached hulls in the combined plot
        if k not in results and os.path.exists(f"error_vs_position_{k}.npz"):
            results[k] = dict(np.load(f"error_vs_position_{k}.npz"))
    if results:
        save_plots(results)


if __name__ == "__main__":
    main()
