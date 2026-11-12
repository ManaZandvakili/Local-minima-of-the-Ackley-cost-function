"""Phase 4 check: the trained RL policy against the tuned PID, same layouts, same
gusty winds (identical seeds), same allocation; plus training curves and videos.

Usage (from the airship/ folder):
    python scripts/evaluate_rl.py --rl results/rl/<run> --pid results/pid/<run> --workers 4

Outputs: <rl run>/eval/ with eval.csv, report.md, training_curves.png and videos.
"""
import argparse
import csv
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.closed_loop import Scenario, run_episode                 # noqa: E402
from airship.config import load_config                                 # noqa: E402
from airship.controllers.pid import PIDGains                           # noqa: E402
from airship.controllers.rl_policy import RLController                 # noqa: E402
from airship.controllers.tuning import make_pid                        # noqa: E402
from airship.dynamics import DynamicsParams                            # noqa: E402
from airship.feasibility import axis_capacities                        # noqa: E402
from airship.layouts import load                                        # noqa: E402

DIRECTIONS = (0.0, 120.0, 240.0)


def scenarios(base_seed=500):
    """Gusty station keeping: 3 m/s mean wind from 3 directions, turbulence and gusts,
    1.2 m initial offset, 45 deg heading change at 30 s."""
    out = []
    for i, d in enumerate(DIRECTIONS):
        out.append(Scenario.from_config({
            "name": f"turbulent_dir{int(d)}", "duration_s": 60.0,
            "wind": {"mean_speed_m_s": 3.0, "mean_dir_deg": d, "turbulence_sigma_m_s": [0.5, 0.5, 0.25],
                     "gust_rate_per_min": 2.0},
            "wind_seed": base_seed + i, "init_pos": [1.0, -0.6, 0.3], "yaw_step_time_s": 30.0,
            "yaw_step_to_deg": 45.0}))
    return out


CKPT = None   # set from --checkpoint


def eval_layout(args):
    label, layout_path, gains, rl_dir, held_out, ckpt = args
    params = DynamicsParams.from_config(load_config("sphere"))
    lay = load(ROOT / layout_path)
    caps = axis_capacities(lay.M(), lay.u_min, lay.u_max)
    rl = RLController.load(rl_dir, caps, ckpt)
    rows = []
    for sc in scenarios():
        for name, ctrl in (("PID", make_pid(PIDGains(**gains), params)), ("RL", rl)):
            r = run_episode(lay, params, ctrl, sc)
            rows.append({"layout": label, "n_props": lay.n_props, "held_out": held_out, "controller": name,
                         "scenario": sc.name, **{k: r[k] for k in ("rms_pos_m", "max_pos_m", "rms_yaw_rad",
                                                                   "effort_N2s", "saturation_pct", "survived")}})
    return rows


def group(label):
    for key in ("maxN11_sphere_random_N9", "sphere_random_N15", "sphere_random_N22"):
        if label.startswith(key):
            return key
    return label


def training_curves(rl_dir, out_png):
    cur = json.load(open(Path(rl_dir) / "curriculum.json"))["log"]
    steps = np.array([c["steps"] for c in cur]) / 1e6
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    ax[0].plot(steps, [c["rms_pos_m"] for c in cur], color="#1f4e9c")
    ax[0].set(xlabel="training steps [millions]", ylabel="episode RMS position error [m]", yscale="log",
              title="Tracking error during training")
    ax[1].plot(steps, [c["survival"] for c in cur], color="#1f4e9c")
    ax[1].set(xlabel="training steps [millions]", ylabel="fraction of episodes survived", ylim=(-0.05, 1.05),
              title="Survival")
    names = [c["stage"] for c in cur]
    for a in ax[:2]:
        for i in range(1, len(names)):
            if names[i] != names[i - 1]:
                a.axvline(steps[i], color="#999", lw=0.8, ls="--")
                a.text(steps[i], a.get_ylim()[1], f" {names[i]}", fontsize=8, va="top", rotation=90)
    try:
        from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
        tb = sorted((Path(rl_dir) / "tb").glob("*/events.*"))[-1]
        ea = EventAccumulator(str(tb)); ea.Reload()
        for tag in [t for t in ea.Tags()["scalars"] if t.startswith("reward_terms/")]:
            ev = ea.Scalars(tag)
            ax[2].plot([e.step / 1e6 for e in ev], [e.value for e in ev], label=tag.split("/")[1])
        ax[2].legend(fontsize=8, ncol=2, frameon=False)
        ax[2].set(xlabel="training steps [millions]", ylabel="mean per-step value", title="Reward terms")
    except Exception as e:                                   # tensorboard missing or no events yet
        ax[2].set_title(f"reward terms unavailable ({type(e).__name__})")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out_png, dpi=140); plt.close(fig)


def videos(rl_dir, pid_dir, out_dir):
    from airship.viz.animate import animate
    params = DynamicsParams.from_config(load_config("sphere"))
    sc = scenarios()[1]
    gains = lambda lab: PIDGains(**json.load(open(Path(pid_dir) / "gains" / f"{lab}.json"))["gains"])
    path = lambda lab: ROOT / json.load(open(Path(pid_dir) / "gains" / f"{lab}.json"))["path"]

    def fly(lab, ctrl_name, scenario):
        lay = load(path(lab))
        ctrl = (make_pid(gains(lab), params) if ctrl_name == "PID"
                else RLController.load(rl_dir, axis_capacities(lay.M(), lay.u_min, lay.u_max), CKPT))
        return run_episode(lay, params, ctrl, scenario, log=True)["log"], lay

    made = []
    a, la = fly("sphere_knee", "PID", sc); b, lb = fly("sphere_knee", "RL", sc)
    made.append(animate([(a, la, "PID - optimized layout (15)"), (b, lb, "RL - optimized layout (15)")],
                        Path(out_dir) / "pid_vs_rl_knee.mp4"))
    a, la = fly("maxN11_sphere_knee", "RL", sc); b, lb = fly("maxN11_sphere_random_N9_00", "RL", sc)
    made.append(animate([(a, la, "RL - optimized layout (9)"), (b, lb, "RL - random layout (9)")],
                        Path(out_dir) / "optimized_vs_random_N9_rl.mp4"))
    fsc = Scenario.from_config({"name": "failure", "duration_s": 45.0,
                                "wind": {"mean_speed_m_s": 2.0, "mean_dir_deg": 60.0,
                                         "turbulence_sigma_m_s": [0.3, 0.3, 0.15]},
                                "wind_seed": 777, "fail_prop": 0, "fail_time_s": 15.0})
    a, la = fly("sphere_knee", "PID", fsc); b, lb = fly("sphere_knee", "RL", fsc)
    made.append(animate([(a, la, "PID - motor 0 fails at 15 s"), (b, lb, "RL - motor 0 fails at 15 s")],
                        Path(out_dir) / "motor_failure_knee.mp4"))
    return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rl", required=True)
    ap.add_argument("--pid", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--no-videos", action="store_true")
    ap.add_argument("--checkpoint", type=int, help="evaluate the checkpoint saved after this many steps")
    ap.add_argument("--only", nargs="*", help="only layouts whose label starts with one of these")
    a = ap.parse_args()
    global CKPT
    CKPT = a.checkpoint
    rl_dir, pid_dir = Path(a.rl).resolve(), Path(a.pid).resolve()
    out = rl_dir / ("eval" if a.checkpoint is None else f"eval_ckpt{a.checkpoint}")
    out.mkdir(exist_ok=True)
    import pickle
    pool = pickle.load(open(rl_dir / "pool.pkl", "rb"))
    held = {L.name for L in pool.holdout}                   # "<run folder>:<layout name>"
    jobs = []
    for f in sorted((pid_dir / "gains").glob("*.json")):
        g = json.load(open(f))
        run_folder = Path(g["path"]).parent.name               # e.g. phase1_maxN11_20261112_024218
        name = Path(g["path"]).stem
        jobs.append((g["layout"], g["path"], g["gains"], str(rl_dir), f"{run_folder}:{name}" in held, a.checkpoint))
    if a.only:
        jobs = [j for j in jobs if any(j[0].startswith(o) for o in a.only)]
    with Pool(a.workers) as p:
        rows = [r for rs in p.map(eval_layout, jobs) for r in rs]
    with open(out / "eval.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # per layout: median over the three wind directions, then per group
    per = {}
    for r in rows:
        per.setdefault((r["layout"], r["controller"]), []).append(r)
    summ = {}
    for (lab, ctrl), rs in per.items():
        summ[(lab, ctrl)] = {"rms": np.median([x["rms_pos_m"] for x in rs]),
                             "effort": np.median([x["effort_N2s"] for x in rs]),
                             "yaw": np.median([x["rms_yaw_rad"] for x in rs]),
                             "surv": np.mean([x["survived"] for x in rs]), "held": rs[0]["held_out"],
                             "n": rs[0]["n_props"]}
    L = ["# Phase 4 evaluation: RL vs PID\n",
         f"RL run `{rl_dir.name}`, PID gains `{pid_dir.name}`. Gusty station keeping: 3 m/s mean wind from "
         f"{len(DIRECTIONS)} directions, turbulence 0.5 m/s, gusts, 45° heading change; identical wind seeds "
         "for both controllers. Values: median over the wind directions (per layout), then median and IQR "
         "over the layouts in a group.\n",
         "| Layout group | N | PID RMS pos [m] | RL RMS pos [m] | PID effort [N²s] | RL effort [N²s] | "
         "PID survival | RL survival |", "|---|---|---|---|---|---|---|---|"]
    groups = {}
    for (lab, ctrl), s in summ.items():
        groups.setdefault(group(lab), {}).setdefault(ctrl, []).append(s)

    def fmt(v, f):
        return f"{np.median(v):{f}}" if len(v) == 1 else \
            f"{np.median(v):{f}} [{np.percentile(v, 25):{f}}-{np.percentile(v, 75):{f}}]"
    for gname, d in groups.items():
        P, R = d["PID"], d["RL"]
        L.append(f"| {gname}{' (' + str(len(P)) + ')' if len(P) > 1 else ''} | {P[0]['n']} | "
                 f"{fmt([s['rms'] for s in P], '.3f')} | {fmt([s['rms'] for s in R], '.3f')} | "
                 f"{fmt([s['effort'] for s in P], '.0f')} | {fmt([s['effort'] for s in R], '.0f')} | "
                 f"{np.mean([s['surv'] for s in P]):.2f} | {np.mean([s['surv'] for s in R]):.2f} |")
    held_in = [k for k, s in summ.items() if not s["held"]]
    pid_med = np.median([summ[(l, "PID")]["rms"] for l, c in held_in if c == "PID"])
    rl_med = np.median([summ[(l, "RL")]["rms"] for l, c in held_in if c == "RL"])
    L.append(f"\n**Brief acceptance (held-in layouts, turbulent scenario):** median RMS position error "
             f"PID {pid_med:.3f} m, RL {rl_med:.3f} m -> RL {'beats' if rl_med < pid_med else 'does not beat'} PID.\n")
    ho = [k for k, s in summ.items() if s["held"]]
    if ho:
        L.append(f"Held-out layouts ({len(ho) // 2}): median RMS PID "
                 f"{np.median([summ[(l, 'PID')]['rms'] for l, c in ho if c == 'PID']):.3f} m, RL "
                 f"{np.median([summ[(l, 'RL')]['rms'] for l, c in ho if c == 'RL']):.3f} m.\n")
    training_curves(rl_dir, out / "training_curves.png")
    L.append("Training curves: `training_curves.png`.\n")
    if not a.no_videos:
        for v in videos(rl_dir, pid_dir, out):
            L.append(f"Video: `{Path(v).name}`")
    text = "\n".join(L) + "\n"
    (out / "report.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
