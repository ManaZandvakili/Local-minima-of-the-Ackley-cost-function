"""Phase 1: run the push-only-aware layout optimizer on the sphere.

For every seed: 20 PSO trials -> Pareto front (authority, effort, count) ->
knee and paper-rule pick. Then, over all seeds: the overall knee and pick,
20 matched-N random feasible layouts, full capacity checks, layout JSON
exports and a markdown report.

Usage (from the airship/ folder):
    python scripts/run_optimizer.py                 # all seeds in configs/optimizer.yaml
    python scripts/run_optimizer.py --seeds 0 1 --workers 2
    python scripts/run_optimizer.py --summarize results/runs/phase1_<timestamp>

Outputs go to a new timestamped folder (never overwritten):
    results/runs/phase1_<timestamp>/     config, per-seed results, summary.csv, report.md
    results/layouts/phase1_<timestamp>/  layout JSON files
"""
import argparse
import csv
import json
import pickle
import sys
import time
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.config import hull_from_config, load_config                   # noqa: E402
from airship.feasibility import AXIS_LABELS, Requirements, assess            # noqa: E402
from airship.layouts import Layout, load, save                               # noqa: E402
from airship.optimizer.objectives import OBJECTIVES, PAPER_KEYS, GeometricRules  # noqa: E402
from airship.optimizer.pareto import knee_index, paper_rule_order, pareto_filter, pareto_points  # noqa: E402
from airship.optimizer.pso import PSOConfig, run_trials                      # noqa: E402
from airship.optimizer.random_layouts import sample_matched                  # noqa: E402


def setup(opt):
    hull_cfg = load_config(opt["hull_config"])
    hull = hull_from_config(hull_cfg)
    R = hull.semi_axes_m[0] if hull.shape == "sphere" else hull_cfg["hull"]["diameter_m"] / 2
    rules = GeometricRules(length_m=R, **opt["geometry_rules"])
    f = dict(opt["feasibility"])
    n_sampled, sample_seed = f.pop("n_sampled_directions"), f.pop("sample_seed")
    req = Requirements.from_config(f)
    return dict(hull_cfg=hull_cfg, hull=hull, R=R, rules=rules, req=req,
                u_max=float(hull_cfg["propulsion"]["u_max_N"]),
                cands=hull.candidate_points(opt["K_candidates"]),
                cfg=PSOConfig.from_config(opt["pso"]), n_sampled=n_sampled, sample_seed=sample_seed)


def run_seed(args):
    opt, seed, run_dir = args
    s = setup(opt)
    t0 = time.time()
    trials = run_trials(s["hull"], s["cands"], s["cfg"], s["req"], s["u_max"], s["rules"],
                        opt["trials_per_seed"], seed, log=lambda m: print(m, flush=True))
    rows = [{"trial": r["trial"], "weights": r["weights"].tolist(), "score": r["eval"].score,
             "positions_m": r["eval"].positions_m, "directions": r["eval"].directions,
             "objectives": r["eval"].objectives, "axis_capacity": r["eval"].axis_capacity,
             "shortfall": r["eval"].shortfall, "gate": r["eval"].gate,
             "excess_props": r["eval"].excess_props,
             "convergence": r["convergence"], "n_evals": r["n_evals"]} for r in trials]
    with open(Path(run_dir) / f"seed{seed}.pkl", "wb") as fh:
        pickle.dump({"seed": seed, "trials": rows, "seconds": time.time() - t0}, fh)
    return seed, time.time() - t0


def make_layout(name, p, d, s, source, objectives, seed, notes):
    return Layout(name=name, shape=s["hull"].shape, semi_axes_m=s["hull"].semi_axes_m,
                  positions_m=p, directions=d, u_max_N=s["u_max"], u_min_N=0.0, source=source,
                  objectives={PAPER_KEYS[k]: objectives[k] for k in OBJECTIVES} if objectives else {},
                  seed=seed, notes=notes)


def full_check(s, p, d):
    from airship.wrench import build_M
    return assess(build_M(p, d), 0.0, s["u_max"], s["req"], s["R"], s["n_sampled"], s["sample_seed"])


def summarize(run_dir, opt):
    run_dir = Path(run_dir)
    s = setup(opt)
    lay_dir = ROOT / "results" / "layouts" / run_dir.name
    lay_dir.mkdir(parents=True, exist_ok=True)
    sel = opt["selection"]
    seeds = sorted(pickle.load(open(p, "rb"))["seed"] for p in run_dir.glob("seed*.pkl"))
    per_seed, pool, csv_rows = [], [], []
    random_cache = {}

    def matched(n):
        if n not in random_cache:
            random_cache[n] = sample_matched(n, s["hull"], s["cands"], s["req"], s["u_max"], s["rules"],
                                             n_target=sel["random_matched"], seed=sel["random_seed_offset"] + n)
        return random_cache[n]

    for seed in seeds:
        data = pickle.load(open(run_dir / f"seed{seed}.pkl", "rb"))
        trials = [t for t in data["trials"]
                  if not t["gate"] and t["shortfall"] == 0.0 and t.get("excess_props", 0) == 0]
        objs = [t["objectives"] for t in trials]
        pidx = pareto_filter(pareto_points(objs))
        knee = pidx[knee_index(pareto_points(objs)[pidx])]
        pick = None
        for k in paper_rule_order(objs, pidx):
            if len(matched(len(trials[k]["positions_m"]))[0]) >= sel["random_min"]:
                pick = k
                break
        for k in pidx:
            t = trials[k]
            name = f"sphere_s{seed}_t{t['trial']:02d}"
            fr = full_check(s, t["positions_m"], t["directions"])
            tags = (["knee"] if k == knee else []) + (["paper_pick"] if k == pick else [])
            lay = make_layout(name, t["positions_m"], t["directions"], s, "optimizer", t["objectives"], seed,
                              {"run": run_dir.name, "trial": t["trial"], "weights": t["weights"], "tags": tags,
                               "authority": t["objectives"]["authority"],
                               "planar_authority": t["objectives"]["planar_authority"],
                               "feasibility": fr.as_dict()})
            save(lay, lay_dir / f"{name}.json")
            pool.append((seed, k, t, fr, lay))
            csv_rows.append({"seed": seed, "trial": t["trial"], "N": len(t["positions_m"]), "tags": " ".join(tags),
                             **{PAPER_KEYS[o]: round(t["objectives"][o], 4) for o in OBJECTIVES},
                             "authority": round(t["objectives"]["authority"], 4),
                             "capacity_min_N": round(fr.capacity_min_N, 3), "feasible": fr.feasible})
        per_seed.append({"seed": seed, "n_feasible_trials": len(trials), "n_trials": len(data["trials"]),
                         "pareto": len(pidx), "knee_N": len(trials[knee]["positions_m"]),
                         "pick_N": len(trials[pick]["positions_m"]) if pick is not None else None,
                         "seconds": data["seconds"]})

    # overall knee and paper-rule pick over the union of all seeds' Pareto layouts
    objs = [p[2]["objectives"] for p in pool]
    gidx = pareto_filter(pareto_points(objs))
    g_knee = pool[gidx[knee_index(pareto_points(objs)[gidx])]]
    g_pick = None
    for k in paper_rule_order(objs, gidx):
        if len(matched(len(pool[k][2]["positions_m"]))[0]) >= sel["random_min"]:
            g_pick = pool[k]
            break
    chosen = {}
    for label, entry in (("sphere_knee", g_knee), ("sphere_paper_pick", g_pick)):
        if entry is None:
            continue
        lay = entry[4]
        out = Layout.from_dict({**lay.to_dict(), "name": label,
                                "notes": {**lay.notes, "same_as": lay.name}})
        save(out, lay_dir / f"{label}.json")
        chosen[label] = out

    # matched-N random layouts for each chosen layout
    rand_info = {}
    for label, lay in chosen.items():
        n = lay.n_props
        if n in rand_info:
            rand_info[n]["matched_to"].append(label)
            continue
        samples, attempts = matched(n)
        rcaps = []
        for i, (p, d) in enumerate(samples):
            fr = full_check(s, p, d)
            rcaps.append(fr.capacity_min_N)
            save(make_layout(f"sphere_random_N{n}_{i:02d}", p, d, s, "random", None,
                             sel["random_seed_offset"] + n,
                             {"run": run_dir.name, "matched_to": label, "feasibility": fr.as_dict()}),
                 lay_dir / f"sphere_random_N{n}_{i:02d}.json")
        rand_info[n] = {"N": n, "matched_to": [label], "count": len(samples), "attempts": attempts,
                        "capacity_min_N_median": float(np.median(rcaps)) if rcaps else None}

    with open(run_dir / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(csv_rows[0]))
        w.writeheader(); w.writerows(csv_rows)

    refs = {}
    for f in sorted((ROOT / "results" / "layouts").glob("*.json")):
        L = load(f)
        fr = assess(L.M(), L.u_min, L.u_max, s["req"], s["R"], s["n_sampled"], s["sample_seed"])
        refs[L.name] = (L, fr)

    report = write_report(run_dir, opt, s, per_seed, csv_rows, chosen, rand_info, refs)
    (run_dir / "report.md").write_text(report)
    json.dump({"per_seed": per_seed, "random": list(rand_info.values())}, open(run_dir / "summary.json", "w"), indent=2)
    print(report)


def write_report(run_dir, opt, s, per_seed, rows, chosen, rand_info, refs):
    L = []
    L.append(f"# Phase 1 optimizer run `{run_dir.name}`\n")
    L.append(f"Sphere, K = {opt['K_candidates']}, {len(per_seed)} seeds x {opt['trials_per_seed']} trials, "
             f"push-only thrust 0-{s['u_max']:g} N per propeller"
             + (f", at most {opt['pso']['max_props']} propellers" if opt['pso'].get('max_props') else "") + ".\n")
    L.append("## Per seed\n")
    L.append("| Seed | Feasible trials | Pareto layouts | Knee N | Paper-pick N | Minutes |")
    L.append("|---|---|---|---|---|---|")
    for p in per_seed:
        L.append(f"| {p['seed']} | {p['n_feasible_trials']}/{p['n_trials']} | {p['pareto']} | {p['knee_N']} | "
                 f"{p['pick_N']} | {p['seconds'] / 60:.1f} |")
    L.append("\n## Pareto layouts (all seeds)\n")
    L.append("| Seed | Trial | N | Count | Effort | Pos. sym. | Dir. sym. | Geometric | Authority | "
             "Capacity min [N] | Tags |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        L.append(f"| {r['seed']} | {r['trial']} | {r['N']} | {r['J1_count']:.3f} | {r['J2_effort']:.3f} | "
                 f"{r['J3_position_symmetry']:.3f} | {r['J4_direction_symmetry']:.3f} | {r['J5_geometric']:.3f} | "
                 f"{r['authority']:.3f} | {r['capacity_min_N']:.2f} | {r['tags']} |")
    Ns = [r["N"] for r in rows]
    L.append(f"\nN over all Pareto layouts: min {min(Ns)}, median {int(np.median(Ns))}, max {max(Ns)}.\n")
    L.append("## Chosen layouts (over all seeds)\n")
    for k, lay in chosen.items():
        fr = lay.notes["feasibility"]
        L.append(f"- **{k}**: N = {lay.n_props} (= {lay.notes['same_as']}), capacity min "
                 f"{fr['capacity_min_N']:.2f} N, objectives "
                 + ", ".join(f"{kk} {v:.3f}" for kk, v in lay.objectives.items()))
    for r in rand_info.values():
        L.append(f"\nMatched random layouts for {', '.join(r['matched_to'])}: {r['count']} with N = {r['N']} "
                 f"({r['attempts']} attempts), median capacity min {r['capacity_min_N_median']:.2f} N.")
    L.append("")
    L.append("## Reference layouts under the new test\n")
    L.append("| Layout | N | Thrust range [N] | Feasible (all 6 DOF) | Capacity min [N] | Axes with zero capacity |")
    L.append("|---|---|---|---|---|---|")
    for name, (lay, fr) in refs.items():
        zero = [a for a, c in zip(AXIS_LABELS, fr.axis_capacity) if c <= 1e-9]
        rng = f"{lay.u_min.min():g} to {lay.u_max.max():g}"
        L.append(f"| {name} | {lay.n_props} | {rng} | {'yes' if fr.feasible else 'no'} | "
                 f"{fr.capacity_min_N:.2f} | {', '.join(zero) or 'none'} |")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="optimizer")
    ap.add_argument("--seeds", type=int, nargs="*")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--summarize", help="only (re)build the summary of an existing run folder")
    ap.add_argument("--max-props", type=int, help="cap the number of propellers (graded penalty above it)")
    a = ap.parse_args()
    opt = load_config(a.config)
    if a.summarize:
        summarize(a.summarize, opt)
        return
    if a.seeds:
        opt["seeds"] = a.seeds
    tag = ""
    if a.max_props:
        opt["pso"]["max_props"] = a.max_props
        tag = f"_maxN{a.max_props}"
    run_dir = ROOT / "results" / "runs" / f"phase1{tag}_{datetime.now():%Y%m%d_%H%M%S}"
    run_dir.mkdir(parents=True)
    yaml.safe_dump({"optimizer": opt, "hull": load_config(opt["hull_config"])}, open(run_dir / "config.yaml", "w"))
    print(f"run folder: {run_dir}", flush=True)
    with Pool(min(a.workers, len(opt["seeds"]))) as pool:
        for seed, sec in pool.imap_unordered(run_seed, [(opt, s, str(run_dir)) for s in opt["seeds"]]):
            print(f"seed {seed} done in {sec / 60:.1f} min", flush=True)
    summarize(run_dir, opt)


if __name__ == "__main__":
    main()
