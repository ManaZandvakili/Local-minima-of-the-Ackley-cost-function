"""Export the simulator's numbers to MATLAB (.mat files in matlab/data/).

1. airship_params.mat - the physical model: mass matrices, weight, buoyancy,
   centre of gravity, drag, damping, motor time constant, time step. These are
   the exact numbers the Python simulator uses (configs/sphere.yaml + mass budget).
2. layouts/<name>.mat - propeller layouts: positions, directions, the 6xN
   thrust-to-wrench matrix A and the thrust limits.
3. reference_flight.mat - one 60 s gusty flight of the reversible knee layout
   flown by the tuned PID in Python: the initial state, every thrust command
   (10 Hz), every wind sample (100 Hz) and the resulting states (10 Hz). MATLAB
   and Simulink replay the same commands and wind and must land on the same states.
4. controllers/ - for the closed-loop model: the tuned PID gains of every layout
   (pid_<layout>.mat) and the trained residual and plain learning networks
   (<kind>_<case>.mat: weights, observation normalisation, scales).
5. validation/ - closed-loop flights in Python (classic, residual, plain on the
   8-propeller reversible layout; gusty station keeping and the 50 m circle) with
   their wind, so the MATLAB loop and the Simulink model can be checked against them.

Usage (from the airship/ folder):
    python scripts/export_matlab.py
"""
import csv
import json
import pickle
import re
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.allocation import Allocator                         # noqa: E402
from airship.config import load_config                           # noqa: E402
from airship.controllers.pid import PIDGains                     # noqa: E402
from airship.controllers.tuning import make_pid                  # noqa: E402
from airship.dynamics import DynamicsParams, Simulator, quat_from_euler  # noqa: E402
from airship.feasibility import axis_capacities                  # noqa: E402
from airship.layouts import load                                 # noqa: E402
from airship.path_following import CirclePath                    # noqa: E402
from airship.wind import WindConfig, WindModel                   # noqa: E402

OUT = ROOT / "matlab" / "data"
L = ROOT / "results" / "layouts"
LAYOUTS = {                                     # MATLAB name: layout file
    "knee14_reversible": L / "phase1_bidir_20261117_145131" / "sphere_knee.json",
    "knee8_reversible_capped": L / "phase1_bidir_maxN11_20261117_151234" / "sphere_knee.json",
    "knee15_oneway": L / "phase1_20261112_015323" / "sphere_knee.json",
    "knee9_oneway_capped": L / "phase1_maxN11_20261112_024218" / "sphere_knee.json",
    "paper16_reversible": L / "sphere_paper_best.json",
    "paper16_oneway": L / "oneway" / "sphere_paper_best.json",
    "starlab6_reversible": L / "sphere_starlab_6thruster.json",
}
PID_RUN = {"reversible": ROOT / "results" / "pid" / "20261117_153423",
           "oneway": ROOT / "results" / "pid" / "20261117_203431"}
RL_RUN = {("residual", "reversible"): ROOT / "results" / "rl" / "20261117_155720_bidir_residual_seed0",
          ("plain", "reversible"): ROOT / "results" / "rl" / "20261117_155720_bidir_plain_seed0",
          ("residual", "oneway"): ROOT / "results" / "rl" / "20261117_205532_oneway_residual_seed0",
          ("plain", "oneway"): ROOT / "results" / "rl" / "20261117_205532_oneway_plain_seed0"}
# MATLAB layout name: (propeller type, label of its tuned PID gains)
GAIN_LABEL = {"knee14_reversible": ("reversible", "sphere_knee"),
              "knee8_reversible_capped": ("reversible", "maxN11_sphere_knee"),
              "knee15_oneway": ("oneway", "sphere_knee"),
              "knee9_oneway_capped": ("oneway", "maxN11_sphere_knee"),
              "paper16_reversible": ("reversible", "sphere_paper_best"),
              "paper16_oneway": ("oneway", "sphere_paper_best"),
              "starlab6_reversible": ("reversible", "sphere_starlab_6thruster")}
# typical random layouts (median of the 20, classic controller, gusty station keeping)
TYPICAL = {"random8_reversible_typical": ("reversible", "maxN11_sphere_random_N8_"),
           "random9_oneway_typical": ("oneway", "maxN11_sphere_random_N9_")}
GAINS = PID_RUN["reversible"] / "gains" / "sphere_knee.json"


def physical(params, sim):
    return {
        "Minv": sim.Minv, "Mrb": sim.Mrb, "Ma": sim.Ma,
        "weight": sim.weight, "buoy": sim.buoy, "r_g": params.r_g.reshape(3, 1),
        "k_drag": params.k_drag, "d_lin": params.d_lin.reshape(3, 1), "d_quad": params.d_quad.reshape(3, 1),
        "tau_m": params.tau_m, "dt": sim.dt, "n_sub": float(sim.n_sub), "g": 9.81,
        "mass_kg": params.mass_kg, "buoyancy_kg": params.buoyancy_kg, "added_mass_kg": params.added_mass_kg,
        "inertia_kgm2": params.I_o, "added_inertia_kgm2": params.J_a,
        "tau_dist": np.zeros((3, 1)),
        "notes": "Exported from airship/scripts/export_matlab.py; SI units; body frame z up, origin at the hull centre",
    }


def layout_struct(name, lay, path):
    return {"name": name, "source_file": str(path.relative_to(ROOT)), "n_props": float(lay.n_props),
            "positions_m": lay.positions_m, "directions": lay.directions, "A": lay.M(),
            "u_lo": lay.u_min.reshape(-1, 1).astype(float), "u_hi": lay.u_max.reshape(-1, 1).astype(float),
            "caps": axis_capacities(lay.M(), lay.u_min, lay.u_max).reshape(-1, 1)}


def typical_random(case, prefix):
    """Label of the random layout whose classic-controller gusty error is the median of the 20."""
    rows = csv.DictReader(open(RL_RUN[("residual", case)] / "eval" / "eval.csv"))
    per = {}
    for r in rows:
        if r["controller"] == "PID" and r["layout"].startswith(prefix):
            per.setdefault(r["layout"], []).append(float(r["rms_pos_m"]))
    med = sorted((float(np.median(v)), k) for k, v in per.items())
    return med[len(med) // 2][1]


def pid_struct(gfile, params):
    g = json.load(open(gfile))
    pid = make_pid(PIDGains(**g["gains"]), params)
    return {"gains": g["gains"], "source_file": str(Path(gfile).relative_to(ROOT)),
            "kp": pid.kp.reshape(-1, 1), "ki": pid.ki.reshape(-1, 1), "kd": pid.kd.reshape(-1, 1),
            "kp_rp": pid.g.kp_rp, "kd_rp": pid.g.kd_rp, "kp_yaw": pid.g.kp_yaw, "ki_yaw": pid.g.ki_yaw,
            "kd_yaw": pid.g.kd_yaw, "mass_eff": pid.m, "inertia": pid.I.reshape(-1, 1), "dt": pid.dt,
            "i_max_pos": pid.i_max_pos, "i_max_yaw": pid.i_max_yaw}


def network_struct(run_dir):
    """Weights of the deterministic policy (mean action), normalisation and scales."""
    import yaml
    from stable_baselines3 import PPO
    cfg = yaml.safe_load(open(run_dir / "config.yaml"))["train_rl"]
    model = PPO.load(run_dir / "model", device="cpu")
    vn = pickle.load(open(run_dir / "vecnormalize.pkl", "rb"))
    pn = model.policy.mlp_extractor.policy_net
    an = model.policy.action_net
    w = lambda m: m.weight.detach().numpy().astype(np.float64)
    b = lambda m: m.bias.detach().numpy().astype(np.float64).reshape(-1, 1)
    e, r = cfg["env"], cfg.get("residual", {}) or {}
    return {"W1": w(pn[0]), "b1": b(pn[0]), "W2": w(pn[2]), "b2": b(pn[2]), "W3": w(an), "b3": b(an),
            "obs_mean": vn.obs_rms.mean.reshape(-1, 1), "obs_var": vn.obs_rms.var.reshape(-1, 1),
            "clip_obs": float(vn.clip_obs), "epsilon": float(vn.epsilon),
            "action_scale": np.asarray(e["action_scale"], float).reshape(-1, 1),
            "residual": float(bool(r.get("enabled", False))),
            "residual_scale": np.asarray(r.get("scale", [0] * 6), float).reshape(-1, 1),
            "pos_scale_m": e["pos_scale_m"], "vel_scale_m_s": e["vel_scale_m_s"],
            "rate_scale_rad_s": e["rate_scale_rad_s"], "integral_clip_m_s": e["integral_clip_m_s"],
            "source_run": str(run_dir.relative_to(ROOT))}


def validation_flight(params, lay, gfile, kind, case, scenario):
    """Closed loop as in closed_loop.run_episode / path_following.run_path_episode, recording the wind."""
    from airship.allocation import Allocator as Al
    from airship.controllers.rl_policy import RLController
    g = json.load(open(gfile))
    caps = axis_capacities(lay.M(), lay.u_min, lay.u_max)
    ctrl = make_pid(PIDGains(**g["gains"]), params)
    if kind != "pid":
        ctrl = RLController.load(RL_RUN[(kind, case)], caps, None, base_pid=ctrl)
    if scenario == "station":
        wcfg = WindConfig(mean_speed_m_s=3.0, mean_dir_deg=0.0, turbulence_sigma_m_s=(0.5, 0.5, 0.25),
                          gust_rate_per_min=2.0)
        sim = Simulator(params, lay, WindModel(wcfg, 500, 65.0))
        x0 = sim.reset(position=(1.0, -0.6, 0.3))
        K = 600
        ref = lambda t: (np.zeros(3), np.zeros(3), np.zeros(3), np.deg2rad(45.0) if t >= 30.0 else 0.0, 0.0)
    else:
        path = CirclePath()
        wcfg = WindConfig(mean_speed_m_s=3.0, mean_dir_deg=0.0, turbulence_sigma_m_s=(0.5, 0.5, 0.25),
                          gust_rate_per_min=2.0)
        sim = Simulator(params, lay, WindModel(wcfg, 7, path.duration_s + 5))
        x0 = sim.reset()
        K = int(round(path.duration_s * 10))
        ref = path.ref
    alloc = Al(lay.M(), lay.u_min, lay.u_max)
    ctrl.reset()
    X, U, TAU, W, REF = (np.zeros((K, sim.x.size)), np.zeros((K, lay.n_props)), np.zeros((K, 6)), [],
                         np.zeros((K, 11)))
    for k in range(K):
        t = k * 0.1
        p_ref, v_ref, a_ref, yaw_ref, yaw_rate = ref(t)
        tau = ctrl.act(sim.x, p_ref, yaw_ref, v_ref=v_ref, a_ref=a_ref, yaw_rate_ref=yaw_rate)
        res = alloc.solve(tau)
        ctrl.notify(res)
        w = sim.wind.sample(sim.n_sub, sim.dt)
        X[k] = sim.step(res.u, winds=w)
        U[k], TAU[k] = res.u, tau
        REF[k] = np.concatenate([p_ref, v_ref, a_ref, [yaw_ref, yaw_rate]])
        W.append(w)
    return {"x0": x0.reshape(-1, 1), "X": X, "U": U, "tau": TAU, "ref": REF}, np.vstack(W)


def reference_flight(params, lay):  # noqa: D401
    """60 s gusty flight with the tuned PID; records everything needed to replay it."""
    wcfg = WindConfig(mean_speed_m_s=3.0, mean_dir_deg=120.0, turbulence_sigma_m_s=(0.5, 0.5, 0.25),
                      gust_rate_per_min=2.0)
    sim = Simulator(params, lay, WindModel(wcfg, seed=101, duration_s=70.0))
    x0 = sim.reset(position=(1.0, -0.6, 0.3), rpy=(0.0, 0.0, 0.0))
    g = json.load(open(GAINS))
    pid = make_pid(PIDGains(**g["gains"]), params)
    alloc = Allocator(lay.M(), lay.u_min, lay.u_max)
    pid.reset()
    K = 600
    cmd, winds, xs = np.zeros((K, lay.n_props)), [], np.zeros((K, sim.x.size))
    for k in range(K):
        yaw_ref = np.deg2rad(45.0) if k * 0.1 >= 30.0 else 0.0
        res = alloc.solve(pid.act(sim.x, np.zeros(3), yaw_ref))
        pid.notify(res)
        w = sim.wind.sample(sim.n_sub, sim.dt)
        cmd[k], xs[k] = res.u, sim.step(res.u, winds=w)
        winds.append(w)
    return {"x0": x0.reshape(-1, 1), "cmd": cmd, "wind": np.vstack(winds), "x_python": xs,
            "t_cmd": np.arange(K) * 0.1, "t_x": (np.arange(K) + 1) * 0.1,
            "layout": "knee14_reversible",
            "description": "Reversible 14-propeller knee layout, 3 m/s wind toward 120 deg, turbulence 0.5 m/s, "
                           "2 gusts/min, seed 101, tuned PID, heading step to 45 deg at 30 s. cmd(k,:) is held "
                           "for 0.1 s; wind(j,:) is held for 0.01 s; x_python(k,:) is the state after command k."}


def main():
    (OUT / "layouts").mkdir(parents=True, exist_ok=True)
    params = DynamicsParams.from_config(load_config("sphere"))
    lays = {n: load(p) for n, p in LAYOUTS.items()}
    sim = Simulator(params, lays["knee14_reversible"], WindModel(WindConfig(), seed=0))
    savemat(OUT / "airship_params.mat", {"phys": physical(params, sim)})
    for n, lay in lays.items():
        savemat(OUT / "layouts" / f"{n}.mat", {"layout": layout_struct(n, lay, LAYOUTS[n])})
        print(f"layout {n}: {lay.n_props} propellers, thrust {lay.u_min.min():+.0f} to {lay.u_max.max():.0f} N")
    ref = reference_flight(params, lays["knee14_reversible"])
    savemat(OUT / "reference_flight.mat", {"ref": ref})
    print(f"reference flight: {len(ref['cmd'])} commands, {len(ref['wind'])} wind samples, "
          f"final position {ref['x_python'][-1, :3].round(3)}")

    # closed loop: typical random layouts, PID gains of every layout, the trained networks
    gain_files = {n: PID_RUN[c] / "gains" / f"{lab}.json" for n, (c, lab) in GAIN_LABEL.items()}
    for n, (case, prefix) in TYPICAL.items():
        lab = typical_random(case, prefix)
        gf = PID_RUN[case] / "gains" / f"{lab}.json"
        lp = ROOT / json.load(open(gf))["path"]
        lays[n] = load(lp)
        savemat(OUT / "layouts" / f"{n}.mat", {"layout": layout_struct(n, lays[n], lp)})
        gain_files[n] = gf
        print(f"layout {n}: {lab} ({lays[n].n_props} propellers)")
    (OUT / "controllers").mkdir(exist_ok=True)
    for n, gf in gain_files.items():
        savemat(OUT / "controllers" / f"pid_{n}.mat", {"pid": pid_struct(gf, params)})
    for (kind, case), run in RL_RUN.items():
        savemat(OUT / "controllers" / f"{kind}_{case}.mat", {"net": network_struct(run)}, do_compression=True)
        print(f"network {kind}_{case} <- {run.name}")

    # validation flights (Python) on the 8-propeller reversible layout
    (OUT / "validation").mkdir(exist_ok=True)
    n = "knee8_reversible_capped"
    for scen in ("station", "circle"):
        wind = None
        for kind in ("pid", "residual", "plain"):
            rec, w = validation_flight(params, lays[n], gain_files[n], kind, "reversible", scen)
            wind = w if wind is None else wind
            assert np.array_equal(w, wind)
            err = rec["X"][:, :3] - (0 if scen == "station" else np.vstack([rec["ref"][1:, :3],
                                                                         CirclePath().ref(len(rec["X"]) * 0.1)[0]]))
            if scen == "circle":                      # keep the files small: positions only
                rec = {"x0": rec["x0"], "X": rec["X"][:, :3]}
            rec.update({"layout": n, "controller": kind, "scenario": scen})
            savemat(OUT / "validation" / f"python_{kind}_{scen}.mat", {"py": rec}, do_compression=True)
            print(f"validation {kind} {scen}: RMS error {np.sqrt(np.mean(np.sum(err ** 2, 1))):.3f} m")
        savemat(OUT / "validation" / f"wind_{scen}.mat", {"wind": wind}, do_compression=True)


if __name__ == "__main__":
    main()
