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

Usage (from the airship/ folder):
    python scripts/export_matlab.py
"""
import json
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
from airship.layouts import load                                 # noqa: E402
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
GAINS = ROOT / "results" / "pid" / "20261117_153423" / "gains" / "sphere_knee.json"


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
            "u_lo": lay.u_min.reshape(-1, 1).astype(float), "u_hi": lay.u_max.reshape(-1, 1).astype(float)}


def reference_flight(params, lay):
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


if __name__ == "__main__":
    main()
