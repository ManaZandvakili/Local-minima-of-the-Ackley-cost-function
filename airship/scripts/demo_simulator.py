"""Phase 2 demo: the 6-DOF sphere with no control input.

Left:   released from a 10 deg roll in still air (the gondola pendulum).
Middle/right: drifting with zero thrust in 3 m/s wind with turbulence and gusts.

Usage (from the airship/ folder):
    python scripts/demo_simulator.py      -> results/figures/phase2_demo.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.config import load_config                                    # noqa: E402
from airship.dynamics import DynamicsParams, Simulator, euler_from_quat, quat_to_rot  # noqa: E402
from airship.layouts import load                                            # noqa: E402
from airship.wind import WindConfig, WindModel                              # noqa: E402

CFG = load_config("sphere")
LAYOUT = load(ROOT / "results" / "layouts" / "sphere_paper_best.json")


def fly(wind_cfg, seconds, rpy=(0, 0, 0), seed=1):
    sim = Simulator(DynamicsParams.from_config(CFG), LAYOUT, WindModel(wind_cfg, seed))
    sim.reset(rpy=rpy)
    t, xs, winds = [], [], []
    for k in range(int(seconds * 10)):
        xs.append(sim.step(np.zeros(sim.n_props)))
        winds.append(sim.last_wind)
        t.append(sim.t)
    return np.array(t), np.array(xs), np.array(winds)


def main():
    out = ROOT / "results" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    t, x, _ = fly(WindConfig(turbulence_sigma_m_s=(0, 0, 0)), 40, rpy=(np.deg2rad(10), 0, 0))
    roll = np.degrees([euler_from_quat(q)[0] for q in x[:, 3:7]])
    ax[0].plot(t, roll, color="#1f4e9c")
    ax[0].set(xlabel="time [s]", ylabel="roll [deg]", title="Released from 10° roll, still air")

    w = WindConfig(mean_speed_m_s=3.0, mean_dir_deg=0.0, gust_rate_per_min=3.0)
    t, x, winds = fly(w, 120)
    ax[1].plot(x[:, 0], x[:, 1], color="#1f4e9c")
    ax[1].plot(0, 0, "ko")
    ax[1].set(xlabel="x [m]", ylabel="y [m]", title="Zero thrust, 3 m/s wind + turbulence + gusts")
    ax[1].axis("equal")
    v = np.array([quat_to_rot(q) @ vb for q, vb in zip(x[:, 3:7], x[:, 7:10])])
    ax[2].plot(t, winds[:, 0], color="#b9c0c7", label="wind x")
    ax[2].plot(t, v[:, 0], color="#1f4e9c", label="airship x-velocity")
    rp = np.degrees([euler_from_quat(q)[:2] for q in x[:, 3:7]])
    ax2 = ax[2].twinx()
    ax2.plot(t, rp[:, 1], color="#e07b00", lw=0.8, label="pitch")
    ax2.set_ylabel("pitch [deg]", color="#e07b00")
    ax[2].set(xlabel="time [s]", ylabel="speed [m/s]", title="Speed follows the wind; gusts tilt it")
    ax[2].legend(loc="upper left", frameon=False)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "phase2_demo.png", dpi=150)
    print(f"saved {out / 'phase2_demo.png'}")


if __name__ == "__main__":
    main()
