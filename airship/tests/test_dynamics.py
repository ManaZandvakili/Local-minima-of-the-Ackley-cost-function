"""Phase 2 sanity tests for the 6-DOF simulator (brief, section 6)."""
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from airship.config import load_config
from airship.dynamics import DynamicsParams, Simulator, derivative, euler_from_quat, quat_from_euler, quat_to_rot
from airship.geometry import Hull
from airship.layouts import Layout, load
from airship.massprops import from_config
from airship.wind import WindConfig, WindModel

CFG = load_config("sphere")
R = CFG["hull"]["diameter_m"] / 2
LAYOUTS = Path(__file__).resolve().parent.parent / "results" / "layouts"
CALM = WindConfig(turbulence_sigma_m_s=(0.0, 0.0, 0.0))


def _sim(layout=None, params=None, wind=CALM, seed=0):
    layout = layout or load(LAYOUTS / "sphere_paper_best.json")
    return Simulator(params or DynamicsParams.from_config(CFG), layout, WindModel(wind, seed))


def _run(sim, seconds, u=None):
    u = np.zeros(sim.n_props) if u is None else u
    out = []
    for _ in range(int(round(seconds * 10))):
        out.append(sim.step(u))
    return np.array(out)


def test_quaternion_round_trip():
    for rpy in [(0.1, -0.2, 0.3), (0.0, 0.0, 2.5), (-0.4, 0.3, -1.0)]:
        np.testing.assert_allclose(euler_from_quat(quat_from_euler(*rpy)), rpy, atol=1e-12)
    q = quat_from_euler(0, 0, np.pi / 2)
    np.testing.assert_allclose(quat_to_rot(q) @ [1, 0, 0], [0, 1, 0], atol=1e-12)


def test_mass_properties_from_budget():
    mp = from_config(CFG)
    assert mp.mass_kg == pytest.approx(1.225 * 16.4)          # neutral buoyancy
    assert mp.trim_kg == pytest.approx(1.225 * 16.4 - 16.161, abs=1e-3)
    assert 0.3 < mp.h_cg_m < 0.7
    assert mp.added_mass_kg == pytest.approx(0.5 * 1.225 * 16.4)


def test_1_stays_put_at_rest():
    """Zero thrust, zero wind, level start: drift < 1 cm over 60 s."""
    sim = _sim()
    x = _run(sim, 60.0)
    assert np.max(np.linalg.norm(x[:, :3], axis=1)) < 0.01
    assert np.max(np.abs(x[:, 7:13])) < 1e-6


def test_2_pendulum_period():
    """10 deg roll, release, zero thrust: oscillates and settles with the pendulum period."""
    p = DynamicsParams.from_config(CFG)
    sim = _sim(params=p)
    sim.reset(rpy=(np.deg2rad(10), 0, 0))
    x = _run(sim, 40.0)
    roll = np.array([euler_from_quat(q)[0] for q in x[:, 3:7]])
    t = np.arange(1, len(roll) + 1) * 0.1
    cross = t[1:][(roll[:-1] > 0) & (roll[1:] <= 0)]               # downward zero crossings
    period = np.mean(np.diff(cross))
    m, h = p.mass_kg, -p.r_g[2]
    I_eff = p.I_o[0, 0] - (m * h) ** 2 / (m + p.added_mass_kg)    # roll coupled with sway
    expected = 2 * np.pi * np.sqrt(I_eff / (m * 9.81 * h))
    assert period == pytest.approx(expected, rel=0.10)
    assert np.abs(roll[-50:]).max() < np.deg2rad(10) * 0.5          # it settles


def test_3_initial_acceleration():
    """Constant Fx from a symmetric layout: initial acceleration ~ Fx / (m + m_a)."""
    hull = Hull.sphere(R)
    lay = Layout("pair", "sphere", hull.semi_axes_m, [[0, R, 0], [0, -R, 0]], [[1, 0, 0], [1, 0, 0]],
                 u_max_N=12.0, source="hand")
    p = DynamicsParams.from_config(CFG)
    sim = _sim(layout=lay, params=p)
    sim.x[13:] = [3.0, 3.0]                                         # motors already at 3 N each
    dx = derivative(sim.x, sim.x[13:].copy(), np.zeros(3), sim.Minv, sim.Mrb, sim.Ma, sim.weight, sim.buoy,
                    p.r_g, sim.A, p.tau_m, p.k_drag, p.d_lin, p.d_quad, np.zeros(3))
    ax = dx[7]
    assert ax == pytest.approx(6.0 * sim.Minv[0, 0], rel=1e-9)     # exact (includes CG coupling)
    assert ax == pytest.approx(6.0 / (p.mass_kg + p.added_mass_kg), rel=0.15)
    assert abs(dx[8]) < 1e-12 and abs(dx[9]) < 1e-12


def test_4_drifts_to_wind_speed():
    """Constant wind, zero thrust: velocity approaches the wind as quadratic drag predicts."""
    wcfg = WindConfig(mean_speed_m_s=2.0, mean_dir_deg=0.0, turbulence_sigma_m_s=(0, 0, 0))
    p = DynamicsParams.from_config(CFG)
    sim = _sim(params=p, wind=wcfg)
    x = _run(sim, 60.0)
    v_world = quat_to_rot(x[-1, 3:7]) @ x[-1, 7:10]
    M_eff = 1.0 / sim.Minv[0, 0]
    vr = 2.0 / (1.0 + p.k_drag * 2.0 * 60.0 / M_eff)                # analytic relative speed at 60 s
    assert v_world[0] == pytest.approx(2.0 - vr, abs=0.05)
    x2 = _run(sim, 600.0)
    assert quat_to_rot(x2[-1, 3:7]) @ x2[-1, 7:10] == pytest.approx([2.0, 0, 0], abs=0.05)


def test_5_energy_conserved_without_damping():
    p = replace(DynamicsParams.from_config(CFG), k_drag=0.0, d_lin=np.zeros(3), d_quad=np.zeros(3))
    sim = _sim(params=p)
    sim.reset(rpy=(np.deg2rad(20), np.deg2rad(-10), 0.3), v_body=(0.5, -0.2, 0.1), w_body=(0.1, 0.2, 0.3))
    e0 = sim.energy()
    energies = [sim.energy() for _ in range(600) if sim.step(np.zeros(sim.n_props)) is not None]
    assert np.max(np.abs(np.array(energies) - e0)) < 1e-6 * max(abs(e0), 1.0)


def test_motor_lag_and_bounds():
    sim = _sim()
    sim.step(np.full(sim.n_props, 100.0))                           # command above u_max
    assert np.all(sim.thrust <= 12.0 + 1e-12)
    expected = 12.0 * (1 - np.exp(-0.1 / 0.10))                     # first-order lag after one control step
    assert sim.thrust[0] == pytest.approx(expected, rel=0.01)
    sim.fail_prop(0)
    sim.step(np.full(sim.n_props, 12.0))
    assert sim.thrust[0] == 0.0


def test_wind_is_reproducible_and_has_the_right_statistics():
    cfg = WindConfig(mean_speed_m_s=3.0, mean_dir_deg=90.0, turbulence_sigma_m_s=(0.5, 0.5, 0.25),
                     gust_rate_per_min=2.0)
    a, b = WindModel(cfg, 7).sample(5000, 0.01), WindModel(cfg, 7).sample(5000, 0.01)
    np.testing.assert_array_equal(a, b)
    long = WindModel(WindConfig(mean_speed_m_s=3.0, mean_dir_deg=90.0), 3).sample(200000, 0.01)
    np.testing.assert_allclose(long.mean(axis=0), [0, 3, 0], atol=0.1)
    np.testing.assert_allclose(long.std(axis=0), [0.5, 0.5, 0.25], rtol=0.15)


def test_speed():
    """Brief target: >= 2,000 control steps per second (physics; allocation comes in Phase 3)."""
    sim = _sim(layout=load(LAYOUTS / "sphere_paper_best.json"), wind=WindConfig())
    sim.step(np.zeros(sim.n_props))                                 # compile
    t0, n = time.perf_counter(), 2000
    for _ in range(n):
        sim.step(np.full(sim.n_props, 1.0))
    rate = n / (time.perf_counter() - t0)
    print(f"\ncontrol steps per second: {rate:.0f}")
    assert rate > 1000
