from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import lsq_linear

from airship.allocation import DEFAULT_W, Allocator
from airship.closed_loop import Scenario, run_episode
from airship.config import load_config
from airship.controllers.pid import PIDController, PIDGains, attitude_error
from airship.controllers.tuning import make_pid
from airship.dynamics import DynamicsParams, quat_from_euler
from airship.feasibility import axis_capacities
from airship.layouts import load
from airship.wind import WindConfig

LAYOUTS = Path(__file__).resolve().parent.parent / "results" / "layouts"
LAY = load(LAYOUTS / "sphere_paper_best.json")
PARAMS = DynamicsParams.from_config(load_config("sphere"))


def test_feasible_wrench_is_reproduced():
    al = Allocator(LAY.M(), 0.0, 12.0)
    caps = axis_capacities(LAY.M(), 0.0, 12.0)
    for k, tau in enumerate(np.eye(6)):
        for s in (1, -1):
            target = s * tau * 0.5 * caps[2 * k + (s < 0)]          # half the capacity along each axis
            r = al.solve(target)
            assert np.linalg.norm(LAY.M() @ r.u - target) < 1e-6
            assert np.all(r.u >= -1e-12) and np.all(r.u <= 12 + 1e-12)


def test_infeasible_wrench_gives_bounded_best_effort():
    al = Allocator(LAY.M(), 0.0, 12.0)
    r = al.solve([500.0, 0, 0, 0, 0, 0])
    assert np.all(r.u >= -1e-12) and np.all(r.u <= 12 + 1e-12)
    assert np.linalg.norm(r.gap) > 1.0
    assert r.saturated.any()


def test_matches_scipy_bvls():
    W = np.array(DEFAULT_W)
    M, n, eps = LAY.M(), LAY.n_props, 1e-8
    A = np.vstack([W[:, None] * M, np.sqrt(eps) * np.eye(n)])
    al = Allocator(M, 0.0, 12.0, eps=eps)
    f = lambda u, t: np.sum((W * (M @ u - t)) ** 2) + eps * u @ u
    for t in np.random.default_rng(0).normal(scale=10, size=(100, 6)):
        ref = lsq_linear(A, np.r_[W * t, np.zeros(n)], bounds=(0, 12), method="bvls").x
        assert f(al.solve(t).u, t) <= f(ref, t) + 1e-9


def test_failed_prop_and_reversible_bounds():
    al = Allocator(LAY.M(), 0.0, 12.0)
    hi = np.full(LAY.n_props, 12.0); hi[3] = 0.0
    al.set_bounds(0.0, hi)
    r = al.solve([2.0, -1.0, 0.5, 0.0, 0.1, 0.3])
    assert r.u[3] == 0.0
    star = load(LAYOUTS / "sphere_starlab_6thruster.json")
    al2 = Allocator(star.M(), star.u_min, star.u_max)
    r2 = al2.solve([-3.0, 2.0, 1.0, 0.0, 0.0, -2.0])              # needs reverse thrust
    assert np.linalg.norm(r2.gap[[0, 1, 2, 5]]) < 1e-6 and (r2.u < 0).any()


def test_attitude_error_signs():
    np.testing.assert_allclose(attitude_error(np.array([1.0, 0, 0, 0]), 0.1), [0, 0, 0.1], atol=1e-12)
    e = attitude_error(quat_from_euler(0.2, 0.0, 0.0), 0.0)
    assert e[0] == pytest.approx(-0.2, abs=1e-9)                     # rolled +0.2 -> error -0.2


def test_pid_zero_at_reference():
    c = PIDController(PIDGains(), 30.0, [39.0, 39.0, 21.0])
    x = np.zeros(13 + LAY.n_props); x[3] = 1.0
    np.testing.assert_allclose(c.act(x, np.zeros(3), 0.0), np.zeros(6), atol=1e-12)


def test_pid_holds_station_in_steady_wind():
    sc = Scenario(duration_s=60, wind=WindConfig(mean_speed_m_s=2.0, turbulence_sigma_m_s=(0, 0, 0)))
    r = run_episode(load(LAYOUTS / "sphere_paper_best.json"), PARAMS, make_pid(PIDGains(), PARAMS), sc)
    assert r["survived"] and r["rms_pos_m"] < 0.5


def test_motor_failure_recovery_is_reported():
    sc = Scenario(duration_s=40, wind=WindConfig(mean_speed_m_s=1.0, turbulence_sigma_m_s=(0, 0, 0)),
                  fail_prop=0, fail_time_s=10.0)
    r = run_episode(LAY, PARAMS, make_pid(PIDGains(), PARAMS), sc)
    assert r["survived"] and "recovery_s" in r
