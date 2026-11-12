import numpy as np
import pytest

from airship.feasibility import (AXIS_LABELS, Requirements, assess, axis_capacities, capacity,
                                 sample_directions, shortfall)
from airship.geometry import Hull, unit
from airship.wrench import build_M

R = 1.575


def test_two_opposed_props_at_centre():
    # two props at the origin pushing +x and -x: +-Fx up to u_max, nothing else
    M = build_M([[0, 0, 0], [0, 0, 0]], [[1, 0, 0], [-1, 0, 0]])
    caps = dict(zip(AXIS_LABELS, axis_capacities(M, 0.0, 2.0)))
    assert caps["+Fx"] == pytest.approx(2.0) and caps["-Fx"] == pytest.approx(2.0)
    assert all(caps[k] == 0.0 for k in AXIS_LABELS if k[1:] != "Fx")


def test_reversible_thrust():
    M = build_M([[0, 0, 0]], [[1, 0, 0]])
    assert capacity(M, [-1, 0, 0, 0, 0, 0], 0.0, 3.0) == 0.0           # push-only: cannot pull
    assert capacity(M, [-1, 0, 0, 0, 0, 0], -3.0, 3.0) == pytest.approx(3.0)


def test_four_tangential_push_only_props_yaw_one_way():
    az = np.deg2rad([0, 90, 180, 270])
    r = np.column_stack([np.cos(az), np.sin(az), np.zeros(4)]) * R
    d = np.cross([0, 0, 1.0], r / R)
    caps = dict(zip(AXIS_LABELS, axis_capacities(build_M(r, d), 0.0, 12.0)))
    assert caps["+Mz"] == pytest.approx(4 * 12 * R)
    assert caps["-Mz"] == 0.0 and caps["+Fx"] == 0.0 and caps["-Fx"] == 0.0


def test_bounds_must_include_zero():
    M = build_M([[0, 0, 0]], [[1, 0, 0]])
    with pytest.raises(ValueError):
        capacity(M, [1, 0, 0, 0, 0, 0], 0.5, 2.0)


def test_axis_capacities_scale_with_u_max():
    hull = Hull.sphere(R)
    p = hull.candidate_points(30)[:20]
    d = unit(np.random.default_rng(4).normal(size=(20, 3)))
    M = build_M(p, d)
    np.testing.assert_allclose(axis_capacities(M, 0, 12.0), 4.8 * axis_capacities(M, 0, 2.5), rtol=1e-6)


def test_all_axes_reachable_implies_every_direction():
    # positive capacity on all 12 signed axes => every sampled direction reachable
    hull = Hull.sphere(R)
    rng = np.random.default_rng(11)
    for _ in range(20):
        p = hull.candidate_points(30)[rng.choice(27, 14, replace=False)]
        d = unit(rng.normal(size=(14, 3)))
        M = build_M(p, d)
        if np.all(axis_capacities(M, 0, 12) > 1e-9):
            for e in sample_directions(30, seed=1):
                assert capacity(M, e, 0, 12) > 1e-9


def test_fewer_than_seven_push_only_props_never_feasible():
    hull = Hull.sphere(R)
    rng = np.random.default_rng(2)
    req = Requirements()
    for _ in range(50):
        p = hull.candidate_points(30)[rng.choice(27, 6, replace=False)]
        d = unit(rng.normal(size=(6, 3)))
        assert shortfall(axis_capacities(build_M(p, d), 0, 12), req) > 0


def test_requirements_and_shortfall():
    req = Requirements(min_force_N=2.0, min_yaw_Nm=1.0, min_roll_pitch_Nm=0.5)
    m = req.minima()
    np.testing.assert_allclose(m, [2, 2, 2, 2, 2, 2, .5, .5, .5, .5, 1, 1])
    caps = m * 2
    assert shortfall(caps, req) == 0.0
    caps[0] = 1.0                       # +Fx at half the minimum
    assert shortfall(caps, req) == pytest.approx(0.5)
    caps[1] = 0.0                       # -Fx missing entirely
    assert shortfall(caps, req) == pytest.approx(1.5)
    planar = Requirements(required_dofs=("Fx", "Fy", "Mz"))
    assert np.isnan(planar.minima()[4:10]).all()


def test_assess_paper_sphere_layout():
    from airship.layouts import load
    from pathlib import Path
    lay = load(Path(__file__).resolve().parent.parent / "results" / "layouts" / "sphere_paper_best.json")
    r = assess(lay.M(), 0.0, 12.0, Requirements(), R, n_sampled=20)
    assert r.feasible and r.capacity_min_N > 0 and r.shortfall == 0.0
