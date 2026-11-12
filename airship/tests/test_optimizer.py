from pathlib import Path

import numpy as np
import pytest

from airship.feasibility import Requirements, axis_capacities, shortfall
from airship.geometry import Hull
from airship.layouts import load
from airship.optimizer.objectives import GeometricRules, evaluate_objectives
from airship.optimizer.pareto import knee_index, pareto_filter
from airship.optimizer.pso import LayoutPSO, PSOConfig, run_trials
from airship.optimizer.random_layouts import sample_matched
from airship.wrench import build_M
from conftest import R_REF

LAYOUTS = Path(__file__).resolve().parent.parent / "results" / "layouts"
SHAPES = {"sphere": ("sphere", None), "cigar": ("ellipsoid", (1.0, 1.0, 2.5)),
          "oblate": ("ellipsoid", (2.0, 2.0, 0.5))}


@pytest.mark.parametrize("key", list(SHAPES))
def test_objectives_match_final(final, key):
    """The ported objectives equal Final.py's on the paper layouts."""
    shape, params = SHAPES[key]
    lay = load(LAYOUTS / f"{key}_paper_best.json")
    obj = evaluate_objectives(lay.positions_m, lay.directions, lay.hull, 30,
                              (2, 2, 0.3, 0.3, 0.3, 2), GeometricRules())
    pn, d = lay.positions_m / R_REF, lay.directions
    met = final.controllability_metrics(final.build_M_matrix(pn, d), final.InnerPSOConfig().dof_weights, pn, d)
    assert obj["count"] == pytest.approx(len(pn) / 15)
    assert obj["effort"] == pytest.approx(met["effort"], rel=1e-6)
    assert obj["pos_sym"] == pytest.approx(met["sym_pos"], rel=1e-9)
    assert obj["dir_sym"] == pytest.approx(met["sym_dir"], abs=1e-9)
    assert obj["authority"] == pytest.approx(met["sv_sum"], rel=1e-6)
    assert obj["geometric"] == pytest.approx(
        final.geometric_constraints_metric(pn, d, shape=shape, shape_params=params), rel=1e-9, abs=1e-12)


def _small_run(seed):
    hull = Hull.sphere(R_REF)
    return run_trials(hull, hull.candidate_points(30), PSOConfig(swarm_size=10, iters=5), Requirements(),
                      12.0, GeometricRules(), num_trials=2, seed=seed, log=lambda m: None)


def test_pso_is_deterministic():
    a, b = _small_run(3), _small_run(3)
    for ta, tb in zip(a, b):
        assert ta["eval"].score == tb["eval"].score
        np.testing.assert_array_equal(ta["eval"].positions_m, tb["eval"].positions_m)
        np.testing.assert_array_equal(ta["weights"], tb["weights"])


def test_pso_feasible_results_have_seven_or_more_props():
    for t in _small_run(5) + _small_run(6):
        ev = t["eval"]
        if not ev.gate and ev.shortfall == 0:
            assert ev.positions_m.shape[0] >= 7


def test_decode_points_thrust_outward():
    hull = Hull.sphere(R_REF)
    pso = LayoutPSO(hull, hull.candidate_points(30), np.ones(5) / 5, PSOConfig(), Requirements(),
                    12.0, GeometricRules(), seed=0)
    mask, p, d = pso.decode(pso.X[0])
    assert np.all(np.sum(d * hull.normals(p), axis=1) >= 0)


def test_pareto_filter_matches_final(final):
    pts = np.random.default_rng(0).normal(size=(40, 3))
    assert pareto_filter(pts).tolist() == final.pareto_filter(pts)[1].tolist()


def test_knee_index():
    pts = np.array([[0.0, 1.0], [0.4, 0.4], [1.0, 0.0]])
    assert knee_index(pts) == 1


def test_matched_random_layouts_are_feasible():
    hull = Hull.sphere(R_REF)
    rules, req = GeometricRules(), Requirements()
    samples, attempts = sample_matched(14, hull, hull.candidate_points(30), req, 12.0, rules,
                                       n_target=3, seed=1)
    assert len(samples) == 3 and attempts >= 3
    for p, d in samples:
        assert p.shape == (14, 3)
        assert not hull.gondola_position_violation(p).any()
        dist = np.linalg.norm((p[:, None] - p[None]) / R_REF, axis=2) + np.eye(14) * 9
        assert dist.min() >= rules.d_min(hull)
        assert shortfall(axis_capacities(build_M(p, d), 0, 12.0), req) == 0
