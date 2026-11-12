import numpy as np
import pytest

from airship.geometry import Hull, gondola_direction_violation, unit
from conftest import FINAL_SHAPES, R_REF


def _hull(key):
    shape, _, abc = FINAL_SHAPES[key]
    return Hull(shape, tuple(np.multiply(abc, R_REF)))


@pytest.mark.parametrize("key", list(FINAL_SHAPES))
def test_points_on_surface_and_normals(key):
    hull = _hull(key)
    p = hull.candidate_points(30)
    np.testing.assert_allclose(hull.surface_distance(p), 0, atol=1e-12)
    np.testing.assert_allclose(hull.project(p), p, atol=1e-12)
    n = hull.normals(p)
    np.testing.assert_allclose(np.linalg.norm(n, axis=1), 1, atol=1e-12)
    assert np.all(np.sum(n * p, axis=1) > 0)          # outward


@pytest.mark.parametrize("key", list(FINAL_SHAPES))
def test_normals_and_projection_match_final(final, key):
    shape, params, _ = FINAL_SHAPES[key]
    hull = _hull(key)
    pn = final.fibonacci_sphere(30) if shape == "sphere" else final.fibonacci_ellipsoid(30, *params)
    np.testing.assert_allclose(hull.normals(pn * R_REF), final.surface_normals(pn, shape, params), atol=1e-12)
    off = pn * 1.3                                     # points off the surface
    np.testing.assert_allclose(hull.project(off * R_REF) / R_REF,
                               final.project_onto_surface(off, shape, params), atol=1e-12)


@pytest.mark.parametrize("key", list(FINAL_SHAPES))
def test_gondola_gates_match_final(final, key):
    shape, params, _ = FINAL_SHAPES[key]
    hull = _hull(key)
    pn = final.fibonacci_sphere(30) if shape == "sphere" else final.fibonacci_ellipsoid(30, *params)
    ours = hull.gondola_position_violation(pn * R_REF)
    theirs = [final.gondola_position_hard_check(pn[i:i + 1], shape, params) for i in range(30)]
    assert ours.tolist() == theirs

    d = unit(np.random.default_rng(1).normal(size=(200, 3)))
    theirs_d = [final.gondola_direction_hard_check(d[i:i + 1]) for i in range(200)]
    assert gondola_direction_violation(d).tolist() == theirs_d


def test_sphere_keepout_rejects_bottom_three():
    # paper: 3 of the 30 sphere candidates (i = 27, 28, 29) fall in the 35 deg cone
    hull = Hull.sphere(R_REF)
    assert np.flatnonzero(hull.gondola_position_violation(hull.candidate_points(30))).tolist() == [27, 28, 29]


def test_invalid_hulls():
    with pytest.raises(ValueError):
        Hull("sphere", (1.0, 1.0, 2.0))
    with pytest.raises(ValueError):
        Hull("ellipsoid", (1.0, -1.0, 2.0))
    with pytest.raises(ValueError):
        Hull("cube", (1.0, 1.0, 1.0))
