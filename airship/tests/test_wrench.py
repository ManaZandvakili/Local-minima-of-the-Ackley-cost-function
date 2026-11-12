import numpy as np
import pytest

from airship.geometry import Hull, unit
from airship.wrench import build_M, scale_torque_rows
from conftest import FINAL_SHAPES, R_REF


def _old_M_phys(final, pos_norm, dirs):
    """How Final.py's simulator builds the physical matrix:
    build_M_matrix on normalised positions, then torque rows * R."""
    M = final.build_M_matrix(pos_norm, dirs).copy()
    M[3:] *= final.AIRSHIP_RADIUS_M
    return M


@pytest.mark.parametrize("key", list(FINAL_SHAPES))
def test_M_matches_final_py(final, key):
    """Phase 0 acceptance: old (normalised x R) and new (metres) M are identical."""
    shape, params, abc = FINAL_SHAPES[key]
    pos_norm = (final.fibonacci_sphere(30) if shape == "sphere"
                else final.fibonacci_ellipsoid(30, *params))
    rng = np.random.default_rng(0)
    dirs = unit(rng.normal(size=pos_norm.shape))

    hull = Hull(shape, tuple(np.multiply(abc, R_REF)))
    pos_m = hull.candidate_points(30)
    np.testing.assert_allclose(pos_m, pos_norm * R_REF, rtol=0, atol=1e-12)

    np.testing.assert_allclose(build_M(pos_m, dirs), _old_M_phys(final, pos_norm, dirs),
                               rtol=0, atol=1e-12)


def test_M_matches_final_py_on_paper_layouts(final):
    """Same check on the exported paper layouts (if they have been exported)."""
    from pathlib import Path
    from airship.layouts import load
    paths = sorted((Path(__file__).resolve().parent.parent / "results" / "layouts").glob("*_paper_best.json"))
    if not paths:
        pytest.skip("run scripts/export_reference_layouts.py first")
    for p in paths:
        lay = load(p)
        np.testing.assert_allclose(lay.M(), _old_M_phys(final, lay.positions_m / R_REF, lay.directions),
                                   rtol=0, atol=1e-12)


def test_torque_units():
    # one propeller at +x on a 1.575 m sphere, pushing +y: Fy = 1 N, Mz = +1.575 N*m
    M = build_M([[1.575, 0, 0]], [[0, 1, 0]])
    np.testing.assert_allclose(M[:, 0], [0, 1, 0, 0, 0, 1.575])


def test_scale_torque_rows():
    M = build_M([[1.575, 0, 0]], [[0, 1, 0]])
    Ms = scale_torque_rows(M, 1.575)
    np.testing.assert_allclose(Ms[:, 0], [0, 1, 0, 0, 0, 1.0])
    np.testing.assert_allclose(M[5, 0], 1.575)      # original untouched


def test_bad_shapes_raise():
    with pytest.raises(ValueError):
        build_M(np.zeros((3, 3)), np.zeros((2, 3)))
