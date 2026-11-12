import numpy as np
import pytest

from airship.geometry import Hull, unit
from airship.layouts import Layout, load, perturb, remove_prop, save

R = 1.575


def _layout(n=12, seed=0, u_max=2.5):
    hull = Hull.sphere(R)
    p = hull.candidate_points(30)[:n]
    d = unit(np.random.default_rng(seed).normal(size=(n, 3)))
    return Layout("test", "sphere", hull.semi_axes_m, p, d, u_max_N=u_max,
                  source="optimizer", objectives={"J1_count": n / 15}, seed=seed)


def test_save_load_roundtrip(tmp_path):
    lay = _layout()
    back = load(save(lay, tmp_path / "x.json"))
    np.testing.assert_array_equal(back.positions_m, lay.positions_m)
    np.testing.assert_array_equal(back.directions, lay.directions)
    assert back.to_dict() == lay.to_dict()


def test_per_prop_u_max_roundtrip(tmp_path):
    lay = _layout(n=5, u_max=[2.5, 2.5, 0.0, 2.5, 1.0])
    back = load(save(lay, tmp_path / "x.json"))
    np.testing.assert_array_equal(back.u_max, [2.5, 2.5, 0.0, 2.5, 1.0])


def test_remove_prop():
    lay = _layout()
    f = remove_prop(lay, 3)
    assert f.n_props == lay.n_props - 1
    assert f.source == "failure" and f.notes == {"parent": "test", "removed_index": 3}
    np.testing.assert_array_equal(f.M(), np.delete(lay.M(), 3, axis=1))
    with pytest.raises(IndexError):
        remove_prop(lay, 12)


def test_remove_prop_per_prop_u_max():
    lay = _layout(n=4, u_max=[1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(remove_prop(lay, 1).u_max, [1.0, 3.0, 4.0])


def test_perturb_stays_valid_and_is_reproducible():
    lay = _layout()
    a = perturb(lay, 0.02, 2.0, np.random.default_rng(7))
    b = perturb(lay, 0.02, 2.0, np.random.default_rng(7))
    np.testing.assert_array_equal(a.positions_m, b.positions_m)
    np.testing.assert_array_equal(a.directions, b.directions)
    np.testing.assert_allclose(a.hull.surface_distance(a.positions_m), 0, atol=1e-9)
    np.testing.assert_allclose(np.linalg.norm(a.directions, axis=1), 1, atol=1e-12)
    assert a.source == "perturbed" and a.notes["parent"] == "test"


def test_perturb_magnitudes():
    lay = _layout(n=27)
    rng = np.random.default_rng(3)
    shifts, tilts = [], []
    for _ in range(40):
        a = perturb(lay, 0.02, 2.0, rng)
        shifts.append(np.linalg.norm(a.positions_m - lay.positions_m, axis=1))
        tilts.append(np.degrees(np.arccos(np.clip(np.sum(a.directions * lay.directions, axis=1), -1, 1))))
    # tangential 2-D Gaussian: mean shift = sigma*sqrt(pi/2) ~ 0.025 m; mean |tilt| = sigma*sqrt(2/pi) ~ 1.6 deg
    assert 0.020 < np.mean(shifts) < 0.030
    assert 1.4 < np.mean(tilts) < 1.8


def test_zero_perturbation_is_identity():
    lay = _layout()
    a = perturb(lay, 0.0, 0.0, np.random.default_rng(0))
    np.testing.assert_allclose(a.positions_m, lay.positions_m, atol=1e-12)
    np.testing.assert_allclose(a.directions, lay.directions, atol=1e-12)


def test_validation():
    hull = Hull.sphere(R)
    p = hull.candidate_points(4)
    d = unit(np.ones((4, 3)))
    with pytest.raises(ValueError, match="off the hull"):
        Layout("x", "sphere", hull.semi_axes_m, p * 1.01, d)
    with pytest.raises(ValueError, match="unit vectors"):
        Layout("x", "sphere", hull.semi_axes_m, p, d * 2)
    with pytest.raises(ValueError, match="source"):
        Layout("x", "sphere", hull.semi_axes_m, p, d, source="magic")
    with pytest.raises(ValueError, match="u_max_N"):
        Layout("x", "sphere", hull.semi_axes_m, p, d, u_max_N=[1, 2])
