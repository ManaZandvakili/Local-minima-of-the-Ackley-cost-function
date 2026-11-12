"""Export the reference layouts to the Phase 0 JSON format.

1. The paper's best layouts (sphere 16, cigar 13, oblate 16 props), from the
   best_config_*.npz files that clean_best_configs.py produced by re-running
   Final.py with the paper's cost function (direction symmetry on).
2. The as-built STARLab 6-thruster set (forward, side and lift pairs, reversible).

Usage (from the airship/ folder):
    python scripts/export_reference_layouts.py
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")

from airship.config import hull_from_config, load_config      # noqa: E402
from airship.geometry import Hull, gondola_direction_violation  # noqa: E402
from airship.layouts import Layout, save                        # noqa: E402

OUT = ROOT / "results" / "layouts"
NPZ_DIR = ROOT.parent                     # best_config_*.npz live at the repo root
CFG = load_config("sphere")
R = CFG["hull"]["diameter_m"] / 2.0
U_MAX = CFG["propulsion"]["u_max_N"]
SHAPES = {                                 # key: (shape, Final.py normalised params)
    "sphere": ("sphere", None),
    "cigar": ("ellipsoid", (1.0, 1.0, 2.5)),
    "oblate": ("ellipsoid", (2.0, 2.0, 0.5)),
}


def load_final():
    spec = importlib.util.spec_from_file_location("final_reference", ROOT / "reference" / "Final.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def paper_objectives(final, pos_norm, dirs, shape, params):
    """The paper's five objectives, computed exactly as Final.py does (normalised units)."""
    np.random.seed(0)                      # Final.py adds tiny random noise before its SVD
    cfg = final.InnerPSOConfig()
    met = final.controllability_metrics(final.build_M_matrix(pos_norm, dirs), cfg.dof_weights, pos_norm, dirs)
    return {
        "J1_count": len(pos_norm) / 15.0,
        "J2_effort": met["effort"],
        "J3_position_symmetry": met["sym_pos"],
        "J4_direction_symmetry": met["sym_dir"],
        "J5_geometric": final.geometric_constraints_metric(pos_norm, dirs, shape=shape, shape_params=params),
    }


def export_paper_layouts(final):
    for key, (shape, params) in SHAPES.items():
        npz = NPZ_DIR / f"best_config_{key}.npz"
        if not npz.exists():
            print(f"[skip] {npz.name} not found")
            continue
        d = np.load(npz)
        pos_norm, dirs = d["positions"], d["directions"]
        abc = np.array(params if params else (1.0, 1.0, 1.0)) * R
        hull = Hull(shape, tuple(abc))
        lay = Layout(
            name=f"{key}_paper_best", shape=shape, semi_axes_m=tuple(abc),
            positions_m=pos_norm * R, directions=dirs, u_max_N=U_MAX, source="optimizer",
            objectives=paper_objectives(final, pos_norm, dirs, shape, params), seed=None,
            notes={
                "origin": "Final.py PSO with the paper's cost (J6 on): best-authority Pareto layout, "
                          "20 trials, seeds t*17+3 (clean_best_configs.py)",
                "objectives_units": "Final.py normalised units (hull radius = 1)",
                "one_way_feasibility": "not checked yet (Phase 1)",
            },
        )
        assert not hull.gondola_position_violation(lay.positions_m).any()
        assert not gondola_direction_violation(lay.directions).any()
        print(f"{save(lay, OUT / f'{lay.name}.json').relative_to(ROOT)}  N={lay.n_props}")


GONDOLA_MOTOR_THRUST_N = 6.0   # guess: 107 g motor + propeller (gondola pair), no data
GONDOLA_DROP_M = 0.25          # assumption: gondola motor height below the hull bottom
GONDOLA_HALF_SPACING_M = 0.25  # assumption: lift motors 0.25 m either side of centre, along x


def export_starlab_6thruster():
    """The as-built STARLab thruster set (project page, drawings, mass budget).

    - Forward pair: hull motors at +-y on the equator, thrust along x.
    - Side pair: hull motors at +-x on the equator, thrust along y
      (sideways together, turning when opposed).
    - Lift pair: two smaller motors on the gondola, thrust along z.
    All motors are reversible (Mana: yes). Body +x is taken as forward; hull
    motors are placed on the hull surface (arm offset ignored).
    """
    hull = hull_from_config(CFG)
    z_g = -(R + GONDOLA_DROP_M)
    pos = np.array([[0, R, 0], [0, -R, 0],                  # forward pair
                    [R, 0, 0], [-R, 0, 0],                  # side pair
                    [GONDOLA_HALF_SPACING_M, 0, z_g], [-GONDOLA_HALF_SPACING_M, 0, z_g]])  # lift pair
    dirs = np.array([[1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]], dtype=float)
    u_max = np.array([U_MAX] * 4 + [GONDOLA_MOTOR_THRUST_N] * 2)
    lay = Layout(
        name="sphere_starlab_6thruster", shape="sphere", semi_axes_m=hull.semi_axes_m,
        positions_m=pos, directions=dirs, u_max_N=u_max, u_min_N=-u_max, allow_off_hull=True,
        source="hand",
        notes={
            "origin": "ACCESS project page (forward, side and lift pairs), SFOC drawing DWG 2, mass budget",
            "assumed": "body +x forward; hull motors on the hull surface (arms ignored); gondola motors "
                       f"{GONDOLA_DROP_M} m below the hull bottom, +-{GONDOLA_HALF_SPACING_M} m along x; "
                       f"hull motors {U_MAX} N (estimated), gondola motors {GONDOLA_MOTOR_THRUST_N} N (guess); "
                       "reverse thrust taken equal to forward thrust",
            "reversible": "yes (Mana)",
        },
    )
    old = OUT / "sphere_starlab_4prop.json"
    if old.exists():
        old.unlink()                                        # superseded by the 6-thruster layout
    print(f"{save(lay, OUT / f'{lay.name}.json').relative_to(ROOT)}  N={lay.n_props}")


if __name__ == "__main__":
    export_paper_layouts(load_final())
    export_starlab_6thruster()
