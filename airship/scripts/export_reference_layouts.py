"""Export the reference layouts to the Phase 0 JSON format.

1. The paper's best layouts (sphere 16, cigar 13, oblate 16 props), from the
   best_config_*.npz files that clean_best_configs.py produced by re-running
   Final.py with the paper's cost function (direction symmetry on).
2. The current STARLab 4-motor layout, read from the SFOC drawings.

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


def export_starlab_4prop():
    """Four equatorial BLDC motors at 90 deg spacing (SFOC drawing DWG 2, sheets 1-2).

    The drawings show each propeller disc facing along the hull surface, i.e. a
    horizontal thrust axis tangent to the equator. The push direction (sign)
    of each motor is not stated, so all four are assumed to push
    counter-clockwise seen from above; the motor arms are ignored and the
    motors are placed on the hull surface.
    """
    hull = hull_from_config(CFG)
    az = np.deg2rad([0, 90, 180, 270])
    r_hat = np.column_stack([np.cos(az), np.sin(az), np.zeros(4)])
    dirs = np.cross([0.0, 0.0, 1.0], r_hat)                       # tangential, CCW
    lay = Layout(
        name="sphere_starlab_4prop", shape="sphere", semi_axes_m=hull.semi_axes_m,
        positions_m=r_hat * R, directions=dirs, u_max_N=U_MAX, source="hand",
        notes={
            "origin": "SFOC application, drawing DWG 2 (four station-keeping BLDC motors on the equator)",
            "assumed": "push direction of each motor (all CCW seen from above); motor arm offset ignored; "
                       "u_max estimated (Kavan C3536-1250kv, 11-inch prop; see configs/sphere.yaml)",
            "controllable_dofs": "Fx, Fy, Mz only (all thrust is horizontal)",
        },
    )
    print(f"{save(lay, OUT / f'{lay.name}.json').relative_to(ROOT)}  N={lay.n_props}")


if __name__ == "__main__":
    export_paper_layouts(load_final())
    export_starlab_4prop()
