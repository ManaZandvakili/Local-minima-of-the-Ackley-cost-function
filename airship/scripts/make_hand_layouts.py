"""Hand-designed symmetric reference layouts (brief, section 8): what an engineer
might pick without an optimizer. Each passes the push-only feasibility test.

    sphere_hand_cube8:   8 props at the cube vertices, thrust tangent to the hull,
                         tilted alternately up and down
    sphere_hand_rings12: two rings of 6 (above and below the equator, staggered 30 deg),
                         thrust tangent, alternately clockwise / counter-clockwise, tilted

Usage (from the airship/ folder):
    python scripts/make_hand_layouts.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from airship.config import hull_from_config, load_config      # noqa: E402
from airship.feasibility import Requirements, assess          # noqa: E402
from airship.geometry import gondola_direction_violation, unit  # noqa: E402
from airship.layouts import Layout, save                      # noqa: E402

CFG = load_config("sphere")
HULL = hull_from_config(CFG)
R = HULL.semi_axes_m[0]
U_MAX = CFG["propulsion"]["u_max_N"]


def tangent_dirs(p_hat, ccw, tilt_rad):
    """Unit thrust directions tangent to the hull: horizontal tangent (+-) rotated
    by tilt about the outward normal."""
    z = np.array([0.0, 0.0, 1.0])
    out = []
    for p, s, t in zip(p_hat, ccw, tilt_rad):
        h = unit(np.cross(z, p)) * s                 # horizontal tangent
        v = unit(np.cross(p, h))                     # tangent pointing "up" along the hull
        out.append(unit(np.cos(t) * h + np.sin(t) * v))
    return np.array(out)


def cube8(tilt_deg):
    p = unit(np.array([[x, y, z] for z in (1, -1) for x in (1, -1) for y in (1, -1)], dtype=float))
    # sign pattern found by a small random search (seed 0, 3000 tries) for push-only feasibility
    ccw = np.array([-1, 1, 1, -1, -1, 1, 1, -1])
    tilt = np.deg2rad(tilt_deg) * np.array([1, 1, 1, -1, -1, -1, -1, 1])
    return p * R, tangent_dirs(p, ccw, tilt)


def rings12(tilt_deg, z_up=0.45, z_lo=-0.35):
    pts, ccw, tilt = [], [], []
    for k, (zz, off) in enumerate(((z_up, 0.0), (z_lo, 30.0))):
        r = np.sqrt(1 - zz**2)
        for j in range(6):
            a = np.deg2rad(off + 60 * j)
            pts.append([r * np.cos(a), r * np.sin(a), zz])
            ccw.append(1 if j % 2 == 0 else -1)
            tilt.append(np.deg2rad(tilt_deg) * (1 if (j + k) % 2 == 0 else -1))
    p = unit(np.array(pts))
    return p * R, tangent_dirs(p, np.array(ccw), np.array(tilt))


def first_feasible(make, tilts, name, note):
    req = Requirements()
    for t in tilts:
        p, d = make(t)
        if HULL.gondola_position_violation(p).any() or gondola_direction_violation(d).any():
            continue
        L = Layout(name, "sphere", HULL.semi_axes_m, p, d, u_max_N=U_MAX, source="hand",
                   notes={"design": note, "tilt_deg": float(t)})
        fr = assess(L.M(), 0.0, U_MAX, req, R)
        if fr.feasible:
            L = Layout.from_dict({**L.to_dict(), "notes": {**L.notes, "feasibility": fr.as_dict()}})
            print(f"{save(L, ROOT / 'results' / 'layouts' / f'{name}.json').relative_to(ROOT)}  "
                  f"N={L.n_props} tilt={t} deg capacity min {fr.capacity_min_N:.2f} N")
            return L
    print(f"[warn] no feasible {name} for tilts {list(tilts)}")


if __name__ == "__main__":
    first_feasible(cube8, [35, 30, 40, 45, 25], "sphere_hand_cube8",
                   "8 props at the cube vertices, tangential thrust tilted alternately up/down")
    first_feasible(rings12, [30, 40, 45, 35, 25, 20, 50], "sphere_hand_rings12",
                   "two staggered rings of 6, tangential thrust alternating CW/CCW and tilted")
