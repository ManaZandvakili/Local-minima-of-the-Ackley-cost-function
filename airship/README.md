# airship: optimized propeller layouts and closed-loop control (sphere)

Journal-version code, following `RL_SPHERE_BRIEF.md`. Everything is in physical
units (metres, newtons, N·m). The conference-paper code is frozen in
`reference/Final.py` and is never modified; functions ported from it say so.

## Setup

```bash
pip install numpy scipy matplotlib pyyaml pytest
cd airship
python -m pytest -q                              # all tests
```

## Status

| Phase | What | Status |
|---|---|---|
| 0 | Physical-unit wrench matrix, layout JSON format, layout helpers | **done** |
| 1 | One-way-thrust feasibility, reliable optimizer, multi-seed Pareto export | next |
| 2 | 6-DOF sphere simulator | — |
| 3 | Bounded allocation and PID baseline | — |
| 4 | RL environment and training | — |
| 5 | Experiments and paper outputs | — |

## Phase 0: what's here

| File | Purpose |
|---|---|
| `airship/geometry.py` | `Hull` (sphere/ellipsoid in metres): Fibonacci candidate points, normals, projection, gondola keep-out gates |
| `airship/wrench.py` | `build_M(positions_m, directions)`: 6×N wrench matrix (force rows in N/N, torque rows in N·m/N); `scale_torque_rows` for unit-consistent conditioning |
| `airship/layouts.py` | `Layout` + JSON `load` / `save`, `remove_prop` (motor failure), `perturb` (manufacturing tolerance) |
| `airship/config.py` | loads `configs/*.yaml` |
| `configs/sphere.yaml` | hull, mass, propulsion, operating limits; each value tagged confirmed / derived / assumption with its source |
| `scripts/export_reference_layouts.py` | writes the paper's best layouts and the current STARLab 4-motor layout to `results/layouts/` |
| `tests/` | geometry, wrench and layout tests, including the Phase 0 acceptance test against `reference/Final.py` |

Export the reference layouts:

```bash
python scripts/export_reference_layouts.py
```

### Layout JSON format

```json
{
  "name": "sphere_paper_best",
  "shape": "sphere", "semi_axes_m": [1.575, 1.575, 1.575],
  "positions_m": [[x, y, z], ...],
  "directions": [[dx, dy, dz], ...],
  "u_max_N": 12.0,
  "source": "optimizer | random | hand | perturbed | failure",
  "objectives": {"J1_count": 1.067, "J2_effort": 1.207, ...},
  "seed": null,
  "notes": {"origin": "..."}
}
```

Body frame: origin at the hull centre (centre of buoyancy), z up, gondola at −z.
Propellers must lie on the hull surface; thrust directions must be unit vectors.
`u_max_N` may be one number or one value per propeller (a failed motor is 0).

### Units: how this relates to Final.py

`Final.py` stores every position divided by the sphere radius R = 1.575 m
(so its cigar is (1, 1, 2.5) and its oblate (2, 2, 0.5)). Here
`positions_m = positions_normalised × 1.575`, and the matrix built from metres
equals `Final.py`'s matrix with its torque rows multiplied by R, which is what
its simulator does. `tests/test_wrench.py` checks this to 1e-12 for all three hulls.

Open questions and assumptions are tracked in `results/REPORT.md`.
