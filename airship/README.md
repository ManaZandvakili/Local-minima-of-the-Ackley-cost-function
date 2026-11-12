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
| 1 | One-way-thrust feasibility, reliable optimizer, multi-seed Pareto export | **done** |
| 2 | 6-DOF sphere simulator | **done** |
| 3 | Bounded allocation and PID baseline | **done** |
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
| `scripts/export_reference_layouts.py` | writes the paper's best layouts and the STARLab 6-thruster set to `results/layouts/` |
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

## Phase 1: what's here

| File | Purpose |
|---|---|
| `airship/feasibility.py` | bounded-thrust capacity LP; 12 signed-axis capacities (the RL capability vector); `Requirements` and a graded `shortfall`; `assess` adds 50 mixed directions |
| `airship/optimizer/objectives.py` | the paper's five objectives, authority and conditioning, ported from Final.py (identical values, checked in tests) |
| `airship/optimizer/pso.py` | the PSO, ported; rank-6 gate replaced by the capacity test, no top-6 fallback, seeded, no SVD noise |
| `airship/optimizer/pareto.py` | Pareto filter (Final.py rule), knee point, paper pick order |
| `airship/optimizer/random_layouts.py` | matched-N random feasible layouts with the new test |
| `configs/optimizer.yaml` | seeds, PSO settings, capacity requirements, geometric constants |
| `scripts/run_optimizer.py` | runs all seeds in parallel, writes a timestamped run folder, layouts and a report |

```bash
python scripts/run_optimizer.py --workers 4          # ~15-20 min for 5 seeds on 4 cores
python scripts/run_optimizer.py --summarize results/runs/phase1_<timestamp>
```

Layout files also carry `u_min_N` (0 = push-only, negative = reversible) and
`allow_off_hull` (motors on arms or on the gondola).

## Phase 2: what's here

| File | Purpose |
|---|---|
| `airship/massprops.py` | mass, CG and inertia from the mass budget in `configs/sphere.yaml` (neutral-buoyancy trim added automatically) |
| `airship/wind.py` | mean wind + Ornstein-Uhlenbeck turbulence + 1-cosine gusts, all from one seed |
| `airship/dynamics.py` | 6-DOF Fossen-style model (added mass, Coriolis, drag, buoyancy vs weight at the CG, motor lag), RK4 at 100 Hz, numba-compiled; `Simulator` steps at 10 Hz |
| `tests/test_dynamics.py` | the brief's five sanity checks, motor lag, wind statistics and speed |
| `scripts/demo_simulator.py` | figure: pendulum release and zero-thrust drift in gusty wind |

```bash
pip install numba
python scripts/demo_simulator.py          # results/figures/phase2_demo.png
```

## Phase 3: what's here

| File | Purpose |
|---|---|
| `airship/allocation.py` | best-effort bounded allocation: min ‖W(Mu − τ)‖² + ε‖u‖², u_min ≤ u ≤ u_max; exact active-set solver (numba, warm-started, ~16 µs); returns thrusts, achieved wrench, gap, saturation flags |
| `airship/controllers/pid.py` | cascaded PID: world-frame position PID with anti-windup, quaternion PD on roll/pitch, PID on yaw; outputs a body-frame wrench |
| `airship/controllers/tuning.py` | per-layout PSO tuning with Final.py's procedure; identical bounds, cost, scenario and budget for every layout |
| `airship/closed_loop.py` | `Scenario` and `run_episode`: controller → allocation → simulator, with RMS/max position error, yaw error, effort, saturation, wrench gap, survival and recovery time |
| `configs/pid.yaml` | allocation weights, tuning bounds/cost, tuning and acceptance scenarios |
| `scripts/tune_pid.py` | tunes all Phase 3-5 layouts in parallel, then flies the acceptance (2 m/s steady) and a turbulent check |

```bash
python scripts/tune_pid.py --run results/layouts/phase1_20261112_015323 --workers 4   # ~15 min
```

Open questions and assumptions are tracked in `results/REPORT.md`.
