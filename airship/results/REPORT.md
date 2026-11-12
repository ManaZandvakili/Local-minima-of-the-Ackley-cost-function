# Report: journal version (sphere)

Plain-language log of results, decisions, assumptions and open questions.

## Decisions from Mana

| # | Question | Answer |
|---|---|---|
| 1 | Hull size | 3.15 m diameter |
| 2 | Mass | mass budget: 16.16 kg with helium (13.29 kg without), about 3.9 kg of spare lift trimmed with lead bars and the ballonet; modelled as neutrally buoyant, about 20 kg in total |
| 3 | Objective set | the paper's five: count, effort, position symmetry, direction symmetry, geometric |
| 4 | Path-following radius | 50 m |
| 5 | Thrust | hull motors: Kavan C3536 (1250 kV in the quotation, 1000 kV in the mass budget), 11-inch props, 4S battery. No published thrust figure: **12 N per motor estimated** |
| 6 | STARLab reference layout | six reversible thrusters: forward pair and side pair on the hull, lift pair on the gondola |
| 7 | Compute for RL | the cloud, for now |

Base code: `Final.py` (the paper version), frozen in `reference/`.

## Phase 0 results

- New package `airship/`: geometry, wrench matrix and layouts in physical units.
- **Acceptance test passes:** the wrench matrix built from metres equals Final.py's
  (normalised positions, torque rows × 1.575 m) to within 1e-12, for all three hulls.
- Reference layouts in `results/layouts/`: the paper's best sphere (16), cigar (13)
  and oblate (16) layouts, and the STARLab 6-thruster set.

## Phase 1 results (run `phase1_20261112_015323`)

Full report: `results/runs/phase1_20261112_015323/report.md`.
Layouts: `results/layouts/phase1_20261112_015323/` (61 Pareto layouts, the two
chosen layouts, 40 matched random layouts).

**What changed in the optimizer.** The rank-6 test is replaced by a push-only
capacity test: for each of the 12 directions (±Fx, ±Fy, ±Fz, ±roll, ±pitch, ±yaw),
the largest force or moment the layout can produce with every propeller between
0 and 12 N. A layout is accepted only if it reaches at least 1 N (forces), 1 N·m
(yaw) and 0.25 N·m (roll, pitch) in every direction. The "force 6 propellers"
fallback is removed, the random SVD noise is removed, and every random number
comes from a seed. The five objectives are unchanged and give exactly the same
values as Final.py (tested).

**Runs.** 5 seeds × 20 PSO trials, about 9-10 minutes per seed.
- 96 of 100 trials found a feasible layout. The other 4 never found a layout that
  passes the gondola keep-out rules (an all-or-nothing rejection inherited from Final.py).
- Every feasible layout has at least 7 propellers, as push-only control requires.

| Seed | Feasible trials | Pareto layouts | Knee N | Paper-rule pick N |
|---|---|---|---|---|
| 0 | 20/20 | 9 | 15 | 20 |
| 1 | 17/20 | 9 | 15 | 19 |
| 2 | 20/20 | 14 | 12 | 19 |
| 3 | 20/20 | 14 | 16 | 21 |
| 4 | 19/20 | 13 | 15 | 22 |

Spread across seeds: the knee point lands at 12-16 propellers (15 in three of five
seeds). Pareto layouts range from 7 to 23 propellers (median 13), because each
trial weights the five objectives differently.

**Chosen layouts (over all seeds):**

| Layout | N | Effort | Pos. sym. | Dir. sym. | Geometric | Capacity min |
|---|---|---|---|---|---|---|
| sphere_knee (best compromise) | 15 | 1.096 | 0.102 | 0.056 | 0.818 | 6.2 N |
| sphere_paper_pick (paper's rule: highest planar authority) | 22 | 0.760 | 0.106 | 0.075 | 4.787 | 19.7 N |
| sphere_paper_best (conference paper, for comparison) | 16 | 1.207 | 0.119 | 0.029 | 0.000 | 10.0 N |

**Against matched random layouts** (20 each, same N, all passing the same tests):
- The knee layout (N = 15) beats **all 20** random layouts on effort (1.10 vs a
  median of 2.11), direction symmetry, geometric suitability and authority, and
  18 of 20 on position symmetry.
- Its worst-direction capacity (6.2 N) is about the same as the random median
  (6.1 N). Capacity is a constraint here, not an objective, so the optimizer only
  makes sure it is high enough; it does not maximise it.
- The paper-rule pick (N = 22) beats all 20 of its random layouts on effort,
  direction symmetry, geometric suitability and authority.

**Reference layouts under the new test:**

| Layout | Feasible (all 6 DOF) | Capacity min | Missing directions |
|---|---|---|---|
| sphere_paper_best (16) | yes | 10.0 N | none |
| cigar_paper_best (13) | yes | 15.4 N | none |
| oblate_paper_best (16) | yes | 24.6 N | none |
| sphere_starlab_6thruster (reversible) | no | 0 | roll (±Mx) |

So the conference-paper layouts still hold up with push-only propellers. The
current STARLab thruster set cannot control roll at all (expected: no thruster
produces a roll moment; the gondola pendulum is what keeps it level).

### Decisions needed before Phase 2

1. **Which optimized layout to use as "the" optimized layout.** The paper's rule
   (highest planar authority) now picks the layout with the most propellers (22),
   because authority grows with the number of propellers. **Recommendation: use the
   knee layout (15 propellers)** as the main optimized layout, and keep the 22-prop
   pick and the minimum-N Pareto layout (7 props) as extra comparison points.
2. **Spacing rule.** The geometric objective still uses the paper's 7-inch
   propeller for the spacing rule (2.5 diameters = 0.44 m). The real propellers are
   11-inch (0.70 m). Kept as in the paper for now; switching to 11-inch would make
   the optimizer keep propellers further apart.
3. **The 4 rejected trials.** Making the gondola rejection graded (still never
   accepted, but with a gradient towards passing) would let those trials recover.
   It changes the search slightly compared with the paper.

## Phase 2 results: the 6-DOF simulator

**Mass model (from the mass budget).** Components 16.16 kg, plus 3.93 kg of trim
ballast on the lower hull for neutral buoyancy: 20.09 kg in total. The centre of
gravity is **0.50 m below the hull centre** (the gondola, blower and ballast below,
the Starlink dummy and GPS on top). Inertia about the hull centre: 39.0 kg·m² in
roll and pitch, 20.8 kg·m² in yaw. Added mass 10.0 kg on each axis.

**Sanity checks (all pass):**

| Check | Result |
|---|---|
| 1. At rest, no wind, no thrust, 60 s | no drift (0 m) |
| 2. Released from 10° roll | swing period 3.79 s vs 3.78 s predicted (0.2 %); settles |
| 3. Push from a symmetric motor pair | acceleration equals the model's exact value; 9 % above Fx/(m + m_a) because the low CG couples surge with pitch |
| 4. Constant 2 m/s wind, no thrust | 1.77 m/s after 60 s vs 1.79 m/s from the drag formula; reaches 2 m/s |
| 5. No damping, no wind, no thrust | energy conserved to 1e-8 J over 60 s |

Speed: about **6,100 control steps per second** on one core (brief target 2,000;
allocation is added in Phase 3). Demo figure: `results/figures/phase2_demo.png`.

**Simplifications** (documented in `airship/dynamics.py`): the wind is held constant
within each 0.01 s physics step; the force from the air itself accelerating
(Froude-Krylov) is not modelled; helium is counted as a solid sphere in the inertia
(as the brief suggests; real gas lags behind, so true roll/yaw inertia is somewhat lower).

### Practical issue found: motor mass vs lift

Each hull motor with its mount and ESC weighs about **0.55 kg** (mass budget: 2.08 kg
for four mounts with motors and propellers, 0.13 kg for four ESCs). The airship has
3.93 kg of spare lift. So with these motors it can carry **at most about 11 hull
motors** (the 4 it has plus 7), with no ballast left.
- The 15-propeller optimized layout would be about **2.1 kg too heavy**.
- The 22-propeller layout would be about **9.9 kg too heavy**.

The simulator currently uses the **same mass for every layout**, so comparisons
isolate the effect of where the propellers are and which way they point (matched
random layouts have the same N anyway). Options for the real airship: smaller
motors, a larger hull, or limiting the optimizer to N <= 11.

## Phase 3 results: allocation and PID baseline

**Allocation.** Weighted bounded least squares (priorities Fx, Fy, yaw > Fz > roll,
pitch), solved exactly by a small active-set solver: same answers as scipy's
`lsq_linear` (to 1e-14), 16 µs per step instead of 600 µs. Reproduces any achievable
wrench to 1e-6; for unachievable ones it returns the closest achievable wrench and
reports the gap. Failed motors are handled by setting their bounds to zero.

**PID.** Position PID in the world frame (with anti-windup), quaternion PD on roll and
pitch, PID on yaw; output is a wrench into the same allocation the RL agent will use.
Tuned **per layout** with the same PSO procedure as Final.py (12 particles, 30
iterations), identical bounds, cost (5·RMS position + 1·RMS yaw + 0.001·∫Σu²dt) and
tuning scenario (1.5 m/s wind with light turbulence, 1 m initial offset, 45° heading
change). The full loop runs at about 3,400 control steps per second.

**Acceptance (brief: RMS < 0.5 m in 2 m/s steady wind):** passed by every layout;
the 15-propeller knee layout holds **0.086 m**.

**First comparison: one gusty run, 3 m/s mean wind, turbulence and gusts, one seed**
(full statistics come in Phase 5; run `results/pid/20261112_030317`):

| Layout | N | RMS position | Max position | Effort [N²s] | Saturation | Survived |
|---|---|---|---|---|---|---|
| Optimized knee | 15 | **0.30 m** | 0.71 m | 24,186 | 8.5 % | yes |
| Random, median of 20 | 15 | 0.47 m | 1.59 m | 17,967 | 8.1 % | 17/20 |
| Optimized paper-rule pick | 22 | **0.23 m** | 0.52 m | **6,245** | 0.2 % | yes |
| Random, median of 20 | 22 | 0.28 m | 0.72 m | 17,020 | 2.0 % | 20/20 |
| Optimized, capped (buildable) | 9 | **0.29 m** | 1.56 m | 16,832 | 6.6 % | yes |
| Random, median of 20 | 9 | 3.67 m | 8.84 m | 20,328 | 23.8 % | **12/20** |
| Minimum-N Pareto layout | 7 | 0.79 m | 3.64 m | 22,914 | 24.1 % | yes |
| Conference paper layout | 16 | 0.26 m | 0.68 m | 19,256 | 3.7 % | yes |
| STARLab 6 thrusters (reversible) | 6 | 0.28 m | 1.01 m | 10,783 | 1.9 % | yes |

What this shows so far:
- In **steady** wind every layout tracks about the same (when the motors can deliver
  what the controller asks, the layout only changes the effort). The differences
  appear in **gusty** wind, when motors start hitting their limits.
- The optimized layout beats the typical random layout with the same number of
  propellers at every N. The gap is largest with few propellers: at **N = 9**, the
  optimized layout holds 0.29 m while the typical random layout drifts 3.7 m and
  8 of 20 random layouts lose control.
- The knee layout's effort is higher than the random median in this single run.
  Effort depends strongly on each layout's tuned gains (the tuning cost weighs effort
  lightly), so this needs the multi-seed comparison in Phase 5 before drawing a
  conclusion.
- The **STARLab layout does well in station keeping** with only 6 thrusters, because
  its thrusters are reversible and point straight along x, y and z. It cannot control
  roll at all (the pendulum keeps it level), and it uses a different kind of motor
  (reversible), so it is not a like-for-like comparison with push-only layouts.

## Open questions

**B. Maximum thrust per propeller.** No published thrust figure exists for the
Kavan C3536 motors. Motors of this class with an 11-inch propeller give roughly
12-13 N at full throttle, so 12 N is used. The quotation lists 1250 kV, the mass
budget 1000 kV; with the 4S battery, 1000 kV is the better match for an 11-inch
propeller. A thrust-stand measurement would settle it. In Phase 1 the value only
scales the capacity numbers; from Phase 2 on it matters a lot.

**C. STARLab thrusters.** Modelled as: forward pair at ±y on the equator pushing
along x, side pair at ±x pushing along y, lift pair on the gondola (0.25 m below
the hull, ±0.25 m along x) pushing along z, all reversible. Gondola motor thrust is
unknown and guessed at 6 N (107 g motor and propeller). Positions of the gondola
motors and which way is "forward" are assumptions.

**D. Centre of gravity.** Computed from the mass budget: 0.50 m below the hull centre.
Assumed: gondola 0.25 m below the hull bottom, trim ballast on the lower hull at
z = -1.45 m, the 1 kg "mass dummy" and GPS on top (Starlink). Please check these
positions; they set the pendulum stability.

**E. Wind range for training (proposal).** The SFOC operating limits are 5.6 m/s
sustained and 8.3 m/s gusts, so the proposal is to train on mean wind 0-5.6 m/s with
gusts up to 8.3 m/s (the brief suggested 0-4 m/s).

**F. Path radius vs operating area.** The 50 m path keeps the comparison with the
paper, but the SFOC limits outdoor tests to a 20 m radius, so it can only be done
in simulation. Station keeping, the primary task, fits inside the limits.

## Assumptions in use

All values are tagged in `configs/sphere.yaml` and `configs/optimizer.yaml`:
air density 1.225 kg/m³, neutral-buoyancy trim (component positions above),
rotational damping (linear 2, quadratic 5), OU turbulence 0.5/0.5/0.25 m/s with 2 s
correlation time, 12 N hull-motor thrust (estimated),
6 N gondola-motor thrust (guess), 0.10 s motor time constant, capacity minima
(1 N, 1 N·m yaw, 0.25 N·m roll/pitch), the 7-inch spacing rule, and the STARLab
thruster geometry above.

Sources: SFOC application (Oct 2026 draft), Aero Drum official quotation (Aug 2025),
ACCESS project page and mass budget, Kavan C3536 product listings.
