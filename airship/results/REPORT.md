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

**D. Centre of gravity.** The gondola and electronics weigh 2.79 kg and hang below
the hull; with the lead bars this sets how far the centre of gravity sits below the
hull centre (the pendulum stability). Phase 2 will compute it from the mass budget;
please check the gondola's hanging distance.

**E. Wind range for training (proposal).** The SFOC operating limits are 5.6 m/s
sustained and 8.3 m/s gusts, so the proposal is to train on mean wind 0-5.6 m/s with
gusts up to 8.3 m/s (the brief suggested 0-4 m/s).

**F. Path radius vs operating area.** The 50 m path keeps the comparison with the
paper, but the SFOC limits outdoor tests to a 20 m radius, so it can only be done
in simulation. Station keeping, the primary task, fits inside the limits.

## Assumptions in use

All values are tagged in `configs/sphere.yaml` and `configs/optimizer.yaml`:
air density 1.225 kg/m³, neutral-buoyancy trim, 12 N hull-motor thrust (estimated),
6 N gondola-motor thrust (guess), 0.10 s motor time constant, capacity minima
(1 N, 1 N·m yaw, 0.25 N·m roll/pitch), the 7-inch spacing rule, and the STARLab
thruster geometry above.

Sources: SFOC application (Oct 2026 draft), Aero Drum official quotation (Aug 2025),
ACCESS project page and mass budget, Kavan C3536 product listings.
