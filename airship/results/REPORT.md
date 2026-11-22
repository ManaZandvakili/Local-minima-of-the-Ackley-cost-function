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

## Phase 4 results: RL training (seed 0)

Run `results/rl/20261112_045357_seed0` (evaluation in its `eval/` folder).

**Setup.** One PPO policy (2×256 tanh, Stable-Baselines3) for every layout: it outputs a
wrench (±50 N, ±50 N·m) that goes through the same allocation as the PID. Trained
5 million steps (59 minutes on 4 cores) on 600 layouts (88 Pareto, 60 random, 3
hand-designed and reference layouts, 295 perturbed copies, 150 one-motor-failure
variants); 21 layouts held out. Curriculum: calm air with the knee layout -> calm air,
all layouts -> steady wind -> turbulence, gusts, mid-flight failures and randomised
physics (reached at 434k steps; every stage passed its target).

Two fixes were needed before it learned (first attempt is not reported as a result):
the per-step penalties are now bounded so that surviving always pays more than crashing
(otherwise the agent learned to end episodes early), and the initial exploration noise
was reduced (random ±30 N pushes flung the airship out of the area).

**Training curves** (`eval/training_curves.png`): learning plateaus after about 1 million
steps in the final stage at about 1.1 m RMS and 88 % survival (over all training layouts,
including random, perturbed and failed ones).

**RL vs PID** (gusty station keeping: 3 m/s mean wind from 3 directions, turbulence,
gusts, 45° heading change; identical wind seeds; median over directions, then over layouts):

| Layout group | N | PID RMS | RL RMS | PID effort | RL effort | PID survival | RL survival |
|---|---|---|---|---|---|---|---|
| Optimized knee | 15 | **0.36 m** | 0.41 m | 13,180 | 25,048 | 1.00 | 1.00 |
| Random (20) | 15 | 0.51 m | 0.84 m | 22,271 | 22,912 | 0.85 | 0.87 |
| Optimized paper-rule pick | 22 | **0.29 m** | 0.50 m | 7,245 | 10,423 | 1.00 | 1.00 |
| Random (20) | 22 | 0.36 m | 0.52 m | 10,528 | 16,389 | 1.00 | 1.00 |
| Optimized capped | 9 | **1.02 m** | 3.51 m | 20,837 | 17,751 | 1.00 | 0.67 |
| Random (20) | 9 | 4.69 m | 4.97 m | 15,437 | 17,232 | 0.32 | 0.30 |
| Conference paper layout | 16 | 0.32 m | 0.41 m | 12,304 | 19,691 | 1.00 | 1.00 |
| STARLab 6 thrusters (reversible) | 6 | 0.37 m | 0.62 m | 10,870 | 13,409 | 1.00 | 1.00 |

**Brief acceptance: not met.** On held-in layouts the median RMS position error is
0.49 m with PID and 0.81 m with RL; on held-out layouts 0.49 m (PID) and 1.50 m (RL).
The RL policy is close to PID on the optimized layouts (0.41 vs 0.36 m on the knee)
but uses more effort, and it generalises worse to layouts it has not seen.

**The layout effect holds for both controllers.** With PID and with RL, the optimized
layout beats the random median at the same N (N = 15: 0.36 vs 0.51 m with PID, 0.41 vs
0.84 m with RL; N = 9: 1.02 vs 4.69 m with PID, and 8 in 10 random 9-prop layouts lose
control).

**Proposed next step (brief, section 8):** residual RL, where the policy adds a
correction to the PID wrench (wrench = PID + scale · action). It starts from the PID's
performance and only has to learn improvements, so it usually trains much faster.
Then 5 training seeds for whichever version is used.

Videos (`eval/`): PID vs RL on the knee layout in gusty wind; optimized vs random
9-prop layouts flown by RL; a motor failing mid-flight (PID vs RL).

## Redo with reversible propellers, and residual RL

**Why.** The motors on the STARLab airship spin both ways, so each propeller can push
or pull. Everything above assumed push-only propellers (thrust from 0 to 12 N). The
redo allows -12 to +12 N (`reverse_thrust: true` in `configs/sphere.yaml`; reverse
thrust is assumed as strong as forward thrust, see open question G). All layouts were
re-optimised, the PID re-tuned, and two RL versions trained. The push-only results
above are kept for comparison.

### Optimizer (runs `phase1_bidir_20261117_145131`, capped `phase1_bidir_maxN11_20261117_151234`)

With reversible thrust a layout needs only 6 propellers (push-only needed 7), and every
layout gets about twice the authority.

| Layout | N | Capacity min (weakest direction) | Matched random median (20) |
|---|---|---|---|
| Knee (no cap) | 14 | 49.0 N | 31.8 N |
| Paper-rule pick (no cap) | 23 | 77.0 N | 70.9 N |
| Knee (N ≤ 11) | 8 | 24.0 N | 8.5 N |
| Paper-rule pick (N ≤ 11) | 11 | 34.3 N | 24.4 N |

Feasible trials: 96/100 (no cap). Pareto layouts: N 6-23 (median 12) without the cap,
6-11 (median 9) with it. Under reversible thrust the conference-paper layout reaches
48.3 N, the two hand-designed layouts 21.3 N (cube, 8) and 44.1 N (rings, 12). The
STARLab 6-thruster set still has no roll authority: it relies on the low centre of
gravity, which is enough in practice.

### PID baseline (run `results/pid/20261117_153423`, 66 layouts)

Same tuning procedure as before. Gusty check (3 m/s, turbulence 0.5 m/s, gusts):

| Layout | N | RMS pos | Effort [N²s] | Saturation |
|---|---|---|---|---|
| Optimized knee | 14 | **0.23 m** | 4,553 | 0 % |
| Random (median of 20) | 14 | 0.29 m | 7,274 | 0.2 % |
| Optimized paper-rule pick | 23 | 0.20 m | 2,974 | 0 % |
| Random (median of 20) | 23 | 0.20 m | 3,771 | 0 % |
| Optimized capped knee | 8 | **0.26 m** | 8,088 | 0 % |
| Random (median of 20) | 8 | 0.27 m | 21,815 | 15 % |
| STARLab 6 thrusters | 6 | 0.28 m | 10,783 | 1.9 % |

All 66 layouts pass. Reversible thrust makes every layout much easier to fly (the
push-only capped knee, 9 props, had 1.02 m in the evaluation scenario; the reversible
capped knee, 8 props, has 0.33 m). The optimized
layouts still win, mainly on effort: at N = 8 the random layouts need 2.7× the
effort and saturate 15 % of the time. At N = 23 the optimized and random layouts are
about equal, as expected: with that many propellers any layout has plenty of authority.

### RL: residual vs plain (seed 0, 5 million steps each, about 75 minutes)

- **Plain RL** (`results/rl/20261117_155720_bidir_plain_seed0`): the policy outputs
  the whole wrench (±50 N, ±50 N·m), as before.
- **Residual RL** (`results/rl/20261117_155720_bidir_residual_seed0`): wrench =
  tuned PID + 15 × action (±15 N, ±15 N·m correction). The PID uses each layout's
  tuned gains (the median gains for layouts that were not tuned, mostly Pareto
  layouts other than the knee). The policy also sees the PID's wrench.

Same layout pool for both (660 training layouts, 21 held out), same curriculum, rewards
and PPO settings. Both versions passed every curriculum stage on its target
(full stage by about 310k steps) and finished training with 100 % survival.

**Evaluation** (gusty station keeping: 3 m/s mean wind from 3 directions, turbulence
0.5 m/s, gusts, 45° heading change, identical wind seeds; RMS position error, median over
directions, then over layouts):

| Layout group | N | PID | Residual RL | Plain RL |
|---|---|---|---|---|
| Optimized knee | 14 | 0.292 m | **0.235 m** | 0.264 m |
| Random (20) | 14 | 0.364 m | 0.256 m | **0.232 m** |
| Optimized paper-rule pick | 23 | 0.254 m | **0.227 m** | 0.333 m |
| Random (20) | 23 | 0.256 m | **0.219 m** | 0.315 m |
| Optimized capped knee | 8 | 0.330 m | 0.252 m | **0.238 m** |
| Random (20) | 8 | 0.335 m | **0.299 m** | 0.361 m (1 crash in 60) |
| Conference paper layout | 16 | 0.263 m | **0.224 m** | 0.260 m |
| STARLab 6 thrusters | 6 | 0.368 m | **0.308 m** | 0.358 m |
| **Held-in layouts (median)** | | 0.318 m | **0.252 m** | 0.295 m |
| **Held-out layouts (14, median)** | | 0.283 m | **0.235 m** | 0.319 m |

Effort: residual RL uses about 5 % more than PID; plain RL 5-25 % more.

**Motor failure** (each propeller in turn fails at 15 s, 2 m/s wind; `eval/failures.md`
in the residual run): every controller survives every single failure on the 14-prop
knee, the 8-prop capped knee and the STARLab set. RMS position error after the failure
(median over which motor fails):

| Layout | PID | Residual RL | Plain RL |
|---|---|---|---|
| Optimized knee (14) | 0.070 m | 0.047 m | 0.037 m |
| Optimized capped knee (8) | 0.071 m | 0.052 m | 0.031 m |
| STARLab 6 thrusters | 0.079 m | 0.091 m | 0.289 m |

**What this shows.**

1. **Residual RL beats the PID on every layout group**, by 10-30 %, including the 14
   layouts it never saw in training (0.235 vs 0.283 m). **The brief's acceptance test
   is met** (held-in 0.252 vs 0.318 m).
2. **Plain RL is not reliable.** It is best on some layouts and worse than PID on
   others (the 23-prop layouts, the held-out layouts), and it lost one random 8-prop
   layout. Residual RL keeps the PID's reliability and adds a consistent improvement.
   Residual RL is the recommended controller.
3. **The layout effect is smaller with reversible thrust.** With the PID and with
   residual RL, the optimized knee still beats the random median at the same N
   (PID 0.29 vs 0.36 m, residual RL 0.235 vs 0.256 m), and it needs much less effort
   at small N. But the gap is smaller than with push-only propellers, where random
   9-prop layouts often lost control. With plain RL the random 14-prop layouts were
   even slightly better than the knee, so the layout comparison should be made with
   the PID or residual RL.
4. These are seed-0 results. Phase 5 should repeat residual RL with 5 seeds before the
   numbers go in the paper.

Videos (`eval/` folder of each RL run): PID vs RL on the knee layout in gusty wind,
optimized vs random 8-prop layouts flown by RL, and a motor failing mid-flight.

## Case study: one-way vs reversible propellers

The same full pipeline was run for both propeller types, so the two cases can be
compared directly:

| Step | One-way (push only, 0 to 12 N) | Reversible (-12 to +12 N) |
|---|---|---|
| Optimizer | `phase1_20261112_015323`, capped `phase1_maxN11_20261112_024218` | `phase1_bidir_20261117_145131`, capped `phase1_bidir_maxN11_20261117_151234` |
| Reference layouts | `results/layouts/oneway/` (STARLab set left out: its 6 thrusters need reverse thrust) | `results/layouts/` |
| PID tuning | `results/pid/20261117_203431` (65 layouts) | `results/pid/20261117_153423` (66 layouts) |
| Residual RL | `results/rl/20261117_205532_oneway_residual_seed0` | `results/rl/20261117_155720_bidir_residual_seed0` |
| Plain RL | `results/rl/20261117_205532_oneway_plain_seed0` | `results/rl/20261117_155720_bidir_plain_seed0` |

Identical settings in both cases: optimizer, PID tuning procedure, RL curriculum,
rewards and PPO budget (5 million steps, seed 0), and the same gusty evaluation
(3 m/s mean wind from 3 directions, turbulence 0.5 m/s, gusts, 45° heading change).
The one-way PID numbers reproduce the earlier Phase 3-4 evaluation exactly.

Figures: `results/figures/case_study_rms.png` (RMS position error, log scale) and
`results/figures/case_study_effort.png` (control effort), made by
`scripts/case_study_figure.py`.

### Optimized vs random layouts, same number of propellers

RMS position error (optimized layout vs median of 20 matched random layouts) and the
share of random layouts that kept control in all three wind directions:

| Case | N | PID: optimized | PID: random | Residual RL: optimized | Residual RL: random | Random kept control |
|---|---|---|---|---|---|---|
| One-way | 9 | **1.02 m** | 4.69 m | **0.91 m** | 5.05 m | 32 % |
| One-way | 15 | **0.36 m** | 0.51 m | **0.31 m** | 0.46 m | 85 % |
| One-way | 22 | **0.29 m** | 0.36 m | 0.33 m | 0.34 m | 100 % |
| Reversible | 8 | 0.33 m | 0.34 m | **0.25 m** | 0.30 m | 100 % |
| Reversible | 14 | **0.29 m** | 0.36 m | **0.24 m** | 0.26 m | 100 % |
| Reversible | 23 | 0.25 m | 0.26 m | 0.23 m | 0.22 m | 100 % |

Control effort with the PID [N²s], optimized vs random median:

| Case | N | Optimized | Random | Random saturation (tuning check) |
|---|---|---|---|---|
| One-way | 9 | 20,837 | 15,437* | 24 % |
| One-way | 15 | 13,180 | 22,271 | 8 % |
| One-way | 22 | 7,245 | 10,528 | 2 % |
| Reversible | 8 | 7,516 | 25,021 | 15 % |
| Reversible | 14 | 4,363 | 8,508 | 0.2 % |
| Reversible | 23 | 2,907 | 3,786 | 0 % |

\* Lower only because 68 % of these episodes ended early in a crash, so less thrust was
used in total; it is not an efficiency gain.

### Controllers in each case

Median RMS position error over the layouts:

| Case | Layouts | PID | Residual RL | Plain RL |
|---|---|---|---|---|
| One-way | held in | 0.499 m | **0.442 m** | 0.762 m |
| One-way | held out (14) | **0.488 m** | 0.591 m | 1.025 m |
| Reversible | held in | 0.318 m | **0.252 m** | 0.295 m |
| Reversible | held out (14) | 0.283 m | **0.235 m** | 0.319 m |

Motor failure (each propeller in turn fails at 15 s, 2 m/s wind): survivals, then
RMS position error after the failure (median over which motor fails):

| Case | Layout | PID | Residual RL | Plain RL |
|---|---|---|---|---|
| One-way | knee, 15 | 15/15, 0.089 m | 15/15, 0.059 m | 15/15, 0.125 m |
| One-way | capped knee, 9 | 8/9, 0.088 m | 7/9, 0.175 m | 8/9, 0.428 m |
| Reversible | knee, 14 | 14/14, 0.070 m | 14/14, 0.047 m | 14/14, 0.037 m |
| Reversible | capped knee, 8 | 8/8, 0.071 m | 8/8, 0.052 m | 8/8, 0.031 m |

### What the case study shows

1. **With one-way propellers, the layout decides whether the airship can be
   controlled at all.** At 9 propellers the optimized layout holds about 1 m in gusty
   wind, while two in three random layouts lose control, with the PID and with RL
   alike. No controller fixes a layout that cannot push in every direction. At 15
   propellers 15 % of random layouts still lose control; only at 22 is every random
   layout controllable.
2. **With reversible propellers every layout can be controlled**, and the layout
   mostly decides how efficiently: the optimized layouts reach similar or better
   accuracy with up to 3.3× less effort and no saturation. With many propellers
   (23) any layout is good enough.
3. **The controller matters more when the hardware is good.** Residual RL improves on
   the PID by about 20 % with reversible propellers, on seen and unseen layouts. With
   one-way propellers the gain is smaller (about 10 % on training layouts), it does not
   carry over to unseen layouts (0.59 vs 0.49 m), and on the 9-prop capped knee it
   survived fewer motor failures than the PID (7/9 vs 8/9). Plain RL is worse than PID
   in the one-way case.
4. **The real airship sits between the two cases.** Real propellers usually give
   50-80 % of their thrust in reverse. The closer they get to the reversible case, the
   less the layout matters for controllability and the more it matters for efficiency.

All numbers are seed 0; Phase 5 should repeat the RL runs with 5 seeds.

## Path following: the conference paper's circle

Script `scripts/path_following.py`; results `results/path/20261122_053908_reversible/` and
`results/path/20261122_060340_oneway/` (report.md, eval.csv, figures `path_sphere_knee.png`,
`path_maxN11_sphere_knee.png`). Videos `results/videos/path_*.mp4`.

Test: the paper's reference (Final.py): 50 m circle at 0.5 m/s, 650 s (one lap), nose along the path, seed 7.
Wind 1.5 m/s (the paper's) and 3 m/s, each toward 0/120/240 deg, with turbulence 0.5 m/s and 2 gusts per minute.
The controllers are the station-keeping ones, used as they are: the classic controller gets the
moving reference (position, velocity, acceleration, heading, heading rate), and the learning controllers see their
error relative to the moving target (no retraining). "Completed" = error under 3 m for the whole lap (the paper's
rule). Values below: median over the 3 directions (and over the 20 random layouts).

**3 m/s wind, classic / residual / plain, RMS error [m], completed share**

| Case | Layout | Classic | Residual | Plain | Completed |
|---|---|---|---|---|---|
| reversible | optimized knee, 14 | 0.236 | 0.160 | 0.200 | 100 / 100 / 100% |
| reversible | random 14 (20) | 0.323 | 0.204 | 0.228 | 98 / 100 / 100% |
| reversible | optimized capped, 8 | 0.289 | 0.185 | 0.226 | 100 / 100 / 100% |
| reversible | random 8 (20) | 0.511 | 0.499 | 0.957 | 37 / 42 / 12% |
| one-way | optimized knee, 15 | 0.346 | 0.293 | 0.538 | 100 / 100 / 33% |
| one-way | random 15 (20) | 1.747 | 1.712 | 2.086 | 5 / 7 / 3% |
| one-way | optimized capped, 9 | 1.446 | 1.614 | 2.454 | 0 / 0 / 0% |
| one-way | random 9 (20) | 4.780 | 4.826 | 5.283 | 0 / 0 / 0% |
| reversible | paper best, 16 | 0.200 | 0.138 | 0.213 | 100 / 100 / 100% |
| one-way | paper best, 16 | 0.288 | 0.223 | 0.497 | 100 / 100 / 67% |
| reversible | STARLab 6-thruster | 2.625 | 2.210 | 2.048 | 0 / 0 / 0% |

(For flights that lost the path, the RMS is over the part flown, so it understates them; read the completed share.)

At the paper's 1.5 m/s every optimized layout completes with every controller (0.07 to 0.19 m). The one-way random
layouts already struggle there: 9 propellers 27% completed, 15 propellers 88%.

Findings:
- The station-keeping result carries over to path following. With one-way propellers, the layout decides whether
  the lap can be flown: the optimized 15-propeller layout completes every 3 m/s lap with the classic controller and
  residual learning, while the random 15-propeller layouts complete 5%. With reversible propellers everything
  flies at 14 propellers, but at 8 propellers only the optimized layout completes every lap (random: 37%).
- In 3 m/s wind, residual learning is the best controller on 5 of the 6 optimized layouts, with no retraining for
  a moving target (for example 0.236 to 0.160 m on the reversible knee); the exception is the one-way 9-propeller
  layout, which no controller can fly. At 1.5 m/s plain learning is slightly ahead on two reversible layouts, but
  in 3 m/s it is the least reliable (one-way knee 33% completed).
- Nine one-way propellers are too few for 3 m/s on the circle: even the optimized layout falls more than 3 m
  behind (but it survives longer than random layouts: 500 s vs about 156 s in the 0 deg video).
- The STARLab 6-thruster layout cannot fly the circle in 3 m/s wind (2 to 2.6 m error); it completes at 1.5 m/s.

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

**G. Reverse thrust.** The redo assumes a propeller pulls as hard as it pushes
(-12 to +12 N). Real propellers are pitched for one direction and usually give less
thrust in reverse (often 50-80 %). `reverse_thrust_frac` in `configs/sphere.yaml`
sets this; a thrust-stand test in both directions would settle it.

## Assumptions in use

All values are tagged in `configs/sphere.yaml` and `configs/optimizer.yaml`:
air density 1.225 kg/m³, neutral-buoyancy trim (component positions above),
rotational damping (linear 2, quadratic 5), OU turbulence 0.5/0.5/0.25 m/s with 2 s
correlation time, 12 N hull-motor thrust (estimated; reversible, equal in reverse since the redo),
6 N gondola-motor thrust (guess), 0.10 s motor time constant, capacity minima
(1 N, 1 N·m yaw, 0.25 N·m roll/pitch), the 7-inch spacing rule, and the STARLab
thruster geometry above.

Sources: SFOC application (Oct 2026 draft), Aero Drum official quotation (Aug 2025),
ACCESS project page and mass budget, Kavan C3536 product listings.
