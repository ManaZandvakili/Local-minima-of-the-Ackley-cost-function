# Report: journal version (sphere)

Plain-language log of results, decisions, assumptions and open questions.

## Decisions from Mana (Section 11)

| # | Question | Answer |
|---|---|---|
| 1 | Hull size | 3.15 m diameter |
| 2 | Mass | from the SFOC application: final operating mass 16 kg, gross lift 20 kg (see open question A) |
| 3 | Objective set | the paper's five: count, effort, position symmetry, direction symmetry, geometric |
| 4 | Path-following radius | 50 m |
| 5 | Thrust | from the Aero Drum quotation: 4 x Kavan C3536-1250kv, 11-inch props, 50 A ESCs. No published thrust figure, so **12 N per motor is estimated** (open question B) |
| 6 | STARLab 4-motor layout | include as a reference (exported; see open question C) |
| 7 | Compute for RL | the cloud, for now |

Base code: `Final.py` (the paper version), frozen in `reference/`.

## Phase 0 results

- New package `airship/` with geometry, wrench matrix and layouts in physical units.
- **Acceptance test passes:** the wrench matrix built from metres equals Final.py's
  (normalised positions, torque rows × 1.575 m) to within 1e-12, for the sphere,
  cigar and oblate. 26 tests pass.
- Reference layouts exported to `results/layouts/`:

| Layout | N | rank of M | Count | Effort | Pos. sym. | Dir. sym. | Geometric |
|---|---|---|---|---|---|---|---|
| sphere_paper_best | 16 | 6 | 1.067 | 1.207 | 0.119 | 0.029 | 0.000 |
| cigar_paper_best | 13 | 6 | 0.867 | 1.393 | 0.407 | 0.046 | 0.124 |
| oblate_paper_best | 16 | 6 | 1.067 | 1.097 | 0.050 | 0.000 | 2.633 |
| sphere_starlab_4prop | 4 | 3 | — | — | — | — | — |

Objective values are in Final.py's normalised units, so they match the paper.
Whether each layout can really fly with push-only propellers is checked in Phase 1.

## Open questions for Mana

**A. Weight vs buoyancy.** The SFOC application lists a final operating mass of
16 kg and a gross lift of 20 kg. The 20 kg equals the mass of air the hull
displaces (1.225 kg/m³ × 16.4 m³ = 20.1 kg). The hull holds about 2.8 kg of helium.
- If the 16 kg *excludes* helium, the airship weighs about 18.8 kg against 20.1 kg of
  buoyancy, i.e. about 1.3 kg of net lift (it pulls up on its tethers).
- If it *includes* helium, the net lift is about 4 kg.
- Arctic air at −15 °C is denser (≈ 1.37 kg/m³), raising buoyancy to about 22.4 kg.

The simulator (Phase 2) models free flight, so for now it will assume the airship is
**trimmed to neutral buoyancy** (ballonet and ballast). Please confirm, or say if
you want net lift and a tether modelled. The section "3.2 Configuration
discrepancies" in the SFOC is still empty, and one line says "net positive
buoyancy 20 kg", which contradicts the other values.

**B. Maximum thrust per propeller.** The Aero Drum quotation specifies 4 x Kavan
C3536-1250kv motors, Kavan 50 A ESCs and 11-inch propellers. Kavan publishes no
static thrust for this motor: it is rated 420 W maximum (30 s, 4S), with an 11x5
propeller recommended on 3S. Motors of this class with an 11-inch propeller give
roughly 1.2-1.3 kg (12-13 N) at full throttle, so **12 N is used, as an estimate**.
That is about 5 times the 2.5 N in Final.py (a 7-inch propeller). A thrust-stand
measurement would settle it. On 4S with an 11-inch propeller the motor would be
past its rating, so the battery voltage matters too (the SFOC lists 4S Li-ion in
one place and 6S LiPo in another).

Effect: in Phase 1 the maximum thrust only scales the capacity numbers (which
layouts pass does not change). From Phase 2 on it matters a lot: with 12 N per
propeller there is far more authority against wind than with 2.5 N.

**C. STARLab 4-motor layout.** The drawings show four motors on the equator at 90°
spacing, each pushing horizontally along the hull. Assumed: all four push
counter-clockwise seen from above, motors placed on the hull surface (arms ignored).
- With push-only propellers this layout can only yaw one way, and cannot produce a
  pure push in any direction (every push also yaws it). Checked with a small LP:
  maximum +Mz 75.6 N·m (with 12 N motors), but 0 for −Mz, +Fx and −Fx.
- With reversible ESCs it can control Fx, Fy and yaw, but never Fz, roll or pitch.

The project page describes the flight software differently: "four thrusters: two
side thrusters for sideways movement and turning, a forward pair, and a lift pair".
That adds up to six, not four, and a "lift pair" would mean some motors push
vertically, unlike the SFOC drawing. The firmware is built on ArduSub, which
normally assumes reversible thrusters, while the quoted Kavan 50 A ESCs are
ordinary one-direction airplane ESCs unless reconfigured.

Please confirm, for each of the 4 motors: where it sits, which way it pushes, and
whether its ESC can reverse.

**D. Mass budget.** Still needed: gondola mass and how far the centre of gravity
sits below the hull centre (assumed 0.30 m); both set the pendulum stability in
Phase 2. The project page links a mass budget, but that page is blocked from this
cloud session; pasting its numbers here is enough. Known so far: ballonet
0.28 m³, 3 x 1 kg lead balance bars plus 6 balance pockets (Aero Drum quotation).

**E. Wind range for training (proposal).** The brief suggests 0–4 m/s. The SFOC
operating limits are 5.6 m/s sustained and 8.3 m/s gusts, so the proposal is to
train on mean wind 0–5.6 m/s with gusts up to 8.3 m/s.

**F. Path radius vs operating area.** The 50 m path keeps the comparison with the
paper, but the SFOC limits outdoor tests to a 20 m radius, so it can only be done
in simulation. Station keeping, the primary task, fits inside the limits.

## Assumptions in use

All values are tagged in `configs/sphere.yaml`. Current assumptions: air density
1.225 kg/m³, neutral-buoyancy trim, CG 0.30 m below the hull centre,
maximum thrust (12 N, estimated), 0.10 s motor time constant, and the STARLab
motor push directions.

Sources: SFOC application (Oct 2026 draft), Aero Drum official quotation (Aug 2025),
ACCESS project page text, Kavan C3536-1250 product listings.
