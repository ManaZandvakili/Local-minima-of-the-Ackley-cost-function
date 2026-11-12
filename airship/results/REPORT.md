# Report: journal version (sphere)

Plain-language log of results, decisions, assumptions and open questions.

## Decisions from Mana (Section 11)

| # | Question | Answer |
|---|---|---|
| 1 | Hull size | 3.15 m diameter |
| 2 | Mass | from the SFOC application: final operating mass 16 kg, gross lift 20 kg (see open question A) |
| 3 | Objective set | the paper's five: count, effort, position symmetry, direction symmetry, geometric |
| 4 | Path-following radius | 50 m |
| 5 | Thrust | "from the SFOC application"; it has no thrust data, so 2.5 N from Final.py is kept (open question B) |
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

**B. Maximum thrust per propeller.** The SFOC application gives no motor, propeller
or thrust data. 2.5 N per propeller (Kavan 7-inch, from Final.py) is kept. Please
send the motor and propeller model, or thrust-stand data if you have it.

**C. STARLab 4-motor layout.** The drawings show four motors on the equator at 90°
spacing, each pushing horizontally along the hull. Assumed: all four push
counter-clockwise seen from above, motors placed on the hull surface (arms ignored).
- With push-only propellers this layout can only yaw one way, and cannot produce a
  pure push in any direction (every push also yaws it). Checked with a small LP:
  maximum +Mz 15.75 N·m, but 0 for −Mz, +Fx and −Fx.
- With reversible ESCs it can control Fx, Fy and yaw, but never Fz, roll or pitch.

Please confirm the push direction of each motor and whether the ESCs are reversible.

**D. Other values not in the SFOC application:** gondola mass, and how far the
centre of gravity sits below the hull centre (assumed 0.30 m). Both set the
pendulum stability in Phase 2.

**E. Wind range for training (proposal).** The brief suggests 0–4 m/s. The SFOC
operating limits are 5.6 m/s sustained and 8.3 m/s gusts, so the proposal is to
train on mean wind 0–5.6 m/s with gusts up to 8.3 m/s.

**F. Path radius vs operating area.** The 50 m path keeps the comparison with the
paper, but the SFOC limits outdoor tests to a 20 m radius, so it can only be done
in simulation. Station keeping, the primary task, fits inside the limits.

## Assumptions in use

All values are tagged in `configs/sphere.yaml`. Current assumptions: air density
1.225 kg/m³, neutral-buoyancy trim, CG 0.30 m below the hull centre, 2.5 N
maximum thrust, 0.10 s motor time constant, and the STARLab motor push directions.
