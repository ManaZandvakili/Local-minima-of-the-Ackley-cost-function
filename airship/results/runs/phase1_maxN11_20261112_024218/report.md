# Phase 1 optimizer run `phase1_maxN11_20261112_024218`

Sphere, K = 30, 5 seeds x 20 trials, push-only thrust 0-12 N per propeller, at most 11 propellers.

## Per seed

| Seed | Feasible trials | Pareto layouts | Knee N | Paper-pick N | Minutes |
|---|---|---|---|---|---|
| 0 | 20/20 | 6 | 9 | 11 | 10.3 |
| 1 | 17/20 | 6 | 9 | 11 | 9.4 |
| 2 | 20/20 | 4 | 10 | 11 | 10.1 |
| 3 | 20/20 | 7 | 9 | 11 | 9.7 |
| 4 | 19/20 | 6 | 9 | 11 | 9.2 |

## Pareto layouts (all seeds)

| Seed | Trial | N | Count | Effort | Pos. sym. | Dir. sym. | Geometric | Authority | Capacity min [N] | Tags |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 9 | 0.600 | 1.591 | 0.103 | 0.269 | 0.134 | 9.514 | 6.11 |  |
| 0 | 1 | 11 | 0.733 | 1.463 | 0.099 | 0.167 | 0.685 | 10.732 | 7.54 | paper_pick |
| 0 | 8 | 10 | 0.667 | 2.649 | 0.103 | 0.001 | 0.000 | 9.681 | 7.36 |  |
| 0 | 10 | 8 | 0.533 | 2.553 | 0.018 | 0.031 | 0.000 | 8.851 | 6.23 |  |
| 0 | 12 | 7 | 0.467 | 3.060 | 0.239 | 0.001 | 0.000 | 8.048 | 6.06 |  |
| 0 | 14 | 9 | 0.600 | 1.823 | 0.152 | 0.135 | 0.000 | 9.680 | 10.52 | knee |
| 1 | 1 | 10 | 0.667 | 1.827 | 0.150 | 0.007 | 0.000 | 10.105 | 6.83 |  |
| 1 | 10 | 10 | 0.667 | 1.656 | 0.047 | 0.061 | 0.000 | 9.905 | 3.06 |  |
| 1 | 11 | 9 | 0.600 | 1.803 | 0.132 | 0.022 | 0.236 | 9.625 | 6.76 | knee |
| 1 | 12 | 11 | 0.733 | 1.550 | 0.216 | 0.116 | 0.170 | 10.297 | 4.94 | paper_pick |
| 1 | 14 | 11 | 0.733 | 1.618 | 0.153 | 0.136 | 0.189 | 10.679 | 8.74 |  |
| 1 | 16 | 8 | 0.533 | 2.273 | 0.074 | 0.315 | 0.345 | 8.973 | 1.90 |  |
| 2 | 2 | 10 | 0.667 | 1.680 | 0.071 | 0.069 | 0.000 | 10.348 | 3.54 | knee |
| 2 | 6 | 9 | 0.600 | 3.820 | 0.209 | 0.000 | 0.000 | 9.001 | 3.27 |  |
| 2 | 11 | 7 | 0.467 | 2.605 | 0.142 | 0.008 | 0.004 | 8.228 | 6.81 |  |
| 2 | 19 | 11 | 0.733 | 1.438 | 0.026 | 0.208 | 0.231 | 10.931 | 2.03 | paper_pick |
| 3 | 0 | 11 | 0.733 | 1.590 | 0.096 | 0.142 | 0.154 | 10.828 | 9.09 |  |
| 3 | 1 | 9 | 0.600 | 1.892 | 0.185 | 0.119 | 0.000 | 9.395 | 10.48 |  |
| 3 | 6 | 7 | 0.467 | 3.058 | 0.089 | 0.130 | 0.000 | 8.101 | 3.19 |  |
| 3 | 13 | 11 | 0.733 | 1.552 | 0.117 | 0.076 | 0.862 | 10.602 | 10.47 |  |
| 3 | 15 | 11 | 0.733 | 1.528 | 0.140 | 0.029 | 0.000 | 10.273 | 0.59 | paper_pick |
| 3 | 17 | 10 | 0.667 | 1.757 | 0.147 | 0.279 | 0.000 | 10.240 | 2.98 |  |
| 3 | 19 | 9 | 0.600 | 2.003 | 0.039 | 0.151 | 0.000 | 9.657 | 6.90 | knee |
| 4 | 0 | 9 | 0.600 | 1.783 | 0.276 | 0.187 | 0.000 | 9.403 | 6.28 | knee |
| 4 | 6 | 7 | 0.467 | 3.042 | 0.279 | 0.104 | 0.015 | 7.929 | 1.69 |  |
| 4 | 7 | 8 | 0.533 | 2.203 | 0.171 | 0.012 | 0.000 | 9.258 | 1.00 |  |
| 4 | 8 | 11 | 0.733 | 1.509 | 0.224 | 0.105 | 0.251 | 10.711 | 9.20 | paper_pick |
| 4 | 11 | 10 | 0.667 | 1.875 | 0.105 | 0.004 | 0.000 | 9.833 | 3.31 |  |
| 4 | 13 | 10 | 0.667 | 1.786 | 0.127 | 0.001 | 0.412 | 9.745 | 3.14 |  |

N over all Pareto layouts: min 7, median 10, max 11.

## Chosen layouts (over all seeds)

- **sphere_knee**: N = 9 (= sphere_s0_t00), capacity min 6.11 N, objectives J1_count 0.600, J2_effort 1.591, J3_position_symmetry 0.103, J4_direction_symmetry 0.269, J5_geometric 0.134
- **sphere_paper_pick**: N = 11 (= sphere_s2_t19), capacity min 2.03 N, objectives J1_count 0.733, J2_effort 1.438, J3_position_symmetry 0.026, J4_direction_symmetry 0.208, J5_geometric 0.231

Matched random layouts for sphere_knee: 20 with N = 9 (633 attempts), median capacity min 1.49 N.

Matched random layouts for sphere_paper_pick: 20 with N = 11 (166 attempts), median capacity min 3.25 N.

## Reference layouts under the new test

| Layout | N | Thrust range [N] | Feasible (all 6 DOF) | Capacity min [N] | Axes with zero capacity |
|---|---|---|---|---|---|
| cigar_paper_best | 13 | 0 to 12 | yes | 15.41 | none |
| oblate_paper_best | 16 | 0 to 12 | yes | 24.60 | none |
| sphere_paper_best | 16 | 0 to 12 | yes | 10.01 | none |
| sphere_starlab_6thruster | 6 | -12 to 12 | no | 0.00 | +Mx, -Mx |
