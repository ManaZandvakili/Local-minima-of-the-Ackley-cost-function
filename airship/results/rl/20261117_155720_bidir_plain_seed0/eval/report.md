# Phase 4 evaluation: RL vs PID

RL run `20261117_155720_bidir_plain_seed0`, PID gains `20261117_153423`. Gusty station keeping: 3 m/s mean wind from 3 directions, turbulence 0.5 m/s, gusts, 45° heading change; identical wind seeds for both controllers. Values: median over the wind directions (per layout), then median and IQR over the layouts in a group.

| Layout group | N | PID RMS pos [m] | RL RMS pos [m] | PID effort [N²s] | RL effort [N²s] | PID survival | RL survival |
|---|---|---|---|---|---|---|---|
| maxN11_sphere_knee | 8 | 0.330 | 0.238 | 7516 | 8537 | 1.00 | 1.00 |
| maxN11_sphere_random_N8 (20) | 8 | 0.335 [0.303-0.412] | 0.361 [0.298-0.604] | 25021 [19362-28213] | 25422 [22727-27505] | 1.00 | 0.95 |
| sphere_knee | 14 | 0.292 | 0.264 | 4363 | 4793 | 1.00 | 1.00 |
| sphere_paper_best | 16 | 0.263 | 0.260 | 4223 | 4480 | 1.00 | 1.00 |
| sphere_paper_pick | 23 | 0.254 | 0.333 | 2907 | 3595 | 1.00 | 1.00 |
| sphere_random_N14 (20) | 14 | 0.364 [0.318-0.396] | 0.232 [0.223-0.238] | 8508 [6786-10826] | 9177 [7323-11513] | 1.00 | 1.00 |
| sphere_random_N23 (20) | 23 | 0.256 [0.250-0.269] | 0.315 [0.306-0.322] | 3786 [3488-4190] | 4232 [3851-4597] | 1.00 | 1.00 |
| sphere_s1_t16 | 6 | 0.321 | 0.370 | 10547 | 13061 | 1.00 | 1.00 |
| sphere_starlab_6thruster | 6 | 0.368 | 0.358 | 10870 | 14566 | 1.00 | 1.00 |

**Brief acceptance (held-in layouts, turbulent scenario):** median RMS position error PID 0.318 m, RL 0.295 m -> RL beats PID.

Held-out layouts (14): median RMS PID 0.283 m, RL 0.319 m.

Training curves: `training_curves.png`.

Video: `pid_vs_rl_knee.mp4`
Video: `optimized_vs_random_N8_rl.mp4`
Video: `motor_failure_knee.mp4`
