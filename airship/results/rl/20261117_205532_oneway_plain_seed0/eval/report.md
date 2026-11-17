# Phase 4 evaluation: RL vs PID

RL run `20261117_205532_oneway_plain_seed0`, PID gains `20261117_203431`. Gusty station keeping: 3 m/s mean wind from 3 directions, turbulence 0.5 m/s, gusts, 45° heading change; identical wind seeds for both controllers. Values: median over the wind directions (per layout), then median and IQR over the layouts in a group.

| Layout group | N | PID RMS pos [m] | RL RMS pos [m] | PID effort [N²s] | RL effort [N²s] | PID survival | RL survival |
|---|---|---|---|---|---|---|---|
| maxN11_sphere_knee | 9 | 1.016 | 2.970 | 20837 | 18492 | 1.00 | 0.67 |
| maxN11_sphere_random_N9 (20) | 9 | 4.688 [3.712-5.158] | 5.271 [4.763-5.407] | 15437 [7126-25321] | 8676 [6179-22451] | 0.32 | 0.28 |
| sphere_knee | 15 | 0.362 | 0.460 | 13180 | 19990 | 1.00 | 1.00 |
| sphere_paper_best | 16 | 0.324 | 0.475 | 12304 | 20087 | 1.00 | 1.00 |
| sphere_paper_pick | 22 | 0.287 | 0.508 | 7245 | 8424 | 1.00 | 1.00 |
| sphere_random_N15 (20) | 15 | 0.507 [0.372-1.171] | 0.754 [0.599-3.008] | 22271 [16974-24050] | 24013 [20008-28712] | 0.85 | 0.88 |
| sphere_random_N22 (20) | 22 | 0.356 [0.300-0.408] | 0.627 [0.576-0.671] | 10528 [9608-15633] | 15841 [11139-24373] | 1.00 | 1.00 |
| sphere_s0_t12 | 7 | 2.629 | 3.736 | 25055 | 24590 | 1.00 | 0.67 |

**Brief acceptance (held-in layouts, turbulent scenario):** median RMS position error PID 0.499 m, RL 0.762 m -> RL does not beat PID.

Held-out layouts (14): median RMS PID 0.488 m, RL 1.025 m.

Training curves: `training_curves.png`.

Video: `pid_vs_rl_knee.mp4`
Video: `optimized_vs_random_N9_rl.mp4`
Video: `motor_failure_knee.mp4`
