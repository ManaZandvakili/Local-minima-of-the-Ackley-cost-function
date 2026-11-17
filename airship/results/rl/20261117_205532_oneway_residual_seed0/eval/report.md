# Phase 4 evaluation: RL vs PID

RL run `20261117_205532_oneway_residual_seed0`, PID gains `20261117_203431`. Gusty station keeping: 3 m/s mean wind from 3 directions, turbulence 0.5 m/s, gusts, 45° heading change; identical wind seeds for both controllers. Values: median over the wind directions (per layout), then median and IQR over the layouts in a group.

| Layout group | N | PID RMS pos [m] | RL RMS pos [m] | PID effort [N²s] | RL effort [N²s] | PID survival | RL survival |
|---|---|---|---|---|---|---|---|
| maxN11_sphere_knee | 9 | 1.016 | 0.909 | 20837 | 20700 | 1.00 | 1.00 |
| maxN11_sphere_random_N9 (20) | 9 | 4.688 [3.712-5.158] | 5.045 [3.965-5.207] | 15437 [7126-25321] | 16465 [6930-24968] | 0.32 | 0.33 |
| sphere_knee | 15 | 0.362 | 0.307 | 13180 | 11385 | 1.00 | 1.00 |
| sphere_paper_best | 16 | 0.324 | 0.278 | 12304 | 13333 | 1.00 | 1.00 |
| sphere_paper_pick | 22 | 0.287 | 0.329 | 7245 | 9142 | 1.00 | 1.00 |
| sphere_random_N15 (20) | 15 | 0.507 [0.372-1.171] | 0.455 [0.332-1.005] | 22271 [16974-24050] | 21713 [18909-25326] | 0.85 | 0.85 |
| sphere_random_N22 (20) | 22 | 0.356 [0.300-0.408] | 0.335 [0.306-0.341] | 10528 [9608-15633] | 11783 [9447-15960] | 1.00 | 1.00 |
| sphere_s0_t12 | 7 | 2.629 | 2.616 | 25055 | 24254 | 1.00 | 1.00 |

**Brief acceptance (held-in layouts, turbulent scenario):** median RMS position error PID 0.499 m, RL 0.442 m -> RL beats PID.

Held-out layouts (14): median RMS PID 0.488 m, RL 0.591 m.

Training curves: `training_curves.png`.

Video: `pid_vs_rl_knee.mp4`
Video: `optimized_vs_random_N9_rl.mp4`
Video: `motor_failure_knee.mp4`
