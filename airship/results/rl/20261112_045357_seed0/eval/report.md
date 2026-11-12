# Phase 4 evaluation: RL vs PID

RL run `20261112_045357_seed0`, PID gains `20261112_030317`. Gusty station keeping: 3 m/s mean wind from 3 directions, turbulence 0.5 m/s, gusts, 45° heading change; identical wind seeds for both controllers. Values: median over the wind directions (per layout), then median and IQR over the layouts in a group.

| Layout group | N | PID RMS pos [m] | RL RMS pos [m] | PID effort [N²s] | RL effort [N²s] | PID survival | RL survival |
|---|---|---|---|---|---|---|---|
| maxN11_sphere_knee | 9 | 1.016 | 3.506 | 20837 | 17751 | 1.00 | 0.67 |
| maxN11_sphere_random_N9 (20) | 9 | 4.688 [3.712-5.158] | 4.973 [3.780-5.282] | 15437 [7126-25321] | 17232 [6263-25787] | 0.32 | 0.30 |
| sphere_knee | 15 | 0.362 | 0.407 | 13180 | 25048 | 1.00 | 1.00 |
| sphere_paper_best | 16 | 0.324 | 0.409 | 12304 | 19691 | 1.00 | 1.00 |
| sphere_paper_pick | 22 | 0.287 | 0.498 | 7245 | 10423 | 1.00 | 1.00 |
| sphere_random_N15 (20) | 15 | 0.507 [0.372-1.171] | 0.841 [0.570-2.195] | 22271 [16974-24050] | 22912 [18262-31235] | 0.85 | 0.87 |
| sphere_random_N22 (20) | 22 | 0.356 [0.300-0.408] | 0.524 [0.459-0.621] | 10528 [9608-15633] | 16389 [14051-21405] | 1.00 | 1.00 |
| sphere_s0_t12 | 7 | 2.629 | 3.608 | 25055 | 26934 | 1.00 | 1.00 |
| sphere_starlab_6thruster | 6 | 0.368 | 0.617 | 10870 | 13409 | 1.00 | 1.00 |

**Brief acceptance (held-in layouts, turbulent scenario):** median RMS position error PID 0.489 m, RL 0.813 m -> RL does not beat PID.

Held-out layouts (14): median RMS PID 0.488 m, RL 1.495 m.

Training curves: `training_curves.png`.

Video: `pid_vs_rl_knee.mp4`
Video: `optimized_vs_random_N9_rl.mp4`
Video: `motor_failure_knee.mp4`
