# Motor-failure comparison

Each propeller in turn fails at 15 s (2 m/s wind, light turbulence, 45 s). PID gains `20261117_203431`. RMS position error after the failure: median over the failed propeller, survivors only.

| Layout | N | PID survived | 20261117_205532_oneway_residual_seed0 survived | 20261117_205532_oneway_plain_seed0 survived | PID RMS after [m] | 20261117_205532_oneway_residual_seed0 RMS after [m] | 20261117_205532_oneway_plain_seed0 RMS after [m] |
|---|---|---|---|---|---|---|---|
| sphere_knee | 15 | 15/15 | 15/15 | 15/15 | 0.089 | 0.059 | 0.125 |
| maxN11_sphere_knee | 9 | 8/9 | 7/9 | 8/9 | 0.088 | 0.175 | 0.428 |
