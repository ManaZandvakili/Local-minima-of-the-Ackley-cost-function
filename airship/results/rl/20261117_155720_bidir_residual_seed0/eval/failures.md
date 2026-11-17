# Motor-failure comparison

Each propeller in turn fails at 15 s (2 m/s wind, light turbulence, 45 s). PID gains `20261117_153423`. RMS position error after the failure: median over the failed propeller, survivors only.

| Layout | N | PID survived | 20261117_155720_bidir_residual_seed0 survived | 20261117_155720_bidir_plain_seed0 survived | PID RMS after [m] | 20261117_155720_bidir_residual_seed0 RMS after [m] | 20261117_155720_bidir_plain_seed0 RMS after [m] |
|---|---|---|---|---|---|---|---|
| sphere_knee | 14 | 14/14 | 14/14 | 14/14 | 0.070 | 0.047 | 0.037 |
| maxN11_sphere_knee | 8 | 8/8 | 8/8 | 8/8 | 0.071 | 0.052 | 0.031 |
| sphere_starlab_6thruster | 6 | 6/6 | 6/6 | 6/6 | 0.079 | 0.091 | 0.289 |
