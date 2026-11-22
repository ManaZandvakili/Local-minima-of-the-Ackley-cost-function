# Path following: reversible propellers

The paper's circle: radius 50 m, 0.5 m/s, 650 s (about one lap), nose along the path. Wind from 3 directions with turbulence (0.5 m/s) and gusts. Error = distance from the reference point; values = median over the 3 directions, then median over the 20 random layouts. "Completed" = error stayed under 3 m for the whole lap (the paper's rule), share of flights.

PID gains `20261117_153423`; learning runs `20261117_155720_bidir_residual_seed0`, `20261117_155720_bidir_plain_seed0`.


## Wind 1.5 m/s

| Layout | N | PID | Residual | Plain | Completed (PID / residual / plain) |
|---|---|---|---|---|---|
| maxN11_sphere_knee (optimized) | 8 | 0.116 m | 0.076 m | 0.073 m | 100% / 100% / 100% |
| random, 8 propellers (20) | 8 | 0.115 m [0.10-0.13] | 0.088 m [0.08-0.11] | 0.128 m [0.12-0.17] | 95% / 95% / 95% |
| sphere_knee (optimized) | 14 | 0.100 m | 0.078 m | 0.069 m | 100% / 100% / 100% |
| random, 14 propellers (20) | 14 | 0.121 m [0.11-0.13] | 0.082 m [0.08-0.09] | 0.077 m [0.07-0.09] | 100% / 100% / 100% |
| sphere_paper_pick (optimized) | 23 | 0.087 m | 0.073 m | 0.113 m | 100% / 100% / 100% |
| random, 23 propellers (20) | 23 | 0.088 m [0.09-0.09] | 0.070 m [0.07-0.07] | 0.113 m [0.11-0.12] | 100% / 100% / 100% |
| sphere_paper_best | 16 | 0.095 m | 0.069 m | 0.069 m | 100% / 100% / 100% |
| sphere_starlab_6thruster | 6 | 0.116 m | 0.093 m | 0.684 m | 100% / 100% / 100% |

## Wind 3.0 m/s

| Layout | N | PID | Residual | Plain | Completed (PID / residual / plain) |
|---|---|---|---|---|---|
| maxN11_sphere_knee (optimized) | 8 | 0.289 m | 0.185 m | 0.226 m | 100% / 100% / 100% |
| random, 8 propellers (20) | 8 | 0.511 m [0.36-0.78] | 0.499 m [0.31-0.99] | 0.957 m [0.67-1.94] | 37% / 42% / 12% |
| sphere_knee (optimized) | 14 | 0.236 m | 0.160 m | 0.200 m | 100% / 100% / 100% |
| random, 14 propellers (20) | 14 | 0.323 m [0.27-0.37] | 0.204 m [0.18-0.22] | 0.228 m [0.22-0.24] | 98% / 100% / 100% |
| sphere_paper_pick (optimized) | 23 | 0.191 m | 0.140 m | 0.239 m | 100% / 100% / 100% |
| random, 23 propellers (20) | 23 | 0.193 m [0.18-0.21] | 0.136 m [0.13-0.14] | 0.252 m [0.24-0.26] | 100% / 100% / 100% |
| sphere_paper_best | 16 | 0.200 m | 0.138 m | 0.213 m | 100% / 100% / 100% |
| sphere_starlab_6thruster | 6 | 2.625 m | 2.210 m | 2.048 m | 0% / 0% / 0% |
