# Path following: one-way propellers

The paper's circle: radius 50 m, 0.5 m/s, 650 s (about one lap), nose along the path. Wind from 3 directions with turbulence (0.5 m/s) and gusts. Error = distance from the reference point; values = median over the 3 directions, then median over the 20 random layouts. "Completed" = error stayed under 3 m for the whole lap (the paper's rule), share of flights.

PID gains `20261117_203431`; learning runs `20261117_205532_oneway_residual_seed0`, `20261117_205532_oneway_plain_seed0`.


## Wind 1.5 m/s

| Layout | N | PID | Residual | Plain | Completed (PID / residual / plain) |
|---|---|---|---|---|---|
| maxN11_sphere_knee (optimized) | 9 | 0.127 m | 0.119 m | 0.192 m | 100% / 100% / 100% |
| random, 9 propellers (20) | 9 | 0.835 m [0.40-1.27] | 0.820 m [0.39-1.47] | 0.981 m [0.58-1.95] | 27% / 30% / 22% |
| sphere_knee (optimized) | 15 | 0.118 m | 0.101 m | 0.144 m | 100% / 100% / 100% |
| random, 15 propellers (20) | 15 | 0.175 m [0.12-0.25] | 0.152 m [0.12-0.25] | 0.282 m [0.18-0.52] | 88% / 88% / 77% |
| sphere_paper_pick (optimized) | 22 | 0.106 m | 0.152 m | 0.154 m | 100% / 100% / 100% |
| random, 22 propellers (20) | 22 | 0.123 m [0.11-0.13] | 0.113 m [0.11-0.13] | 0.247 m [0.17-0.42] | 100% / 100% / 100% |
| sphere_paper_best | 16 | 0.119 m | 0.096 m | 0.136 m | 100% / 100% / 100% |

## Wind 3.0 m/s

| Layout | N | PID | Residual | Plain | Completed (PID / residual / plain) |
|---|---|---|---|---|---|
| maxN11_sphere_knee (optimized) | 9 | 1.446 m | 1.614 m | 2.454 m | 0% / 0% / 0% |
| random, 9 propellers (20) | 9 | 4.780 m [3.81-5.49] | 4.826 m [4.19-5.44] | 5.283 m [5.10-5.51] | 0% / 0% / 0% |
| sphere_knee (optimized) | 15 | 0.346 m | 0.293 m | 0.538 m | 100% / 100% / 33% |
| random, 15 propellers (20) | 15 | 1.747 m [0.96-2.32] | 1.712 m [0.91-2.82] | 2.086 m [1.66-2.74] | 5% / 7% / 3% |
| sphere_paper_pick (optimized) | 22 | 0.233 m | 0.222 m | 0.430 m | 100% / 100% / 100% |
| random, 22 propellers (20) | 22 | 0.387 m [0.30-0.55] | 0.303 m [0.24-0.48] | 0.730 m [0.59-0.89] | 65% / 63% / 45% |
| sphere_paper_best | 16 | 0.288 m | 0.223 m | 0.497 m | 100% / 100% / 67% |
