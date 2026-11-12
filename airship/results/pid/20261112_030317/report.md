# PID tuning `20261112_030317`

Same bounds, cost, tuning scenario (1.5 m/s wind + light turbulence, 1 m offset, 45° heading step) and PSO budget for every layout.

| Layout | N | Tuning cost | 2 m/s steady: RMS pos [m] | effort [N²s] | 3 m/s turbulent: RMS pos [m] | max pos [m] | effort [N²s] | saturation % | survived |
|---|---|---|---|---|---|---|---|---|---|
| sphere_knee | 15 | 2.650 | 0.086 | 2299 | 0.301 | 0.71 | 24186 | 8.5 | True |
| sphere_paper_pick | 22 | 1.861 | 0.072 | 1082 | 0.230 | 0.52 | 6245 | 0.2 | True |
| sphere_s0_t12 | 7 | 4.724 | 0.050 | 9540 | 0.786 | 3.64 | 22914 | 24.1 | True |
| sphere_random_N15 (median [IQR] of 20) | 15 | 3.123 [2.633-5.520] | 0.084 [0.064-0.100] | 3115 [2147-5901] | 0.468 [0.274-1.222] | 1.59 [0.85-3.41] | 17967 [14722-22368] | 8.1 [4.2-14.5] | 17/20 |
| sphere_random_N22 (median [IQR] of 20) | 22 | 2.192 [1.872-2.455] | 0.085 [0.072-0.094] | 1497 [1104-1968] | 0.281 [0.235-0.324] | 0.72 [0.59-0.81] | 17020 [8774-20035] | 2.0 [0.1-3.4] | 20/20 |
| sphere_paper_best | 16 | 2.386 | 0.081 | 1729 | 0.261 | 0.68 | 19256 | 3.7 | True |
| sphere_starlab_6thruster | 6 | 2.224 | 0.084 | 1762 | 0.283 | 1.01 | 10783 | 1.9 | True |
| maxN11_sphere_knee | 9 | 5.275 | 0.059 | 10516 | 0.286 | 1.56 | 16832 | 6.6 | True |
| maxN11_sphere_random_N9 (median [IQR] of 20) | 9 | 9.455 [6.156-17.597] | 0.149 [0.079-0.312] | 14802 [8276-23098] | 3.674 [1.118-5.239] | 8.84 [4.16-10.03] | 20328 [11637-24279] | 23.8 [15.7-31.9] | 12/20 |
