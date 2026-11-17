# PID tuning `20261117_153423`

Same bounds, cost, tuning scenario (1.5 m/s wind + light turbulence, 1 m offset, 45° heading step) and PSO budget for every layout.

| Layout | N | Tuning cost | 2 m/s steady: RMS pos [m] | effort [N²s] | 3 m/s turbulent: RMS pos [m] | max pos [m] | effort [N²s] | saturation % | survived |
|---|---|---|---|---|---|---|---|---|---|
| sphere_knee | 14 | 1.592 | 0.073 | 667 | 0.231 | 0.52 | 4553 | 0.0 | True |
| sphere_paper_pick | 23 | 1.410 | 0.064 | 448 | 0.201 | 0.44 | 2974 | 0.0 | True |
| sphere_s1_t16 | 6 | 2.143 | 0.072 | 1342 | 0.289 | 1.18 | 10762 | 4.6 | True |
| sphere_random_N14 (median [IQR] of 20) | 14 | 1.909 [1.762-2.270] | 0.089 [0.079-0.095] | 1187 [930-1903] | 0.286 [0.250-0.308] | 0.67 [0.58-0.74] | 7274 [6535-8884] | 0.2 [0.0-0.3] | 20/20 |
| sphere_random_N23 (median [IQR] of 20) | 23 | 1.489 [1.468-1.514] | 0.065 [0.063-0.068] | 528 [479-602] | 0.201 [0.197-0.211] | 0.44 [0.43-0.47] | 3771 [3278-4110] | 0.0 [0.0-0.0] | 20/20 |
| sphere_paper_best | 16 | 1.554 | 0.066 | 596 | 0.211 | 0.47 | 4369 | 0.0 | True |
| sphere_starlab_6thruster | 6 | 2.224 | 0.084 | 1762 | 0.283 | 1.01 | 10783 | 1.9 | True |
| maxN11_sphere_knee | 8 | 1.907 | 0.081 | 1137 | 0.258 | 0.60 | 8088 | 0.0 | True |
| maxN11_sphere_random_N8 (median [IQR] of 20) | 8 | 3.253 [2.788-4.306] | 0.076 [0.067-0.085] | 4823 [2595-6534] | 0.271 [0.243-0.339] | 0.79 [0.65-1.25] | 21815 [16059-26359] | 15.0 [7.5-23.4] | 20/20 |
