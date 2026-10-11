# Final combination grid (spec128-129)

Settings: E1 241600 (CROP-seq vs bulk), E2 67200 (bulk across lines), E3 1680 x 160 Cm-KO pairs. Corrected sample set: clean H2052 controls (H28_NTC_R1-R3), H2052_NTC_R3 dropped, RUNX2 and SNAI2 added, RUNX1 excluded.

## E1 CROP-seq vs bulk

| target | settings | % p<0.05 | best z | grid-wide p (target) | grid-wide p (all targets) |
|---|---|---|---|---|---|
| H2052 | 26560 | 3.8 | 2.24 | 0.82 |   1 |
| MSTO | 26880 | 0.1 | 2.36 | 0.87 | 0.99 |
| H2452 | 26880 | 5.8 | 2.96 | 0.16 | 0.74 |
| H28 | 26880 | 4.7 | 2.87 | 0.3 | 0.82 |
| pooled 4 lines | 26880 | 2.7 | 2.88 | 0.23 | 0.81 |
| pooled w/o H2052 | 26880 | 1.5 | 2.58 | 0.44 | 0.96 |
| pooled w/o MSTO | 26880 | 4.6 | 2.54 | 0.47 | 0.97 |
| pooled w/o H2452 | 26880 | 2.1 | 2.41 | 0.67 | 0.99 |
| pooled w/o H28 | 26880 | 2.8 | 2.51 | 0.56 | 0.98 |

TF x target with the own KO ranked 1-2 in at least half of the settings: TEAD4 in H2052 (76%); PITX1 in H2052 (67%); TEAD4 in MSTO (73%); TEAD4 in H2452 (99%); PITX1 in H2452 (68%); TEAD4 in H28 (56%); PITX1 in H28 (84%); BPTF in H28 (50%); TEAD4 in pooled 4 lines (96%); PITX1 in pooled 4 lines (73%); TEAD4 in pooled w/o H2052 (98%); PITX1 in pooled w/o H2052 (65%); TEAD4 in pooled w/o MSTO (91%); PITX1 in pooled w/o MSTO (82%); TEAD4 in pooled w/o H2452 (82%); PITX1 in pooled w/o H2452 (69%); TEAD4 in pooled w/o H28 (96%); PITX1 in pooled w/o H28 (58%)

## E2 bulk across lines

| line set | settings | % p<0.05 | % p<0.001 | median z | max z |
|---|---|---|---|---|---|
| all | 13440 | 31.4 | 8.4 | 0.99 | 5.62 |
| no_H2052 | 13440 | 18.1 | 1.9 | 0.70 | 4.81 |
| no_MSTO | 13440 | 25.6 | 6.9 | 0.85 | 5.30 |
| no_H2452 | 13440 | 27.8 | 1.7 | 1.08 | 4.94 |
| no_H28 | 13440 | 16.8 | 1.0 | 0.55 | 3.92 |

KOs with per-KO p < 0.05 in at least half of the settings (all 4 lines): 

## E3 Cm lower in KO than NTC

| pair | in text | screen NTC-KO | screen p | % settings supported lower | H2052 | MSTO | H2452 | H28 |
|---|---|---|---|---|---|---|---|---|
| Cm20 TEAD4 | yes | +0.82 | 2e-24 | 85 | +2.03 | -0.91 | +2.00 | +1.55 |
| Cm8 TWIST1 | yes | +0.17 | 1e-04 | 83 | +2.20 | +1.21 | -0.64 | +2.47 |
| Cm1 SOX9 | yes | +0.39 | 0.0044 | 42 | +1.72 | -1.26 | +0.83 | +1.98 |
| Cm10 TCF3 | yes | +0.36 | 2e-06 | 1 | +0.75 | -0.10 | +0.83 | -1.10 |
| Cm9 TWIST1 | yes | +0.31 | 2e-13 | 0 | -0.03 | -1.14 | -2.16 | +2.10 |
| Cm9 PITX1 | yes | +0.08 | 0.016 | 0 | +1.11 | -0.85 | -2.14 | +1.74 |
| Cm16 PITX1 | yes | +0.74 | 4e-90 | 0 | +0.93 | -0.54 | -1.12 | +0.85 |
| Cm8 TEAD4 |  | +0.30 | 3e-04 | 99 | +1.94 | +0.86 | +1.24 | +1.67 |
| Cm13 TEAD4 |  | +0.26 | 0.0092 | 78 | +1.83 | +1.03 | +0.12 | +0.62 |
| Cm13 TWIST1 |  | +0.12 | 0.0089 | 59 | +2.14 | +0.44 | -1.26 | +1.42 |
| Cm16 SOX9 |  | +0.58 | 3e-05 | 38 | +0.02 | +0.17 | +0.88 | +1.07 |
| Cm16 TEAD2 |  | +0.88 | 4e-05 | 30 | +0.36 | +0.82 | -0.22 | +0.52 |
| Cm18 TEAD4 |  | +0.46 | 1e-08 | 29 | +1.36 | -2.29 | +1.96 | +0.83 |
| Cm1 TEAD2 |  | +0.51 | 0.011 | 27 | +1.18 | -0.94 | +1.15 | +0.68 |
| Cm6 TEAD4 |  | +0.38 | 1e-05 | 23 | -1.15 | +1.19 | +1.28 | +0.74 |
| Cm12 TEAD2 |  | +0.66 | 0.0013 | 22 | +0.40 | +1.26 | +0.33 | +0.54 |
