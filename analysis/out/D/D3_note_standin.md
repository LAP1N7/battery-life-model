# D3 note (input: STAND-IN prototype capacity_fit_AB.csv (EFC based, NOT team A's file): `..\..\prototype\results\mohtat\capacity_fit_AB.csv`, columns A, B)

Computed, not interpreted. Criterion proposed in D0_spec: rho >= 0.5 and p < 0.05 overall, and pooled-within-temperature rho >= 0.3 with the same sign. n is small (<= 21; 7 per temperature group).

Cells in A's file: 21/21 (missing: ). Shares at t_ref = last check-up hours with SOH >= 0.80; mode shares/expansion at 10 % measured loss ('matched_how' column; cells that never reach it use their last check-up).
A,B rounding: |A_file - A_recomputed| max 1.0e-02, |B_file - B_recomputed| max 1.4e-03; lin_share file vs recomputed Spearman 0.88 (cells with a negative A or B in the unrounded fit: 5).
Reference caveat: modes of H cells are relative to the first averaged check-up (47-294 Ah), not BOL (7 cells flagged); expansion artifact cells (flagged in D2): 03, 09, 13, 16, 21.

## A,B as in file

| test | subset | n | rho | p |
|---|---|---|---|---|
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | all | 21 | 0.86 | 0.000 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within C | 7 | -0.36 | 0.432 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within H | 7 | 0.96 | 0.000 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within R | 7 | 0.82 | 0.025 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | pooled within T (rank-demeaned) | 21 | 0.69 | 0.001 |
| H-D2a lin share ~ LAM share | all | 21 | 0.86 | 0.000 |
| H-D2a lin share ~ LAM share | within C | 7 | -0.36 | 0.432 |
| H-D2a lin share ~ LAM share | within H | 7 | 0.96 | 0.000 |
| H-D2a lin share ~ LAM share | within R | 7 | 0.82 | 0.025 |
| H-D2a lin share ~ LAM share | pooled within T (rank-demeaned) | 21 | 0.69 | 0.001 |
| secondary: lin share ~ LAM_NE share | all | 21 | 0.75 | 0.000 |
| secondary: lin share ~ LAM_NE share | within C | 7 | -0.46 | 0.294 |
| secondary: lin share ~ LAM_NE share | within H | 7 | 0.93 | 0.003 |
| secondary: lin share ~ LAM_NE share | within R | 7 | 0.82 | 0.025 |
| secondary: lin share ~ LAM_NE share | pooled within T (rank-demeaned) | 21 | 0.51 | 0.027 |
| secondary: lin share ~ LAM_PE share | all | 21 | -0.33 | 0.140 |
| secondary: lin share ~ LAM_PE share | within C | 7 | -0.39 | 0.383 |
| secondary: lin share ~ LAM_PE share | within H | 7 | -0.70 | 0.077 |
| secondary: lin share ~ LAM_PE share | within R | 7 | -0.33 | 0.465 |
| secondary: lin share ~ LAM_PE share | pooled within T (rank-demeaned) | 21 | -0.42 | 0.070 |
| H-D2b irreversible expansion (clean) ~ lin share | all | 21 | -0.47 | 0.033 |
| H-D2b irreversible expansion (clean) ~ lin share | within C | 7 | 0.61 | 0.148 |
| H-D2b irreversible expansion (clean) ~ lin share | within H | 7 | -0.74 | 0.057 |
| H-D2b irreversible expansion (clean) ~ lin share | within R | 7 | -0.19 | 0.691 |
| H-D2b irreversible expansion (clean) ~ lin share | pooled within T (rank-demeaned) | 21 | -0.16 | 0.503 |
| secondary: expansion dX_irr_r1 ~ lin share | all | 20 | -0.21 | 0.366 |
| secondary: expansion dX_irr_r1 ~ lin share | within C | 7 | 0.64 | 0.119 |
| secondary: expansion dX_irr_r1 ~ lin share | within H | 6 | -0.88 | 0.021 |
| secondary: expansion dX_irr_r1 ~ lin share | within R | 7 | -0.37 | 0.413 |
| secondary: expansion dX_irr_r1 ~ lin share | pooled within T (rank-demeaned) | 20 | -0.31 | 0.213 |
