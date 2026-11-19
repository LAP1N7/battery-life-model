# D3 note (input: P1_pred (team A): `out\P1_pred.csv`, columns A_full, B_full)

Computed, not interpreted. Criterion proposed in D0_spec: rho >= 0.5 and p < 0.05 overall, and pooled-within-temperature rho >= 0.3 with the same sign. n is small (<= 21; 7 per temperature group).

Cells in A's file: 17/21 (missing: 06, 09, 11, 17). Shares at t_ref = last check-up hours with SOH >= 0.80; mode shares/expansion at 10 % measured loss ('matched_how' column; cells that never reach it use their last check-up).
A,B rounding: |A_file - A_recomputed| max 4.9e-05, |B_file - B_recomputed| max 4.9e-05; lin_share file vs recomputed Spearman 0.93 (cells with a negative A or B in the unrounded fit: 5).
Reference caveat: modes of H cells are relative to the first averaged check-up (47-294 Ah), not BOL (7 cells flagged); expansion artifact cells (flagged in D2): 03, 09, 13, 16, 21.

## A,B as in file

| test | subset | n | rho | p |
|---|---|---|---|---|
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | all | 17 | 0.69 | 0.002 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within C | 5 | -0.56 | 0.322 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within H | 5 | 0.80 | 0.104 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within R | 7 | 0.86 | 0.012 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | pooled within T (rank-demeaned) | 17 | 0.69 | 0.004 |
| H-D2a lin share ~ LAM share | all | 17 | 0.69 | 0.002 |
| H-D2a lin share ~ LAM share | within C | 5 | -0.56 | 0.322 |
| H-D2a lin share ~ LAM share | within H | 5 | 0.80 | 0.104 |
| H-D2a lin share ~ LAM share | within R | 7 | 0.86 | 0.012 |
| H-D2a lin share ~ LAM share | pooled within T (rank-demeaned) | 17 | 0.69 | 0.004 |
| secondary: lin share ~ LAM_NE share | all | 17 | 0.81 | 0.000 |
| secondary: lin share ~ LAM_NE share | within C | 5 | -0.56 | 0.322 |
| secondary: lin share ~ LAM_NE share | within H | 5 | 0.60 | 0.285 |
| secondary: lin share ~ LAM_NE share | within R | 7 | 0.83 | 0.021 |
| secondary: lin share ~ LAM_NE share | pooled within T (rank-demeaned) | 17 | 0.42 | 0.117 |
| secondary: lin share ~ LAM_PE share | all | 17 | -0.59 | 0.013 |
| secondary: lin share ~ LAM_PE share | within C | 5 | -0.56 | 0.322 |
| secondary: lin share ~ LAM_PE share | within H | 5 | -0.20 | 0.747 |
| secondary: lin share ~ LAM_PE share | within R | 7 | -0.27 | 0.558 |
| secondary: lin share ~ LAM_PE share | pooled within T (rank-demeaned) | 17 | -0.44 | 0.103 |
| H-D2b irreversible expansion (clean) ~ lin share | all | 17 | -0.23 | 0.373 |
| H-D2b irreversible expansion (clean) ~ lin share | within C | 5 | 0.21 | 0.741 |
| H-D2b irreversible expansion (clean) ~ lin share | within H | 5 | -0.60 | 0.285 |
| H-D2b irreversible expansion (clean) ~ lin share | within R | 7 | -0.31 | 0.504 |
| H-D2b irreversible expansion (clean) ~ lin share | pooled within T (rank-demeaned) | 17 | -0.18 | 0.523 |
| secondary: expansion dX_irr_r1 ~ lin share | all | 16 | -0.05 | 0.866 |
| secondary: expansion dX_irr_r1 ~ lin share | within C | 5 | 0.15 | 0.805 |
| secondary: expansion dX_irr_r1 ~ lin share | within H | 4 | -0.80 | 0.200 |
| secondary: expansion dX_irr_r1 ~ lin share | within R | 7 | -0.58 | 0.175 |
| secondary: expansion dX_irr_r1 ~ lin share | pooled within T (rank-demeaned) | 16 | -0.42 | 0.135 |

## A,B recomputed unrounded (same cells as file)

| test | subset | n | rho | p |
|---|---|---|---|---|
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | all | 17 | 0.78 | 0.000 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within C | 5 | 0.60 | 0.285 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within H | 5 | 0.80 | 0.104 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within R | 7 | 0.86 | 0.014 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | pooled within T (rank-demeaned) | 17 | 0.83 | 0.000 |
| H-D2a lin share ~ LAM share | all | 17 | 0.78 | 0.000 |
| H-D2a lin share ~ LAM share | within C | 5 | 0.60 | 0.285 |
| H-D2a lin share ~ LAM share | within H | 5 | 0.80 | 0.104 |
| H-D2a lin share ~ LAM share | within R | 7 | 0.86 | 0.014 |
| H-D2a lin share ~ LAM share | pooled within T (rank-demeaned) | 17 | 0.83 | 0.000 |
| secondary: lin share ~ LAM_NE share | all | 17 | 0.91 | 0.000 |
| secondary: lin share ~ LAM_NE share | within C | 5 | 0.30 | 0.624 |
| secondary: lin share ~ LAM_NE share | within H | 5 | 0.60 | 0.285 |
| secondary: lin share ~ LAM_NE share | within R | 7 | 0.89 | 0.007 |
| secondary: lin share ~ LAM_NE share | pooled within T (rank-demeaned) | 17 | 0.70 | 0.004 |
| secondary: lin share ~ LAM_PE share | all | 17 | -0.50 | 0.040 |
| secondary: lin share ~ LAM_PE share | within C | 5 | 0.60 | 0.285 |
| secondary: lin share ~ LAM_PE share | within H | 5 | -0.20 | 0.747 |
| secondary: lin share ~ LAM_PE share | within R | 7 | -0.39 | 0.383 |
| secondary: lin share ~ LAM_PE share | pooled within T (rank-demeaned) | 17 | -0.30 | 0.285 |
| H-D2b irreversible expansion (clean) ~ lin share | all | 17 | -0.32 | 0.216 |
| H-D2b irreversible expansion (clean) ~ lin share | within C | 5 | -0.10 | 0.873 |
| H-D2b irreversible expansion (clean) ~ lin share | within H | 5 | -0.60 | 0.285 |
| H-D2b irreversible expansion (clean) ~ lin share | within R | 7 | -0.43 | 0.337 |
| H-D2b irreversible expansion (clean) ~ lin share | pooled within T (rank-demeaned) | 17 | -0.26 | 0.345 |
| secondary: expansion dX_irr_r1 ~ lin share | all | 16 | -0.06 | 0.812 |
| secondary: expansion dX_irr_r1 ~ lin share | within C | 5 | 0.60 | 0.285 |
| secondary: expansion dX_irr_r1 ~ lin share | within H | 4 | -0.80 | 0.200 |
| secondary: expansion dX_irr_r1 ~ lin share | within R | 7 | -0.64 | 0.119 |
| secondary: expansion dX_irr_r1 ~ lin share | pooled within T (rank-demeaned) | 16 | -0.35 | 0.213 |

## A,B recomputed unrounded, all 21 cells

| test | subset | n | rho | p |
|---|---|---|---|---|
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | all | 21 | 0.86 | 0.000 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within C | 7 | 0.57 | 0.180 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within H | 7 | 0.93 | 0.003 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | within R | 7 | 0.86 | 0.014 |
| H-D1 sqrt share ~ LLI share (== H-D2a mirrored) | pooled within T (rank-demeaned) | 21 | 0.81 | 0.000 |
| H-D2a lin share ~ LAM share | all | 21 | 0.86 | 0.000 |
| H-D2a lin share ~ LAM share | within C | 7 | 0.57 | 0.180 |
| H-D2a lin share ~ LAM share | within H | 7 | 0.93 | 0.003 |
| H-D2a lin share ~ LAM share | within R | 7 | 0.86 | 0.014 |
| H-D2a lin share ~ LAM share | pooled within T (rank-demeaned) | 21 | 0.81 | 0.000 |
| secondary: lin share ~ LAM_NE share | all | 21 | 0.93 | 0.000 |
| secondary: lin share ~ LAM_NE share | within C | 7 | 0.43 | 0.337 |
| secondary: lin share ~ LAM_NE share | within H | 7 | 0.86 | 0.014 |
| secondary: lin share ~ LAM_NE share | within R | 7 | 0.89 | 0.007 |
| secondary: lin share ~ LAM_NE share | pooled within T (rank-demeaned) | 21 | 0.70 | 0.001 |
| secondary: lin share ~ LAM_PE share | all | 21 | -0.55 | 0.009 |
| secondary: lin share ~ LAM_PE share | within C | 7 | 0.50 | 0.253 |
| secondary: lin share ~ LAM_PE share | within H | 7 | -0.71 | 0.071 |
| secondary: lin share ~ LAM_PE share | within R | 7 | -0.39 | 0.383 |
| secondary: lin share ~ LAM_PE share | pooled within T (rank-demeaned) | 21 | -0.27 | 0.269 |
| H-D2b irreversible expansion (clean) ~ lin share | all | 21 | -0.42 | 0.058 |
| H-D2b irreversible expansion (clean) ~ lin share | within C | 7 | -0.21 | 0.645 |
| H-D2b irreversible expansion (clean) ~ lin share | within H | 7 | -0.64 | 0.119 |
| H-D2b irreversible expansion (clean) ~ lin share | within R | 7 | -0.43 | 0.337 |
| H-D2b irreversible expansion (clean) ~ lin share | pooled within T (rank-demeaned) | 21 | -0.33 | 0.167 |
| secondary: expansion dX_irr_r1 ~ lin share | all | 20 | 0.02 | 0.925 |
| secondary: expansion dX_irr_r1 ~ lin share | within C | 7 | 0.11 | 0.819 |
| secondary: expansion dX_irr_r1 ~ lin share | within H | 6 | -0.60 | 0.208 |
| secondary: expansion dX_irr_r1 ~ lin share | within R | 7 | -0.64 | 0.119 |
| secondary: expansion dX_irr_r1 ~ lin share | pooled within T (rank-demeaned) | 20 | -0.28 | 0.261 |
