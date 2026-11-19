Write ONE MATLAB script `model_dev_main.m` (local functions at the end of the same file). Output ONLY MATLAB code, no markdown fences, no commentary. Comments may be Korean. MATLAB R2023a+, base MATLAB only (no Statistics/Optimization/Curve Fitting toolboxes: use `\`, `fminsearch`, implement bootstrap and percentile CI yourself).

## Context
Mohtat 2020 data: 21 NMC/graphite pouch cells, 7 cycling protocols × 3 ambient temperatures (T: 'C' −5 °C, 'R' 25 °C, 'H' 45 °C), one cell per protocol×T. Each row = one check-up: capacity loss L (%, column loss), degradation modes LLI, LAM_PE, LAM_NE (%), resistance increases dR0, dRarc (mOhm) relative to the cell's reference check-up. Hot cells' reference is not BOL (ref_not_bol true; first usable check-up at ~47–300 Ah), cold/room cells' reference is BOL (ah = 0). Earlier per-temperature model M4: ΔR = s_cell·L^p_T (p_T shared within T). Pre-registered analyses D1–D3 below; implement exactly, no extra criteria.

## Inputs
- JOINED = '..\mohtat_modes\mohtat_modes_joined.csv'. Columns: cell, T, T_C, chg (numeric charge C-rate), dis (text), dod, ah, days, efc, loss, LLI, LAM_PE, LAM_NE, ref_not_bol, dR0, dRarc. Head attached. Rows ordered by cell then ah; the reference row has loss = 0 and dR = 0.
- TEFF = 'cell_temp.csv' with columns cell, T_eff_mean, T_eff_chg_mean (°C, measured cell surface temperature during ageing cycles). Head attached. Use T_eff_mean as T_eff.
- Protocol key = (chg, dis, dod). Each key has one cell per T.
- Colours: C [0 0.45 0.70], R [0 0.62 0.45], H [0.84 0.37 0].

## D1 Reference-point sensitivity
For each C and R cell, find the H cell with the same protocol key; A_ref = that H cell's reference ah (ah of its row with loss = 0). New reference for the C/R cell = its first check-up with ah ≥ A_ref (if none, drop the cell and log it). Recompute for rows after the new reference: L' = 100·(1 − (1 − loss/100)/(1 − loss_ref/100)) (capacity loss relative to new reference), LLI', LAM_PE', LAM_NE' = value − value_ref (percentage points), dR0' = dR0 − dR0_ref, dRarc' = dRarc − dRarc_ref, days' = days − days_ref. H cells unchanged.
For both versions (orig, shifted): fit M4 per T and target (ΔR = s_cell·L^p_T on rows with L > 0, finite dR; closed-form s_cell for fixed p, fminsearch on p, p in [0.2, 4]) and E2 per T and mode (log m = log A_cell + z_T·log days, rows m > 0.05 and days > 0, fixed effects by cell, OLS). Cell bootstrap (1000, rng(1), resample cells within T) for p_T and z_T; also the bootstrap distribution of differences p_C − p_H and p_R − p_H (and z) by pairing bootstrap draws (independent resamples per T). Report point, 95 % percentile CI, and whether the CI of the difference excludes 0. Save `out\D1_ref_sensitivity.csv` (version, target/mode, T, estimate, ci_lo, ci_hi, n_cells, n_rows) and `out\D1_differences.csv` (version, quantity, pair, diff, ci_lo, ci_hi, excludes0).

## D2 Structured models for s
Use orig version data (rows L > 0, finite dR). For target in {dR0, dRarc}, fit on the original scale by least squares over all 21 cells jointly:
 S0: s per cell, p per T (= M4; params 21+3)
 S1: log s = a_T (3), p per T (3)
 S2: log s = a_T + α·log(chg), p per T
 S3: log s = a0 + β·T_amb + α·log(chg), p = b0 + b1·T_amb  (T_amb = T_C)
 S4: as S3 with T_eff instead of T_amb (both in log s and p)
Use fminsearch on all nonlinear params with several starts (from S1/M4 solutions), MaxFunEvals large enough. Report params, RSS, n, k, AIC = n·log(RSS/n)+2k, R² (centred, pooled). LOCO for S1–S4: for each cell, refit without it, predict its rows from its T/chg/T_eff, collect residuals; report LOCO RMSE (pooled over all held-out rows) and median per-cell pred/meas at its last check-up. Cell bootstrap (500, rng(2)) CIs for α, β, b1 in S2–S4. Save `out\D2_models.csv` (target, model, params as text, k, AIC, R2, LOCO_RMSE, ratio_median) and `out\D2_params_ci.csv`. Print the pre-registered checks: for each target, LOCO RMSE of S2, S3, S4 relative to S1 (flag ≤ 0.90), and S4 vs S3.

## D3 Accelerated-test transfer (hot → room/cold)
For both versions (orig, shifted): fit on H cells only: (a) m = k·L through origin for LLI, LAM_PE, LAM_NE; (b) ΔR = s·L^p pooled (single s, p) for dR0, dRarc. For every C and R cell, interpolate measured m and ΔR at L = 5 and 10 % (cummax L; skip if not reached), predict with H fits, ratio = pred/meas (skip meas ≤ 0). Report per T × quantity × level: n, median ratio, IQR, and flag "reproduced" if median within [0.8, 1.25]. Save `out\D3_transfer.csv` and `out\D3_transfer_cells.csv`.

## Figures (exportgraphics 150 dpi, white background)
- `out\fig_D1_ref.png`: 2×1: p_T (R0, arc) and z_T (3 modes) with CIs, orig vs shifted side by side per T.
- `out\fig_D2_models.png`: 1×2 (R0, arc): LOCO RMSE and AIC per model S1–S4 (bars), plus predicted vs measured scatter for the best-LOCO model coloured by T.
- `out\fig_D3_transfer.png`: 1×2 (L=5, L=10): median pred/meas (log scale) per quantity for C and R with IQR bars, band 0.8–1.25 shaded.

Sectioned code (%% 0 설정, D1, D2, D3, 그림), try/catch per section with log to `out\log.txt`, create `out\`. No hard-coded results. Runtime under three minutes.
