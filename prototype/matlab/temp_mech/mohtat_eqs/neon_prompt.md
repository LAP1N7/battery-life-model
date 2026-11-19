Write ONE MATLAB script `mohtat_eqs_main.m` (local functions at the end of the same file). Output ONLY MATLAB code, no markdown fences, no commentary. Comments may be Korean. MATLAB R2023a+, base MATLAB only (no Statistics/Optimization/Curve Fitting toolboxes: use `\`, `fminsearch`, and implement bootstrap yourself).

## Context
Mohtat 2020 data: 21 NMC/graphite 5 Ah pouch cells aged at −5 °C (T='C'), 25 °C ('R'), 45 °C ('H'), 7 cycling protocols per temperature, one cell per protocol. Each row = one check-up with capacity loss L (%), DVA degradation modes LLI, LAM_PE, LAM_NE (%), and EIS resistance increases dR0, dRarc (mOhm) relative to the reference check-up. Goal: express per-temperature behaviour as compact equations with fitted constants and CIs, to be written up like a results section ("at temperature T, this cell follows ΔR = ...").

## Input
- JOINED = '..\mohtat_modes\mohtat_modes_joined.csv'. Columns: cell, T, T_C, chg (charge C-rate, numeric), dis (text), dod, ah, days, efc, loss, LLI, LAM_PE, LAM_NE, ref_not_bol, dR0, dRarc. Head attached. Use rows with loss > 0 and days > 0. For resistance equations additionally require finite dR0 and dRarc and dR > 0 where logs are taken. Ignore cells with fewer than 3 usable rows for per-cell fits (still use them in pooled fits).
- Colours: C [0 0.45 0.70], R [0 0.62 0.45], H [0.84 0.37 0].

## E1  Mode partition vs capacity loss (per T)
For each T and each mode m ∈ {LLI, LAM_PE, LAM_NE}: (a) m = k·L through origin (least squares), (b) m = k·L^q (nonlinear LS with fminsearch on log-parameters, start from (a) and q=1). Report k, q, R² (centred, 1−SSres/SStot), RMSE, n rows, n cells, and AIC (n·log(RSS/n)+2·npar) for (a) vs (b).

## E2  Mode time-shape (per T)
For each T and mode: log m = log A_cell + z_T·log days (cell-specific intercepts = fixed effects, one shared exponent z_T per temperature), rows with m > 0.05. Report z_T, the median and range of A_cell (back-transformed, units %/day^z), R² on the original scale. Also the same with efc instead of days (z_T_efc).

## E3  Resistance vs capacity loss (per T) — main result
For each T and target ∈ {dR0, dRarc}, compare four models:
 M1 pooled linear: ΔR = s·L
 M2 pooled power: ΔR = s·L^p
 M3 cell-specific linear: ΔR = s_cell·L
 M4 cell-specific slope + shared exponent: ΔR = s_cell·L^p_T
Fit by least squares on the original scale (fminsearch for p; for fixed p the s's are closed-form). Report parameters, R² (centred, pooled over all rows of that T), RMSE, AIC, n. For M4 report p_T and median/range of s_cell. Cell-cluster bootstrap (1000 resamples of cells with replacement within T, rng(1)) for s (M1), s and p (M2), p_T (M4): 95 % percentile CI.
Also print a variance decomposition: R² of M1, gain from cell-specific slopes (M3 − M1), gain from curvature (M4 − M3).

## E4  Protocol dependence of the per-cell slope (per T)
Using s_cell from M3 (both targets) and the cell's chg (charge C-rate) and dod: print a table cell, T, chg, dis, dod, s_cell_R0, s_cell_arc, ratio s_arc/s_R0. Spearman (ranks with average ties + Pearson, implement yourself) between s_cell and chg within each T and pooled; report n. No inference beyond that.

## E5  Composite equation check
For each T: combine E2 (LLI time law, use the median A_cell and z_T) with the E1(a) relation L ≈ LLI/k_LLI and E3 M2 to get ΔR(t) = s·(A·t^z / k_LLI)^p. For each cell with ≥ 3 rows, compute predicted vs measured dR0 and dRarc over its check-ups; report median ratio pred/meas and R² per T. This shows how well a single per-temperature equation describes an individual cell (cells differ by protocol, so expect spread; just report it).

## Outputs
Create `out\`. Save `out\E1_modes_vs_loss.csv`, `out\E2_mode_time.csv`, `out\E3_resistance_models.csv` (one row per T × target × model with all parameters and CIs as columns), `out\E3_cell_slopes.csv` (cell, T, target, s_cell M3, s_cell M4), `out\E4_protocol.csv`, `out\E5_composite.csv`. Print compact tables for each. Figures (exportgraphics 150 dpi, white background):
- `out\fig_E3_fits.png`: 2×3 tiles (rows dR0, dRarc; cols C, R, H): points ΔR vs L coloured by cell (thin lines per cell), pooled M1 line dashed black, M2 curve solid black, title with s, p and R² of M1/M2/M4.
- `out\fig_E1_modes.png`: 1×3 tiles (LLI, LAM_PE, LAM_NE) vs L, points coloured by temperature, M1 line per T.
Sectioned code (%% 0 설정, E1 … E5, 그림). try/catch with a log to `out\log.txt`. No hard-coded results. Runtime under a minute.
