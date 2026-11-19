Write ONE MATLAB script file `kit_models_main.m` (local functions at the end of the same file). Output ONLY MATLAB code, no markdown fences, no commentary. Comments may be Korean or English. MATLAB R2023a+, base MATLAB plus Optimization Toolbox (`lsqnonlin`) — if `lsqnonlin` is unavailable (check with `exist('lsqnonlin','file')`), fall back to `fminsearch` on the sum of squares with parameter transforms for bounds. No other toolboxes.

## Context
Battery resistance ageing study (KIT dataset, 228 cells). Two EIS-derived components at SOC 50 %, measured at 25 °C: R0 (high-frequency intercept, mOhm) and R_arc (mid-frequency arc width, mOhm). We test three things:
 (A) Arc decomposition: does the drop of the arc's instantaneous activation energy with ageing come from a NEW resistance with a lower Ea added to the original charge-transfer resistance?
 (B) Leave-one-temperature-out (LOTO) model comparison: which degradation formula predicts a held-out ageing temperature best?
 (C) Real-use check: models trained on constant-condition cells predict WLTP drive-profile cells.

## Inputs (same folder; define paths at top)
1. `kit_rt_soc50.csv`: columns cell, day, efc, soh, R0, R_arc, T_age (°C: 0/10/25/40), age_type (1 calendar, 2 constant cycling, 3 WLTP profile), chg_rate, dis_rate, cal_soc (calendar storage SOC %, 0 for others), v_min, v_max, group. One row per check-up spectrum (25 °C, SOC 50).
2. `..\kit_optemp\out\optemp_pairs.csv`: paired spectra of the same check-up measured at the ageing temperature (op) and at ~25 °C (rt). Columns used: cell, soc, T_meas_op, T_meas_rt, R0_op, R0_rt, R_arc_op, R_arc_rt, soh, day, efc, T_age, age_type. Use soc == 50 only, and rows with finite positive R_arc_op, R_arc_rt.

Constants: Rg = 8.314e-3 kJ/mol/K, Tref = 298.15 K. Temperatures in K = °C + 273.15. Define at top: I_WLTP = 1.0 (assumed effective charge C-rate for WLTP cells, metadata has 0) and SOC_CYC = 50 (effective SOC for cycling/WLTP cells).

## Preprocessing (both A and B)
- Per cell sort by day. Baseline = mean of the first two rows with soh >= 0.99 (skip cell if none). ΔR0 = R0 − R0_bol, ΔRarc = R_arc − Rarc_bol. Shift day and efc so the baseline row is 0 (day − day_bol, efc − efc_bol; calendar cells have efc 0 → force 0).
- Keep rows up to the first time soh < 0.80 (inclusive) to avoid end-of-life knees.
- Stressors per row: TK = T_age+273.15; x = (1/TK − 1/Tref)/Rg; SOC = cal_soc for calendar, SOC_CYC otherwise; I = chg_rate for cycling, I_WLTP for WLTP, 0 for calendar.

## (A) Arc decomposition with the paired data
Model for a cell with baseline arc R_b (its R_arc_rt at the first pair with soh >= 0.99; if none skip the cell) and measured increment ΔR = R_arc_rt − R_b at 25 °C:
  R_arc_op_pred = R_b · exp(E_ct·(1/T_op − 1/T_rt)/Rg) + max(ΔR,0) · exp(E_new·(1/T_op − 1/T_rt)/Rg)
(T_op, T_rt in K from T_meas_op, T_meas_rt). Residual = log(R_arc_op) − log(R_arc_op_pred).
Fit variants by least squares over all pairs (T_age in {0,10,40}):
  A0: single E (E_new = E_ct), one parameter — the "sensitivity is constant" null model.
  A1: E_ct and one shared E_new for all ageing temperatures (2 params).
  A2: E_ct shared, separate E_new for T_age 0, 10, 40 (4 params).
Bounds E in [−50, 150] kJ/mol. Report for each: parameters with 95 % CI from the Jacobian (if lsqnonlin) or from a cell-level bootstrap (200 resamples of cells, always do the bootstrap for A1 and A2), RSS, n, AIC = n·log(RSS/n) + 2k, BIC. Also report the same three variants for R0 (replace R_arc by R0).
Also, per T_age group, a "fresh-only" check: E fitted on pairs with soh >= 0.98 only (single E), to see whether fresh Ea differs by measurement window.
Write `out\A_arc_decomposition.csv` (variant, component, param, value, ci_lo, ci_hi, RSS, n, k, AIC, BIC) and print a compact table. Figure `out\fig_A_decomp.png`: apparent Ea (from the pairs, as in the source data: Ea = Rg·log(R_op/R_rt)/(1/T_op−1/T_rt)) vs SOH, binned medians per T_age (dots), and the A1 and A2 model-implied apparent Ea curves per T_age (lines), arc panel and R0 panel. The model-implied apparent Ea uses the same formula applied to R_arc_op_pred.

## (B) LOTO model comparison (age_type 1 and 2 cells only for training/testing)
Candidate formulas for ΔR (fit separately for ΔR0 and ΔRarc, residual in mOhm, each cell weighted 1/sqrt(n_rows of that cell)):
  M1 (classic accelerated test, one process): ΔR = exp(a − Ea·x) · day^z. Params a, Ea, z (z in [0.1, 2]).
  M2 (BLAST-Lite-type structure: time power + cycle power, each Arrhenius; mirrors BLAST-Lite NMC resistance law r1·q1·t^r2 + r3·q3·EFC^r4): ΔR = exp(a1 − E1·x + s·(SOC−50)/50) · day^z + exp(a2 − E2·x) · max(I,1e-3)^m · efc^w. Params a1,E1,s,z,a2,E2,m,w (z,w in [0.1,2]; E in [−150,150]; s in [−6,6]; m in [−1,3]).
  M3 (current model ④: same as M2 with z = 0.5 and w = 1 fixed).
  M4 (M3 plus a low-temperature branch on the cycle term): cycle term multiplied by [exp(−E2·x) + exp(c + Ec·x)] instead of exp(−E2·x), Ec in [0,150], c in [−10,10].
For LOTO: for each held-out T in {0,10,25,40}, fit on the other three temperatures (calendar + cycling cells), predict all rows of held-out cells. Metrics per (model, component, held-out T): RMSE (mOhm), median of pred/meas at each cell's last row (only cells with meas > 0.05 mOhm at last row), fraction of cells over-predicted. Also an in-sample fit on all four temperatures (RMSE, AIC, parameters). Use multistart (at least 6 random starts within bounds, fixed rng(1)) and keep the best.
Write `out\B_loto.csv`, `out\B_insample_params.csv`; print a compact table: rows model, columns held-out T, cells "RMSE | median ratio" for each component. Figure `out\fig_B_loto.png`: grouped bars of median pred/meas ratio (log y axis, line at 1) by held-out T, one panel per component, bars per model.

## (C) WLTP check
Fit M1–M4 in-sample on age_type 1+2 cells (all temperatures), predict age_type 3 (WLTP) cells. Metrics per model and component: RMSE, median pred/meas at last row, per T_age. Write `out\C_wltp.csv`, print table. Note in a printed line that I_WLTP is an assumption; rerun M2–M4 predictions with I_WLTP = 0.5 and 2.0 and report the median ratio range as a sensitivity line.

## Style / robustness
Sectioned code (%% 0 설정, 1 전처리, 2 (A) 아크 분해, 3 (B) LOTO, 4 (C) WLTP, 5 그림). Create `out\`. Print progress. No globals. Guard every fit with try/catch and log failures to `out\fit_log.txt`. Total runtime should be under ~15 minutes; if the multistart is slow, reduce starts but keep >= 4. Save figures with exportgraphics (150 dpi, white background). Never hard-code results.
