Write ONE MATLAB script `temp_eqs_main.m` (local functions at the end of the same file). Output ONLY MATLAB code, no markdown fences, no commentary. Comments may be Korean. MATLAB R2023a+, base MATLAB only (no Statistics/Optimization toolboxes: use `\`, `fminsearch`, and implement bootstrap and Spearman yourself; Spearman = Pearson on average-tie ranks).

## Context
Two public ageing datasets with EIS check-ups. For each check-up we have capacity loss L (%) and increases of ohmic resistance R0 and of the mid-frequency arc R_arc relative to the cell's reference check-up. Earlier finding (Mohtat only): at equal capacity loss, cold-aged cells grow the arc more than R0, hot-aged cells grow R0 more. Three tasks: (K) re-test this "cold = arc-type, hot = R0-type" pattern on KIT cells aged at 0/10/25/40 °C; (S) turn the per-temperature resistance-per-loss slopes s_T into a function of temperature s(T); (P) within-temperature test whether modes beyond LLI explain resistance at equal loss (Mohtat).

## Inputs
1. KIT = '..\kit_models\kit_rt_soc50.csv'. Columns: cell, day, efc, soh (fraction), R0, R_arc (mOhm, measured at 25 °C, SOC 50), T_age (0/10/25/40), age_type (1 calendar, 2 constant cycling, 3 WLTP), chg_rate, dis_rate, cal_soc, v_min, v_max, group. Head attached. Reference = the cell's first check-up (smallest day). L = 100·(soh_ref − soh)/soh_ref. dR0 = R0 − R0_ref, dRarc = R_arc − R_arc_ref. Protocol key for pairing = (age_type, chg_rate, dis_rate, v_min, v_max, cal_soc).
2. MOH = '..\mohtat_modes\mohtat_modes_joined.csv'. Columns: cell, T ('C' −5 °C, 'R' 25 °C, 'H' 45 °C), T_C, chg, dis, dod, ah, days, efc, loss (L, %), LLI, LAM_PE, LAM_NE (%), ref_not_bol, dR0, dRarc (mOhm, already relative to reference). Head attached.
3. MREF = 'mohtat_ref_resistance.csv': cell, T, R0_ref, Rarc_ref (mOhm) for normalisation.
Colours by temperature: −5 [0 0.30 0.60], 0 [0 0.45 0.70], 10 [0.35 0.70 0.90], 25 [0 0.62 0.45], 40 [0.90 0.60 0], 45 [0.84 0.37 0].

## (K) KIT equal-loss comparison
- Cells with age_type 2 and 3 (cycling). Per cell, make L monotone with cummax; interpolate dR0, dRarc at L = 5, 10, 20 % (skip levels the cell never reaches). ratio = dRarc/dR0 (NaN if dR0 <= 0).
- Per T_age × level: n, medians of dR0, dRarc, ratio, and median of dR0/R0_ref and dRarc/R_arc_ref (relative). Save `out\K_equal_loss_cells.csv`, `out\K_equal_loss_groups.csv`.
- Paired test: for each protocol key present at both 0 °C and 25 °C (and separately 10 vs 25, 40 vs 25), compare ratio at the same level; report count of pairs where colder ratio > 25 °C ratio, n pairs, and median difference. Also a sign test p-value (binomial, implement with nchoosek sums). Save `out\K_paired.csv`.
- Also calendar cells (age_type 1) at L = 5 % only, same table (may be few).

## (S) Slope as function of temperature
- Per cell slope through origin s_cell = Σ(L·ΔR)/Σ(L²) on rows with L > 0, for dR0 and dRarc, and normalised slopes s_cell/R_ref (unit: fraction of initial resistance per % loss, ×100 to print as %/%). KIT: age_type 2 cells, require ≥ 3 rows with L > 0 and max L ≥ 5. Mohtat: all cells with ≥ 3 usable rows (finite dR, loss > 0); R_ref from MREF.
- Per dataset × temperature: median s, IQR, n. Save `out\S_cell_slopes.csv`, `out\S_group.csv`.
- Fit on cell-level normalised slopes (log scale, least squares on log s, cells with s > 0), with x = 1000/T_K:
  A1 single Arrhenius: log s = a + b·(x − x_ref), x_ref = 1000/298.15. Report apparent Ea = b·R·1000 J/mol in kJ/mol (positive Ea = slope grows with temperature).
  A2 two-branch: s = exp(a1 − Eh·(x − x_ref)·1000/R) + exp(a2 + Ec·(x − x_ref)·1000/R) with Eh, Ec >= 0 (fminsearch on log-parameters, several starts), i.e. one branch rising with T and one rising as T falls.
  Fit per dataset (KIT 4 temperatures; Mohtat 3 temperatures — A2 with 4 parameters on 3 temperatures is not identifiable: still fit A1 and print a note) and jointly (both datasets, dataset-specific a / a1,a2 offsets, shared Ea / Eh, Ec). Compare by AIC. Cell bootstrap (1000, rng(1), resample cells within dataset × temperature) CIs for Ea, Eh, Ec. Do this for R0, arc, and for the ratio s_arc/s_R0 (A1 only).
- Save `out\S_arrhenius.csv` (dataset, target, model, params, CIs, AIC, n). Print a compact table and print the final equation strings, e.g. "s_R0(T) = s_ref·exp[−Ea/R·(1/T − 1/298.15)], s_ref = ..., Ea = ... kJ/mol".

## (P) Mohtat within-temperature mode test at equal loss
- From MOH, per cell interpolate LLI, LAM_PE, LAM_NE, dR0, dRarc at L = 5 and 10 % (cummax L). Within each T (and pooled with T as fixed effect = demean within T): partial Spearman of dR0 and dRarc with LAM_PE and with LAM_NE controlling for LLI (regress ranks on LLI rank, correlate residuals). Report rho, n. Also within T the correlation matrix of LLI, LAM_PE, LAM_NE across check-ups (all rows, not interpolated) and the condition number of [LLI LAM_PE LAM_NE]. Save `out\P_mode_partial.csv`, `out\P_collinearity.csv`. Print with a line: n per T is 7 cells at most, treat as exploratory.

## Figures (exportgraphics 150 dpi, white background)
- `out\fig_K_equal_loss.png`: 1×3 (L = 5, 10, 20 %): ratio dRarc/dR0 per cell vs T_age (jittered points + median marker), KIT cycling cells, plus Mohtat 5/10 % medians as diamonds at −5/25/45 °C.
- `out\fig_S_arrhenius.png`: 1×3 (R0, arc, ratio): log normalised slope vs 1000/T, cell points coloured by temperature, KIT circles, Mohtat triangles, A1 lines per dataset, A2 curve for KIT dashed.

Sectioned code (%% 0 설정, K, S, P, 그림). try/catch per section, log to `out\log.txt`. Create `out\`. No hard-coded results. Runtime under two minutes.
