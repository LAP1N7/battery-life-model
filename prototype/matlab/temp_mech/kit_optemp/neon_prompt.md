Write ONE MATLAB script file `kit_optemp_main.m` (plus local functions at the end of the same file). Output ONLY the MATLAB code, no markdown fences, no commentary. Comments in the code may be Korean or English. Target MATLAB R2023a or newer, base MATLAB only (no toolboxes except Statistics and Machine Learning Toolbox for `corr(...,'Type','Spearman')`; if you use it, guard with a fallback that ranks manually).

## Goal
KIT battery dataset: at each check-up the same cell was measured by EIS twice: once at its ageing (operating) temperature (`is_rt == 0`) and once at room temperature ~25 °C (`is_rt == 1`). We separate:
- ageing effect = change of a resistance component measured at 25 °C over life;
- instantaneous temperature effect = apparent activation energy of each resistance component, from the pair (op-T spectrum, RT spectrum) of the same check-up.
We want, per component, the instantaneous Ea and whether it changes with ageing (SOH) and with ageing temperature.

## Inputs (paths relative to the script folder; define them at the top as variables)
- RAW_DIR = '..\..\..\data\kit_radar\eis_addendum\cell_eisv2_fixed'   files `cell_eisv2_<cell>.csv`, semicolon-separated, header row. Columns used: timestamp_s, sd_block_id, is_rt, soc_nom, valid, t_avg_degC, freq_Hz, z_amp_mOhm, z_ph_deg. (Other columns exist; read only these with detectImportOptions + SelectedVariableNames, Delimiter ';'.) A sample header + first rows is attached as raw_sample_head.txt.
- `kit_meta.csv` (attached head): cell, T_age (0/10/25/40), age_type (1 calendar, 2 cycling, 3 WLTP profile), chg_rate, dis_rate, cal_soc, group, v_min, v_max.
- `kit_soh.csv`: cell, t (unix s), day, efc, soh — SOH known at RT spectrum times; interpolate linearly in t per cell, no extrapolation beyond 3 days outside the range (NaN then).
- `eis_features.m` (attached; it is on the path at ..\eis_features.m — add `addpath('..')`). Call `d = eis_features(f, re, mim)` with f [Hz], re = Re Z, mim = −Im Z (capacitive positive). It returns R0, R_arc, f_valley etc.

## Spectrum construction
- Keep rows with valid == 1 and freq_Hz not NaN. One spectrum = one sd_block_id. Require >= 25 frequency points.
- Re = z_amp_mOhm .* cos(deg2rad(z_ph_deg)); mim = z_amp_mOhm .* sin(deg2rad(z_ph_deg)). IMPORTANT: in this file the sign convention already gives capacitive = positive for this expression, so pass it directly as mim (no extra minus sign). Add a sanity check: for each spectrum the median of mim over 1–100 Hz must be > 0, otherwise flag and skip.
- Per spectrum store: cell, block, is_rt, soc (soc_nom), t (first timestamp_s), T_meas (mean t_avg_degC), R0, R_arc (from eis_features), and a low-frequency proxy R_lf = Re at the lowest frequency minus (R0 + R_arc) when R_arc is valid (diffusion tail, NaN otherwise), and Re_lowf = Re at lowest frequency (total polarization proxy). Units mOhm.

## Pairing
- Use only cells with T_age in {0, 10, 40} (cells at 25 °C have no separate op-T measurement; verify by counting is_rt==0 spectra per T_age and printing the counts).
- For each op-T spectrum (is_rt==0) find the RT spectrum (is_rt==1) of the same cell and same soc closest in time; accept if |Δt| < 2 days. Keep pairs with |T_meas_op − T_meas_rt| >= 8 K.
- For each pair and each component X in {R0, R_arc, R_lf, Re_lowf}: Ea_X = 8.314e-3 * log(X_op / X_rt) / (1/(T_op+273.15) − 1/(T_rt+273.15)) in kJ/mol (positive = resistance higher when colder). Skip if either value <= 0 or NaN.
- Attach soh (interpolated at the RT spectrum time), day, efc, and the meta columns.

## Outputs (write to subfolder `out\`, create it)
1. `optemp_pairs.csv`: one row per pair with all fields above.
2. `optemp_summary.csv`: for each component × T_age × age_type(1/2/3) × SOH bin (>=0.98, 0.95–0.98, 0.90–0.95, <0.90) × soc: n_pairs, n_cells, median Ea, 25th and 75th percentile.
3. `optemp_trend.csv`: for each component × T_age × soc: Spearman rho of Ea vs (1 − soh) using cell-level aggregation to avoid pseudo-replication: for each cell take the slope of Ea vs (1−soh) by least squares (cells with >= 4 pairs), then report median slope, fraction of cells with positive slope, n_cells, and a sign-test p-value (binomial, two-sided).
4. Print to the Command Window: pair counts per T_age; the SOC-50 median Ea per component and T_age at SOH >= 0.98 (fresh) and at SOH < 0.95 (aged).
5. Figures (save PNG with exportgraphics, 150 dpi, white background):
   - `fig_optemp_Ea_vs_soh.png`: 2×2 tiles, one per component, SOC 50 only, x = SOH (reverse axis so ageing goes right), y = Ea [kJ/mol], points coloured by T_age (0 °C blue [0 0.45 0.70], 10 °C light blue [0.34 0.71 0.91], 40 °C orange [0.84 0.37 0]), plus per-T_age binned medians as lines.
   - `fig_optemp_nyquist_example.png`: one 10 °C cycling cell and one 40 °C cycling cell (pick the cell with the most pairs in each group), SOC 50: Nyquist (Re vs −Im, axis equal) of op-T and RT spectra at the first pair and at the pair closest to SOH 0.90 (or the last pair). Solid = RT, dashed = op-T. Title with cell id, T and SOH.

## Robustness / style
- Loop over cells with a progress print every 20 cells; total runtime should be a few minutes. Preallocate or collect rows in cell arrays then build tables once.
- Wrap per-file reading in try/catch; log skipped files and reasons into `out\skipped.txt`.
- No global variables. Clear, sectioned code (%% headers): 0 설정, 1 스펙트럼 추출, 2 짝짓기, 3 Ea 계산, 4 요약/추세, 5 그림.
- Do not hard-code results; everything computed from the files.
