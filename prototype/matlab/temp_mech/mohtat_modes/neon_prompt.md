Write ONE MATLAB script `mohtat_modes_main.m` (local functions at the end of the same file). Output ONLY MATLAB code, no markdown fences, no commentary. Comments may be Korean. MATLAB R2023a+, base MATLAB only (no Signal/Statistics/Optimization toolboxes: implement smoothing, peak tracking, nonnegative least squares via base `lsqnonneg`, and bootstrap yourself).

## Context
Mohtat 2020 dataset: 21 NMC/graphite 5 Ah pouch cells, aged at −5 °C (T='C'), 25 °C ('R'), 45 °C ('H'); one cell per condition; 7 cycling protocols per temperature. At each check-up there is a slow C/20 charge/discharge curve. We already have per-check-up degradation modes from half-cell DVA (LLI, LAM_PE, LAM_NE in %) and EIS resistance changes (dR0, dRarc in mOhm). Three analyses:
 (1) ICA (dQ/dV) of the C/20 charge curve: how the graphite/NMC peaks change with ageing at each temperature, and whether ICA features move with the DVA modes.
 (2) Mode kinetics: each mode m(t) = a·sqrt(days) + b·EFC per cell; which shape dominates at each temperature.
 (3) Mode → resistance link: per temperature, ΔR = c_LLI·LLI + c_PE·LAM_PE + c_NE·LAM_NE, with cell-bootstrap CIs.

## Inputs (define at top)
- JOINED = 'mohtat_modes_joined.csv' (same folder). Columns: cell (1–21), T ('C','R','H'), T_C, chg (charge C-rate), dis, dod, ah (cumulative Ah, identifies the check-up), days (since reference check-up), efc (equivalent full cycles since reference), loss (% capacity loss vs reference), LLI, LAM_PE, LAM_NE (%), ref_not_bol (true for hot cells whose reference is the first check-up after BOL), dR0, dRarc (mOhm vs reference, NaN where no EIS). One row per check-up. Head attached.
- RAW = '..\..\..\Mothat\data' with subfolders '01'…'21', each containing 'OCV_wExpansion.csv' (comma-separated, one header row, 9 numeric columns in this order: time_s, current_mA, voltage_V, expansion, temperature_C, Q, capacity_Ah, cycle, total_Ah_throughput). Read with readmatrix (NumHeaderLines 1). Files can be large; read one at a time, keep only needed columns.
- Colours: C [0 0.45 0.70], R [0 0.62 0.45], H [0.84 0.37 0].

## (1) ICA
- In each OCV file group rows by the 9th column (total_Ah_throughput, identifies the check-up; values are repeated within a check-up — group with uniquetol, tolerance 1e-6 relative). Charge segment = rows with 100 < current_mA < 400; require >= 500 rows. Q = capacity_Ah (column 7), V = voltage (column 3). Sort by Q, remove non-increasing V jitter by keeping a monotone envelope (cummax of V), then interpolate Q on a uniform V grid from 3.45 to 4.18 V, step 0.002 V; dQ/dV by central difference; smooth with a Gaussian kernel of sigma = 4 grid points (implement with conv), dQ/dV in Ah/V.
- Peak tracking per cell: on the reference check-up curve (the check-up whose total_Ah is closest to the smallest `ah` of that cell in JOINED), find local maxima with prominence >= 5 % of the curve max inside V windows W1 = [3.50, 3.70], W2 = [3.70, 3.85], W3 = [3.85, 4.10] (one peak per window: the highest). For every later check-up record, per window: peak height (max of dQ/dV inside the window), peak position (V at the max), window area (integral of dQ/dV over the window = capacity charged in that window, Ah). Also the total charge capacity.
- Normalise: height_rel = height/height_ref, area_rel = area/area_ref, shift_mV = 1000·(pos − pos_ref).
- Join ICA features to JOINED by cell and nearest ah (|Δah| <= 5). Save `out\ica_features.csv`.
- Spearman (implement via ranks + Pearson; ties by average rank) between each ICA feature (3 windows × {height_rel, area_rel, shift_mV}) and each DVA mode (LLI, LAM_PE, LAM_NE), computed (a) pooled over all check-ups with loss > 0 and (b) within each temperature group. Also partial version controlling for loss: correlate residuals after linear regression of both variables on loss. Save `out\ica_vs_modes.csv` (feature, mode, scope, rho, rho_partial, n).
- Figure `out\fig_ica_curves.png`: 1×3 tiles for cells 2 (C), 1 (R), 3 (H): dQ/dV vs V at the reference check-up (dashed) and at the check-up closest to 10 % loss (solid), with the three windows shaded lightly; title with cell, T, loss.
- Figure `out\fig_ica_peaks_vs_loss.png`: 3×3 tiles (rows = windows, cols = height_rel / area_rel / shift_mV) vs loss, points coloured by temperature.

## (2) Mode kinetics
- For each cell and each mode m ∈ {LLI, LAM_PE, LAM_NE}: fit m = a·sqrt(days) + b·efc with a, b >= 0 using lsqnonneg on rows with days >= 0. Report a, b, RMSE, R², sqrt-share at the last check-up = a·sqrt(days_end)/(a·sqrt(days_end) + b·efc_end) (NaN if both 0), and a shape-only alternative: fit m = A·days^z (log-linear on rows with m > 0.05 and days > 0) and report z.
- IMPORTANT note to print and to write in a comment: within one constant-protocol cell efc is proportional to days, so 'a' and 'b' are distinguished only by curve shape (sqrt vs linear), not by time vs cycling. Call them sqrt-shape and linear-shape terms.
- Group summaries per temperature × mode: median and range of sqrt-share and z, n cells. Save `out\mode_kinetics_cells.csv`, `out\mode_kinetics_groups.csv`.
- Figure `out\fig_mode_trajectories.png`: 1×3 tiles (LLI, LAM_PE, LAM_NE) vs days (log x is NOT needed; use linear), each cell a thin line coloured by temperature, plus the fitted curve dashed.

## (3) Mode → resistance link
- For each temperature group separately and for ΔR ∈ {dR0, dRarc}: rows with finite ΔR and loss > 0. Fit ΔR = c_LLI·LLI + c_PE·LAM_PE + c_NE·LAM_NE (no intercept; unconstrained least squares with `\`). Report coefficients (mOhm per %p), R², n rows, n cells, and the condition number of the design matrix and pairwise correlations among the three modes in that group (collinearity diagnostic).
- Cell-cluster bootstrap: 1000 resamples of cells (with replacement) within the group, rng(1); 95 % percentile CI for each coefficient; flag 'clear' if the CI excludes 0.
- Also a pooled model across all temperatures with temperature-specific coefficients and a model with shared coefficients; compare by AIC = n·log(RSS/n)+2k.
- Also report, per temperature, the simple sensitivities dΔR/dloss (slope through origin of ΔR on loss) with bootstrap CI — "mOhm per % capacity loss".
- Save `out\mode_resistance_link.csv`, `out\resistance_per_loss.csv`. Print compact tables.
- Figure `out\fig_link_coeffs.png`: coefficient ± CI bars per mode, grouped by temperature, panels dR0 and dRarc, horizontal zero line.

## Output / style
- Create `out\`. Sectioned code (%% 0 설정, 1 ICA, 2 모드 진행식, 3 모드-저항 연결, 4 그림). Print progress and the key tables. try/catch around each cell's file reading; log to `out\log.txt`. exportgraphics 150 dpi, white background. No hard-coded results. Runtime should be a few minutes.
