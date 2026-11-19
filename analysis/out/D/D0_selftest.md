# D0 self-test: parameter recovery (synthetic, known truth)

OCP pair for fitting: `mohtat_peyman`; BOL truth Cp=5.6 Ah, Cn=6.2 Ah, n_Li=5.26 Ah. Each case = BOL + one aged check-up with independent noise; modes are recovered relative to the fitted BOL (as in D1). Error = recovered - true, in percentage points (pp). 'random' cases: LLI 0-25 %, LAM_PE 0-15 %, LAM_NE 0-25 % (uniform, independent); plus 5 structured patterns x 4 repeats.

Nuisance (scale 1): resistance 5-20 mohm x 0.25 A with low-SOC rise, charge/discharge asymmetry +-30 %, hysteresis 4-12 mV with SOC modulation, white 0.4 mV + three slow sinusoids ~1 mV, charge start voltage 3.0-3.25 V. 'nuisance_x2' doubles all; 'ocp_mismatch' (x1, x3, x6) generates with a graphite OCP stretched by 1 % (3 %, 6 %) in stoichiometry and a 5 mV (15, 30 mV) NMC ripple but fits with the nominal OCP (this is the structural-error test; real D1 fits have 2-7 mV RMSE); 'chg_only_fallback' fits the charge branch only with a fitted offset.

## Recovery error (pp), random cases and all cases

| index | scenario | subset | n | fit_rmse_median_mV | LLI_bias | LLI_rmse | LLI_p95 | LAM_PE_bias | LAM_PE_rmse | LAM_PE_p95 | LAM_NE_bias | LAM_NE_rmse | LAM_NE_p95 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | clean | all | 60 | 0.0 | -0.0 | 0.0 | 0.0 | -0.0 | 0.0 | 0.0 | -0.0 | 0.0 | 0.0 |
| 1 | clean | random only | 40 | 0.0 | -0.0 | 0.0 | 0.0 | -0.0 | 0.0 | 0.0 | -0.0 | 0.0 | 0.0 |
| 2 | nuisance | all | 60 | 0.54 | -0.01 | 0.08 | 0.14 | -0.01 | 0.36 | 0.84 | -0.02 | 0.35 | 0.63 |
| 3 | nuisance | random only | 40 | 0.56 | -0.01 | 0.08 | 0.13 | 0.01 | 0.38 | 0.86 | 0.01 | 0.3 | 0.56 |
| 4 | nuisance_x2 | all | 60 | 0.99 | -0.02 | 0.15 | 0.33 | -0.03 | 0.62 | 1.19 | -0.07 | 0.65 | 1.2 |
| 5 | nuisance_x2 | random only | 40 | 0.99 | -0.01 | 0.13 | 0.23 | -0.06 | 0.61 | 1.19 | -0.12 | 0.53 | 1.01 |
| 6 | ocp_mismatch | all | 60 | 2.61 | 0.08 | 0.33 | 0.64 | 0.42 | 0.75 | 1.43 | 0.13 | 0.96 | 1.96 |
| 7 | ocp_mismatch | random only | 40 | 2.61 | 0.13 | 0.25 | 0.47 | 0.41 | 0.69 | 1.43 | 0.27 | 0.74 | 1.28 |
| 8 | ocp_mismatch_x3 | all | 60 | 7.52 | -0.29 | 1.96 | 4.52 | 1.3 | 2.08 | 3.34 | -1.06 | 5.73 | 13.67 |
| 9 | ocp_mismatch_x3 | random only | 40 | 7.45 | -0.31 | 2.09 | 4.73 | 1.4 | 2.14 | 3.5 | -1.18 | 6.36 | 16.21 |
| 10 | ocp_mismatch_x6 | all | 60 | 14.67 | 10.62 | 13.67 | 18.98 | -2.03 | 3.58 | 7.9 | 14.69 | 21.76 | 32.23 |
| 11 | ocp_mismatch_x6 | random only | 40 | 14.72 | 9.52 | 12.87 | 17.79 | -1.6 | 3.32 | 7.39 | 13.05 | 20.79 | 30.49 |
| 12 | chg_only_fallback | all | 60 | 0.83 | -0.04 | 0.23 | 0.42 | 0.1 | 0.83 | 1.54 | -0.15 | 0.55 | 0.96 |
| 13 | chg_only_fallback | random only | 40 | 0.83 | -0.05 | 0.26 | 0.42 | 0.09 | 0.81 | 1.74 | -0.08 | 0.56 | 1.29 |

Slope of recovered vs true (1 = perfect, <1 = shrinkage / poor identifiability), random cases:

| index | scenario | LLI_slope | LAM_PE_slope | LAM_NE_slope |
|---|---|---|---|---|
| 1 | clean | 1.0 | 1.0 | 1.0 |
| 3 | nuisance | 1.0 | 1.03 | 1.0 |
| 5 | nuisance_x2 | 1.0 | 1.01 | 1.02 |
| 7 | ocp_mismatch | 1.02 | 0.92 | 0.97 |
| 9 | ocp_mismatch_x3 | 1.22 | 0.76 | 0.85 |
| 11 | ocp_mismatch_x6 | 1.81 | 0.62 | 0.81 |
| 13 | chg_only_fallback | 1.0 | 1.01 | 0.99 |

## Identifiable? (random cases, rule: RMSE < 1.5 pp)

- clean: LLI identifiable, LAM_PE identifiable, LAM_NE identifiable
- nuisance: LLI identifiable, LAM_PE identifiable, LAM_NE identifiable
- nuisance_x2: LLI identifiable, LAM_PE identifiable, LAM_NE identifiable
- ocp_mismatch: LLI identifiable, LAM_PE identifiable, LAM_NE identifiable
- ocp_mismatch_x3: LLI POOR, LAM_PE POOR, LAM_NE POOR
- ocp_mismatch_x6: LLI POOR, LAM_PE POOR, LAM_NE POOR
- chg_only_fallback: LLI identifiable, LAM_PE identifiable, LAM_NE identifiable

## Structured patterns, scenario `nuisance` (mean recovered vs truth, pp)

| kind | true_LLI | LLI | true_PE | LAM_PE | true_NE | LAM_NE |
|---|---|---|---|---|---|---|
| LAM_NE_only | 0.0 | 0.0 | 0.0 | -0.2 | 10.0 | 9.8 |
| LAM_PE_only | 0.0 | 0.1 | 10.0 | 10.0 | 0.0 | 0.1 |
| LLI+LAM_NE | 10.0 | 9.9 | 0.0 | 0.3 | 10.0 | 9.8 |
| LLI_only | 10.0 | 9.9 | 0.0 | -0.0 | 0.0 | -0.1 |
| all | 10.0 | 10.0 | 5.0 | 4.7 | 8.0 | 8.0 |

## Collinearity and the LAM_PE question

error correlations (nuisance, random): LLI-LAM_NE +0.20, LAM_PE-LLI +0.29, LAM_PE-LAM_NE -0.09. A strong negative/positive value means the fit trades one mode against the other.
- LLI-limited (true LLI > LAM_PE + 5 pp) = False: n=21, LAM_PE rmse 0.30 pp, LLI rmse 0.08, LAM_NE rmse 0.30
- LLI-limited (true LLI > LAM_PE + 5 pp) = True: n=19, LAM_PE rmse 0.45 pp, LLI rmse 0.08, LAM_NE rmse 0.30

## Capacity-from-modes identity check

capacity ratio implied by TRUE modes vs ratio in the synthetic data (3.0 V loaded discharge end): median |diff| 0.01 pp, max 0.08 pp (differences come from the loaded 3.0 V cut-off vs the OCV limit; this is the floor of completion criterion C2).