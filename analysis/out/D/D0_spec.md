# Track C / D0 spec: degradation-mode pipeline (DVA half-cell fit + irreversible expansion)

Status: written before the final code; location `Mothatnalysis\out\D\` (code `Mothatnalysis\stepD*.py`, track C owns only `stepD*` files). Items marked "(fixed after sensitivity run)" are the only free choices that were
allowed to be tuned on the data; they are listed in section 9 with the reason. Nothing else was tuned toward the pass criteria.

## 0. Honesty note on blindness (read first)

An exploratory Mohtat analysis already exists (`Desktop\켄\prototype\mohtat_crosscheck.py`, results in
`prototype\results\mohtat\summary.md`). It looked at the same 21 cells, the same C/20 curves, the same OCPs, and it reported
sqrt-t share vs LLI share (rho 0.75 overall / 0.60 within temperature), linear share vs LAM share (0.35 / 0.20), and
irreversible expansion behaviour. Therefore the H-D1 / H-D2 criteria below are **NOT blind for Mohtat**. Mohtat results
are a regression/consistency check of the new, portable code and a sanity check against the exploratory numbers. The blind
test of H-D1/H-D2 has to be run on a **new dataset** (different cells, ideally different lab) with the criteria frozen here.
Differences between this pipeline and the prototype (charge+discharge average vs charge-only fit, fixed window, top-aligned
capacity axis, multistart) are intentional, so numbers will not match exactly; section 8 lists the expected differences.

## 1. Purpose

Per check-up of a cell, estimate from one slow (C/20-ish) charge/discharge curve pair:

* observed pattern: the shape of V(Q) (and dV/dQ) inside a fixed voltage window,
* degradation mode: LLI, LAM_PE, LAM_NE from a two-electrode half-cell fit,
* mechanism: **inferred only** (e.g. LLI-dominant is consistent with SEI growth or plating, but the fit cannot tell which;
  LAM_NE could be graphite cracking, isolation, or a fitting artifact of OCP mismatch). Mechanism wording never appears
  in a result column; only in discussion text, and always as "consistent with".

Second quantity: irreversible thickness change at the same SOC and same rest condition (D2), to be joined with team A's
sqrt-t (A) and linear (B) capacity-fade components (D3) for the hypotheses H-D1 / H-D2.

The tool must be reusable: dataset-agnostic core (`stepD_core.py`), dataset adapter (`stepD_loader_*.py`) producing a standard long
table, pluggable OCP pair. See `D_portability.md`.

## 2. Model and definitions

V_cell(Q) = U_p(y0 - Q/Cp) - U_n(x0 + Q/Cn)

* Cp, Cn: electrode capacities [Ah] (Cp = NMC, Cn = graphite), y0, x0: stoichiometries at the **empty state**, Q counted from
  the empty state (empty = model OCV equal to the window minimum, 3.0 V for Mohtat).
* Internally the fit uses the equivalent top-of-charge parametrisation (xT, yT at the full state), because the top end is the
  common anchor of the charge and discharge branch (section 3). y0, x0 are then reported at the empty state by root finding.
* n_Li = y0*Cp + x0*Cn (Ah of cyclable Li in the two electrodes). n_Li is invariant to where along the charge it is evaluated.
* LAM_PE = 1 - Cp/Cp_ref, LAM_NE = 1 - Cn/Cn_ref, LLI = 1 - n_Li/n_Li_ref. ref = reference check-up of the same cell (BOL where
  an averaged curve exists, see 3.4).
* Capacity from modes: given (Cp, Cn, n_Li) and the OCP pair, the cell capacity between the voltage limits is fully determined
  (solve Li-conserving states at V = Vmax and V = Vmin). Used for the completion check 2 (section 6).
* SOH (column `soh`) = measured C/20 discharge capacity / that of the reference check-up.

## 3. Inputs

3.1 Standard long table (adapter output; `stepD_loaders.py` validates): cell, checkup_idx, ah_throughput, direction
(charge|discharge), Q_Ah (Ah since the start of that branch, increasing), V, T, optional thickness_um. One charge and one
discharge branch per check-up, quasi-OCV rate (C/20 for Mohtat, 0.25 A). The charge branch includes the CV tail if present.
3.2 Chemistry-matched OCP pair from the registry (Mohtat: PyBaMM Mohtat2020 Peyman-MPM fits, section 7/C1).
3.3 Fixed cell voltage window [3.0, 4.2] V for the cycler limits. Fit window on the averaged curve: fixed voltage range
[Vfit_lo, Vfit_hi] (default 3.30-4.18 V, fixed after sensitivity run, section 9), identical for all check-ups and cells.
3.4 Curve: charge and discharge branches are aligned at the **top of charge** (s = capacity below the top state; the top is a
shared state for both branches, whereas the bottom is not: charge starts after a 1C discharge or a 0-3 h rest, at 3.02-3.3 V),
and the voltage of the two branches is averaged at equal s over their common s range. The average removes the +/- I*R and the
symmetric part of hysteresis. Fallback if one branch is missing or truncated (e.g. high-temperature BOL has no C/20 discharge;
cell 09 last charge stops at 3.71 V): single-branch fit with a fitted constant polarisation offset Voff (sign constrained);
flagged `mode_used = chg|dis`. Single-branch results are not mixed silently with averaged results: the reference is an averaged
fit wherever available and the flag marks any check-up whose curve type differs from the reference.
Reference for high-temperature cells: those cells have no BOL C/20 discharge, so the first check-up with both branches
(47 Ah ... 294 Ah throughput depending on the cell) is the reference and is flagged `ref_not_bol`. Sensitivities reported:
(a) BOL charge-only reference (all check-ups refit charge-only, `chg_bolref`), (b) same-batch room-temperature BOL reference
(absolute Cp, Cn, n_Li of the R cell of the same condition group, i.e. cells 01, 04, ... ; `avg_roommate_bol`; applied to the C and H cells), (c) reference = first averaged check-up
after BOL for every cell (`avg_ref_after_bol`), (d) fit window lower edge 3.0 V (`avg_window_3`). All in `D1_sensitivity.csv`.

## 4. Outputs (all in `out\D\`)

| file | content |
|---|---|
| `D0_spec.md` | this file |
| `D0_selftest.csv`, `D0_selftest.md`, `D0_selftest.png` | parameter-recovery test (synthetic cells with known Cp, Cn, x0, y0) |
| `D1_dva.csv` | cell, ah, soh, Cp, Cn, y0, x0, LLI, LAM_PE, LAM_NE, rmse_mV + n_starts_converged, multistart spread, Jacobian sd, flags, pass flags |
| `D1_sensitivity.csv`, `D1_summary.md`, `D1_fit_cell*.png` | alternative references, summary, 3 example fit figures |
| `D2_expansion.csv` | irreversible thickness per check-up, SOC/rest rule, artifact flags |
| `D3_join.csv`, `D3_correlations.csv`, `D3_note.md` | per-cell join and Spearman tests (input = team A `out\P1_pred.csv`; `_standin` files = prototype stand-in test) |
| `D_portability.md` | what another dataset must provide, known limits |

## 5. Decision criteria

Per check-up completion criteria (reported as pass rate per check-up and per cell; reference check-up excluded from the rate):

* C1 (fit): fit RMSE on the averaged curve inside the fit window < 10 mV.
* C2 (capacity): |capacity loss implied by (LLI, LAM_PE, LAM_NE) - measured C/20 discharge capacity loss| < 1 percentage point, both relative to the reference check-up.
* OCP fitness: BOL (or reference) RMSE is the fitness test of the OCP pair for that cell; C1 failing at BOL means the OCP pair is not usable for that cell/dataset.

Identifiability (self-test and per check-up): a mode is called identifiable if its recovery error is below 1.5 pp (self-test, realistic noise)
and the Jacobian-based sd is below 1.5 pp; otherwise reported as "poorly identifiable" and not interpreted.

Hypotheses (compute, do not over-interpret; definitions fixed here):

* shares: sqrt share = 1 - lin share, lin share = B*x/(A*sqrt(x)+B*x) at the cell's end-of-test x (team A provides A_full, B_full;
  input path and column names are CLI arguments); LLI share = LLI/(LLI+LAM_PE+LAM_NE); LAM share = 1 - LLI share, both at the
  10 % measured-loss matched point (linear interpolation between check-ups; fallback: last check-up, flagged).
* H-D1: Spearman rho(sqrt share, LLI share) >= 0.5 and p < 0.05.
* H-D2a: Spearman rho(lin share, LAM share) >= 0.5 and p < 0.05.
* H-D2b: Spearman rho(irreversible expansion at the matched point, lin share) >= 0.5 and p < 0.05.
* Temperature confounds (R/C/H groups differ in both modes and curve shape): every test is also reported within each
  temperature group (n = 7) and pooled within group (group-demeaned ranks). A hypothesis counts as supported only if
  the criterion holds overall **and** the within-group pooled value keeps the same sign with rho >= 0.3 (soft guard).
  With n = 21 (7 per group) these are screening numbers, not tests.

## 6. Expected problems (known before running)

1. LAM_PE poorly identifiable when Li inventory (anode side) limits the capacity: the cathode capacity only shows in the
   slope of the high-voltage end; prototype saw LAM_PE rank correlation 0.2-0.5 across fit variants.
2. LLI and LAM_NE are collinear when capacity is LLI-limited.
3. OCP mismatch (NMC111 vs 532 label, fit to this cell family only) shows up as spurious LAM/LLI; BOL RMSE is the guard.
4. BOL first charge of a fresh cell includes SEI formation charge; the BOL charge branch is then not the same Li state as the
   discharge branch at the bottom (top alignment limits but does not remove this).
5. High-temperature cells: no BOL C/20 discharge, first averaged check-up at 47-294 Ah, charge branch starts at 3.17-3.3 V
   after a 3 h rest (common s range shorter, window lower edge matters).
6. Thickness sensor: baseline drop (cells 09, 21: -80...-100 um after BOL), step events (03, 21), frozen/NaN signal (13 last),
   all flagged in D2, never silently removed.
7. Graphite hysteresis and kinetic asymmetry between charge and discharge leave a residual that the average does not cancel.
8. n = 21 cells with one cell per condition: correlations carry no replicate uncertainty.
9. Python target is 3.14; the local interpreter is 3.13.5. Only numpy/pandas/scipy/scikit-learn/matplotlib features common to both are used; not verified on 3.14.

## 7. C1: reference curves (OCP source and chemistry claim)

* OCP source: PyBaMM `Mohtat2020` parameter set, functions `graphite_ocp_PeymanMPM` and `NMC_ocp_PeymanMPM`
  (file `.../lithium_ion/Mohtat2020.py`, last commit 43fc19759669a6b6f4a58bbc3107e496afa3da66 on repo `main`; the path given
  in the task, `src/pybamm/...develop`, now returns 404, the file lives at `packages/pybamm/src/pybamm/input/parameters/lithium_ion/Mohtat2020.py`).
  Re-implemented in `stepD_ocp_mohtat.py`; coefficients re-checked against the fetched file.
* Chemistry statement in that file: `get_parameter_values` docstring says "Parameters for a graphite/NMC532 pouch cell from the
  paper Mohtat2020 and references therein". The OCP docstrings cite only "Peyman MPM manuscript (to be submitted)".
  The earlier prototype notes the dataset paper's cell table as NMC111, so the label differs (532 in PyBaMM vs 111 in the paper);
  both are NMC with similar OCP shape, and the OCP is a fit to this cell family, so the mismatch matters mainly for absolute stoichiometries.
  Fitness is judged by BOL RMSE (criterion C1), not by the label.

## 8. Expected differences vs the prototype

Prototype primary fit = charge branch only with fitted R, per-cell DE global search at BOL, bottom-aligned Q axis, window 1-99.5 % of Q,
modes in sensitivity variants. Here: averaged branches, top alignment, fixed voltage window, multistart with convergence counting,
generic bounds. Expect similar LLI, differences in LAM_NE/LAM_PE of 1-3 pp (prototype's own method spread), different H-cell reference.

## 9. Change log of free choices and deviations (filled in after the runs)

Free choices tuned on data (the only ones): (1) fit window lower edge 3.30 V, upper 4.18 V. Reason: H-cell charge branches start at 3.17-3.3 V after a rest, so a lower edge below
about 3.3 V would change the covered range between check-ups; checked on cells 01 and 03 (window 3.0 / 3.3 / 3.45 V, upper 4.10 / 4.18 V): LLI changes < 1 pp, LAM_NE 1-3 pp, LAM_PE 1-2 pp; later reported for all cells as the
sensitivity `avg_window_3` in `D1_sensitivity.csv`. (2) rest window 400-580 s for the thickness (D2), chosen because the shortest top rest in the data is about 600 s; checked: difference to the end-of-rest value median 0.1 um.
Nothing was tuned toward C1/C2.

Deviations / additions after writing the spec:
* Files live in `Mothatnalysis` with `stepD_*` helper names (coordinator instruction); data via `analysis/../data`.
* C/20 branch selection uses the same current windows as `common.load_c20` (charge 100-400 mA, discharge -400..-100 mA, > 500 samples); capacity = cycler `Capacity [Ah]` of the discharge, so `ah` and `soh` match the team's `out/c20.csv`
  (soh = capacity / first check-up with both branches; max |soh diff| 6.7e-5, same 241 check-ups). `loss_meas` for criterion C2 is relative to the mode reference instead.
* Sensitivity `avg_ref_after_bol` added after seeing that criterion C2 fails mostly through a constant +1.2 pp offset created at the BOL reference of fresh cells; both rates are reported, the BOL-reference rate is the pre-specified one.
* Top-rest for D2 is only accepted if the preceding charge reached the top (V > 4.15 V); otherwise `flag_no_top_state`.
* H-D1 and H-D2a are the same rank statistic (complementary shares) and are reported as a mirrored pair; secondary variants with LAM_NE-only and LAM_PE-only shares added.
* Self-test scenarios `ocp_mismatch_x3` and `ocp_mismatch_x6` added after the first self-test showed no identifiability problem under measurement nuisance alone.
