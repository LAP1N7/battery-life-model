# D_portability: running the mechanism tool (D1/D2) on another dataset

Files (all prefixed `stepD`, owned by track C): `stepD_core.py` (dataset-agnostic: OCP registry, curve averaging, bounded multistart fit,
modes vs reference, capacity-from-modes, identifiability), `stepD_pipeline.py` (fit all cells, modes, pass flags), `stepD_loaders.py` (loader registry +
validation), `stepD_loader_mohtat.py` and `stepD_ocp_mohtat.py` (the only Mohtat-specific code), `stepD1_dva.py`, `stepD2_expansion.py`,
`stepD3_compare_prep.py`, `stepD0_selftest.py`.

## 1. What a dataset must provide

1. **A quasi-OCV curve pair per check-up**: one CC charge branch and one CC discharge branch at about C/20 (C/10 is the practical upper limit;
   the averaging cancels the symmetric I*R and hysteresis, not kinetic limits at higher rates). Both branches should span (almost) the full voltage window and
   the charge branch should end at the upper limit (include the CV tail if the protocol has one, so that the discharge starts from the same state).
   If only one branch exists the pipeline falls back to a single-branch fit with a fitted constant polarisation offset (flagged `mode_used = chg|dis`);
   results are then not interchangeable with averaged ones.
2. **A fixed voltage window** (`vmin`, `vmax`, e.g. 3.0 / 4.2 V). The fit uses a fixed sub-window on the averaged curve (`--fit-vlo/--fit-vhi`, default 3.30-4.18 V),
   identical for all check-ups and cells; choose it so that every check-up covers it (the lower edge is limited by where the shorter branch starts).
3. **Chemistry-matched half-cell OCPs** `Up(y)` and `Un(x)` registered as `OCPPair(name, Up, Un, x_range, y_range, source, chemistry)`
   (see `stepD_ocp_mohtat.py`; external file via `--ocp-module my_ocp.py --ocp name`). The BOL RMSE on the reference check-up is the fitness test:
   if it is above 10 mV the OCP pair (or the data) is not suitable, and no mode value should be interpreted.
4. **A reference check-up per cell** (BOL, or the first check-up with an averaged curve; flagged when not BOL).
5. Optional for D2: thickness/strain signal and a rest at the top state (table `rest`, columns cell, checkup_idx, ah_throughput, t_rest_s, thickness_um).

## 2. Standard long table (adapter output)

`load(data_dir, cells) -> dict(curves, rest, meta, info)`

| table | columns |
|---|---|
| curves | cell, checkup_idx, ah_throughput, direction (`charge` or `discharge`), Q_Ah (Ah since the start of that branch, increasing), V, T; optional thickness_um, cap_cycler_Ah |
| rest (optional) | cell, checkup_idx, ah_throughput, phase (`top_rest`), t_rest_s (0 = start of the rest), V, thickness_um, T |
| meta | cell plus free columns; recommended `T_group`, `room_mate` (same-batch room-temperature cell used for the sensitivity reference) |
| info | vmin, vmax, ocp (registry name), c_nom_ah, dataset |

Use: `python -I -B -X utf8 stepD1_dva.py --dataset mohtat` or, for another dataset,
`python -I -B -X utf8 stepD1_dva.py --dataset mylab --loader-spec path\to\my_loader.py:load --data-dir D:\x --ocp my_pair --ocp-module path\to\my_ocp.py`.
`--loader-spec` / `--ocp-module` execute the given Python file: use only files you wrote or reviewed, never files inside a downloaded data folder.
Run with `python -I -B` (isolated; the step scripts add their own folder to `sys.path`).

Procedure to bring up a new dataset: (1) write the loader, run `stepD0_selftest.py --ocp my_pair --ocp-module ...` to see what is identifiable with the new OCP
shapes; (2) run D1 on one cell and look at the fit figure and the BOL RMSE; (3) only then run all cells.

## 3. Known limits (what makes results unreliable)

* **Blended / silicon-containing anodes (e.g. NMC / graphite+SiO such as the KIT cells)**: a graphite-only `Un` cannot describe the curve; a blended anode OCP is needed
  (graphite + Si/SiO weighted by a fitted fraction), and Si/SiO shows large voltage hysteresis and OCP drift with ageing and cycling history.
  The charge/discharge average does not remove this (it is not symmetric about the OCV), so it appears as spurious LLI/LAM_NE. Treat results as invalid unless a blended OCP is
  registered and the BOL RMSE passes.
* **No slow check-up curve** (only C/2 or faster, or dynamic cycles): not usable; polarisation then dominates the shape. Pulses or EIS alone do not give mode information here.
* **OCP mismatch** (different NMC grade, e.g. 111 vs 532 vs 811; different graphite): shifts stoichiometries and leaks into LAM_PE/LAM_NE. The Mohtat OCP is a fit to that cell family itself, so Mohtat results
  are optimistic compared with a dataset that uses a literature OCP.
* **LAM_PE is poorly identifiable when Li inventory limits the capacity** and LLI/LAM_NE are collinear (see `D0_selftest.md` for numbers under realistic nuisance).
* **Reference**: LLI/LAM are relative to the reference check-up. A fresh-cell BOL curve contains formation effects (in Mohtat the BOL reference gives a systematic +1.2 pp excess
  of implied vs measured capacity loss); a missing BOL (high-temperature Mohtat cells) moves the zero point. Cross-cell comparison of absolute mode values needs the same reference convention.
* **Window coverage**: if the charge branch starts after a long rest (Mohtat H cells start at 3.17-3.3 V) the common range of the average shrinks; keep `fit_vlo` above that.
* **Temperature of the check-up**: the OCP pair has no temperature dependence; check-ups should be at one temperature (Mohtat: all at about 25 C).
* **Thickness (D2)**: needs a top-state rest of at least 580 s (configurable constants in `stepD2_expansion.py`); a pouch-cell stack sensor has baseline jumps, so flags are produced rather than silent corrections.
* Python 3.14 target; developed on 3.13.5 with numpy 2.3, scipy 1.16, pandas 2.3, matplotlib 3.10 (3.14 not available locally, not tested).
