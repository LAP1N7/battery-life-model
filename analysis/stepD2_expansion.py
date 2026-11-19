# -*- coding: utf-8 -*-
"""
stepD2_expansion.py - irreversible thickness change at the SAME SOC and the SAME rest condition per check-up.

Rule (fixed in D0_spec.md):
  SOC state   : top of charge (cell charged to 4.2 V: C/20 CC, plus CV tail if the protocol has one) = 100 % SOC by the voltage limit.
  rest rule   : thickness = median over rest time 400-580 s after the current dropped below 20 mA (the shortest rest in the data is ~600 s;
                longer rests are cut to the same window, so every check-up is evaluated at the same rest time).
  irreversible: dX_irr = X_top(k) - X_top(reference), reference = first check-up (BOL). Because the sensor baseline can jump after BOL,
                dX_irr_r1 uses the first check-up after BOL as the reference, and dX_irr_clean additionally removes step events.
  The bottom state (start of the C/20 charge) is reported as a diagnostic only: its rest history differs between check-ups (0 s to 3 h).

Artifact flags (never silently removed; clean columns are additional): sensor_nan (thickness missing > 20 % in the charge), frozen (reversible swing < 5 um),
baseline_drop (first increment < -20 um), step_event (increment > 80 um and > 4x the median positive increment of the cell).

    python -I -B -X utf8 stepD2_expansion.py [--dataset mohtat] [--cells ...]
"""
import os
import sys
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import pandas as pd

import stepD_loaders as dl

OUT_DEFAULT = os.path.join(HERE, "out", "D")
T_LO, T_HI = 400.0, 580.0          # rest-time window [s]
SWING_MIN = 5.0                    # um
STEP_ABS, STEP_REL = 80.0, 4.0     # step-event thresholds
DROP_MIN = -20.0
NAN_MAX = 0.20


def per_checkup(rest, curves):
    rows = []
    for (cell, k), g in rest.groupby(["cell", "checkup_idx"], sort=True):
        g = g.sort_values("t_rest_s")
        t = g["t_rest_s"].values
        x = g["thickness_um"].values
        m = (t >= T_LO) & (t <= T_HI) & np.isfinite(x)
        rec = dict(cell=cell, checkup_idx=int(k), ah=float(g["ah_throughput"].iloc[0]), rest_len_s=float(t.max()))
        rec["X_top_um"] = float(np.median(x[m])) if m.sum() >= 5 else np.nan
        e = (t >= t.max() - 60) & np.isfinite(x)
        rec["X_top_restend_um"] = float(np.median(x[e])) if e.sum() >= 3 else np.nan
        rec["flag_short_rest"] = bool(t.max() < T_HI)
        rows.append(rec)
    D = pd.DataFrame(rows)
    base = curves.groupby(["cell", "checkup_idx"], as_index=False)["ah_throughput"].max().rename(columns={"ah_throughput": "ah"})
    D = base.merge(D.drop(columns="ah"), on=["cell", "checkup_idx"], how="left")
    # bottom state and charge-branch NaN fraction from the charge branch
    ch = curves[curves["direction"] == "charge"]
    b = []
    for (cell, k), g in ch.groupby(["cell", "checkup_idx"], sort=True):
        x = g["thickness_um"].values
        fin = np.isfinite(x)
        xb = float(np.median(x[:20][np.isfinite(x[:20])])) if np.isfinite(x[:20]).any() else np.nan
        b.append(dict(cell=cell, checkup_idx=int(k), X_bottom_um=xb, nan_frac_charge=float(1 - fin.mean())))
    return D.merge(pd.DataFrame(b), on=["cell", "checkup_idx"], how="outer")


def flag_and_reference(D):
    out = []
    for cell, g in D.groupby("cell", sort=False):
        g = g.sort_values("checkup_idx").reset_index(drop=True).copy()
        g["X_swing_um"] = g["X_top_um"] - g["X_bottom_um"]
        g["flag_sensor_nan"] = g["nan_frac_charge"] > NAN_MAX
        g["flag_no_top_state"] = g["X_top_um"].isna() & ~g["flag_sensor_nan"]      # charge did not reach the top / no rest (not a sensor fault)
        g["flag_frozen"] = (g["X_swing_um"] < SWING_MIN) & ~g["flag_no_top_state"]
        invalid = g["flag_sensor_nan"] | g["flag_frozen"] | g["flag_no_top_state"]
        X = g["X_top_um"].where(~invalid)
        g["X_valid_um"] = X
        g["dX_irr_um"] = X - X.iloc[0]
        inc = X.diff()
        pos = inc[inc > 0]
        med = float(pos.median()) if len(pos) else 0.0
        g["flag_baseline_drop"] = False
        if len(g) > 1 and np.isfinite(inc.iloc[1]) and inc.iloc[1] < DROP_MIN:
            g.loc[1:, "flag_baseline_drop"] = True        # all later check-ups sit on the shifted baseline
        g["flag_step_event"] = (inc > STEP_ABS) & (inc > STEP_REL * max(med, 1e-9))
        g["flag_step_event"] = g["flag_step_event"].fillna(False)
        ref1 = X.iloc[1] if len(g) > 1 else np.nan
        g["dX_irr_r1_um"] = X - ref1
        # clean: drop step-event increments (replace by median increment) and start from the first post-BOL check-up when the baseline dropped
        incc = inc.copy()
        incc[g["flag_step_event"]] = med
        base = X.iloc[1] if g["flag_baseline_drop"].any() else X.iloc[0]
        clean = incc.fillna(0.0)
        if g["flag_baseline_drop"].any():
            clean.iloc[1] = 0.0
        clean = clean.cumsum()
        clean = clean - (clean.iloc[1] if g["flag_baseline_drop"].any() else 0.0)
        g["dX_irr_clean_um"] = clean.where(~invalid)
        g["flag_artifact_cell"] = bool(g["flag_sensor_nan"].any() or g["flag_frozen"].any() or g["flag_baseline_drop"].any() or g["flag_step_event"].any())
        g["artifact_kinds"] = ",".join(k for k, f in (("sensor_nan_or_frozen", (g["flag_sensor_nan"] | g["flag_frozen"]).any()), ("baseline_drop", g["flag_baseline_drop"].any()),
                                                    ("step_event", g["flag_step_event"].any())) if f)
        out.append(g)
    return pd.concat(out, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default="mohtat")
    ap.add_argument("--loader-spec", default=None)
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--cells", default=None)
    ap.add_argument("--out-dir", default=OUT_DEFAULT)
    a = ap.parse_args()
    ds = dl.load_dataset(a.loader_spec or a.dataset, a.data_dir, a.cells.split(",") if a.cells else None)
    if ds["rest"] is None or not len(ds["rest"]) or "thickness_um" not in ds["curves"].columns:
        raise SystemExit("dataset provides no thickness / top-rest table; D2 cannot run")
    D = flag_and_reference(per_checkup(ds["rest"], ds["curves"]))
    meta = ds["meta"]
    if "T_group" in meta.columns:
        D = D.merge(meta[["cell", "T_group"]], on="cell", how="left")
    d1 = os.path.join(a.out_dir, "D1_dva.csv")
    if os.path.exists(d1):
        S = pd.read_csv(d1, dtype={"cell": str})[["cell", "ah", "soh", "LLI", "LAM_PE", "LAM_NE"]]
        S["ahk"] = S["ah"].round(2)
        D["ahk"] = D["ah"].round(2)
        D = D.merge(S.drop(columns="ah"), on=["cell", "ahk"], how="left").drop(columns="ahk")
    D["soc_rule"] = "top of charge (4.2 V CC/CV), 100 % by voltage limit"
    D["rest_rule"] = f"median of rest time {T_LO:.0f}-{T_HI:.0f} s after I<20 mA"
    cols = ["cell", "T_group", "checkup_idx", "ah", "soh", "X_top_um", "dX_irr_um", "dX_irr_r1_um", "dX_irr_clean_um", "X_bottom_um", "X_swing_um",
            "rest_len_s", "X_top_restend_um", "nan_frac_charge", "flag_short_rest", "flag_no_top_state", "flag_sensor_nan", "flag_frozen", "flag_baseline_drop", "flag_step_event",
            "flag_artifact_cell", "artifact_kinds", "LLI", "LAM_PE", "LAM_NE", "soc_rule", "rest_rule"]
    D = D[[c for c in cols if c in D.columns]]
    os.makedirs(a.out_dir, exist_ok=True)
    D.to_csv(os.path.join(a.out_dir, "D2_expansion.csv"), index=False, float_format="%.5g")
    # ---- diagnostics ----
    L = ["# D2 summary\n", f"rule: top of charge, thickness median at rest time {T_LO:.0f}-{T_HI:.0f} s; reference = first check-up (BOL).\n"]
    long_rest = D[D["rest_len_s"] > 3000]
    diff = (long_rest["X_top_um"] - long_rest["X_top_restend_um"]).abs()
    L.append(f"- rest-rule sensitivity: {len(long_rest)} check-ups have a ~1 h rest; |X(400-580 s) - X(end of rest)|: median {diff.median():.2f} um, 95th pct {diff.quantile(.95):.2f}, max {diff.max():.2f} um "
             f"(irreversible changes are 20-500 um, so the rest rule is not a limiting factor).")
    L.append(f"- short rests (< {T_HI:.0f} s): {int(D['flag_short_rest'].sum())} check-ups (they are excluded from X_top and show NaN).")
    art = D.groupby("cell")[["flag_sensor_nan", "flag_frozen", "flag_baseline_drop", "flag_step_event"]].any()
    art["kinds"] = D.groupby("cell")["artifact_kinds"].first()
    art = art[art["kinds"] != ""]
    L.append("- artifact cells flagged automatically: " + (", ".join(f"{c} ({r['kinds']})" for c, r in art.iterrows()) or "none") +
             ". Previously known: 03, 09, 13, 21.")
    ev = D[D["flag_step_event"]][["cell", "checkup_idx", "ah"]]
    for _, r in ev.iterrows():
        L.append(f"  - step event: cell {r['cell']} at check-up {int(r['checkup_idx'])} ({r['ah']:.0f} Ah)")
    # compare with the team's expansion.csv (x_irr based on bottom state)
    f = os.path.join(HERE, "out", "expansion.csv")
    if a.dataset == "mohtat" and os.path.exists(f):
        T = pd.read_csv(f)
        T["cell"] = T["cell"].map(lambda v: f"{int(v):02d}")
        T["ahk"] = T["ah"].round(2)
        M = D.assign(ahk=D["ah"].round(2)).merge(T[["cell", "ahk", "x_irr", "x_dis"]], on=["cell", "ahk"], how="inner")
        M = M[np.isfinite(M["dX_irr_um"]) & np.isfinite(M["x_irr"]) & (M["ah"] > 0)]
        dd = M["dX_irr_um"] - M["x_irr"]
        L.append(f"- vs team out/expansion.csv (x_irr = bottom-state thickness, first 20 charge samples, no rest rule): matched {len(M)} check-ups after BOL, "
                 f"Spearman(top-state dX_irr, x_irr) {M['dX_irr_um'].corr(M['x_irr'], method='spearman'):.3f}; difference median {dd.median():+.1f} um, 95th pct |diff| {dd.abs().quantile(.95):.1f} um. "
                 f"My X_bottom (first 20 samples of the C/20 charge) differs from team x_dis by more than 1 um in {int(((M['X_bottom_um'] - M['x_dis']).abs() > 1).sum())} check-ups "
                 f"(cells {', '.join(sorted(M.loc[(M['X_bottom_um'] - M['x_dis']).abs() > 1, 'cell'].unique()))}; median offset {((M['x_dis'] - M['X_bottom_um']).abs())[(M['X_bottom_um'] - M['x_dis']).abs() > 1].median():.0f} um): the team's selection 100<I<400 mA picks up the CV-taper rows of a preceding "
                 f"1C charge at the top of charge, so for these check-ups x_dis is a TOP-state thickness while for the others it is a bottom-state thickness (mixed states in x_irr).")
    last = D.dropna(subset=["dX_irr_um"]).sort_values("ah").groupby("cell").tail(1)
    if "T_group" in D.columns:
        L.append("\n- last VALID check-up dX_irr_um by temperature group (median [min, max]):")
        for g, d in last.groupby("T_group"):
            v = d["dX_irr_um"].dropna()
            L.append(f"  - {g}: {v.median():.0f} [{v.min():.0f}, {v.max():.0f}] um (n={len(v)})")
    with open(os.path.join(a.out_dir, "D2_summary.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
