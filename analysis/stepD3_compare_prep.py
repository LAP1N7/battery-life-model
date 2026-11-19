# -*- coding: utf-8 -*-
"""
stepD3_compare_prep.py - join team A's sqrt(t) / linear capacity-fade components with track C results (D1 modes, D2 irreversible expansion)
and compute (not over-interpret) the Spearman tests of H-D1 / H-D2.

A's definition (stepP1_prototype.py): Q_loss(t) = A*sqrt(t) + B*t, t = cycling hours (out/checkup_time.csv), loss = 1 - cap/cap_first, A_full/B_full from a least-squares
fit over check-ups with SOH >= 0.80, columns A_full, B_full in out/P1_pred.csv (rounded to 4 decimals there).
Shares are evaluated at t_ref = hours of the last check-up with SOH >= 0.80 (the end of the range A fitted):
    lin_share = B*t_ref / (A*sqrt(t_ref) + B*t_ref),   sqrt_share = 1 - lin_share.
Because P1_pred.csv is rounded to 1e-4 (B is 0.0000-0.0002 for several cells) the same fit is also recomputed here, unrounded, with A's recipe
(columns *_re); both versions are tested.

Mode shares from D1 at the matched 10 % measured capacity loss (linear interpolation):  LLI_share = LLI/(LLI+LAM_PE+LAM_NE), LAM_share = 1 - LLI_share.
NOTE: H-D1 (sqrt share vs LLI share) and H-D2a (linear share vs LAM share) are the SAME rank statistic by construction (both shares are complements),
so they are reported once as a pair of mirrored tests; secondary variants use LAM_NE share and LAM_PE share separately.

    python -I -B -X utf8 stepD3_compare_prep.py --p1 out\\P1_pred.csv [--colA A_full --colB B_full] [--level 0.10]
    stand-in test (not A's file):  --p1 <path>\\capacity_fit_AB.csv --colA A --colB B --label stand-in
"""
import os
import sys
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import pandas as pd
from scipy import stats

import stepD_pipeline as dp

OUT = os.path.join(HERE, "out")
OUT_D = os.path.join(OUT, "D")
TEMP = {c: "RCH"[(c - 1) % 3] for c in range(1, 22)}     # same mapping as common.py


def recompute_AB(cell, c20, hours, soh_min=0.80):
    g = c20[c20["cell"] == cell].merge(hours[hours["cell"] == cell][["cell", "ah", "hours"]], on=["cell", "ah"], how="left").dropna(subset=["hours"]).reset_index(drop=True)
    if len(g) < 3:
        return np.nan, np.nan, np.nan
    g["soh"] = g["cap"] / g["cap"].iloc[0]
    gg = g[g["soh"] >= soh_min]
    if len(gg) < 3:
        return np.nan, np.nan, np.nan
    X = np.c_[np.sqrt(gg["hours"]), gg["hours"]]
    ab, *_ = np.linalg.lstsq(X, 1 - gg["soh"], rcond=None)
    return float(ab[0]), float(ab[1]), float(gg["hours"].iloc[-1])


def lin_share(A, B, t):
    a, b = A * np.sqrt(t), B * t
    tot = a + b
    return np.where(np.abs(tot) > 1e-12, b / tot, np.nan)


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def pooled_within(x, y, grp):
    """Spearman-like correlation after removing group mean of ranks (temperature control). p approximate (dof = n - groups - 1)."""
    d = pd.DataFrame(dict(x=x, y=y, g=grp)).dropna()
    if len(d) < 6 or d["g"].nunique() < 2:
        return np.nan, np.nan, len(d)
    rx, ry = d["x"].rank(), d["y"].rank()
    rx = rx - rx.groupby(d["g"]).transform("mean")
    ry = ry - ry.groupby(d["g"]).transform("mean")
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = len(d) - d["g"].nunique() - 1
    t = r * np.sqrt(dof / max(1e-12, 1 - r * r))
    return r, float(2 * stats.t.sf(abs(t), dof)), len(d)


def run_tests(J, xcol, ycol, name, variant):
    rows = []
    r, p, n = spearman(J[xcol].values, J[ycol].values)
    rows.append(dict(test=name, variant=variant, subset="all", n=n, rho=r, p=p))
    for g, d in J.groupby("T_group"):
        r_, p_, n_ = spearman(d[xcol].values, d[ycol].values)
        rows.append(dict(test=name, variant=variant, subset=f"within {g}", n=n_, rho=r_, p=p_))
    r_, p_, n_ = pooled_within(J[xcol].values, J[ycol].values, J["T_group"].values)
    rows.append(dict(test=name, variant=variant, subset="pooled within T (rank-demeaned)", n=n_, rho=r_, p=p_))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--p1", default=os.path.join(OUT, "P1_pred.csv"), help="team A table with per-cell A and B (CLI argument, default out/P1_pred.csv)")
    ap.add_argument("--colA", default="A_full")
    ap.add_argument("--colB", default="B_full")
    ap.add_argument("--col-share", default=None, help="use this column of --p1 as the linear share directly (stand-in file in EFC units); then only the 'file' variant is tested")
    ap.add_argument("--label", default="P1_pred (team A)")
    ap.add_argument("--d1", default=os.path.join(OUT_D, "D1_dva.csv"))
    ap.add_argument("--d2", default=os.path.join(OUT_D, "D2_expansion.csv"))
    ap.add_argument("--c20", default=os.path.join(OUT, "c20.csv"))
    ap.add_argument("--hours", default=os.path.join(OUT, "checkup_time.csv"))
    ap.add_argument("--level", type=float, default=0.10, help="matched measured capacity loss for mode shares / expansion")
    ap.add_argument("--out-dir", default=OUT_D)
    ap.add_argument("--tag", default="", help="suffix for output file names (use for the stand-in test)")
    a = ap.parse_args()

    P = pd.read_csv(a.p1)
    P["cell"] = P["cell"].map(lambda v: f"{int(v):02d}")
    c20 = pd.read_csv(a.c20)
    hours = pd.read_csv(a.hours)
    D1 = pd.read_csv(a.d1, dtype={"cell": str})
    D2 = pd.read_csv(a.d2, dtype={"cell": str})
    rows = []
    for cell in sorted(D1["cell"].unique()):
        ci = int(cell)
        pr = P[P["cell"] == cell]
        sh_direct = float(pr[a.col_share].iloc[0]) if (a.col_share and len(pr)) else None
        A_p = float(pr[a.colA].iloc[0]) if len(pr) else np.nan
        B_p = float(pr[a.colB].iloc[0]) if len(pr) else np.nan
        A_r, B_r, t_ref = recompute_AB(ci, c20, hours)
        g = D1[(D1["cell"] == cell) & D1["loss_meas"].notna()].sort_values("ah")
        d2 = D2[D2["cell"] == cell].copy()
        d2["ahk"] = d2["ah"].round(2)
        g = g.assign(ahk=g["ah"].round(2)).merge(d2[["ahk", "dX_irr_clean_um", "dX_irr_r1_um", "dX_irr_um", "flag_artifact_cell"]], on="ahk", how="left")
        vals, how = {}, None
        for k in ("LLI", "LAM_PE", "LAM_NE", "dX_irr_clean_um", "dX_irr_r1_um", "dX_irr_um"):
            vals[k], h = dp.matched_value(g, k, a.level)
            if k == "LLI":
                how = h
        tot = vals["LLI"] + vals["LAM_PE"] + vals["LAM_NE"]
        meta_g = D1.loc[D1["cell"] == cell, "T_group"].iloc[0] if "T_group" in D1.columns else TEMP[ci]
        rows.append(dict(
            cell=cell, T_group=meta_g, in_A_file=bool(len(pr)), A=A_p, B=B_p, A_re=A_r, B_re=B_r, t_ref_h=t_ref,
            lin_share=(sh_direct if sh_direct is not None else (float(lin_share(A_p, B_p, t_ref)) if np.isfinite(A_p) and np.isfinite(t_ref) else np.nan)),
            lin_share_re=float(lin_share(A_r, B_r, t_ref)) if np.isfinite(A_r) else np.nan,
            neg_component=bool((np.isfinite(A_r) and (A_r < 0 or B_r < 0))),
            LLI_at=100 * vals["LLI"], LAM_PE_at=100 * vals["LAM_PE"], LAM_NE_at=100 * vals["LAM_NE"], matched_how=how,
            LLI_share=vals["LLI"] / tot, LAM_share=(vals["LAM_PE"] + vals["LAM_NE"]) / tot,
            LAM_NE_share=vals["LAM_NE"] / tot, LAM_PE_share=vals["LAM_PE"] / tot,
            dX_irr_clean_um=vals["dX_irr_clean_um"], dX_irr_r1_um=vals["dX_irr_r1_um"], dX_irr_um=vals["dX_irr_um"],
            flag_ref_not_bol=bool(g["flag_ref_not_bol"].iloc[0]) if len(g) else np.nan,
            flag_expansion_artifact=bool(g["flag_artifact_cell"].eq(True).any()) if len(g) else np.nan,
        ))
    J = pd.DataFrame(rows)
    J["sqrt_share"] = 1 - J["lin_share"]
    J["sqrt_share_re"] = 1 - J["lin_share_re"]
    os.makedirs(a.out_dir, exist_ok=True)
    J.to_csv(os.path.join(a.out_dir, f"D3_join{a.tag}.csv"), index=False, float_format="%.5g")

    res = []
    for variant, sfx, sub in (("A,B as in file", "", J[J["in_A_file"]]), ("A,B recomputed unrounded (same cells as file)", "_re", J[J["in_A_file"]]),
                              ("A,B recomputed unrounded, all 21 cells", "_re", J))[: (1 if a.col_share else 3)]:
        S = sub.copy()
        S["lin"] = S["lin_share" + sfx]
        S["sq"] = 1 - S["lin"]
        res += run_tests(S, "sq", "LLI_share", "H-D1 sqrt share ~ LLI share (== H-D2a mirrored)", variant)
        res += run_tests(S, "lin", "LAM_share", "H-D2a lin share ~ LAM share", variant)
        res += run_tests(S, "lin", "LAM_NE_share", "secondary: lin share ~ LAM_NE share", variant)
        res += run_tests(S, "lin", "LAM_PE_share", "secondary: lin share ~ LAM_PE share", variant)
        res += run_tests(S, "dX_irr_clean_um", "lin", "H-D2b irreversible expansion (clean) ~ lin share", variant)
        res += run_tests(S, "dX_irr_r1_um", "lin", "secondary: expansion dX_irr_r1 ~ lin share", variant)
    R = pd.DataFrame(res)
    R.to_csv(os.path.join(a.out_dir, f"D3_correlations{a.tag}.csv"), index=False, float_format="%.4g")

    # note
    L = [f"# D3 note (input: {a.label}: `{os.path.relpath(a.p1, HERE)}`, columns {a.colA}, {a.colB})\n",
         "Computed, not interpreted. Criterion proposed in D0_spec: rho >= 0.5 and p < 0.05 overall, and pooled-within-temperature rho >= 0.3 with the same sign. n is small (<= 21; 7 per temperature group).\n",
         f"Cells in A's file: {int(J['in_A_file'].sum())}/21 (missing: {', '.join(J.loc[~J['in_A_file'], 'cell'])}). "
         f"Shares at t_ref = last check-up hours with SOH >= 0.80; mode shares/expansion at {100*a.level:.0f} % measured loss ('matched_how' column; cells that never reach it use their last check-up).",
         f"A,B rounding: |A_file - A_recomputed| max {np.nanmax(np.abs(J['A'] - J['A_re'])):.1e}, |B_file - B_recomputed| max {np.nanmax(np.abs(J['B'] - J['B_re'])):.1e}; "
         f"lin_share file vs recomputed Spearman {J['lin_share'].corr(J['lin_share_re'], method='spearman'):.2f} (cells with a negative A or B in the unrounded fit: {int(J['neg_component'].sum())}).",
         f"Reference caveat: modes of H cells are relative to the first averaged check-up (47-294 Ah), not BOL ({int(J['flag_ref_not_bol'].fillna(False).sum())} cells flagged); "
         f"expansion artifact cells (flagged in D2): {', '.join(J.loc[J['flag_expansion_artifact'].fillna(False), 'cell'])}.\n"]
    for variant, d in R.groupby("variant", sort=False):
        L.append(f"## {variant}\n")
        L.append("| test | subset | n | rho | p |\n|---|---|---|---|---|")
        for _, r in d.iterrows():
            L.append(f"| {r['test']} | {r['subset']} | {r['n']} | {r['rho']:.2f} | {r['p']:.3f} |")
        L.append("")
    with open(os.path.join(a.out_dir, f"D3_note{a.tag}.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
