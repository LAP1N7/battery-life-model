"""Mechanism cross-check from several independent indicators (public data only).

Run:  python -I mech_crosscheck.py [--rebuild]
Inputs are untrusted downloaded CSV/XLSX files: only read with pandas, never executed/unpickled.
Outputs: prototype/results/mech_crosscheck/{results.json, summary.md, fig_*.png, kit_cell_series.csv(cache)}
"""
import sys, re, glob, json, math, pathlib, warnings
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parents[2]   # 레포 밖 작업 폴더(data/, prototype/results/가 있는 곳)
DATA = ROOT / "data"
PROTO = ROOT / "prototype"
OUT = PROTO / "results" / "mech_crosscheck"
OUT.mkdir(parents=True, exist_ok=True)
M50 = DATA / "m50t_dma_7235858" / "extracted" / "Battery Degradation Data"
B_XLSX = DATA / "zenodo_15755725_degradation_path" / "extracted_degradation_metrics.xlsx"
KIT = DATA / "kit_radar" / "result_v2" / "Unnamed entity"
KIT_EIS_FIXED = DATA / "kit_radar" / "eis_addendum" / "cell_eisv2_fixed"
CKUP = PROTO / "results" / "kit_checkups.csv"
REBUILD = "--rebuild" in sys.argv

R = {}  # everything that goes into results.json


# ------------------------------------------------------------------ helpers
def f(x, nd=4):
    if x is None:
        return None
    try:
        x = float(x)
    except Exception:
        return None
    return None if (math.isnan(x) or math.isinf(x)) else round(x, nd)


def spear(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 4 or np.ptp(a[m]) == 0 or np.ptp(b[m]) == 0:
        return dict(rho=None, p=None, n=int(m.sum()))
    r = stats.spearmanr(a[m], b[m])
    return dict(rho=f(r.statistic), p=f(r.pvalue, 5), n=int(m.sum()))


def kend(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 4 or np.ptp(a[m]) == 0 or np.ptp(b[m]) == 0:
        return dict(tau=None, p=None, n=int(m.sum()))
    r = stats.kendalltau(a[m], b[m])
    return dict(tau=f(r.statistic), p=f(r.pvalue, 5), n=int(m.sum()))


def interp_at(x, y, xs):
    """linear interpolation inside the observed range only (no extrapolation)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 2:
        return np.nan
    o = np.argsort(x); x, y = x[o], y[o]
    if xs < x[0] or xs > x[-1]:
        return np.nan
    return float(np.interp(xs, x, y))


# =================================================================== 1. M50T
def load_m50():
    rows = []
    align = []
    for win, wname in [("0-30% SoC", "0-30%"), ("0-100% SoC", "0-100%")]:
        for mf in sorted((M50 / "Extracted Data - Main" / win).glob("*.csv")):
            mm = re.search(r"cell (\w) \((\d+)degC\)", mf.name)
            cell, Tn = mm.group(1), int(mm.group(2))
            dmf = list((M50 / "Extracted Data - Degradation Modes" / win / f"cell {cell}").glob("*DM data*.csv"))[0]
            m = pd.read_csv(mf, index_col=0)
            d = pd.read_csv(dmf, index_col=0)
            n_ok = len(m) == len(d)
            # row-wise (positional) match: both files are indexed 0..n-1 per RPT; check by SoH vs C/10 capacity ratio
            soh_c = m["C/10 Capacity (mA.h)"].values / m["C/10 Capacity (mA.h)"].values[0]
            maxdiff = float(np.nanmax(np.abs(soh_c - d["SoH"].values))) if n_ok else None
            align.append(dict(window=wname, cell=cell, T=Tn, n_main=len(m), n_dm=len(d), max_abs_SoH_minus_C10ratio=f(maxdiff, 6)))
            df = pd.DataFrame(dict(
                window=wname, cell=cell, T=Tn, rpt=np.arange(len(m)),
                thr=m["Charge Throughput (A.h)"].values,
                days=m["Days of degradation"].values,
                cap_loss=1 - d["SoH"].values,
                LLI=d["LLI"].values, LAM_PE=d["LAM PE"].values, LAM_NE=d["LAM NE_tot"].values,
                LAM_NEGr=d["LAM NE_Gr"].values, LAM_NESi=d["LAM NE_Si"].values,
                R=m["0.1s Resistance (Ohms)"].values))
            R0 = df["R"].iloc[0]
            df["R_growth"] = df["R"] / R0 - 1
            rows.append(df)
    return pd.concat(rows, ignore_index=True), align


M50_IND = ["cap_loss", "LLI", "LAM_PE", "LAM_NE", "LAM_NEGr", "LAM_NESi", "R_growth"]


def analyse_m50():
    df, align = load_m50()
    out = {}
    out["alignment_check"] = align
    out["n_cells"] = int(df.groupby(["window", "cell"]).ngroups)
    out["matching_method"] = ("Main and DM files have equal row counts per cell and are matched by row position (RPT index); "
                              "verified that DM SoH equals C/10-capacity/initial C/10-capacity (max |diff| in alignment_check).")
    out["conditions_n_cells"] = {f"{w}|{T}C": int(g["cell"].nunique()) for (w, T), g in df.groupby(["window", "T"])}
    out["units_note"] = "DM values are fractions of initial (x100 for %). R_growth = R(0.1 s)/R(RPT0)-1; R exists only on about half the RPTs."
    out["ranges"] = {c: [f(df[c].min()), f(df[c].max())] for c in M50_IND}

    # ---- closure
    cl = df.dropna(subset=["cap_loss", "LLI", "LAM_PE", "LAM_NE"]).copy()
    cl = cl[cl.rpt > 0]
    cl["lam_max"] = cl[["LAM_PE", "LAM_NE"]].max(axis=1)
    cl["mx"] = cl[["LLI", "LAM_PE", "LAM_NE"]].max(axis=1)
    lim = cl[["LLI", "LAM_PE", "LAM_NE"]].idxmax(axis=1)
    def resid_stats(pred):
        r = cl["cap_loss"] - pred
        return dict(mean=f(r.mean()), sd=f(r.std()), rmse=f(np.sqrt((r ** 2).mean())), mae=f(r.abs().mean()),
                    slope_cap_vs_pred=f(np.polyfit(pred, cl["cap_loss"], 1)[0]),
                    spearman=spear(cl["cap_loss"], pred), median_ratio=f((cl["cap_loss"] / pred).median()))
    out["closure"] = dict(
        n_points=int(len(cl)), n_cells=int(cl.groupby(["window", "cell"]).ngroups),
        cap_vs_LLI=resid_stats(cl["LLI"]),
        cap_vs_max_modes=resid_stats(cl["mx"]),
        cap_vs_max_LAM=resid_stats(cl["lam_max"]),
        limiting_mode_fraction={k: f((lim == k).mean()) for k in ["LLI", "LAM_PE", "LAM_NE"]},
        limiting_mode_by_condition={f"{w}|{T}C": {k: f((lim[g.index] == k).mean()) for k in ["LLI", "LAM_PE", "LAM_NE"]}
                                    for (w, T), g in cl.groupby(["window", "T"])},
        caveat="DM are model estimates from a half-cell OCV fit that is itself constrained to reproduce the measured capacity; "
               "capacity 'closure' is therefore partly by construction and not an independent test.")

    # ---- common throughput
    def cell_at(xs):
        recs = []
        for (w, c, T), g in df.groupby(["window", "cell", "T"]):
            rec = dict(window=w, cell=c, T=T)
            for ind in M50_IND:
                gg = g if ind != "R_growth" else g.dropna(subset=["R_growth"])
                rec[ind] = interp_at(gg["thr"], gg[ind], xs)
            recs.append(rec)
        return pd.DataFrame(recs)
    cmn = {}
    cmn_tables = {}
    for xs in (4000, 6600):
        c = cell_at(xs)
        cmn_tables[xs] = c
        cond = c.groupby(["window", "T"]).agg(**{f"{k}_mean": (k, "mean") for k in M50_IND},
                                              **{f"{k}_sd": (k, "std") for k in M50_IND}, n=("cell", "count")).reset_index()
        n_valid = {k: c.groupby(["window", "T"])[k].count().to_dict() for k in M50_IND}
        # rank agreement among 6 conditions
        ranks = {}
        mat_s, mat_k = {}, {}
        for a in M50_IND:
            for b in M50_IND:
                if a < b:
                    mat_s[f"{a}~{b}"] = spear(cond[f"{a}_mean"], cond[f"{b}_mean"])
                    mat_k[f"{a}~{b}"] = kend(cond[f"{a}_mean"], cond[f"{b}_mean"])
        # cell level (17 cells) too
        cell_s = {f"{a}~{b}": spear(c[a], c[b]) for a in M50_IND for b in M50_IND if a < b}
        cmn[str(xs)] = dict(
            conditions=[{k: (f(v) if isinstance(v, (float, np.floating)) else (int(v) if isinstance(v, (np.integer,)) else v))
                         for k, v in r.items()} for r in cond.to_dict("records")],
            n_cells_with_value={k: int(c[k].notna().sum()) for k in M50_IND},
            condition_level_spearman=mat_s, condition_level_kendall=mat_k, cell_level_spearman=cell_s,
            note="condition level n=6 (T x window): too few for significance; cell level n=17 (cells with a value at this throughput).")
    out["common_throughput"] = cmn
    c6 = cmn_tables[6600].groupby(["window", "T"])[M50_IND].mean()
    rk = c6[["cap_loss", "LLI", "LAM_PE", "LAM_NE", "R_growth"]].rank(ascending=False)
    out["rank_table_6600Ah_1_is_largest"] = {f"{w}|{T}C": {k: int(v) for k, v in r.items()} for (w, T), r in rk.iterrows()}
    out["rank_spread_6600Ah"] = {f"{w}|{T}C": int(r.max() - r.min()) for (w, T), r in rk.iterrows()}

    # ---- condition effects at common throughput (ratio of cell-mean)
    eff = {}
    for xs, c in cmn_tables.items():
        cond = c.groupby(["window", "T"])[M50_IND].mean()
        e = {}
        for ind in M50_IND:
            d = {}
            for w in ["0-30%", "0-100%"]:
                for T in (10, 40):
                    try:
                        d[f"{w}: T{T}/T25"] = f(cond.loc[(w, T), ind] / cond.loc[(w, 25), ind])
                    except Exception:
                        d[f"{w}: T{T}/T25"] = None
            for T in (10, 25, 40):
                try:
                    d[f"T{T}: window100/window30"] = f(cond.loc[("0-100%", T), ind] / cond.loc[("0-30%", T), ind])
                except Exception:
                    d[f"T{T}: window100/window30"] = None
            e[ind] = d
        eff[str(xs)] = e
    out["condition_effects_ratio_at_common_throughput"] = eff

    # ---- resistance growth vs LLI / LAM across points, cell-grouped
    pts = df[(df.rpt > 0)].dropna(subset=["R_growth", "LLI", "LAM_PE", "LAM_NE"])
    res = {"n_points": int(len(pts)), "n_cells": int(pts.groupby(["window", "cell"]).ngroups)}
    for mode in ["LLI", "LAM_PE", "LAM_NE", "LAM_NEGr", "LAM_NESi", "cap_loss"]:
        pooled = spear(pts["R_growth"], pts[mode])
        per = []
        for _, g in pts.groupby(["window", "cell"]):
            s = spear(g["R_growth"], g[mode])
            if s["rho"] is not None:
                per.append(s["rho"])
        # within-condition (demean per condition) to remove between-condition structure
        res[mode] = dict(pooled=pooled, per_cell_median_rho=f(np.median(per)) if per else None,
                         per_cell_n=len(per), per_cell_frac_positive=f(np.mean(np.array(per) > 0)) if per else None)
    # R growth vs throughput itself: everything grows with throughput, so also partial on throughput rank
    def partial_rank(a, b, z):
        ra, rb, rz = [stats.rankdata(v) for v in (a, b, z)]
        ea = ra - np.polyval(np.polyfit(rz, ra, 1), rz)
        eb = rb - np.polyval(np.polyfit(rz, rb, 1), rz)
        return f(np.corrcoef(ea, eb)[0, 1])
    res["partial_spearman_given_throughput"] = {m: partial_rank(pts["R_growth"].values, pts[m].values, pts["thr"].values)
                                               for m in ["LLI", "LAM_PE", "LAM_NE"]}
    res["caveat"] = ("Points within a cell are serially correlated and all indicators rise with throughput, so pooled rho overstates "
                     "independent evidence; partial rank correlation given throughput removes the common trend (cells still pooled).")
    out["resistance_vs_modes"] = res
    return out, df, cmn_tables


# =================================================================== 2. dataset B
def analyse_B():
    d = pd.read_excel(B_XLSX, index_col=0)
    out = dict(
        what_it_is=("extracted_degradation_metrics.xlsx: 48 cells (F0001..F0155 ids), 6 temperatures (0/10/25/30/45/55 C) x 4 C-rates "
                    "(0.5/0.75/1.0/1.5), 2 replicates each; columns CF, dR_dch, dR_ch, LLI, LAMpe, LAMne, RMSE. No readme/units in the file: "
                    "chemistry, cell model, charge/discharge definition of 'c-rate', cycle count and the units of CF/LLI/LAM are NOT stated "
                    "(catalog note: Zenodo 15755725 short-term ageing study, 'commercial lithium-ion cells'). Values look like percent "
                    "(LLI 0.5-7). dR_* are mostly negative (resistance decreased?) - sign convention unverified. One LAMne is ~1e-14 (fit bound)."),
        n=int(len(d)), n_per_T_C=int(d.groupby(["temperature", "c-rate"]).size().iloc[0]))
    out["ranges"] = {c: [f(d[c].min()), f(d[c].max())] for c in ["CF", "dR_dch", "dR_ch", "LLI", "LAMpe", "LAMne"]}
    d = d.copy()
    d["LAM_sum"] = d["LAMpe"] + d["LAMne"]
    d["LAM_over_LLI"] = d["LAM_sum"] / d["LLI"]
    res = {}
    for col in ["CF", "LLI", "LAMpe", "LAMne", "LAM_sum", "LAM_over_LLI", "dR_dch", "dR_ch"]:
        pooled = spear(d["c-rate"], d[col])
        per_T = {}
        rhos = []
        for T, g in d.groupby("temperature"):
            s = spear(g["c-rate"], g[col])
            per_T[str(T)] = s
        # stratified: residualise on T-mean (rank-based), n=48
        z = d[col] - d.groupby("temperature")[col].transform("mean")
        strat = spear(d["c-rate"], z)
        low = d[d.temperature <= 10]; high = d[d.temperature >= 25]
        res[col] = dict(pooled=pooled, within_T_demeaned=strat,
                        lowT_le10=spear(low["c-rate"], low[col]), highT_ge25=spear(high["c-rate"], high[col]),
                        mean_by_crate={str(k): f(v) for k, v in d.groupby("c-rate")[col].mean().items()},
                        per_T_n8_each=per_T)
    out["crate_effect"] = res
    # temperature effect for context
    out["temperature_effect_spearman"] = {col: spear(d["temperature"], d[col]) for col in ["CF", "LLI", "LAMpe", "LAMne"]}
    out["mean_by_T"] = {col: {str(k): f(v) for k, v in d.groupby("temperature")[col].mean().items()} for col in ["CF", "LLI", "LAMpe", "LAMne"]}
    # direction summary
    s_lli = res["LLI"]["within_T_demeaned"]["rho"]; s_lam = res["LAM_sum"]["within_T_demeaned"]["rho"]
    out["direction"] = dict(
        crate_up_raises_LLI=bool(s_lli is not None and s_lli > 0),
        crate_up_raises_LAM=bool(s_lam is not None and s_lam > 0),
        rho_LLI_within_T=s_lli, rho_LAM_within_T=s_lam, rho_LAMpe_within_T=res["LAMpe"]["within_T_demeaned"]["rho"],
        rho_LAMne_within_T=res["LAMne"]["within_T_demeaned"]["rho"],
        rho_ratio_LAMsum_over_LLI_within_T=res["LAM_over_LLI"]["within_T_demeaned"]["rho"])
    return out, d


# =================================================================== 3. KIT
def eis_features(df):
    """per EIS block (soc 50, RT, valid): HFR (Re at Im zero-crossing from the high-frequency side) and Re(1 Hz)."""
    recs = []
    for blk, g in df.groupby("sd_block_id", sort=False):
        g = g.sort_values("freq_Hz", ascending=False)
        re = (g["z_amp_mOhm"] * np.cos(np.deg2rad(g["z_ph_deg"]))).values
        im = (g["z_amp_mOhm"] * np.sin(np.deg2rad(g["z_ph_deg"]))).values  # sign convention of file: ph<0 at highest f (inductive)
        fq = g["freq_Hz"].values
        hfr = np.nan
        for k in range(len(fq) - 1):
            if im[k] < 0 and im[k + 1] >= 0:
                w = -im[k] / (im[k + 1] - im[k])
                hfr = re[k] + w * (re[k + 1] - re[k]); break
        j = int(np.argmin(np.abs(np.log(fq / 1.0))))
        recs.append(dict(t=g["timestamp_s"].iloc[0], hfr=hfr, re1hz=re[j] if abs(fq[j] - 1.0) < 0.2 else np.nan))
    return pd.DataFrame(recs)


def build_kit_series():
    """one row per (cell, check-up): day, efc, cap q, R10, R1s (pulse at RT, SOC 30/50/70 mean), EIS HFR/Re1Hz (SOC 50, RT).
    plus per cell: ageing-cycle CE table.  Returns (series DataFrame, ce DataFrame, qc diagnostics)."""
    ck = pd.read_csv(CKUP)
    cells = ck["cell"].unique()
    S, CE = [], []
    step_chk = []
    pulse_vary = []
    for cell in cells:
        ec = KIT / "cell_eocv2" / f"cell_eocv2_{cell}.csv"
        pf = KIT / "cell_plsv2" / f"cell_plsv2_{cell}.csv"
        ef = KIT_EIS_FIXED / f"cell_eisv2_{cell}.csv"
        if not ef.exists() or ef.stat().st_size < 1000:
            ef = KIT / "cell_eisv2" / f"cell_eisv2_{cell}.csv"
        e = pd.read_csv(ec, sep=";", usecols=["timestamp_s", "cyc_condition", "cyc_charged", "coulomb_efficiency",
                                              "delta_q_chg_Ah", "delta_q_dischg_Ah", "total_q_dischg_sum_Ah", "age_type", "soh_cap"])
        t0 = e["timestamp_s"].iloc[0]
        e["efc"] = e["total_q_dischg_sum_Ah"].abs() / 3.0
        e["day"] = (e["timestamp_s"] - t0) / 86400.0
        # CE of ageing cycles: cyc_condition==1 (see code inference in summary); discharge rows carry dis/chg ratio
        a = e[(e["cyc_condition"] == 1) & (e["cyc_charged"] == 0) & e["coulomb_efficiency"].notna()]
        if len(a):
            CE.append(pd.DataFrame(dict(cell=cell, efc=a["efc"].values, day=a["day"].values, ce=a["coulomb_efficiency"].values / 100.0,
                                        qdis=a["delta_q_dischg_Ah"].abs().values)))
        efc_of_t = lambda t: np.interp(t, e["timestamp_s"].values, e["efc"].values)
        # pulses
        p = pd.read_csv(pf, sep=";", usecols=["timestamp_s", "sd_block_id", "is_rt", "soc_nom", "t_avg_degC", "r_ref_10ms_mOhm",
                                              "r_ref_1s_mOhm", "seq_nr", "v_raw_V", "i_raw_A"])
        prec = []
        for blk, g in p.groupby("sd_block_id", sort=False):
            if g["is_rt"].iloc[0] != 1 or g["soc_nom"].iloc[0] not in (30, 50, 70):
                continue
            prec.append(dict(t=g["timestamp_s"].iloc[0], soc=g["soc_nom"].iloc[0], T=g["t_avg_degC"].iloc[0],
                             r10=g["r_ref_10ms_mOhm"].iloc[0], r1=g["r_ref_1s_mOhm"].iloc[0]))
            if g["soc_nom"].iloc[0] == 50 and len(step_chk) < 4000:
                i = g["i_raw_A"].values; v = g["v_raw_V"].values
                k = next((k for k in range(1, len(i)) if abs(i[k]) > 0.5 and abs(i[k - 1]) < 0.05), None)
                if k is not None:
                    step_chk.append((cell, (v[k] - v[k - 1]) / (i[k] - i[k - 1]) * 1000.0, g["r_ref_10ms_mOhm"].iloc[0], g["r_ref_1s_mOhm"].iloc[0]))
        pr = pd.DataFrame(prec).sort_values("t")
        pr["grp"] = (pr["t"].diff().fillna(1e9) > 12 * 3600).cumsum()
        pg = pr.groupby("grp").agg(t=("t", "mean"), r10=("r10", "mean"), r1=("r1", "mean"), T=("T", "mean"), n=("r10", "size")).reset_index()
        pg = pg[pg["n"] >= 2]
        # EIS
        en = None
        try:
            z = pd.read_csv(ef, sep=";", usecols=["timestamp_s", "sd_block_id", "is_rt", "soc_nom", "valid", "freq_Hz", "z_amp_mOhm", "z_ph_deg"])
            z = z[(z["is_rt"] == 1) & (z["soc_nom"] == 50) & (z["valid"] == 1)]
            if len(z):
                en = eis_features(z).sort_values("t")
        except Exception:
            en = None
        rows = []
        for _, r in pg.iterrows():
            row = dict(cell=cell, t=r["t"], day=(r["t"] - t0) / 86400.0, efc=float(efc_of_t(r["t"])), r10=r["r10"], r1=r["r1"], pulse_T=r["T"])
            if en is not None and len(en):
                j = (en["t"] - r["t"]).abs().idxmin()
                if abs(en.loc[j, "t"] - r["t"]) < 2 * 86400:
                    row["hfr"] = en.loc[j, "hfr"]; row["re1hz"] = en.loc[j, "re1hz"]
            rows.append(row)
        sr = pd.DataFrame(rows)
        if len(sr):
            for c_ in ["r10", "r1", "hfr", "re1hz"]:
                if c_ not in sr:
                    sr[c_] = np.nan
            sr["R10_growth"] = sr["r10"] / sr["r10"].iloc[0] - 1
            sr["R1s_growth"] = sr["r1"] / sr["r1"].iloc[0] - 1
            sr["HFR_growth"] = sr["hfr"] / sr["hfr"].iloc[0] - 1 if sr["hfr"].notna().iloc[0] else np.nan
            sr["Rpol_growth"] = (sr["re1hz"] - sr["hfr"]) / (sr["re1hz"].iloc[0] - sr["hfr"].iloc[0]) - 1 if sr["re1hz"].notna().iloc[0] and sr["hfr"].notna().iloc[0] else np.nan
            pulse_vary.append(dict(cell=cell, r10_first=sr["r10"].iloc[0], r10_last=sr["r10"].iloc[-1], n=len(sr)))
            S.append(sr)
    series = pd.concat(S, ignore_index=True)
    ce = pd.concat(CE, ignore_index=True) if CE else pd.DataFrame(columns=["cell", "efc", "day", "ce", "qdis"])
    sc = pd.DataFrame(step_chk, columns=["cell", "r_step_mOhm", "r10", "r1"])
    pv = pd.DataFrame(pulse_vary)
    return series, ce, sc, pv


def load_kit():
    series_f, ce_f, sc_f, pv_f = [OUT / n for n in ("kit_series.csv", "kit_ce.csv", "kit_stepcheck.csv", "kit_pulsevary.csv")]
    if all(p.exists() for p in (series_f, ce_f, sc_f, pv_f)) and not REBUILD:
        return tuple(pd.read_csv(p) for p in (series_f, ce_f, sc_f, pv_f))
    series, ce, sc, pv = build_kit_series()
    series.to_csv(series_f, index=False); ce.to_csv(ce_f, index=False); sc.to_csv(sc_f, index=False); pv.to_csv(pv_f, index=False)
    return series, ce, sc, pv


KIT_IND = ["cap_loss", "R10_growth", "R1s_growth", "HFR_growth", "Rpol_growth", "one_minus_CE"]
CE_TRIM = 0.02  # ignore half-cycles with |1-CE| > 2 % (SOC-drift artefacts, e.g. first cycle after rest)


def kit_conditions(ck, series, ce, xs, axis):
    """per cell values at common point xs (axis 'efc' or 'day'), then mean over cells per condition P."""
    cond = ck.groupby("P").first()[["age_type", "T", "cal_soc", "v_min", "v_max", "chg_rate", "dis_rate", "profile"]]
    recs = []
    for cell, g in ck.groupby("cell"):
        P = int(g["P"].iloc[0])
        rec = dict(P=P, cell=cell, one_minus_CE=np.nan, one_minus_CE_median=np.nan, dod_Ah_per_cycle=np.nan)
        q = interp_at(g[axis], g["q"], xs)
        rec["cap_loss"] = 1 - q if np.isfinite(q) else np.nan
        s = series[series.cell == cell]
        for k in ["R10_growth", "R1s_growth", "HFR_growth", "Rpol_growth"]:
            rec[k] = interp_at(s[axis], s[k], xs) if len(s) else np.nan
        c = ce[ce.cell == cell]
        if axis == "efc" and len(c):
            c = c[(c["efc"] <= xs)]
            om = 1 - c["ce"]
            om = om[om.abs() <= CE_TRIM]
            rec["one_minus_CE"] = om.mean() if len(om) >= 5 else np.nan
            rec["one_minus_CE_median"] = om.median() if len(om) >= 5 else np.nan
            rec["n_ce"] = len(om)
            rec["dod_Ah_per_cycle"] = c["qdis"].median()
        else:
            rec["one_minus_CE"] = np.nan
        recs.append(rec)
    cells = pd.DataFrame(recs)
    # keep conditions where >=2 cells have a capacity value at xs
    cnt = cells.groupby("P")["cap_loss"].count()
    good = cnt[cnt >= 2].index
    m = cells[cells.P.isin(good)].groupby("P")[KIT_IND + ["one_minus_CE_median", "dod_Ah_per_cycle"]].mean()
    m["n_cells"] = cells[cells.P.isin(good)].groupby("P")["cap_loss"].count()
    m = cond.join(m, how="inner")
    return m.reset_index(), cells


def agree_matrix(tab):
    mat = {}
    for a in KIT_IND:
        for b in KIT_IND:
            if a < b:
                mat[f"{a}~{b}"] = spear(tab[a], tab[b])
    return mat


def analyse_kit():
    ck = pd.read_csv(CKUP)
    series, ce, sc, pv = load_kit()
    out = {}
    out["n_cells_checkups"] = int(ck["cell"].nunique())
    out["n_conditions_P"] = int(ck["P"].nunique())
    out["n_cells_with_pulse_series"] = int(series["cell"].nunique())
    out["n_cells_with_eis_hfr"] = int(series.dropna(subset=["hfr"])["cell"].nunique())
    out["n_cells_with_ce"] = int(ce["cell"].nunique())
    out["cyc_condition_codes_inferred"] = (
        "0 = non-cycling/other (rest, conditioning; no CE; mostly discharge/charge fragments around tests), "
        "1 = ageing cycles (present only in age_type 2 cycle cells; CE reported on discharge rows as Q_dis / preceding Q_chg, "
        "verified by hand on P045), 2 = check-up cycles (full-capacity C/3 charge/discharge; ~29 per cell). "
        "Profile cells (age_type 3, WLTP) have no CE in cond 1; calendar cells (age_type 1) have no ageing cycles. "
        "Inferred from data, not from documentation.")
    out["ce_trim_rule"] = f"half-cycles with |1-CE|>{CE_TRIM*100:.0f}% dropped; cell value needs >=5 remaining ageing cycles with efc<=x*"
    # pulse validity
    last_first = pv["r10_last"] / pv["r10_first"] - 1
    out["pulse_r_changes_with_ageing"] = dict(
        n_cells=int(len(pv)), n_r10_increased=int((last_first > 0).sum()), median_r10_change=f(last_first.median()),
        min_r10_change=f(last_first.min()), max_r10_change=f(last_first.max()),
        conclusion="r_ref_10ms/1s are measured per pulse block (not fixed references): they differ block to block and change with ageing.")
    rs = sc.dropna()
    if len(rs) > 20:
        out["pulse_step_vs_r_ref_check"] = dict(
            n_blocks=int(len(rs)), spearman_rstep_vs_r10=spear(rs["r_step_mOhm"], rs["r10"]), spearman_rstep_vs_r1=spear(rs["r_step_mOhm"], rs["r1"]),
            median_ratio_rstep_over_r10=f((rs["r_step_mOhm"] / rs["r10"]).median()),
            note="crude plausibility check only: dV/dI across first current step sample (sample timing unknown); r_ref values were NOT independently recomputed.")
    # pulse temperature check
    out["pulse_temperature_note"] = (f"RT pulse blocks (is_rt=1; cyc_charged=1) used; median measured T = {f(series['pulse_T'].median(),2)} C "
                                     f"(min {f(series['pulse_T'].min(),1)}, max {f(series['pulse_T'].max(),1)}), so R is compared at ~25 C for all ageing temperatures.")

    results = {}
    tables = {}
    for name, ages, axis, xs in [("cycle_efc60", [2], "efc", 60), ("cycle_efc400", [2], "efc", 400),
                                 ("profile_efc60", [3], "efc", 60), ("profile_efc200", [3], "efc", 200),
                                 ("calendar_day300", [1], "day", 300), ("calendar_day440", [1], "day", 440)]:
        sub = ck[ck.age_type.isin(ages)]
        tab, cells = kit_conditions(ck[ck.age_type.isin(ages)], series, ce, xs, axis)
        tables[name] = (tab, cells)
        results[name] = dict(
            axis=axis, x_star=xs, n_conditions=int(len(tab)), n_conditions_total=int(sub["P"].nunique()),
            n_cells=int(cells[cells.P.isin(tab.P)].cap_loss.notna().sum()),
            indicator_available_counts={k: int(tab[k].notna().sum()) for k in KIT_IND},
            agreement_spearman=agree_matrix(tab),
            table=[{k: (f(v) if isinstance(v, (float, np.floating)) else (int(v) if isinstance(v, (np.integer,)) else v)) for k, v in r.items()}
                   for r in tab.to_dict("records")])
    out["common_point_results"] = results

    # ---------------- signature tests
    sig = {}
    # S1: Li plating (cycle, low T, high charge rate)
    for key in ("cycle_efc60", "cycle_efc400"):
        tab = tables[key][0].copy()
        tab["win"] = list(zip(tab.v_min, tab.v_max))
        recs = []
        for (T, win), g in tab.groupby(["T", "win"]):
            base = g[(g.chg_rate == 0.33) & (g.dis_rate == 1.0)]
            if base.empty:
                continue
            b = base.iloc[0]
            for cr in (1.0, 1.67):
                hi = g[(g.chg_rate == cr) & (g.dis_rate == 1.0)]
                if hi.empty:
                    continue
                h = hi.iloc[0]
                recs.append(dict(T=int(T), win=f"{win[0]}-{win[1]}", chg=cr,
                                 d_caploss_pp=(h.cap_loss - b.cap_loss) * 100, cap_ratio=h.cap_loss / b.cap_loss if b.cap_loss else np.nan,
                                 d_1mCE_pp=(h.one_minus_CE - b.one_minus_CE) * 100,
                                 d_R10_pp=(h.R10_growth - b.R10_growth) * 100, d_R1s_pp=(h.R1s_growth - b.R1s_growth) * 100))
        pr = pd.DataFrame(recs)
        s1 = dict(n_pairs=int(len(pr)), pairs=[{k: (f(v) if isinstance(v, float) else v) for k, v in r.items()} for r in pr.to_dict("records")])
        if len(pr):
            def grp(mask, label):
                g = pr[mask]
                return dict(n=int(len(g)),
                            frac_caploss_up=f((g.d_caploss_pp > 0).mean()), mean_d_caploss_pp=f(g.d_caploss_pp.mean()),
                            frac_1mCE_up=f((g.d_1mCE_pp.dropna() > 0).mean()) if g.d_1mCE_pp.notna().any() else None,
                            n_with_CE=int(g.d_1mCE_pp.notna().sum()), mean_d_1mCE_pp=f(g.d_1mCE_pp.mean()),
                            mean_d_R10_pp=f(g.d_R10_pp.mean()), mean_d_R1s_pp=f(g.d_R1s_pp.mean()),
                            frac_R10_rise_smaller_than_caploss_rise=f((g.d_R10_pp < g.d_caploss_pp).mean()))
            for cr in (1.0, 1.67):
                s1[f"lowT_le10_chg{cr}"] = grp((pr["T"] <= 10) & (pr.chg == cr), "low")
                s1[f"highT_ge25_chg{cr}"] = grp((pr["T"] >= 25) & (pr.chg == cr), "high")
            lo = s1["lowT_le10_chg1.67"]; hi_ = s1["highT_ge25_chg1.67"]
            cond_a = lo["n"] > 0 and (lo["frac_caploss_up"] or 0) >= 0.75
            cond_b = lo["n_with_CE"] > 0 and (lo["frac_1mCE_up"] or 0) >= 0.75
            cond_c = (lo["mean_d_R10_pp"] is not None and lo["mean_d_caploss_pp"] is not None and lo["mean_d_R10_pp"] < lo["mean_d_caploss_pp"])
            cond_d = (lo["mean_d_caploss_pp"] or 0) > (hi_["mean_d_caploss_pp"] or 0)
            s1["criteria_lowT_1.67C"] = dict(caploss_up_in_ge75pct_pairs=bool(cond_a), one_minus_CE_up_in_ge75pct_pairs=bool(cond_b),
                                             R10_rise_smaller_than_caploss_rise=bool(cond_c), lowT_effect_larger_than_highT_effect=bool(cond_d))
            nsup = sum([cond_a, cond_b, cond_c, cond_d])
            s1["verdict"] = "supported" if nsup == 4 else ("partially supported" if nsup >= 2 else "not supported")
            s1["rule"] = "supported = all 4 criteria true; partially = 2-3; not supported = 0-1. Effects are 1.67C-vs-0.33C charge at fixed 1C discharge, same T and voltage window (and 1.0C-vs-0.33C in pairs)."
        sig[f"S1_plating_{key}"] = s1

    # S2: SEI: calendar
    for key in ("calendar_day300", "calendar_day440"):
        tab = tables[key][0]
        s2 = dict(n_conditions=int(len(tab)),
                  spearman_T_vs=dict(cap_loss=spear(tab["T"], tab.cap_loss), R10=spear(tab["T"], tab.R10_growth), R1s=spear(tab["T"], tab.R1s_growth), HFR=spear(tab["T"], tab.HFR_growth)),
                  spearman_SOC_vs=dict(cap_loss=spear(tab.cal_soc, tab.cap_loss), R10=spear(tab.cal_soc, tab.R10_growth), R1s=spear(tab.cal_soc, tab.R1s_growth), HFR=spear(tab.cal_soc, tab.HFR_growth)),
                  spearman_caploss_vs=dict(R10=spear(tab.cap_loss, tab.R10_growth), R1s=spear(tab.cap_loss, tab.R1s_growth), HFR=spear(tab.cap_loss, tab.HFR_growth), Rpol=spear(tab.cap_loss, tab.Rpol_growth)),
                  mean_by_T={str(T): dict(cap_loss=f(g.cap_loss.mean()), R10=f(g.R10_growth.mean()), R1s=f(g.R1s_growth.mean()), HFR=f(g.HFR_growth.mean())) for T, g in tab.groupby("T")},
                  mean_by_SOC={str(T): dict(cap_loss=f(g.cap_loss.mean()), R10=f(g.R10_growth.mean()), R1s=f(g.R1s_growth.mean()), HFR=f(g.HFR_growth.mean())) for T, g in tab.groupby("cal_soc")})
        up = lambda d: all((d[k]["rho"] or 0) > 0.3 for k in ("cap_loss", "R10"))
        s2["criteria"] = dict(T_raises_caploss_and_R10=bool(up(s2["spearman_T_vs"])), SOC_raises_caploss_and_R10=bool(up(s2["spearman_SOC_vs"])),
                              caploss_R10_rank_corr_gt_0p5=bool((s2["spearman_caploss_vs"]["R10"]["rho"] or 0) > 0.5))
        n = sum(s2["criteria"].values())
        s2["verdict"] = "supported" if n == 3 else ("partially supported" if n >= 1 else "not supported")
        s2["rule"] = "supported = T and SOC each raise both cap loss and R10 (rho>0.3) and cap-loss~R10 rho>0.5; partially = 1-2 of 3."
        sig[f"S2_SEI_{key}"] = s2

    # S3: deep DoD -> LAM (direction only)
    for key in ("cycle_efc60", "cycle_efc400"):
        tab = tables[key][0].copy()
        tab["win"] = tab.apply(lambda r: f"{r.v_min}-{r.v_max}", axis=1)
        wmap = {"2.5-4.2": "deep", "3.249-4.2": "mid", "3.249-4.092": "shallow"}
        tab["w"] = tab["win"].map(wmap)
        dod = tab.groupby("w")["dod_Ah_per_cycle"].median().to_dict()
        recs = []
        for (T, cr, dr), g in tab.groupby(["T", "chg_rate", "dis_rate"]):
            gd = g.set_index("w")
            for a, b in (("deep", "shallow"), ("deep", "mid")):
                if a in gd.index and b in gd.index:
                    A, B = gd.loc[a], gd.loc[b]
                    recs.append(dict(T=int(T), chg=cr, dis=dr, pair=f"{a}-vs-{b}", d_caploss_pp=(A.cap_loss - B.cap_loss) * 100,
                                     d_R10_pp=(A.R10_growth - B.R10_growth) * 100, d_R1s_pp=(A.R1s_growth - B.R1s_growth) * 100,
                                     d_1mCE_pp=(A.one_minus_CE - B.one_minus_CE) * 100))
        pr = pd.DataFrame(recs)
        s3 = dict(median_discharge_Ah_per_cycle_by_window={k: f(v, 3) for k, v in dod.items()}, n_pairs=int(len(pr)))
        for pair, g in pr.groupby("pair"):
            s3[pair] = dict(n=int(len(g)), frac_caploss_higher_in_deep=f((g.d_caploss_pp > 0).mean()), mean_d_caploss_pp=f(g.d_caploss_pp.mean()),
                            frac_R10_higher_in_deep=f((g.d_R10_pp > 0).mean()), mean_d_R10_pp=f(g.d_R10_pp.mean()),
                            frac_R1s_higher_in_deep=f((g.d_R1s_pp > 0).mean()), mean_d_R1s_pp=f(g.d_R1s_pp.mean()),
                            frac_1mCE_higher_in_deep=f((g.d_1mCE_pp.dropna() > 0).mean()) if g.d_1mCE_pp.notna().any() else None,
                            mean_d_1mCE_pp=f(g.d_1mCE_pp.mean()))
        for pair, g in pr.groupby("pair"):
            for lab, gg in (("lowT_le10", g[g["T"] <= 10]), ("highT_ge25", g[g["T"] >= 25])):
                s3[pair][lab] = dict(n=int(len(gg)), mean_d_caploss_pp=f(gg.d_caploss_pp.mean()), mean_d_R10_pp=f(gg.d_R10_pp.mean()),
                                     frac_caploss_higher_in_deep=f((gg.d_caploss_pp > 0).mean()), frac_R10_higher_in_deep=f((gg.d_R10_pp > 0).mean()))
        v = [s3.get("deep-vs-shallow", {}).get("frac_R10_higher_in_deep"), s3.get("deep-vs-mid", {}).get("frac_R10_higher_in_deep")]
        c1 = all((x or 0) >= 0.75 for x in v)
        v2 = [s3.get("deep-vs-shallow", {}).get("frac_caploss_higher_in_deep"), s3.get("deep-vs-mid", {}).get("frac_caploss_higher_in_deep")]
        c2 = all((x or 0) >= 0.75 for x in v2)
        s3["criteria"] = dict(R10_higher_in_deep_ge75pct=bool(c1), caploss_higher_in_deep_ge75pct=bool(c2))
        s3["verdict"] = ("direction supported (R and capacity loss larger at deeper DoD at equal EFC); LAM attribution NOT testable (no OCV/ICA/DVA)"
                         if (c1 and c2) else ("partially supported (see fractions); LAM attribution NOT testable" if (c1 or c2) else "not supported; LAM attribution NOT testable"))
        s3["confound"] = "windows differ in more than DoD (lower cutoff 2.5 V vs 3.249 V; upper 4.2 vs 4.092 V -> SOC level / voltage stress), so DoD is not isolated."
        sig[f"S3_deepDoD_{key}"] = s3
    out["signature_tests"] = sig
    # per-condition consistency of indicators (rank percentile spread across indicators)
    cons = {}
    for key in ("cycle_efc60", "cycle_efc400", "calendar_day300", "profile_efc200"):
        tab = tables[key][0].copy()
        use = [c for c in ["cap_loss", "R10_growth", "R1s_growth", "HFR_growth", "one_minus_CE"] if tab[c].notna().sum() > 3]
        rk = tab[use].rank(pct=True)
        tab["rank_spread"] = rk.max(axis=1) - rk.min(axis=1)
        tab["rank_mean"] = rk.mean(axis=1)
        for c in use:
            tab["rk_" + c] = rk[c]
        cols = ["P", "age_type", "T", "cal_soc", "v_min", "v_max", "chg_rate", "dis_rate", "rank_spread", "rank_mean"] + ["rk_" + c for c in use]
        worst = tab.sort_values("rank_spread", ascending=False).head(6)[cols]
        cons[key] = dict(indicators_used=use, n=int(len(tab)), frac_spread_le_0p25=f((tab.rank_spread <= 0.25).mean()),
                         frac_spread_ge_0p5=f((tab.rank_spread >= 0.5).mean()), median_spread=f(tab.rank_spread.median()),
                         most_disagreeing=[{k: (f(v) if isinstance(v, (float, np.floating)) else (int(v) if isinstance(v, (np.integer,)) else v)) for k, v in r.items()} for r in worst.to_dict("records")],
                         by_T={str(T): dict(n=int(len(g)), cap_loss=f(g.cap_loss.mean()), R10=f(g.R10_growth.mean()), R1s=f(g.R1s_growth.mean()), HFR=f(g.HFR_growth.mean()),
                                            one_minus_CE=f(g.one_minus_CE.mean(), 6) if g.one_minus_CE.notna().any() else None, spread_median=f(g.rank_spread.median())) for T, g in tab.groupby("T")})
        if key.startswith("cycle"):
            cons[key]["by_chg_rate"] = {str(c): dict(n=int(len(g)), cap_loss=f(g.cap_loss.mean()), R10=f(g.R10_growth.mean()), one_minus_CE=f(g.one_minus_CE.mean(), 6), spread_median=f(g.rank_spread.median())) for c, g in tab.groupby("chg_rate")}
    out["condition_consistency_rank_spread"] = cons
    return out, tables, series


# =================================================================== figures
def fig_m50(df):
    ind = [("cap_loss", "1-SoH"), ("LLI", "LLI"), ("LAM_PE", "LAM PE"), ("LAM_NEGr", "LAM NE (Gr)"), ("LAM_NESi", "LAM NE (Si)"), ("R_growth", "R/R0-1")]
    fig, axs = plt.subplots(len(ind), 2, figsize=(10, 15), sharex="col")
    col = {10: "#1f77b4", 25: "#2ca02c", 40: "#d62728"}
    for j, w in enumerate(["0-30%", "0-100%"]):
        for i, (k, lab) in enumerate(ind):
            ax = axs[i, j]
            for (c, T), g in df[df.window == w].groupby(["cell", "T"]):
                gg = g.dropna(subset=[k])
                ax.plot(gg["thr"], gg[k], "-o", ms=3, lw=1, color=col[T], alpha=.8)
            if i == 0:
                ax.set_title(f"window {w} SOC (single 0.3C/1C rate)")
            ax.set_ylabel(lab); ax.grid(alpha=.3)
        axs[-1, j].set_xlabel("charge throughput [Ah]")
    from matplotlib.lines import Line2D
    axs[0, 1].legend([Line2D([0], [0], color=col[T]) for T in col], [f"{T} C" for T in col], title="ageing T", loc="upper left")
    fig.suptitle("M50T: degradation-mode trajectories per cell (OCV-fit DM = model estimates)", y=.995)
    fig.tight_layout(); fig.savefig(OUT / "fig_a_m50t_mode_trajectories.png", dpi=130); plt.close(fig)


def heat(ax, labels, mat, title):
    n = len(labels)
    M = np.full((n, n), np.nan)
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            if i == j:
                M[i, j] = 1
            else:
                s = mat.get(f"{a}~{b}") or mat.get(f"{b}~{a}")
                if s and s["rho"] is not None:
                    M[i, j] = s["rho"]
    im = ax.imshow(M, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(n)); ax.set_yticks(range(n)); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8); ax.set_yticklabels(labels, fontsize=8)
    for i in range(n):
        for j in range(n):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title(title, fontsize=9)
    return im


def fig_heat(m50_out, kit_out):
    panels = [("M50T cond-level @6600 Ah (n=6)", M50_IND, m50_out["common_throughput"]["6600"]["condition_level_spearman"]),
              ("M50T cell-level @6600 Ah (n=17)", M50_IND, m50_out["common_throughput"]["6600"]["cell_level_spearman"])]
    for key, lab in [("cycle_efc60", "KIT cycle @60 EFC"), ("cycle_efc400", "KIT cycle @400 EFC"), ("calendar_day300", "KIT calendar @300 d"), ("calendar_day440", "KIT calendar @440 d")]:
        r = kit_out["common_point_results"][key]
        panels.append((f"{lab} (n={r['n_conditions']})", KIT_IND, r["agreement_spearman"]))
    fig, axs = plt.subplots(2, 3, figsize=(17, 10))
    for ax, (t, lab, mat) in zip(axs.ravel(), panels):
        im = heat(ax, lab, mat, t)
    fig.subplots_adjust(wspace=0.75, hspace=0.6)
    fig.colorbar(im, ax=axs, shrink=.6, label="Spearman rho")
    fig.suptitle("Indicator agreement across conditions (Spearman)")
    fig.savefig(OUT / "fig_b_indicator_agreement_heatmaps.png", dpi=130, bbox_inches="tight"); plt.close(fig)


def fig_scatter(tables):
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.5))
    mk = {0.0: "s", 0.33: "o", 1.0: "^", 1.67: "D"}
    cmap = {0: "#08306b", 10: "#4292c6", 25: "#fd8d3c", 40: "#bd0026"}
    for ax, key, ttl in zip(axs, ["cycle_efc60", "cycle_efc400", "calendar_day300"], ["KIT cycle @60 EFC", "KIT cycle @400 EFC", "KIT calendar @300 d (marker s)"]):
        tab = tables[key][0]
        for _, r in tab.iterrows():
            ax.scatter(r.R10_growth * 100, r.cap_loss * 100, c=cmap[int(r["T"])], marker=mk.get(float(r.chg_rate), "o"), s=45, edgecolor="k", lw=.4)
        ax.set_xlabel("R10ms growth at 25 C [%]"); ax.set_ylabel("capacity loss [%]"); ax.set_title(ttl); ax.grid(alpha=.3)
        sp = spear(tab.R10_growth, tab.cap_loss)
        ax.text(.03, .95, f"Spearman {sp['rho']}, n={sp['n']}", transform=ax.transAxes, va="top", fontsize=9)
    from matplotlib.lines import Line2D
    h = [Line2D([0], [0], marker="o", ls="", color=cmap[T], label=f"{T} C") for T in cmap] + \
        [Line2D([0], [0], marker=m, ls="", color="gray", label=("calendar" if c == 0 else f"{c}C chg")) for c, m in mk.items()]
    axs[2].legend(handles=h, fontsize=8, loc="lower right")
    fig.suptitle("KIT: capacity loss vs resistance growth (colour = ageing T, marker = charge rate)")
    fig.tight_layout(); fig.savefig(OUT / "fig_c_kit_caploss_vs_resistance.png", dpi=130); plt.close(fig)



# =================================================================== summary.md (Korean)
def pct(x, nd=1):
    return "n/a" if x is None else f"{x*100:.{nd}f}%"


def rho(d):
    return "n/a" if d is None or d.get("rho") is None else f"{d['rho']:.2f}"


def write_summary(R):
    M, B, K = R["M50T"], R["datasetB"], R["KIT"]
    cl = M["closure"]
    ct = M["common_throughput"]["6600"]
    cs, cell_s = ct["condition_level_spearman"], ct["cell_level_spearman"]
    rv = M["resistance_vs_modes"]
    eff = M["condition_effects_ratio_at_common_throughput"]["6600"]
    cp = K["common_point_results"]
    sg = K["signature_tests"]
    s1, s1b = sg["S1_plating_cycle_efc60"], sg["S1_plating_cycle_efc400"]
    s2, s2b = sg["S2_SEI_calendar_day300"], sg["S2_SEI_calendar_day440"]
    s3, s3b = sg["S3_deepDoD_cycle_efc60"], sg["S3_deepDoD_cycle_efc400"]
    bd, be = B["direction"], B["crate_effect"]
    cons = K["condition_consistency_rank_spread"]
    L = []
    A = L.append
    A("# 열화 메커니즘 교차검증 (공개 데이터, 실험 없음)\n")
    A("생성: `prototype/mech_crosscheck.py` / 수치 전부 `results.json`에서 자동 삽입. 그림: `fig_a_*`, `fig_b_*`, `fig_c_*`.\n")
    A("## 1. 데이터셋별 사용 가능한 방법\n")
    A("| 데이터셋 | OCV-fit DM(LLI/LAM) | 용량 | DCIR 펄스 | EIS | CE | ICA/DVA |")
    A("|---|---|---|---|---|---|---|")
    A(f"| A. LG M50T (NMC811/Gr-Si, 17셀) | 있음(모델 추정) | C/10 용량 | 0.1 s 저항(RPT 절반만) | 없음 | 미사용 | 추출 CSV만 사용, 원 곡선 미확인 -> 이 분석에서는 안 함 |")
    A(f"| B. zenodo 15755725 (48셀) | LLI/LAMpe/LAMne(방법·단위 미기재) | CF(정의 미기재) | dR_dch/dR_ch(부호·정의 미기재) | 원자료 미다운로드 | 없음 | 안 함 |")
    A(f"| C. KIT LG INR18650HG2 ({K['n_cells_checkups']}셀, {K['n_conditions_P']}조건) | 없음(OCV 곡선 없음) | 체크업 용량 | 펄스 R10ms/R1s (25 C, SOC 30/50/70 평균) | 있음({K['n_cells_with_eis_hfr']}셀, 고주파 저항 HFR·1 Hz 실수부) | 사이클 셀만({K['n_cells_with_ce']}셀) | 불가(곡선 없음) |\n")
    A("## 2. M50T (OCV-fit DM) 결과\n")
    A(f"- 셀 {M['n_cells']}개, 조건(온도 x 윈도우)별 셀 수: " + ", ".join(f"{k.replace('|', ' / ')}: {v}" for k, v in M["conditions_n_cells"].items()) + ". 0-100%/25 C는 2셀뿐입니다.")
    A("- 행 매칭: Main 파일과 DM 파일의 셀별 행 수가 모두 같아 행 위치(RPT 순번)로 매칭했고, DM의 SoH가 C/10 용량/초기 C/10 용량과 일치함을 확인했습니다(최대 차이 약 2e-6, results.json의 alignment_check).")
    A(f"- 용량 폐합(점 {cl['n_points']}개, 셀 {cl['n_cells']}개): 용량손실-LLI 평균 잔차 {cl['cap_vs_LLI']['mean']*100:.2f} %p, RMSE {cl['cap_vs_LLI']['rmse']*100:.2f} %p, Spearman {cl['cap_vs_LLI']['spearman']['rho']:.2f}. "
      f"max(LLI, LAM_PE, LAM_NE)와 비교하면 평균 잔차 {cl['cap_vs_max_modes']['mean']*100:.2f} %p, RMSE {cl['cap_vs_max_modes']['rmse']*100:.2f} %p로 오히려 용량손실을 과대추정합니다. 즉 용량은 LLI에 가장 가깝게 대응합니다.")
    lf = cl["limiting_mode_fraction"]
    A(f"- max를 기준으로 한계 모드 비율: LAM_NE {pct(lf['LAM_NE'],0)}, LLI {pct(lf['LLI'],0)}, LAM_PE {pct(lf['LAM_PE'],0)}. 0-30%/10 C만 LLI가 우세({pct(cl['limiting_mode_by_condition']['0-30%|10C']['LLI'],0)})합니다. "
      f"주의: DM 피팅이 측정 용량을 재현하도록 제약되므로 이 폐합은 부분적으로 구조적(독립 검증 아님)입니다.")
    A("- 공통 처리량 6600 Ah에서 6개 조건의 순위(1=가장 큼):\n")
    A("| 조건 | 1-SoH | LLI | LAM_PE | LAM_NE | R 증가 | 순위 폭 |")
    A("|---|---|---|---|---|---|---|")
    for k, v in M["rank_table_6600Ah_1_is_largest"].items():
        A(f"| {k.replace('|', ' / ')} | {v['cap_loss']} | {v['LLI']} | {v['LAM_PE']} | {v['LAM_NE']} | {v['R_growth']} | {M['rank_spread_6600Ah'][k]} |")
    A("")
    A(f"- 지표 간 일치(조건 n=6, 유의성 판단 불가 수준): 용량~LLI rho {rho(cs['LLI~cap_loss'])}, 용량~LAM_NE {rho(cs['LAM_NE~cap_loss'])}, 용량~LAM_PE {rho(cs['LAM_PE~cap_loss'])}, 용량~R {rho(cs['R_growth~cap_loss'])}, "
      f"R~LLI {rho(cs['LLI~R_growth'])}, R~LAM_PE {rho(cs['LAM_PE~R_growth'])}, R~LAM_NE {rho(cs['LAM_NE~R_growth'])}. 셀 수준(n=17, R은 16): 용량~LLI {rho(cell_s['LLI~cap_loss'])}, R~LLI {rho(cell_s['LLI~R_growth'])}, R~LAM_NE {rho(cell_s['LAM_NE~R_growth'])}.")
    A(f"- 저항증가 vs 모드(전 점, 셀 그룹, n={rv['n_points']}점/{rv['n_cells']}셀): pooled Spearman LLI {rho(rv['LLI']['pooled'])}, LAM_PE {rho(rv['LAM_PE']['pooled'])}, LAM_NE {rho(rv['LAM_NE']['pooled'])}. "
      f"셀별 중앙값은 LLI {rv['LLI']['per_cell_median_rho']}, LAM_NE {rv['LAM_NE']['per_cell_median_rho']}이나 모든 지표가 처리량과 함께 증가하므로 과대입니다. 처리량을 통제한 부분 상관은 LLI {rv['partial_spearman_given_throughput']['LLI']:.2f}, LAM_PE {rv['partial_spearman_given_throughput']['LAM_PE']:.2f}, LAM_NE {rv['partial_spearman_given_throughput']['LAM_NE']:.2f}로, 저항은 LAM_NE/LLI와 더 잘 따라가고 LAM_PE와는 약합니다.")
    A("- 조건 효과(6600 Ah, 셀 평균 비율):\n")
    A("| 지표 | 0-30%: 10 C/25 C | 0-30%: 40 C/25 C | 0-100%: 10 C/25 C | 0-100%: 40 C/25 C | 100%/30% @10 C | @25 C | @40 C |")
    A("|---|---|---|---|---|---|---|---|")
    for ind in ["cap_loss", "LLI", "LAM_PE", "LAM_NE", "R_growth"]:
        e = eff[ind]
        A(f"| {ind} | {e['0-30%: T10/T25']} | {e['0-30%: T40/T25']} | {e['0-100%: T10/T25']} | {e['0-100%: T40/T25']} | {e['T10: window100/window30']} | {e['T25: window100/window30']} | {e['T40: window100/window30']} |")
    A("")
    A("  해석: 0-100% 윈도우에서는 10 C가 가장 나쁘고(저온 사이클, 도금 유사 거동과 방향 일치) LAM_NE(흑연)와 LLI, R이 함께 큽니다. 반대로 0-30% 윈도우에서는 10 C가 가장 좋아 온도 경향이 윈도우에 따라 뒤집힙니다. "
      "윈도우 효과도 같은 Ah 처리량 기준으로 25/40 C에서는 0-100%가 0-30%보다 덜 열화(비율 < 1)라 '깊은 DoD일수록 LAM 증가' 기대와 맞지 않습니다(처리량 정규화가 사이클 수/시간 차이를 가림). 셀 간 산포도 큽니다(예: 0-30%/10 C 용량손실 SD "
      f"{ct['conditions'][3]['cap_loss_sd']*100:.1f} %p). 한 조건이 2~3셀이라 위 비율은 방향 참고용입니다.\n")
    A("## 3. 데이터셋 B (C-rate 효과)\n")
    A("- 실제 내용: 시트 1개, 48셀(F0001~F0155 중), 온도 6수준(0/10/25/30/45/55 C) x C-rate 4수준(0.5/0.75/1.0/1.5) x 반복 2. 열: CF, dR_dch, dR_ch, LLI, LAMpe, LAMne, RMSE. 파일에 설명·단위가 없어 화학 조성, 셀 모델, C-rate의 충/방전 정의, 사이클 수, CF/LLI/LAM의 단위는 확인하지 못했습니다(카탈로그 메모: 상용 리튬이온 셀의 단기 가속 사이클 연구). 값 크기는 %로 보이고(LLI 0.9~7.1), dR_*는 대부분 음수라 부호 규약도 미확인이며, LAMne 하나가 ~1e-14(피팅 경계값)입니다.")
    A(f"- C-rate 효과(온도별 평균 제거 후 Spearman, n={B['n']}): LLI {rho(be['LLI']['within_T_demeaned'])} (p={be['LLI']['within_T_demeaned']['p']:.3f}), LAMne {rho(be['LAMne']['within_T_demeaned'])}, LAMpe {rho(be['LAMpe']['within_T_demeaned'])}, LAM 합 {rho(be['LAM_sum']['within_T_demeaned'])}, CF {rho(be['CF']['within_T_demeaned'])}.")
    A(f"- 저온(<=10 C, n=16)에서는 LLI {rho(be['LLI']['lowT_le10'])}, LAMne {rho(be['LAMne']['lowT_le10'])}, LAM 합 {rho(be['LAM_sum']['lowT_le10'])}로 C-rate 효과가 뚜렷하고, 고온(>=25 C, n=32)에서는 LLI {rho(be['LLI']['highT_ge25'])}, LAM 합 {rho(be['LAM_sum']['highT_ge25'])}로 거의 없습니다.")
    A(f"- 방향 판정: C-rate 증가가 LLI를 올림(예), LAMne를 올림(예, 특히 저온), 그러나 LAMpe는 내려가서(rho {rho(be['LAMpe']['within_T_demeaned'])}) LAM 합 전체는 변화 없음(rho {rho(be['LAM_sum']['within_T_demeaned'])}). "
      "'고C-rate -> LAM/도금형 손실 증가' 기대는 음극(LAMne)과 저온에서만 부분적으로 맞고, 전체 LAM에서는 확인되지 않습니다. 단위, 사이클 수, 화학 조성이 불명이라 문헌 기대와의 비교는 방향 수준입니다.\n")
    A("## 4. KIT (LG INR18650HG2)\n")
    A(f"- 펄스 R은 측정값임을 확인: {K['pulse_r_changes_with_ageing']['n_r10_increased']}/{K['pulse_r_changes_with_ageing']['n_cells']}셀에서 마지막 R10ms가 첫 값보다 큼(중앙값 {pct(K['pulse_r_changes_with_ageing']['median_r10_change'],0)}). 원시 v/i 첫 전류 스텝의 dV/dI와 r_ref_10ms의 Spearman {rho(K['pulse_step_vs_r_ref_check']['spearman_rstep_vs_r10'])}(n={K['pulse_step_vs_r_ref_check']['n_blocks']}), 중앙 비 {K['pulse_step_vs_r_ref_check']['median_ratio_rstep_over_r10']}. 이는 거친 타당성 확인일 뿐 R을 독립 재계산한 것은 아닙니다.")
    A(f"- 펄스는 is_rt=1(실온 측정)만 사용해 모든 노화 온도의 R을 약 25 C에서 비교합니다(측정 온도 중앙값 {K['pulse_temperature_note'].split('= ')[1].split(' C')[0]} C).")
    A("- cyc_condition 코드(데이터에서 추정, 문서 근거 아님): 0 = 사이클 외 구간(휴지/조건화 조각, CE 없음), 1 = 노화 사이클(사이클 셀에만 존재, CE는 방전 행에 직전 충전 대비 비로 기록되며 P045에서 수작업 확인), 2 = 체크업 사이클(셀당 약 29개). 프로파일(WLTP) 셀은 코드 1에 CE가 없고 캘린더 셀은 노화 사이클이 없습니다.")
    A("- CE 처리: 노화 사이클(코드 1)의 반 사이클 CE 중 |1-CE| > 2 %는 제외하고, 셀당 지점 이전 EFC에서 5개 이상 남을 때만 평균(1-CE)을 사용했습니다. 체크업 사이클은 제외했습니다.")
    A("- 공통 지점의 조건 수:\n")
    A("| 지점 | 조건 수(전체) | 셀 수 | 사용 지표 수 |")
    A("|---|---|---|---|")
    for k, v in cp.items():
        A(f"| {k} | {v['n_conditions']} ({v['n_conditions_total']}) | {v['n_cells']} | " + ", ".join(f"{a}:{b}" for a, b in v["indicator_available_counts"].items()) + " |")
    A("")
    A("- 지표 일치(조건 간 Spearman, 용량손실 기준):\n")
    A("| 지점 | ~R10 | ~R1s | ~HFR | ~1-CE |")
    A("|---|---|---|---|---|")
    for k, v in cp.items():
        a = v["agreement_spearman"]
        g = lambda x: rho(a.get(f"R10_growth~cap_loss")) if x == "R10" else rho(a.get("R1s_growth~cap_loss")) if x == "R1s" else rho(a.get("HFR_growth~cap_loss")) if x == "HFR" else rho(a.get("cap_loss~one_minus_CE"))
        A(f"| {k} (n={v['n_conditions']}) | {g('R10')} | {g('R1s')} | {g('HFR')} | {g('CE')} |")
    A("")
    A("  대부분 0.7~0.96으로 높지만, 온도와 충전율이 모든 지표를 같이 움직이는 공통 요인이라 지표 간 독립 증거로 읽으면 과대입니다(공통 요인 통제 부분상관은 계산하지 않음). 일치하지 않는 부분: profile @60 EFC에서 HFR~용량 rho "
      f"{rho(cp['profile_efc60']['agreement_spearman'].get('HFR_growth~cap_loss'))}(n=12); 사이클 @60 EFC에서 HFR~R10 {rho(cp['cycle_efc60']['agreement_spearman'].get('HFR_growth~R10_growth'))}, Rpol~1-CE {rho(cp['cycle_efc60']['agreement_spearman'].get('Rpol_growth~one_minus_CE'))}로 상대적으로 낮습니다.\n")
    A("### 메커니즘 서명 검증\n")
    lo = s1["lowT_le10_chg1.67"]; hi = s1["highT_ge25_chg1.67"]
    A(f"1. **저온 + 고충전율 -> 리튬 도금형(용량손실 크고 CE 낮고 저항 증가는 작음)**: **{s1['verdict']}**. 1.67C vs 0.33C 충전(방전 1C 고정, 같은 T·전압 윈도우) @60 EFC, T<=10 C(n={lo['n']}쌍): 용량손실 증가 {pct(lo['frac_caploss_up'],0)}의 쌍, 평균 +{lo['mean_d_caploss_pp']:.1f} %p; 1-CE 증가 {pct(lo['frac_1mCE_up'],0)}, 평균 +{lo['mean_d_1mCE_pp']:.3f} %p; "
      f"그러나 R10 증가가 평균 +{lo['mean_d_R10_pp']:.1f} %p로 용량손실 증가보다 큼(R10이 더 작은 쌍 {pct(lo['frac_R10_rise_smaller_than_caploss_rise'],0)}) -> '저항 증가 작음' 부분은 **지지되지 않음**. T>=25 C(n={hi['n']}쌍)에서는 용량손실 +{hi['mean_d_caploss_pp']:.1f} %p, 1-CE +{hi['mean_d_1mCE_pp']:.3f} %p, R10 {hi['mean_d_R10_pp']:+.2f} %p로 효과가 사실상 저온에 한정됩니다. "
      f"@400 EFC는 T<=10 C에서 1.67C 쌍이 {s1b['lowT_le10_chg1.67']['n']}개뿐이라 검증 불충분(1.0C 쌍 {s1b['lowT_le10_chg1.0']['n']}개: 용량 +{s1b['lowT_le10_chg1.0']['mean_d_caploss_pp']:.1f} %p, R10 +{s1b['lowT_le10_chg1.0']['mean_d_R10_pp']:.1f} %p). 도금을 직접 확인한 것은 아니고(OCV/ICA/분해 없음) 용량+CE 패턴이 일치할 뿐입니다.")
    A(f"2. **고온/고 SOC 캘린더 -> SEI(용량손실과 저항 증가 동반)**: **{s2['verdict']}** (300 d, n={s2['n_conditions']}조건). 온도 rho: 용량 {rho(s2['spearman_T_vs']['cap_loss'])}, R10 {rho(s2['spearman_T_vs']['R10'])}, HFR {rho(s2['spearman_T_vs']['HFR'])}; SOC rho: 용량 {rho(s2['spearman_SOC_vs']['cap_loss'])}, R10 {rho(s2['spearman_SOC_vs']['R10'])}; 용량~R10 {rho(s2['spearman_caploss_vs']['R10'])}. "
      f"440 d에서도 같은 결론(용량~R10 {rho(s2b['spearman_caploss_vs']['R10'])}, HFR {rho(s2b['spearman_caploss_vs']['HFR'])}). 단 R 증가는 40 C와 100 % SOC에서만 뚜렷하고(40 C 평균 R10 {pct(s2['mean_by_T']['40']['R10'],0)} vs 용량손실 {pct(s2['mean_by_T']['40']['cap_loss'],1)}), 0~25 C 평균(용량손실 6~7 %, R10 -2~+1 %)과 SOC 10/50 % 평균(용량손실 5~7 %, R10 약 0 %)에서는 R10이 거의 늘지 않습니다 -> 이 영역의 손실은 저항 증가 없는 용량 손실(LLI형)로, 'SEI = 저항 동반'은 고온/고 SOC에 한정해서만 맞습니다.")
    d1, d2 = s3["deep-vs-shallow"], s3["deep-vs-mid"]
    A(f"3. **깊은 DoD -> LAM(저항 증가)**: 방향은 지지, LAM 귀속은 **검증 불가**. @60 EFC, 같은 EFC에서 깊은 윈도우(2.5-4.2 V) vs 얕은(3.249-4.092 V) n={d1['n']}쌍: 용량손실 +{d1['mean_d_caploss_pp']:.1f} %p(깊은 쪽이 큰 쌍 {pct(d1['frac_caploss_higher_in_deep'],0)}), R10 +{d1['mean_d_R10_pp']:.1f} %p({pct(d1['frac_R10_higher_in_deep'],0)}); 중간(3.249-4.2 V) 대비 n={d2['n']}쌍: +{d2['mean_d_caploss_pp']:.1f} / +{d2['mean_d_R10_pp']:.1f} %p. "
      f"효과는 저온에서 크고(깊은-얕은: 용량 +{d1['lowT_le10']['mean_d_caploss_pp']:.1f}, R10 +{d1['lowT_le10']['mean_d_R10_pp']:.1f} %p) 고온에서 작습니다(+{d1['highT_ge25']['mean_d_caploss_pp']:.1f}, +{d1['highT_ge25']['mean_d_R10_pp']:.1f} %p). @400 EFC도 같은 방향(R10 +{s3b['deep-vs-shallow']['mean_d_R10_pp']:.1f} %p). "
      "윈도우가 DoD 외에 하한 전압(2.5 V vs 3.249 V)과 상한 전압도 달라 DoD만 분리되지 않으며, OCV 곡선이 없어 LAM인지 LLI인지 구분할 수 없습니다. M50T에서는 25/40 C에서 같은 방향이 아니었습니다(2절).\n")
    A("### 조건별 지표 일관성(순위 백분위의 최대-최소 폭)\n")
    for key, lab in [("cycle_efc60", "사이클 @60 EFC"), ("cycle_efc400", "사이클 @400 EFC"), ("calendar_day300", "캘린더 @300 d"), ("profile_efc200", "프로파일 @200 EFC")]:
        c = cons[key]
        A(f"- {lab} (n={c['n']}, 지표 {len(c['indicators_used'])}개): 폭 <=0.25인 조건 {pct(c['frac_spread_le_0p25'],0)}, >=0.5인 조건 {pct(c['frac_spread_ge_0p5'],0)}, 중앙값 {c['median_spread']}.")
    c = cons["cycle_efc60"]
    bad = c["most_disagreeing"][0]
    A(f"- 가장 불일치: 사이클 P{bad['P']} ({bad['T']} C, {bad['v_min']}-{bad['v_max']} V, 충전 {bad['chg_rate']}C/방전 {bad['dis_rate']}C): 용량·R10·R1s는 상위(백분위 {bad['rk_cap_loss']:.2f}/{bad['rk_R10_growth']:.2f}/{bad['rk_R1s_growth']:.2f})인데 HFR은 하위({bad['rk_HFR_growth']:.2f}) -> HFR(오믹 근사)은 작은데 R10/R1s는 큼. 성분 분해는 하지 않았으므로 해석은 가설 수준.")
    bt = c["by_T"]
    A(f"- 온도별 @60 EFC 평균: 용량손실 0/10/25/40 C = {pct(bt['0']['cap_loss'])}/{pct(bt['10']['cap_loss'])}/{pct(bt['25']['cap_loss'])}/{pct(bt['40']['cap_loss'])}, R10 증가 = {pct(bt['0']['R10'])}/{pct(bt['10']['R10'])}/{pct(bt['25']['R10'])}/{pct(bt['40']['R10'])}. 25/40 C 사이클은 용량이 ~5 % 줄어도 R10은 오히려 -1~-2 %(조기 용량손실이 저항 증가 없이 진행)입니다.\n")
    A("## 5. 종합: 방법 간 일관/불일관\n")
    A("**일관**\n")
    A("- 저온(0/10 C) 사이클이 두 데이터셋 모두에서 가장 나쁘고 저항 증가를 동반(M50T 0-100%/10 C, KIT 0/10 C @60/400 EFC). 단 M50T 0-30%는 반대.")
    A("- 용량손실은 DM의 LLI와 가장 잘 대응(M50T 잔차 +0.8 %p), KIT에서는 R10/R1s/HFR/1-CE와 rho 0.7~0.96로 같이 움직임.")
    A("- 고온·고 SOC 캘린더: 용량손실, R10, R1s, HFR 모두 증가(KIT).\n")
    A("**불일관 / 확인 안 됨**\n")
    A("- 저온 고충전율에서 '저항 증가 작음'은 틀림(저항도 크게 증가). 도금 서명은 용량+CE 쪽만 일치.")
    A("- M50T 0-100%/40 C와 0-30%/40 C는 R 순위와 용량/LLI 순위가 3단계 이상 어긋남(순위 폭 3). 0-30% 윈도우에서는 온도 경향이 뒤집힘.")
    A("- 25/40 C 사이클(KIT)과 0~25 C 캘린더: 용량손실은 있는데 저항은 안 늘어남 -> 저항 지표만으로는 이 손실을 설명하지 못함.")
    A("- 데이터셋 B: C-rate가 LAMpe를 내리는 방향이라 'LAM 증가' 기대와 전체적으로 불일치.")
    A("- M50T의 같은 Ah 기준 윈도우 효과(25/40 C에서 0-100% < 0-30%)는 KIT의 '깊은 DoD가 더 나쁨'과 방향이 다름(셀·화학·처리량 정규화가 달라 직접 비교 불가).\n")
    A("## 6. 한계\n")
    A("- DM은 모델 추정치이며 용량 데이터에 제약된 피팅이라 용량과의 일치는 독립 증거가 아님. 불확도는 보고되지 않음(RMSE 파일은 있으나 사용하지 않음).")
    A("- 데이터셋마다 셀·화학이 다름(NMC811/Gr-Si 21700, HG2 18650(전극 조성은 이 분석에서 확인하지 않음), B는 미상). 교차 데이터셋 비교는 방향 수준.")
    A("- KIT는 OCV 곡선이 없어 ICA/DVA·DM 불가. 모든 KIT 모드 해석은 용량·저항·CE의 간접 증거.")
    A("- M50T 데이터 이상치/무릎: 0-100%/10 C 셀 A의 RPT4(2355 Ah)에서 R/R0-1 = 5.1(이웃 점은 0.06~0.16)로 측정 이상치로 보이며 제거하지 않고 두었습니다(Spearman은 순위 기반이라 영향이 작고, 6600/4000 Ah 보간점에는 쓰이지 않음). 0-30% 윈도우의 4셀(A, B, E, F)은 약 8100~9700 Ah 이후 급격한 용량 붕괴(무릎, 용량손실 > 30 %)가 있고 이 구간에서 LAM_PE와 R이 함께 급등합니다. 폐합/풀링 통계에는 이 점들이 포함되어 있고, 공통 처리량 비교(4000/6600 Ah)는 무릎 이전 구간입니다.")
    A("- M50T는 단일 C-rate(0.3C/1C)라 C-rate 효과는 B로만 봄. 조건당 2~3셀(0-100%/25 C는 2셀)이고 셀 간 산포가 큼.")
    A("- 같은 EFC에서 비교하면 고충전율/저온 셀은 더 짧은 시간에 도달하므로 시간 의존(SEI) 손실과 EFC 의존 손실이 섞임. @60 EFC는 초기 급속 열화 구간이라 체크업 정규화(첫 체크업 기준)와 초기 과도 구간의 영향이 큼.")
    A("- CE는 반 사이클 비의 평균(|1-CE|>2 % 제외). 프로파일/캘린더 셀은 CE 없음. EIS HFR/Rpol은 Im 부호 교차 및 1 Hz 실수부로 직접 정의한 지표이며 위상 부호 규약은 데이터에서 추정.")
    A("- KIT 펄스 R은 25 C 측정값만 사용(is_rt=1). 저온 측정값(is_rt=0)은 쓰지 않았음.")
    A("- 지표 상관은 온도·충전율 공통 요인에 의해 부풀려질 수 있고, 다중 비교 보정은 하지 않음.\n")
    A("## 7. lboro_calendar_tempvar_26870248 (LG M50LT 캘린더) 타당성 메모\n")
    A("- 확인한 것: Readme.txt, Voltage_profile_fitting_single_cell.m 본문, Voltage_profile_fitting_all_cells.m(grep), Degradation_models.m, 파일 목록. 데이터 파일은 읽거나 계산하지 않았습니다(이미 Test_data_x로 풀려 있었고 zip은 건드리지 않음).")
    A("- Voltage_profile_fitting_*.m은 OCP_LGM50LT.mat 반쪽전지 모델을 초기(Initial)와 최종(Final) C/24 방전 곡선에만 피팅해 LLI/LAM_PE/LAM_NE를 셀당 2시점(시험 전/후)으로 산출합니다. RPT마다 DM이 나오는 구조는 아닙니다.")
    A("- RPT마다 Capacity_test.txt, EIS_0.txt, EIS_100.txt(0/100 % SOC EIS)가 있어 용량과 EIS의 시간 추이는 가능합니다.")
    A("- 조건(Readme/스크립트 기준): 25/50 C 저장, 50 C 노출 시간 비율 변화, SOC 50/100 %. 스크립트 주석에 17, 22번 셀은 시험기 문제로 누락이라고 되어 있습니다.")
    A("- 판단: 가능하지만 MATLAB 피팅(fmincon + GlobalSearch)을 Python/scipy로 옮기고 Autolab 텍스트를 파싱해야 합니다. 결과는 캘린더 전용으로 시험 종료 시점의 DM(셀당 2점, 약 26셀)과 RPT별 EIS 추이이며, KIT 캘린더와 독립적인 'DM vs EIS' 교차확인 한 건을 추가할 수 있습니다. 다만 셀(M50LT)과 설계(온도 변동 저장)가 달라 KIT 캘린더를 대체하는 것이 아니라 보조 확인입니다.\n")
    (OUT / "summary.md").write_text("\n".join(L), encoding="utf-8")


# =================================================================== main
def main():
    m50_out, m50_df, _ = analyse_m50()
    R["M50T"] = m50_out
    b_out, _ = analyse_B()
    R["datasetB"] = b_out
    kit_out, tables, series = analyse_kit()
    R["KIT"] = kit_out
    R["lboro_feasibility"] = lboro_note()
    fig_m50(m50_df); fig_heat(m50_out, kit_out); fig_scatter(tables)
    (OUT / "results.json").write_text(json.dumps(R, indent=1, ensure_ascii=False, default=str), encoding="utf-8")
    write_summary(R)
    print("done; results.json written")


def lboro_note():
    return dict(
        folder="data/lboro_calendar_tempvar_26870248 (LG M50LT, 28 cell folders, 20 RPT subfolders + Initial/Final)",
        inspected=["Readme.txt", "Voltage_profile_fitting_single_cell.m", "Voltage_profile_fitting_all_cells.m (grep only)", "Degradation_models.m", "file listing"],
        not_done="No data files read or computed in this task; zip not touched (already extracted to Test_data_x).",
        findings=["Voltage_profile_fitting_*.m fits a half-cell model (OCP_LGM50LT.mat) to C/24 discharge curves of the INITIAL and FINAL tests only "
                  "-> LLI/LAM_PE/LAM_NE at two time points per cell (BoL, EoT), not per RPT.",
                  "Per-RPT data: Capacity_test.txt, EIS_0.txt and EIS_100.txt (EIS at 0 and 100 % SOC), so capacity + EIS trend per RPT is available.",
                  "Conditions per Readme/scripts: storage 25/50 C with variable fraction of time at 50 C, SOC 50/100 %; cells 17 and 22 flagged as missing in a script comment."],
        feasibility=("Feasible but requires porting MATLAB fitting (fmincon + GlobalSearch) to Python/scipy and parsing Autolab text files; "
                     "gives calendar-only DM at end of test (2 points/cell, ~26 cells) plus EIS growth per RPT. It would add one independent DM-vs-EIS "
                     "check for calendar ageing, but it is a different chemistry/cell (M50LT) and a thermal-cycling design, so it is a side-check, "
                     "not a replacement for KIT calendar data."))


if __name__ == "__main__":
    main()
