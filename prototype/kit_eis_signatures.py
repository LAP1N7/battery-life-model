"""KIT RADAR (LG INR18650HG2): stress condition -> mechanism hypothesis -> EIS/pulse/CE diagnostic signal.

Run:  python -I kit_eis_signatures.py [--rebuild] [--workers N]
Inputs are untrusted CSVs: read with pandas only, never executed.
Reuses mech_crosscheck.py (loaded by path, because `python -I` has no script dir on sys.path):
  kit_series.csv (pulse R per check-up), kit_ce.csv (ageing-cycle CE), kit_checkups.csv (capacity).
Outputs: prototype/results/kit_eis/{spectrum_features.csv, results.json, summary.md, fig_*.png}

Sign convention of the files: Im(file) = -Im(Z_standard) (inductive HF part is negative, the capacitive
arc is positive).  All 'Im' below are file-convention ("-Im Z").
"""
import sys, json, math, pathlib, warnings, importlib.util, itertools
import numpy as np, pandas as pd
from scipy import stats, optimize, signal as spsig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parents[2]   # 레포 밖 작업 폴더(data/, prototype/results/가 있는 곳)
PROTO = ROOT / "prototype"
OUT = PROTO / "results" / "kit_eis"
OUT.mkdir(parents=True, exist_ok=True)
MC = PROTO / "results" / "mech_crosscheck"
EIS_DIR = ROOT / "data" / "kit_radar" / "eis_addendum" / "cell_eisv2_fixed"
KIT = ROOT / "data" / "kit_radar" / "result_v2" / "Unnamed entity"
LIT = ROOT / "data" / "lit_mechanism_matrix.json"
CKUP = PROTO / "results" / "kit_checkups.csv"
REBUILD = "--rebuild" in sys.argv
NW = 12
if "--workers" in sys.argv:
    NW = int(sys.argv[sys.argv.index("--workers") + 1])
PRIMARY_SOC = 50
SOCS = [10, 30, 50, 70, 90]
SOH_PTS = [0.99, 0.98, 0.975, 0.95, 0.90]


def f4(x, nd=4):
    try:
        x = float(x)
    except Exception:
        return None
    return None if (math.isnan(x) or math.isinf(x)) else round(x, nd)


# ============================================================ spectrum analysis
def zmodel(p, w):
    R0, Lu, lR1, lt1, a1, lR2, lt2, a2, lA, n = p
    jw = 1j * w
    return (R0 + jw * Lu * 1e-6 + np.exp(lR1) / (1 + (jw * np.exp(lt1)) ** a1)
            + np.exp(lR2) / (1 + (jw * np.exp(lt2)) ** a2) + np.exp(lA) * jw ** (-n))


LB = np.array([1.0, 0.0, math.log(0.03), math.log(1e-5), 0.55, math.log(0.03), math.log(1e-5), 0.55, -8.0, 0.15])
UB = np.array([200.0, 1000.0, math.log(100.0), math.log(20.0), 1.0, math.log(100.0), math.log(20.0), 1.0, 8.0, 1.0])


def fit_spectrum(f, re, imf, R0_guess):
    """R0 + L + ZARC1 + ZARC2 + CPE tail on all finite points. Z_std = re - 1j*imf."""
    keep = f <= FIT_FMAX  # above ~3 kHz Re(Z) falls below R0 (cable/skin artefact), not modelled
    f, re, imf = f[keep], re[keep], imf[keep]
    w = 2 * math.pi * f
    zd = re - 1j * imf
    wt = 1.0 / np.abs(zd)

    def res(p):
        d = (zmodel(p, w) - zd) * wt
        return np.concatenate([d.real, d.imag])
    best = None
    r0 = min(max(R0_guess if np.isfinite(R0_guess) else float(np.min(re)), 1.5), 190)
    for t1, t2, n0 in itertools.product((2e-3, 1e-2), (0.05, 0.4), (0.4, 0.7)):
        p0 = np.array([r0, 100.0, math.log(1.5), math.log(t1), 0.85, math.log(3.0), math.log(t2), 0.85, math.log(max(re.max() - r0, 1.0) * 0.5), n0])
        p0 = np.clip(p0, LB + 1e-6, UB - 1e-6)
        try:
            s = optimize.least_squares(res, p0, bounds=(LB, UB), x_scale=np.array([5, 50, 1, 1, .3, 1, 1, .3, 1, .3]), max_nfev=300)
        except Exception:
            continue
        if best is None or s.cost < best.cost:
            best = s
    if best is None:
        return None
    p = best.x.copy()
    if p[3] > p[6]:  # order arcs by tau
        p[[2, 3, 4, 5, 6, 7]] = p[[5, 6, 7, 2, 3, 4]]
    rr = math.sqrt(np.mean(best.fun ** 2) * 2)  # rms of |dz|/|z|
    return dict(p=p, rel_rms=rr, ok=bool(best.success))


FIT_FMAX = 3200.0
TAU = np.logspace(-4, 1, 50)  # DRT grid (s)


def drt(f, re, imf, lam=0.03):
    """Tikhonov-regularised non-negative DRT on f<=1.6 kHz (R_inf free, inductive points dropped)."""
    m = f <= 1600
    w = 2 * math.pi * f[m]
    wt = 1.0 / np.hypot(re[m], imf[m])
    x = np.outer(w, TAU)
    Are = 1.0 / (1 + x ** 2)
    Aim = x / (1 + x ** 2)
    A = np.vstack([np.hstack([np.ones((m.sum(), 1)), Are]) * wt[:, None], np.hstack([np.zeros((m.sum(), 1)), Aim]) * wt[:, None]])
    b = np.concatenate([re[m] * wt, imf[m] * wt])
    K = len(TAU)
    D = np.zeros((K - 2, K + 1))
    for i in range(K - 2):
        D[i, 1 + i:1 + i + 3] = [1, -2, 1]
    Aa = np.vstack([A, lam * D])
    ba = np.concatenate([b, np.zeros(K - 2)])
    g, _ = optimize.nnls(Aa, ba, maxiter=5000)
    gam = g[1:]
    fit = A @ g
    rr = math.sqrt(np.mean((fit - b) ** 2) * 2)
    return gam, g[0], rr


def spectrum_features(f, re, imf):
    """model-free features. f descending."""
    out = dict(R0=np.nan, f_apex=np.nan, tau_apex=np.nan, im_apex=np.nan, R_arc_mf=np.nan, Rpol_1Hz=np.nan, lf_angle=np.nan)
    # R0: first Im zero crossing (negative->positive) when scanning down in frequency
    for k in range(len(f) - 1):
        if imf[k] < 0 <= imf[k + 1]:
            wgt = -imf[k] / (imf[k + 1] - imf[k])
            out["R0"] = re[k] + wgt * (re[k + 1] - re[k])
            break
    # apex of the arc: max Im in 1..1000 Hz; parabolic refinement in log f
    band = np.where((f >= 1.0) & (f <= 1000.0))[0]
    if len(band) >= 3 and np.isfinite(out["R0"]):
        k = band[np.argmax(imf[band])]
        fa = f[k]
        if 0 < k < len(f) - 1:
            x = np.log10(f[k - 1:k + 2]); y = imf[k - 1:k + 2]
            try:
                c = np.polyfit(x, y, 2)
                if c[0] < 0:
                    xv = -c[1] / (2 * c[0])
                    if x.min() <= xv <= x.max():
                        fa = 10 ** xv
            except Exception:
                pass
        out.update(f_apex=fa, tau_apex=1 / (2 * math.pi * fa), im_apex=imf[k])
        # R_arc_mf: Re at the local minimum of Im below the apex (f>=0.2 Hz) minus R0; else NaN
        lo = np.where((f < f[k]) & (f >= 0.2))[0]
        if len(lo) >= 2:
            j = lo[np.argmin(imf[lo])]
            if imf[j] < imf[k] * 0.98 and j > k:
                out["R_arc_mf"] = re[j] - out["R0"]
        j1 = int(np.argmin(np.abs(np.log(f / 1.0))))
        if abs(f[j1] - 1.0) < 0.2:
            out["Rpol_1Hz"] = re[j1] - out["R0"]
    lf = f <= 1.0
    if lf.sum() >= 5:
        sl = np.polyfit(re[lf], imf[lf], 1)[0]
        out["lf_angle"] = math.degrees(math.atan(sl))
    return out


def drt_features(gam, rinf):
    g = np.asarray(gam)
    d = dict(drt_rinf=rinf)
    for nm, lo, hi in (("fast", 0, 2e-3), ("mid", 2e-3, 0.1), ("slow", 0.1, 1e9)):
        d["drt_area_" + nm] = float(g[(TAU >= lo) & (TAU < hi)].sum())
    d["drt_area_total"] = float(g.sum())
    pk, pr = spsig.find_peaks(np.concatenate([[0], g, [0]]), prominence=0.05 * max(g.max(), 1e-9))
    pk = pk - 1
    pk = [k for k in pk if 0 <= k < len(g)]
    if pk:
        kk = max(pk, key=lambda k: g[k])
        d["drt_tau_main"] = float(TAU[kk]); d["drt_h_main"] = float(g[kk])
        d["drt_npk"] = len(pk)
    else:
        d["drt_tau_main"] = np.nan; d["drt_h_main"] = np.nan; d["drt_npk"] = 0
    return d


def process_cell(cell):
    """returns list of per-spectrum rows (all valid RT spectra, SOC 10..90)."""
    warnings.filterwarnings("ignore")
    ef = EIS_DIR / f"cell_eisv2_{cell}.csv"
    if not ef.exists() or ef.stat().st_size < 1000:
        return []
    ec = KIT / "cell_eocv2" / f"cell_eocv2_{cell}.csv"
    e = pd.read_csv(ec, sep=";", usecols=["timestamp_s", "total_q_dischg_sum_Ah"])
    t0 = e["timestamp_s"].iloc[0]
    te = e["timestamp_s"].values; efc_e = e["total_q_dischg_sum_Ah"].abs().values / 3.0
    z = pd.read_csv(ef, sep=";", usecols=["timestamp_s", "sd_block_id", "is_rt", "soc_nom", "valid", "t_avg_degC", "freq_Hz", "z_amp_mOhm", "z_ph_deg"])
    nraw = z.groupby("sd_block_id").ngroups
    z = z[(z["is_rt"] == 1) & (z["valid"] == 1) & z["freq_Hz"].notna()]
    rows = []
    for blk, g in z.groupby("sd_block_id", sort=False):
        g = g.dropna(subset=["z_amp_mOhm", "z_ph_deg"]).sort_values("freq_Hz", ascending=False)
        if len(g) < 25:
            continue
        f = g["freq_Hz"].values.astype(float)
        amp = g["z_amp_mOhm"].values; ph = np.deg2rad(g["z_ph_deg"].values)
        re = amp * np.cos(ph); imf = amp * np.sin(ph)
        t = float(g["timestamp_s"].iloc[0])
        row = dict(cell=cell, block=int(blk), soc=int(g["soc_nom"].iloc[0]), t=t, day=(t - t0) / 86400.0,
                   efc=float(np.interp(t, te, efc_e)), T_eis=float(g["t_avg_degC"].mean()), n_freq=len(f))
        row.update(spectrum_features(f, re, imf))
        ft = fit_spectrum(f, re, imf, row["R0"])
        if ft is not None:
            p = ft["p"]
            R1, R2 = math.exp(p[2]), math.exp(p[5])
            row.update(fit_R0=p[0], fit_L_nH=p[1] * 1e3, fit_R1=R1, fit_tau1=math.exp(p[3]), fit_a1=p[4], fit_R2=R2, fit_tau2=math.exp(p[6]),
                       fit_a2=p[7], fit_A=math.exp(p[8]), fit_n=p[9], fit_rel_rms=ft["rel_rms"], fit_conv=int(ft["ok"]),
                       fit_Rarc=R1 + R2, fit_tau_dom=math.exp(p[3]) if R1 >= R2 else math.exp(p[6]))
        if row["soc"] == PRIMARY_SOC or True:
            try:
                gam, rinf, rr = drt(f, re, imf)
                row.update(drt_rel_rms=rr, **drt_features(gam, rinf))
                row["drt_gamma"] = ";".join(f"{v:.4g}" for v in gam)
            except Exception:
                pass
        rows.append(row)
    return rows


def build_features():
    fp = OUT / "spectrum_features.csv"
    if fp.exists() and not REBUILD:
        return pd.read_csv(fp)
    ck = pd.read_csv(CKUP)
    cells = list(ck["cell"].unique())
    from multiprocessing import Pool
    allr = []
    with Pool(NW) as pool:
        for i, r in enumerate(pool.imap_unordered(process_cell, cells, chunksize=2)):
            allr.extend(r)
            if i % 20 == 0:
                print(f"  cells {i+1}/{len(cells)}  spectra {len(allr)}", flush=True)
    df = pd.DataFrame(allr)
    df.to_csv(fp, index=False)
    return df


if __name__ == "__main__" and "--features-only" in sys.argv:
    d = build_features()
    print(len(d))
    sys.exit(0)


# ============================================================ analysis
FIT_TOL = 0.02       # relative rms residual allowed for the equivalent-circuit fit
DRT_TOL = 0.05
CE_TRIM = 0.02
GROUP_ORDER = ["G0", "G1", "G2a", "G2b", "G3", "G5"]
GROUP_LABEL = {
    "G0": "G0 cycle, T<=10C, 0.33C charge (cold, gentle; reference)",
    "G1": "G1 plating-prone: cycle, T<=10C, charge>=1C",
    "G2a": "G2a calendar hot/high-SOC: T>=25C and SOC>=50%",
    "G2b": "G2b calendar mild: rest",
    "G3": "G3 high-T cycling: cycle, T>=25C",
    "G5": "G5 WLTP profile cells (not in the brief; kept as a separate class)",
}
GCOL = {"G0": "#7f7f7f", "G1": "#1f77b4", "G2a": "#d62728", "G2b": "#ff9896", "G3": "#ff7f0e", "G5": "#2ca02c"}
SIGNALS = [  # name, column, label, transform
    ("dR0", "R0", "dR0/R0 (HFR)", "rel"),
    ("dRarc_fit", "fit_Rarc_ok", "dR_arc/R_arc0 (ECM fit R1+R2)", "rel"),
    ("dRarc_mf", "R_arc_mf", "dR_arc/R_arc0 (model-free, Re@Im-min - R0)", "rel"),
    ("dRpol", "Rpol_1Hz", "dRpol(1Hz-R0)/Rpol0", "rel"),
    ("dlogtau_apex", "tau_apex_ok", "dlog(tau apex)", "log"),
    ("dlogtau_fit", "fit_tau_dom_ok", "dlog(tau dominant ECM arc)", "log"),
    ("dLFang", "lf_angle", "d LF angle (deg, <=1 Hz)", "abs"),
    ("dDRT_mid", "drt_area_mid2", "dDRT area 2-100 ms", "rel"),
    ("dDRT_fast", "drt_area_fast2", "dDRT area <2 ms", "rel"),
    ("dDRT_tau", "drt_tau_main2", "dlog(DRT main peak tau)", "log"),
    ("dR10", "R10_growth", "dR(10 ms pulse)/R0", "pulse"),
    ("dR1s", "R1s_growth", "dR(1 s pulse)/R0", "pulse"),
    ("mean_1mCE", "mean_1mCE", "mean (1-CE), ageing cycles", "level"),
]
SIG_NAMES = [s[0] for s in SIGNALS]
SIG_LABEL = {s[0]: s[2] for s in SIGNALS}


def assign_group(r):
    if r.age_type == 1:
        return "G2a" if (r["T"] >= 25 and r.cal_soc >= 50) else "G2b"
    if r.age_type == 3:
        return "G5"
    if r["T"] <= 10:
        return "G1" if r.chg_rate >= 1.0 else "G0"
    return "G3"


def assign_window(r):
    if r.age_type != 2:
        return None
    return {(2.5, 4.2): "deep", (3.249, 4.2): "mid", (3.249, 4.092): "shallow"}.get((round(r.v_min, 3), round(r.v_max, 3)))


def drt_feats2(g):
    """re-derive DRT features from the stored gamma: peak search restricted to tau<=0.5 s and not at the grid edges
    (the low-frequency diffusion tail otherwise always wins at the 10 s edge)."""
    g = np.asarray(g, float)
    d = {}
    for nm, lo, hi in (("fast", 0, 2e-3), ("mid", 2e-3, 0.1), ("slow", 0.1, 0.5), ("tail", 0.5, 1e9)):
        d["drt_area_%s2" % nm] = float(g[(TAU >= lo) & (TAU < hi)].sum())
    k_max = int(np.searchsorted(TAU, 0.5))
    gg = g[:k_max]
    pk, _ = spsig.find_peaks(gg, prominence=0.05 * max(gg.max(), 1e-9)) if gg.max() > 0 else ([], None)
    pk = [k for k in pk if 1 <= k < len(gg) - 1]
    if len(pk):
        kk = max(pk, key=lambda k: gg[k])
        d["drt_tau_main2"] = float(TAU[kk]); d["drt_h_main2"] = float(gg[kk]); d["drt_npk2"] = len(pk)
    else:
        d["drt_tau_main2"] = np.nan; d["drt_h_main2"] = np.nan; d["drt_npk2"] = 0
    return d


def prepare(feat, ck, series, ce):
    info = ck.groupby("cell").first()[["P", "age_type", "T", "cal_soc", "v_min", "v_max", "chg_rate", "dis_rate", "profile"]]
    info["group"] = info.apply(assign_group, axis=1)
    info["window"] = info.apply(assign_window, axis=1)
    feat = feat.copy()
    # SOH from check-up capacity (linear in time, <=3 d outside the range tolerated by clipping)
    soh = np.full(len(feat), np.nan)
    for cell, g in feat.groupby("cell"):
        c = ck[ck.cell == cell].sort_values("day")
        d = g["day"].values
        ok = (d >= c.day.min() - 3) & (d <= c.day.max() + 3)
        v = np.interp(d, c.day.values, c.q.values)
        soh[g.index[ok]] = v[ok[:len(v)]] if False else v[ok]
    feat["soh"] = soh
    # validity masks
    n0 = len(feat)
    ok_fit = (feat.get("fit_conv", 0) == 1) & (feat["fit_rel_rms"] <= FIT_TOL) & (feat["fit_Rarc"] > 0.1) & (feat["fit_Rarc"] < 100)
    feat["fit_L_nH"] = feat["fit_L_nH"] / 1e3  # features CSV stored p[1]*1e3 (unit slip); p[1] is already in nH (1e-6 mOhm s)
    feat["fit_ok"] = ok_fit.fillna(False)
    feat["fit_Rarc_ok"] = np.where(feat.fit_ok, feat["fit_Rarc"], np.nan)
    feat["fit_tau_dom_ok"] = np.where(feat.fit_ok, feat["fit_tau_dom"], np.nan)
    edge = (feat["f_apex"] <= 1.05) | (feat["f_apex"] >= 950)
    feat["apex_interior"] = ~edge & feat["f_apex"].notna()
    feat["tau_apex_ok"] = np.where(feat.apex_interior, feat["tau_apex"], np.nan)
    # DRT re-derivation
    D = []
    for gs, rr in zip(feat["drt_gamma"].values, feat["drt_rel_rms"].values):
        if isinstance(gs, str) and rr <= DRT_TOL:
            D.append(drt_feats2([float(x) for x in gs.split(";")]))
        else:
            D.append({k: np.nan for k in ("drt_area_fast2", "drt_area_mid2", "drt_area_slow2", "drt_area_tail2", "drt_tau_main2", "drt_h_main2", "drt_npk2")})
    feat = pd.concat([feat.reset_index(drop=True), pd.DataFrame(D)], axis=1)
    # pulses (nearest check-up block, within 2 d)
    s = series.sort_values("t")[["cell", "t", "R10_growth", "R1s_growth"]].rename(columns={"t": "t_p"})
    feat = feat.sort_values("t")
    feat = pd.merge_asof(feat, s, left_on="t", right_on="t_p", by="cell", direction="nearest", tolerance=2 * 86400)
    # cumulative mean (1-CE), >=5 trimmed ageing cycles with efc <= spectrum efc
    feat["mean_1mCE"] = np.nan
    cg = {c: g.sort_values("efc") for c, g in ce.groupby("cell")}
    for cell, g in feat.groupby("cell"):
        if cell not in cg:
            continue
        c = cg[cell]
        om = 1 - c["ce"].values
        keep = np.abs(om) <= CE_TRIM
        e_ = c["efc"].values[keep]; o_ = om[keep]
        if len(o_) < 5:
            continue
        cs = np.cumsum(o_)
        idx = np.searchsorted(e_, g["efc"].values, side="right")
        val = np.where(idx >= 5, cs[np.clip(idx - 1, 0, None)] / np.maximum(idx, 1), np.nan)
        feat.loc[g.index, "mean_1mCE"] = val
    feat = feat.sort_values(["cell", "soc", "day"]).reset_index(drop=True)
    feat = feat.join(info[["P", "group", "window", "age_type", "T", "cal_soc", "chg_rate", "dis_rate", "v_min", "v_max"]], on="cell", rsuffix="_c")
    return feat, info


def add_deltas(feat):
    """per (cell, soc): BOL = mean of first two valid values (only if their SOH>=0.99), then deltas."""
    feat = feat.copy()
    for nm, col, lab, kind in SIGNALS:
        if nm != col:
            feat[nm] = np.nan
    for (cell, soc), g in feat.groupby(["cell", "soc"]):
        g = g.sort_values("day")
        for nm, col, lab, kind in SIGNALS:
            if kind in ("pulse", "level"):
                if nm != col:
                    feat.loc[g.index, nm] = g[col].values
                continue
            x = g[col].values.astype(float)
            sh = g["soh"].values
            valid = np.where(np.isfinite(x) & (x > 0 if kind != "abs" else True))[0]
            if len(valid) == 0:
                continue
            first = valid[:2]
            first = [i for i in first if (np.isfinite(sh[i]) and sh[i] >= 0.99) or i == valid[0]]
            base = np.mean(x[first])
            if kind == "rel":
                feat.loc[g.index, nm] = x / base - 1
            elif kind == "log":
                feat.loc[g.index, nm] = np.log(x / base)
            else:
                feat.loc[g.index, nm] = x - base
    return feat


def at_soh(g, col, target, max_gap=0.06):
    """first downward crossing of SOH=target (g sorted by day); linear interpolation of col between the two spectra."""
    sh = g["soh"].values; x = g[col].values
    m = np.isfinite(sh) & np.isfinite(x)
    sh, x = sh[m], x[m]
    if len(sh) < 2:
        return np.nan
    for i in range(1, len(sh)):
        if sh[i - 1] > target >= sh[i] and (sh[i - 1] - sh[i]) <= max_gap:
            w = (sh[i - 1] - target) / (sh[i - 1] - sh[i])
            return x[i - 1] + w * (x[i] - x[i - 1])
    return np.nan


def at_soh_level(g, col, target):
    """slowly varying cumulative quantity (mean 1-CE): value at the first spectrum with SOH <= target (no interpolation, needs finite value)."""
    d = g[np.isfinite(g["soh"].values) & np.isfinite(g[col].values)]
    d = d[d["soh"] <= target]
    return float(d[col].iloc[0]) if len(d) and g["soh"].max() > target else np.nan


def matched_table(feat, soc, targets=SOH_PTS):
    rows = []
    d = feat[feat.soc == soc]
    for cell, g in d.groupby("cell"):
        g = g.sort_values("day")
        r0 = g.iloc[0]
        for tg in targets:
            rec = dict(cell=cell, P=r0["P"], group=r0["group"], window=r0["window"], T=r0["T"], target=tg)
            for nm in SIG_NAMES:
                rec[nm] = at_soh_level(g, nm, tg) if nm == "mean_1mCE" else at_soh(g, nm, tg)
            rows.append(rec)
    return pd.DataFrame(rows)


def cliffs(a, b):
    a = np.asarray(a); b = np.asarray(b)
    if len(a) < 1 or len(b) < 1:
        return np.nan, np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


def cond_level(tab, sig):
    t = tab.dropna(subset=[sig])
    return t.groupby(["P", "group", "target"])[sig].median().reset_index()


def signature_tables(mt):
    """median per group at matched SOH, per %p capacity loss, with n_cells and n_conditions."""
    out = {}
    for tg in SOH_PTS:
        pp = (1 - tg) * 100
        for sig in SIG_NAMES:
            t = mt[mt.target == tg].dropna(subset=[sig])
            for gname in GROUP_ORDER:
                g = t[t.group == gname]
                out[(tg, sig, gname)] = dict(n_cells=int(len(g)), n_cond=int(g.P.nunique()),
                                             median=f4(g[sig].median(), 5) if len(g) else None,
                                             q25=f4(g[sig].quantile(.25), 5) if len(g) else None, q75=f4(g[sig].quantile(.75), 5) if len(g) else None,
                                             per_pp=f4(g[sig].median() / pp, 5) if len(g) else None)
    return out


def separation_tests(mt, min_cond=3):
    """per signal/target: Kruskal-Wallis on condition medians; per group one-vs-rest Cliff's delta (condition level)."""
    res = {}
    for sig in SIG_NAMES:
        for tg in SOH_PTS:
            c = cond_level(mt[mt.target == tg], sig)
            grp = {g_: c[c.group == g_][sig].values for g_ in GROUP_ORDER}
            use = {k: v for k, v in grp.items() if len(v) >= min_cond}
            kw = None
            if len(use) >= 2:
                try:
                    k = stats.kruskal(*use.values()); kw = dict(H=f4(k.statistic, 3), p=f4(k.pvalue, 5), groups=list(use), n={g_: int(len(v)) for g_, v in use.items()})
                except Exception:
                    kw = None
            ovr = {}
            for g_ in GROUP_ORDER:
                a = grp[g_]; b = np.concatenate([grp[x] for x in GROUP_ORDER if x != g_])
                if len(a) >= min_cond and len(b) >= min_cond:
                    d, p = cliffs(a, b)
                    ovr[g_] = dict(delta=f4(d, 3), p=f4(p, 5), n_g=int(len(a)), n_rest=int(len(b)), med_g=f4(np.median(a), 5), med_rest=f4(np.median(b), 5))
                else:
                    ovr[g_] = dict(delta=None, p=None, n_g=int(len(a)), n_rest=int(len(b)))
            res[(sig, tg)] = dict(kw=kw, ovr=ovr)
    return res


def earliest(sep, d_thr=0.5, p_thr=0.01):
    """first SOH point (descending) at which group g separates from the rest (|delta|>=d_thr, p<p_thr) AND is still separated at >=1 later point."""
    out = {}
    for sig in SIG_NAMES:
        for g_ in GROUP_ORDER:
            hits = []
            for tg in SOH_PTS:
                o = sep[(sig, tg)]["ovr"][g_]
                hits.append(o["delta"] is not None and abs(o["delta"]) >= d_thr and o["p"] is not None and o["p"] < p_thr)
            first = None; first_any = None
            for i, h in enumerate(hits):
                if h and first_any is None:
                    first_any = SOH_PTS[i]
                if h and any(hits[i + 1:]):
                    first = SOH_PTS[i]; break
            sign = None
            if first is not None:
                dd = sep[(sig, first)]["ovr"][g_]["delta"]
                sign = "up" if dd > 0 else "down"
            out[(sig, g_)] = dict(earliest_persistent=first, earliest_any=first_any, direction=sign, hits=hits)
    return out


# ------------------------------------------------------------ early-feature classifier
EARLY_FEATS_EIS = ["dR0", "dRarc_fit", "dRarc_mf", "dRpol", "dlogtau_apex", "dlogtau_fit", "dLFang", "dDRT_mid", "dDRT_fast", "dDRT_tau"]
EARLY_FEATS_PULSE = ["dR10", "dR1s"]


def early_cell_table(feat, soc, soh_min=0.975):
    """one row per cell: signals at the last spectrum with SOH>=soh_min (needs >=2 spectra in that range)."""
    rows = []
    d = feat[(feat.soc == soc) & (feat.soh >= soh_min)]
    for cell, g in d.groupby("cell"):
        g = g.sort_values("day")
        if len(g) < 2:
            continue
        # last spectrum, but average the last 2 valid values of every signal to reduce noise
        r = dict(cell=cell, P=g["P"].iloc[0], group=g["group"].iloc[0], soh_last=g["soh"].iloc[-1], n_spec=len(g), day_last=g["day"].iloc[-1])
        for nm in SIG_NAMES:
            v = g[nm].dropna()
            r[nm] = v.iloc[-2:].mean() if len(v) else np.nan
        r["caploss_pp"] = (1 - g["soh"].iloc[-2:].mean()) * 100
        rows.append(r)
    return pd.DataFrame(rows)


def run_logo(X, y, groups, model="logreg", seed=0):
    from sklearn.pipeline import make_pipeline
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.metrics import balanced_accuracy_score, accuracy_score
    steps = [SimpleImputer(strategy="median", keep_empty_features=True), StandardScaler()]
    clf = LogisticRegression(C=0.5, max_iter=500, class_weight="balanced") if model == "logreg" else RandomForestClassifier(
        n_estimators=200, min_samples_leaf=2, class_weight="balanced", random_state=seed, n_jobs=1)
    pipe = make_pipeline(*steps, clf)
    pred = np.empty(len(y), dtype=object)
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        if len(set(y[tr])) < 2:
            pred[te] = y[tr][0]; continue
        pipe.fit(X[tr], y[tr]); pred[te] = pipe.predict(X[te])
    return accuracy_score(y, pred), balanced_accuracy_score(y, pred), pred


def classifier_block(et, feats, label, classes, model="logreg", n_perm=100, seed=1):
    t = et[et.group.isin(classes)].copy()
    t = t[t[feats].notna().any(axis=1)] if feats else t
    X = t[feats].values.astype(float); y = t["group"].values.astype(object); gr = t["P"].values
    if len(set(y)) < 2 or len(t) < 20:
        return dict(label=label, n=int(len(t)))
    acc, bacc, pred = run_logo(X, y, gr, model)
    maj = pd.Series(y).value_counts(normalize=True).max()
    res = dict(label=label, model=model, n_cells=int(len(t)), n_conditions=int(t.P.nunique()), classes={c: int((y == c).sum()) for c in classes},
               features=feats, accuracy=f4(acc, 3), balanced_accuracy=f4(bacc, 3), majority_baseline=f4(maj, 3), chance_balanced=f4(1 / len(set(y)), 3))
    if model == "logreg" and n_perm:
        rng = np.random.default_rng(seed)
        P_ids = np.unique(gr)
        lab_of = {p: y[gr == p][0] for p in P_ids}
        pb = []
        for _ in range(n_perm):
            perm = dict(zip(P_ids, rng.permutation([lab_of[p] for p in P_ids])))
            yp = np.array([perm[p] for p in gr], dtype=object)
            pb.append(run_logo(X, yp, gr, model)[1])
        res["perm_balanced_mean"] = f4(np.mean(pb), 3); res["perm_balanced_p95"] = f4(np.percentile(pb, 95), 3)
        res["perm_p_value"] = f4((np.sum(np.array(pb) >= bacc) + 1) / (n_perm + 1), 3)
        res["note_perm"] = "labels shuffled at condition level (all replicates of a condition share a label), %d permutations" % n_perm
    # per-class recall
    res["recall_by_class"] = {c: f4(np.mean(pred[y == c] == c), 3) for c in classes if (y == c).any()}
    return res


# ------------------------------------------------------------ exponents of capacity loss (literature-expectation check)
def loss_exponents(ck, info):
    rows = []
    for cell, g in ck.groupby("cell"):
        i = info.loc[cell]
        ax = "day" if i.age_type == 1 else "efc"
        g = g[(g[ax] > 1)]
        loss = 1 - g["q"].values
        m = loss > 0.004
        if m.sum() < 5:
            continue
        x = np.log(g[ax].values[m]); yv = np.log(loss[m])
        b = np.polyfit(x, yv, 1)[0]
        rows.append(dict(cell=cell, group=i.group, axis=ax, b=b, n=int(m.sum()), loss_end=float(loss.max())))
    return pd.DataFrame(rows)


# ------------------------------------------------------------ figures
def clip_ylim(ax, vals, lo=2, hi=98, pad=0.1):
    v = vals[np.isfinite(vals)]
    if len(v) > 10:
        a, b = np.percentile(v, [lo, hi]); r = b - a
        ax.set_ylim(a - pad * r, b + pad * r)


def fig_feature_vs_soh(feat, soc, path):
    sigs = ["dR0", "dRarc_fit", "dRarc_mf", "dlogtau_apex", "dLFang", "dR10", "dR1s", "mean_1mCE"]
    fig, axs = plt.subplots(2, 4, figsize=(18, 8.5))
    bins = np.arange(1.0, 0.64, -0.02)
    d = feat[(feat.soc == soc) & feat.soh.notna()]
    for ax, sg in zip(axs.ravel(), sigs):
        for g_ in GROUP_ORDER:
            t = d[(d.group == g_)].dropna(subset=[sg])
            if len(t) < 5:
                continue
            ax.scatter(t.soh, t[sg], s=3, color=GCOL[g_], alpha=0.12, rasterized=True)
            idx = np.digitize(t.soh.values, bins[::-1])
            xs, ms, qa, qb = [], [], [], []
            for k in range(1, len(bins)):
                lo_, hi_ = bins[k], bins[k - 1]
                s_ = t[(t.soh <= hi_) & (t.soh > lo_)]
                if len(s_) >= 6 and s_.P.nunique() >= 3:
                    xs.append((lo_ + hi_) / 2); ms.append(s_[sg].median()); qa.append(s_[sg].quantile(.25)); qb.append(s_[sg].quantile(.75))
            if xs:
                ax.plot(xs, ms, "-o", ms=3, color=GCOL[g_], label=f"{g_} (n={t.cell.nunique()} cells)")
                ax.fill_between(xs, qa, qb, color=GCOL[g_], alpha=0.12, lw=0)
        ax.set_xlim(1.005, 0.68); ax.set_xlabel("check-up SOH (capacity / first check-up)"); ax.set_title(SIG_LABEL[sg], fontsize=9)
        clip_ylim(ax, d[sg].values)
        ax.axhline(0, color="k", lw=0.5)
    axs[0, 0].legend(fontsize=7, loc="best")
    fig.suptitle(f"EIS/pulse/CE signals vs SOH by stress group (EIS at {soc}% SOC, 25C). lines: bin median (bins>=6 spectra, >=3 conditions), band: IQR", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def fig_heatmap(sig_tab, path):
    sigs = [s for s in SIG_NAMES if s != "dDRT_fast"]
    fig, axs = plt.subplots(1, 2, figsize=(17, 4.8))
    for ax, tg in zip(axs, (0.975, 0.95)):
        M = np.full((len(GROUP_ORDER), len(sigs)), np.nan); A = [[""] * len(sigs) for _ in GROUP_ORDER]
        for i, g_ in enumerate(GROUP_ORDER):
            for j, sg in enumerate(sigs):
                e = sig_tab[(tg, sg, g_)]
                if e["per_pp"] is not None and e["n_cond"] >= 2:
                    M[i, j] = e["per_pp"]; A[i][j] = f"{e['per_pp']:.3g}\nn={e['n_cells']}/{e['n_cond']}"
        Z = M / np.nanmax(np.abs(M), axis=0, keepdims=True)
        ax.imshow(Z, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        for i in range(len(GROUP_ORDER)):
            for j in range(len(sigs)):
                ax.text(j, i, A[i][j], ha="center", va="center", fontsize=6.5)
        ax.set_xticks(range(len(sigs))); ax.set_xticklabels([SIG_LABEL[s].split(" (")[0] for s in sigs], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(GROUP_ORDER))); ax.set_yticklabels(GROUP_ORDER)
        ax.set_title(f"median signal change per %-point capacity loss, at SOH {tg} (cells/conditions)", fontsize=9)
    fig.suptitle("Response of each signal per % capacity loss (colour: per-column scaled by max |value|; red = increase)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def fig_early(earl, clf_res, path):
    sigs = SIG_NAMES
    fig, axs = plt.subplots(1, 2, figsize=(17, 6), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axs[0]
    M = np.full((len(sigs), len(GROUP_ORDER)), np.nan); T = [[""] * len(GROUP_ORDER) for _ in sigs]
    for i, sg in enumerate(sigs):
        for j, g_ in enumerate(GROUP_ORDER):
            e = earl[(sg, g_)]
            if e["earliest_persistent"] is not None:
                M[i, j] = (1 - e["earliest_persistent"]) * 100
                T[i][j] = f"{e['earliest_persistent']:g}{'+' if e['direction']=='up' else '-'}"
    ax.imshow(M, cmap="viridis_r", vmin=0, vmax=10, aspect="auto")
    for i in range(len(sigs)):
        for j in range(len(GROUP_ORDER)):
            ax.text(j, i, T[i][j] if T[i][j] else ".", ha="center", va="center", fontsize=8, color="w" if T[i][j] else "0.6")
    ax.set_xticks(range(len(GROUP_ORDER))); ax.set_xticklabels(GROUP_ORDER)
    ax.set_yticks(range(len(sigs))); ax.set_yticklabels([SIG_LABEL[s] for s in sigs], fontsize=7.5)
    ax.set_title("earliest matched-SOH point where group != rest (|Cliff d|>=0.5, p<0.01, persists)\n"
                 "text = SOH, +/- = group higher/lower; '.' = never separated at SOH>=0.90", fontsize=8.5)
    ax = axs[1]
    labs, bacc, base, pk = [], [], [], []
    for r in clf_res:
        if "balanced_accuracy" in r:
            labs.append(f"{r['label']}\n(n={r['n_cells']})"); bacc.append(r["balanced_accuracy"]); base.append(r["chance_balanced"]); pk.append(r.get("perm_balanced_p95", np.nan))
    y = np.arange(len(labs))
    ax.barh(y, bacc, color="#1f77b4", height=0.5, label="balanced acc. (leave-condition-out)")
    ax.plot(base, y, "k|", ms=14, label="chance")
    ax.plot(pk, y, "r|", ms=14, label="permutation 95th pct")
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=7); ax.set_xlim(0, 1); ax.invert_yaxis(); ax.legend(fontsize=7)
    ax.set_title("group classification from early features (SOH>=0.975)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def fig_soc_fit(feat, soc_tab, path):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))
    d0 = feat.sort_values("day").groupby(["cell", "soc"]).first().reset_index()
    for ax, (col, lab) in zip(axs.ravel()[:3], [("R0", "R0 (HFR) at first check-up [mOhm]"), ("R_arc_mf", "R_arc model-free at first check-up [mOhm]"), ("tau_apex_ok", "tau apex [s]")]):
        data = [d0[d0.soc == s][col].dropna().values for s in SOCS]
        ax.boxplot(data, labels=[f"{s}%\nn={len(v)}" for s, v in zip(SOCS, data)], showfliers=False)
        ax.set_title(lab, fontsize=9)
        if col == "tau_apex_ok":
            ax.set_yscale("log")
    ax = axs[1, 0]
    ax.hist(feat[feat.soc == PRIMARY_SOC]["fit_rel_rms"].clip(upper=0.1), bins=60, color="#1f77b4")
    ax.axvline(FIT_TOL, color="r"); ax.set_title("ECM fit: relative rms residual (SOC 50; red = drop threshold)", fontsize=9); ax.set_xlabel("rel. rms")
    # example fits
    ex = feat[(feat.soc == PRIMARY_SOC) & feat.fit_ok].sort_values("day")
    for ax, cond in zip(axs.ravel()[4:], ["first", "aged"]):
        ids = ex.groupby("cell").first().index
        cell = ids[0] if cond == "first" else ex[ex.soh < 0.85].cell.iloc[0]
        row = ex[ex.cell == cell].iloc[0 if cond == "first" else -1]
        zz = pd.read_csv(EIS_DIR / f"cell_eisv2_{cell}.csv", sep=";", usecols=["sd_block_id", "freq_Hz", "z_amp_mOhm", "z_ph_deg"])
        g = zz[zz.sd_block_id == row["block"]].dropna().sort_values("freq_Hz", ascending=False)
        f_ = g.freq_Hz.values; re_ = g.z_amp_mOhm.values * np.cos(np.deg2rad(g.z_ph_deg.values)); im_ = g.z_amp_mOhm.values * np.sin(np.deg2rad(g.z_ph_deg.values))
        p = np.array([row.fit_R0, row.fit_L_nH, math.log(row.fit_R1), math.log(row.fit_tau1), row.fit_a1, math.log(row.fit_R2), math.log(row.fit_tau2), row.fit_a2, math.log(row.fit_A), row.fit_n])
        zm = zmodel(p, 2 * np.pi * f_)
        ax.plot(re_, im_, "o", ms=4, label="data (Im=-ImZ)"); ax.plot(zm.real, -zm.imag, "-", label=f"ECM fit (rel.rms {row.fit_rel_rms:.3f})")
        ax.set_title(f"{cell}, SOH {row.soh:.3f}", fontsize=9); ax.set_xlabel("Re Z [mOhm]"); ax.legend(fontsize=7)
    axs[1, 1].axis("on")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


# ------------------------------------------------------------ misc analyses
def g4_pairs(mt, info):
    """deep vs shallow / deep vs mid windows at the same (T, chg, dis); signal compared at equal SOH (equal capacity loss)."""
    cond = info.groupby("P").first()[["T", "chg_rate", "dis_rate", "window", "age_type"]]
    t = mt.copy()
    t = t.join(cond[["chg_rate", "dis_rate"]], on="P")
    t = t[t.window.notna()]
    out = {}
    for tg in (0.975, 0.95, 0.90):
        for sig in SIG_NAMES:
            c = t[t.target == tg].dropna(subset=[sig]).groupby(["T", "chg_rate", "dis_rate", "window"])[sig].median().unstack("window")
            for other in ("shallow", "mid"):
                if "deep" in c and other in c:
                    d = (c["deep"] - c[other]).dropna()
                    out[(tg, sig, other)] = dict(n_pairs=int(len(d)), median_diff=f4(d.median(), 5) if len(d) else None,
                                                 frac_deep_higher=f4((d > 0).mean(), 3) if len(d) else None,
                                                 p_sign=f4(stats.binomtest(int((d > 0).sum()), len(d), 0.5).pvalue, 4) if len(d) >= 5 else None)
    return out


def soc_sensitivity(feat):
    res = dict(soc_available={int(s): int((feat.soc == s).sum()) for s in sorted(feat.soc.unique())})
    d0 = feat.sort_values("day").groupby(["cell", "soc"]).first().reset_index()
    res["bol_median_by_soc"] = {int(s): dict(n=int(len(g)), R0=f4(g.R0.median(), 3), R_arc_mf=f4(g.R_arc_mf.median(), 3), Rpol_1Hz=f4(g.Rpol_1Hz.median(), 3),
                                              fit_Rarc=f4(g[g.fit_ok].fit_Rarc.median(), 3), tau_apex_ms=f4(g.tau_apex_ok.median() * 1e3, 3),
                                              frac_apex_interior=f4(g.apex_interior.mean(), 3), lf_angle=f4(g.lf_angle.median(), 2),
                                              fit_ok_frac=f4(g.fit_ok.mean(), 3)) for s, g in d0.groupby("soc")}
    tabs = {s: matched_table(feat, s, targets=[0.95]) for s in SOCS}
    res["matched_0.95_median_by_soc_group"] = {}
    for s in SOCS:
        t = tabs[s]
        res["matched_0.95_median_by_soc_group"][int(s)] = {g_: dict(n=int(t[(t.group == g_)].dR0.notna().sum()), dR0=f4(t[t.group == g_].dR0.median(), 4),
                                                                   dRarc_mf=f4(t[t.group == g_].dRarc_mf.median(), 4), dRarc_fit=f4(t[t.group == g_].dRarc_fit.median(), 4),
                                                                   dlogtau=f4(t[t.group == g_].dlogtau_apex.median(), 4)) for g_ in GROUP_ORDER}
    # rank agreement of cell-level matched values (SOH 0.95) between SOC 50 and the other SOCs
    base = tabs[PRIMARY_SOC].set_index("cell")
    agree = {}
    for s in SOCS:
        if s == PRIMARY_SOC:
            continue
        o = tabs[s].set_index("cell")
        agree[int(s)] = {}
        for sg in ("dR0", "dRarc_mf", "dRarc_fit", "dlogtau_apex", "dR10"):
            m = base[sg].notna() & o[sg].reindex(base.index).notna()
            if m.sum() >= 8:
                r = stats.spearmanr(base[sg][m], o[sg].reindex(base.index)[m])
                agree[int(s)][sg] = dict(rho=f4(r.statistic, 3), n=int(m.sum()))
    res["cell_level_rank_agreement_with_soc50_at_soh0.95"] = agree
    return res, tabs


def stability(feat, mt, soc):
    """feature stability: availability, noise (SD of delta in spectra with SOH>=0.99), signal at SOH 0.90 (median over cycle cells, cond-level) and SNR."""
    d = feat[(feat.soc == soc)]
    early = d[(d.soh >= 0.99) & (d.day > 1.5)]
    out = {}
    base_n = len(d)
    for sg in SIG_NAMES:
        e = early[sg].dropna()
        t90 = mt[(mt.target == 0.90)].dropna(subset=[sg])
        t95 = mt[(mt.target == 0.95)].dropna(subset=[sg])
        cyc = t90[t90.group.isin(["G0", "G1", "G3"])]
        mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v))) if len(v) else np.nan
        out[sg] = dict(frac_available=f4(d[sg].notna().mean(), 3), n_early=int(len(e)), noise_sd_early=f4(e.std(), 5) if len(e) > 5 else None,
                       noise_mad_early=f4(mad(e.values), 5) if len(e) > 5 else None,
                       n_cells_at_0p90=int(len(t90)), median_at_0p90_cycle=f4(cyc[sg].median(), 5) if len(cyc) else None,
                       iqr_at_0p90_cycle=f4(cyc[sg].quantile(.75) - cyc[sg].quantile(.25), 5) if len(cyc) else None,
                       median_at_0p95_all=f4(t95[sg].median(), 5) if len(t95) else None, n_cells_at_0p95=int(len(t95)))
        if out[sg]["noise_mad_early"] and out[sg]["median_at_0p90_cycle"] is not None and out[sg]["noise_mad_early"] > 0:
            out[sg]["snr_median_over_noiseMAD_0p90_cycle"] = f4(abs(out[sg]["median_at_0p90_cycle"]) / out[sg]["noise_mad_early"], 2)
    return out


# ------------------------------------------------------------ literature expectations
def lit_check(sig_tab, sep, mt50, expo, g4, early_tab):
    rows = []

    def sg(tg, s, g_, k="median"):
        return sig_tab[(tg, s, g_)][k]
    # plating: CE
    for other in ("G0", "G3"):
        for tg in (0.975, 0.95):
            x = pair_stat(mt50, "mean_1mCE", "G1", other, tg)
            if x:
                rows.append(dict(expectation=f"plating -> cumulative mean (1-CE) higher in G1 than in {other} at SOH {tg}", source="lit matrix #19 (refs 12,22,29)",
                                 measured=f"median G1={x['med_a']} (n={x['n_a']} cond) vs {other}={x['med_b']} (n={x['n_b']}), Cliff d={x['delta']}, p={x['p']}",
                                 verdict=("consistent" if (x["delta"] > 0.3 and x["p"] < 0.05) else ("not consistent" if x["delta"] <= 0 else "weak/n.s."))))
    if len(expo):
        for g_ in GROUP_ORDER:
            e = expo[expo.group == g_]
            if len(e):
                rows.append(dict(expectation=f"capacity-loss time exponent b ({g_}; loss ~ x^b; sqrt-t -> b~0.5, linear -> b~1)", source="to be cited (not in lit matrix)",
                                 measured=f"median b={f4(e.b.median(),2)}, IQR {f4(e.b.quantile(.25),2)}-{f4(e.b.quantile(.75),2)}, n={len(e)} cells (loss>0.4 %, axis={e.axis.iloc[0]})",
                                 verdict=("sqrt-t-like" if abs(e.b.median() - 0.5) < 0.2 else ("linear-like" if abs(e.b.median() - 1) < 0.2 else "neither"))))
    return rows


LIT_MAP = [  # (lit index, comparison, signal, expected direction of A relative to B, note)
    (0, ("G2a", "G2b"), "dR0", "flat_up", "SEI: R0/HFR"), (1, ("G2a", "G2b"), "dRarc_mf", "up", "SEI: mid/high-f arc"),
    (2, ("G2a", "G2b"), "dDRT_fast", "up", "SEI: DRT short-tau peak area"), (3, ("G2a", "G2b"), "dLFang", "flat", "SEI: low-f diffusion"),
    (4, ("G2a", "G2b"), "dR10", "up", "SEI: short pulse"), (5, ("G2a", "G2b"), "dR1s", "up", "SEI: ~1 s pulse (long pulse = 10-30 s in lit)"),
    (13, ("G1", "G0"), "dR0", "ambiguous", "plating: R0"), (14, ("G1", "G0"), "dRarc_mf", "up", "plating: mid-f arc (irreversible growth)"),
    (15, ("G1", "G0"), "dDRT_mid", "ambiguous", "plating: DRT"), (16, ("G1", "G0"), "dLFang", "ambiguous", "plating: low-f"),
    (17, ("G1", "G0"), "dR10", "up", "plating: short pulse"), (18, ("G1", "G0"), "dR1s", "up", "plating: ~1 s pulse"),
    (19, ("G1", "G0"), "mean_1mCE", "up", "plating: CE drop (1-CE up)"),
    (26, ("deep", "shallow"), "dR0", "flat", "LAM_PE: R0"), (27, ("deep", "shallow"), "dRarc_mf", "up", "LAM_PE: mid arc"),
    (30, ("deep", "shallow"), "dR10", "up", "LAM_PE: short pulse"), (31, ("deep", "shallow"), "dR1s", "up", "LAM_PE: long pulse"),
    (32, ("deep", "shallow"), "mean_1mCE", "up", "LAM_PE: CE small decrease (1-CE small up)"),
    (39, ("deep", "shallow"), "dRarc_mf", "up", "LAM_NE: mid arc"), (49, ("deep", "shallow"), "dR0", "up", "contact loss: R0 up"),
]


def pair_stat(mt, sig, ga, gb, tg):
    c = cond_level(mt[mt.target == tg], sig)
    a = c[c.group == ga][sig].values; b = c[c.group == gb][sig].values
    if len(a) < 2 or len(b) < 2:
        return None
    d, p = cliffs(a, b)
    return dict(n_a=int(len(a)), n_b=int(len(b)), med_a=f4(np.median(a), 5), med_b=f4(np.median(b), 5), delta=f4(d, 3), p=f4(p, 4))


def verdict_rule(exp, d, p):
    if exp == "ambiguous":
        return "문헌도 모호: 수치만 보고"
    if d is None:
        return "검정 불가"
    sig = (p is not None and p < 0.05 and abs(d) >= 0.3)
    if d == 0:
        return "차이 없음"
    if exp == "up":
        return "일치(유의)" if (d > 0 and sig) else ("방향만 일치(n.s.)" if d > 0 else ("반대 방향(유의)" if sig else "반대 방향(n.s.)"))
    if exp == "flat_up":
        return "일치(유의한 상승)" if (d > 0 and sig) else ("일치(차이 없음/약한 상승)" if d > -0.3 else "불일치(하락)")
    if exp == "flat":
        return "일치(차이 없음)" if not sig else "불일치(유의한 차이)"
    return "?"


def lit_vs_data(mt, g4):
    try:
        L = json.load(open(LIT, encoding="utf-8"))
    except Exception:
        return []
    rows = []
    for idx, (ga, gb), sig, exp, note in LIT_MAP:
        if idx >= len(L):
            continue
        e = L[idx]
        meas = {}
        if ga == "deep":
            for tg in (0.975, 0.95, 0.90):
                x = g4.get((tg, sig, gb))
                if x and x["n_pairs"]:
                    meas[tg] = dict(n=x["n_pairs"], diff=x["median_diff"], frac=x["frac_deep_higher"], p=x["p_sign"])
            ref = meas.get(0.975) or meas.get(0.95)
            if ref:
                v = " | ".join(f"{tg}: " + verdict_rule(exp, None if m["frac"] is None else (m["frac"] - 0.5) * 2, m["p"]) for tg, m in meas.items())
                txt = "; ".join(f"SOH {tg}: deep-shallow={m['diff']}, deep>shallow in {m['frac']} of {m['n']} pairs (sign p={m['p']})" for tg, m in meas.items())
            else:
                v, txt = "검정 불가", "-"
        else:
            for tg in (0.99, 0.975, 0.95, 0.90):
                x = pair_stat(mt, sig, ga, gb, tg)
                if x:
                    meas[tg] = x
            ref = meas.get(0.95) or meas.get(0.975)
            if ref:
                v = " | ".join(f"{tg}: " + verdict_rule(exp, m["delta"], m["p"]) for tg, m in meas.items())
                txt = "; ".join(f"SOH {tg}: med {ga}={m['med_a']} (n={m['n_a']} cond) vs {gb}={m['med_b']} (n={m['n_b']}), d={m['delta']}, p={m['p']}" for tg, m in meas.items())
            else:
                v, txt = "검정 불가", "-"
        rows.append(dict(lit_idx=idx, mechanism=e["mechanism"], diagnostic=e["diagnostic"], feature=e["feature"], lit_direction=e["direction"], refs=e["refs"],
                         specificity=e["specificity"], compare=f"{ga} vs {gb}", signal=sig, expected=exp, note=note, measured=txt, verdict=v))
    return rows


# ------------------------------------------------------------ report
def md_table(header, rows):
    s = "| " + " | ".join(header) + " |\n|" + "|".join(["---"] * len(header)) + "|\n"
    for r in rows:
        s += "| " + " | ".join("" if v is None else str(v) for v in r) + " |\n"
    return s


def main():
    R = {}
    ck = pd.read_csv(CKUP)
    series = pd.read_csv(MC / "kit_series.csv"); ce = pd.read_csv(MC / "kit_ce.csv")
    raw = build_features()
    feat, info = prepare(raw, ck, series, ce)
    feat = add_deltas(feat)
    feat.drop(columns=["drt_gamma"]).to_csv(OUT / "spectrum_features_enriched.csv", index=False)

    # ---------------- 0. inventory
    inv = dict(n_cells_checkups=int(ck.cell.nunique()), n_cells_with_spectra=int(feat.cell.nunique()), n_spectra_valid_rt=int(len(feat)),
               spectra_by_soc={int(s): int((feat.soc == s).sum()) for s in sorted(feat.soc.unique())},
               cells_with_soc50=int(feat[feat.soc == 50].cell.nunique()),
               spectra_per_cell_soc50=dict(median=f4(feat[feat.soc == 50].groupby("cell").size().median()), min=int(feat[feat.soc == 50].groupby("cell").size().min()),
                                           max=int(feat[feat.soc == 50].groupby("cell").size().max())),
               freq_range_Hz=[f4(0.05), f4(14705.88)], points_per_spectrum=29,
               room_temperature_T_eis=dict(median=f4(feat.T_eis.median(), 2), min=f4(feat.T_eis.min(), 1), max=f4(feat.T_eis.max(), 1)),
               spectra_without_soh=int(feat.soh.isna().sum()),
               group_cells={g_: int((info.group == g_).sum()) for g_ in GROUP_ORDER},
               group_conditions={g_: int(info[info.group == g_].P.nunique()) for g_ in GROUP_ORDER},
               window_cells={str(k): int(v) for k, v in info.window.value_counts().items()})
    R["inventory"] = inv
    qual = {}
    for s in SOCS:
        d = feat[feat.soc == s]
        qual[int(s)] = dict(n=int(len(d)), fit_converged=int((d.fit_conv == 1).sum()), fit_ok=int(d.fit_ok.sum()), fit_dropped=int((~d.fit_ok).sum()),
                            fit_dropped_frac=f4(1 - d.fit_ok.mean(), 4), fit_rel_rms_median_ok=f4(d[d.fit_ok].fit_rel_rms.median(), 4),
                            fit_rel_rms_p95_ok=f4(d[d.fit_ok].fit_rel_rms.quantile(.95), 4),
                            R0_found=int(d.R0.notna().sum()), apex_interior=int(d.apex_interior.sum()), R_arc_mf_found=int(d.R_arc_mf.notna().sum()),
                            lf_angle_found=int(d.lf_angle.notna().sum()), drt_ok=int(d.drt_area_mid2.notna().sum()),
                            drt_peak_found=int(d.drt_tau_main2.notna().sum()))
    R["spectrum_quality"] = qual
    d50 = feat[feat.soc == PRIMARY_SOC]
    R["spectrum_quality_extra"] = dict(
        R0_median_mOhm=f4(d50.R0.median(), 3), R_arc_mf_median=f4(d50.R_arc_mf.median(), 3), fit_Rarc_median=f4(d50[d50.fit_ok].fit_Rarc.median(), 3),
        tau_apex_ms_median=f4(d50.tau_apex_ok.median() * 1e3, 3), lf_angle_median_deg=f4(d50.lf_angle.median(), 1),
        fit_tau1_ms_median=f4(d50[d50.fit_ok].fit_tau1.median() * 1e3, 3), fit_tau2_ms_median=f4(d50[d50.fit_ok].fit_tau2.median() * 1e3, 3),
        fit_params_at_bound=dict(
            alpha1_at_ub=f4((d50[d50.fit_ok].fit_a1 > 0.995).mean(), 3), alpha2_at_ub=f4((d50[d50.fit_ok].fit_a2 > 0.995).mean(), 3),
            R2_at_lb=f4((d50[d50.fit_ok].fit_R2 < 0.035).mean(), 3), R1_at_lb=f4((d50[d50.fit_ok].fit_R1 < 0.035).mean(), 3),
            tau1_at_lb=f4((d50[d50.fit_ok].fit_tau1 < 1.2e-5).mean(), 3), tau2_at_ub=f4((d50[d50.fit_ok].fit_tau2 > 19).mean(), 3),
            L_at_ub=f4((d50[d50.fit_ok].fit_L_nH > 990).mean(), 3)))

    # ---------------- SOC sensitivity
    soc_res, soc_tabs = soc_sensitivity(feat)
    R["soc_sensitivity"] = soc_res

    # ---------------- matched SOH tables (primary SOC)
    mt = soc_tabs_all = matched_table(feat, PRIMARY_SOC)
    mt.to_csv(OUT / "matched_soh_cell_table_soc50.csv", index=False)
    sig_tab = signature_tables(mt)
    sep = separation_tests(mt)
    earl = earliest(sep)
    R["signature_table"] = {f"{tg}|{sig}|{g_}": v for (tg, sig, g_), v in sig_tab.items()}
    R["separation_tests"] = {f"{sig}|{tg}": v for (sig, tg), v in sep.items()}
    R["earliest_separation"] = {f"{sig}|{g_}": v for (sig, g_), v in earl.items()}
    R["reach_counts"] = {str(tg): {g_: dict(n_cells=int(mt[(mt.target == tg) & (mt.group == g_)].dR0.notna().sum()),
                                           n_cells_group=int((info.group == g_).sum()),
                                           n_reached_checkup=int(sum(1 for c in info.index[info.group == g_] if ck[ck.cell == c].q.min() <= tg))) for g_ in GROUP_ORDER} for tg in SOH_PTS}
    R["stability"] = stability(feat, mt, PRIMARY_SOC)
    g4 = g4_pairs(mt, info)
    R["G4_deep_vs_other"] = {f"{tg}|{sig}|{o}": v for (tg, sig, o), v in g4.items()}

    # ---------------- classifier on early features
    et = early_cell_table(feat, PRIMARY_SOC, 0.975)
    et.to_csv(OUT / "early_cell_table_soc50_soh975.csv", index=False)
    cls_all = ["G0", "G1", "G2a", "G2b", "G3", "G5"]
    cyc = ["G0", "G1", "G3"]
    clf = []
    clf.append(classifier_block(et, EARLY_FEATS_EIS, "A EIS only, 6 groups", cls_all))
    clf.append(classifier_block(et, EARLY_FEATS_EIS + EARLY_FEATS_PULSE, "B EIS+pulse, 6 groups", cls_all))
    clf.append(classifier_block(et, ["caploss_pp"], "C capacity loss only, 6 groups", cls_all))
    clf.append(classifier_block(et, EARLY_FEATS_EIS + EARLY_FEATS_PULSE + ["caploss_pp"], "D EIS+pulse+loss, 6 groups", cls_all))
    clf.append(classifier_block(et, EARLY_FEATS_EIS + EARLY_FEATS_PULSE, "E EIS+pulse, cycle G0/G1/G3", cyc))
    clf.append(classifier_block(et, EARLY_FEATS_EIS + EARLY_FEATS_PULSE + ["mean_1mCE"], "F +mean(1-CE), cycle G0/G1/G3", cyc))
    clf.append(classifier_block(et, ["mean_1mCE"], "G mean(1-CE) only, cycle G0/G1/G3", cyc))
    clf.append(classifier_block(et, ["caploss_pp"], "H capacity loss only, cycle G0/G1/G3", cyc))
    clf.append(classifier_block(et, EARLY_FEATS_EIS + EARLY_FEATS_PULSE, "I EIS+pulse, calendar G2a/G2b", ["G2a", "G2b"]))
    clf.append(classifier_block(et, EARLY_FEATS_EIS + EARLY_FEATS_PULSE, "A RF EIS+pulse, 6 groups", cls_all, model="rf", n_perm=0))
    R["early_classifier"] = dict(point="last spectrum(s) with SOH>=0.975 per cell (mean of last 2 valid values), SOC 50, cells with >=2 spectra in range",
                                 n_cells_early=int(len(et)), median_soh_last=f4(et.soh_last.median(), 4), median_caploss_pp=f4(et.caploss_pp.median(), 3),
                                 caploss_pp_by_group={g_: f4(et[et.group == g_].caploss_pp.median(), 3) for g_ in GROUP_ORDER},
                                 results=clf)

    # ---------------- exponents + literature
    expo = loss_exponents(ck, info)
    expo.to_csv(OUT / "loss_exponents.csv", index=False)
    lit = lit_check(sig_tab, sep, mt, expo, g4, et)
    R["literature_check"] = lit
    R["lit_vs_data"] = lit_vs_data(mt, g4)
    if LIT.exists():
        try:
            L = json.load(open(LIT, encoding="utf-8"))
            R["lit_matrix_file"] = dict(exists=True, top_level=list(L.keys())[:30] if isinstance(L, dict) else f"list[{len(L)}]")
        except Exception as ex:
            R["lit_matrix_file"] = dict(exists=True, error=str(ex))
    else:
        R["lit_matrix_file"] = dict(exists=False)

    # ---------------- figures
    fig_feature_vs_soh(feat, PRIMARY_SOC, OUT / "fig_feature_vs_soh.png")
    fig_heatmap(sig_tab, OUT / "fig_signature_heatmap.png")
    fig_early(earl, [c for c in clf if c.get("model") == "logreg"], OUT / "fig_early_separation.png")
    fig_soc_fit(feat, soc_res, OUT / "fig_soc_dependence_and_fit.png")

    with open(OUT / "results.json", "w", encoding="utf-8") as fh:
        json.dump(R, fh, ensure_ascii=False, indent=1, default=lambda o: None if (isinstance(o, float) and math.isnan(o)) else (o.item() if hasattr(o, "item") else str(o)))
    write_summary(R, sig_tab, sep, earl, clf, g4, lit, expo, info)
    print("done")


def write_summary(R, sig_tab, sep, earl, clf, g4, lit, expo, info):
    L = []
    A = L.append
    inv = R["inventory"]; q = R["spectrum_quality"]; ex = R["spectrum_quality_extra"]
    A("# KIT RADAR EIS 시그니처: 스트레스 조건 → 메커니즘 가설 → 진단 신호\n")
    A("자동 생성 (`kit_eis_signatures.py`). 공개 데이터만 사용. 메커니즘 그룹은 **조건 기반 가설 라벨이며 정답(ground truth)이 아님**.\n")
    A("## 1. 데이터와 SOC 선택\n")
    A(f"- 셀 {inv['n_cells_with_spectra']}개, 유효 실온(is_rt=1, valid=1) 스펙트럼 {inv['n_spectra_valid_rt']}개, 스펙트럼당 29점, 0.05 Hz - 14.7 kHz, 측정 온도 중앙값 {inv['room_temperature_T_eis']['median']} C (범위 {inv['room_temperature_T_eis']['min']}-{inv['room_temperature_T_eis']['max']}).")
    A(f"- 존재하는 SOC: {', '.join(str(k)+'% (n='+str(v)+')' for k, v in inv['spectra_by_soc'].items())}. **주 분석은 SOC 50% 하나로 통일** (이전 mech_crosscheck와 동일, 셀당 중앙값 {inv['spectra_per_cell_soc50']['median']}개 스펙트럼). 다른 SOC는 아래 5절에서 민감도로만 비교.")
    A("- 파일의 Im 부호는 표준과 반대(-Im Z): 고주파 유도성 부분이 음수, 용량성 호가 양수. 이하 Im은 파일 부호.")
    A("- 3 kHz 이상에서는 Re(Z)가 R0 아래로 떨어지는 아티팩트(케이블/표피효과 추정)가 있어 ECM 피팅은 f <= 3.2 kHz만 사용 (유도성 L은 자유 파라미터). R0는 Im=0 교차점의 선형 보간.")
    A("- SOH = 체크업 용량 / 첫 체크업 용량을 EIS 시각(day)에 선형 보간 (체크업 범위 밖 3일 초과는 제외, 해당 없는 스펙트럼 %d개)." % inv["spectra_without_soh"])
    A("- 이전 스펙트럼 대비 변화량(Δ)의 BOL = 셀당 처음 2개 유효 스펙트럼 평균(SOH>=0.99일 때). \n")
    A("## 2. 스펙트럼별 특징량과 피팅 품질\n")
    A("특징량: R0(HFR), R_arc(ECM R1+R2 / 모델 프리: Im 극소점의 Re - R0), 아크 정점 주파수 → τ=1/(2π f_apex) (1-1000 Hz 내부 극대, 경계면이면 무효), 저주파(<=1 Hz) 기울기 각도, 자체 구현 Tikhonov-NNLS DRT (피크 위치/면적; 선택 항목).\n")
    A(f"ECM: R0 + L + (R1||CPE1) + (R2||CPE2) + CPE 꼬리(확산), 8개 초기값 중 최소 비용 선택, 가중 1/|Z|. 폐기 기준: 미수렴 또는 상대 rms 잔차 > {FIT_TOL*100:.0f}% 또는 R_arc 범위 밖.\n")
    rows = []
    for s in SOCS:
        e = q[s]
        rows.append([f"{s}%", e["n"], e["fit_ok"], f"{e['fit_dropped']} ({e['fit_dropped_frac']*100:.1f}%)", e["fit_rel_rms_median_ok"], e["fit_rel_rms_p95_ok"], e["R0_found"], e["apex_interior"], e["R_arc_mf_found"], e["lf_angle_found"], e["drt_ok"], e["drt_peak_found"]])
    A(md_table(["SOC", "스펙트럼", "ECM 채택", "ECM 폐기", "rms 중앙", "rms p95", "R0 발견", "아크 정점(내부)", "R_arc 모델프리", "LF 각도", "DRT 채택", "DRT 피크"], rows))
    A(f"\nSOC 50 중앙값: R0 {ex['R0_median_mOhm']} mOhm, R_arc(모델프리) {ex['R_arc_mf_median']}, R_arc(ECM) {ex['fit_Rarc_median']}, τ_apex {ex['tau_apex_ms_median']} ms, LF 각도 {ex['lf_angle_median_deg']} deg; ECM τ1 {ex['fit_tau1_ms_median']} ms, τ2 {ex['fit_tau2_ms_median']} ms. **즉 '아크'는 약 1-3 mOhm(R0의 10-30%)로 작고 모델 파라미터는 경계에 자주 걸림:**")
    A("ECM 채택 스펙트럼 중 경계 파라미터 비율: " + ", ".join(f"{k}={v}" for k, v in ex["fit_params_at_bound"].items()) + ".\n")
    A("## 3. 그룹 정의 (조건 기반 가설)\n")
    rows = [[g_, GROUP_LABEL[g_], inv["group_cells"][g_], inv["group_conditions"][g_]] for g_ in GROUP_ORDER]
    A(md_table(["그룹", "정의", "셀 수", "조건 수"], rows))
    A("\nG4(깊은 창 vs 얕은 창): 2.5-4.2 V(deep), 3.249-4.2 V(mid), 3.249-4.092 V(shallow)를 같은 (T, 충전율, 방전율)에서 짝지어 비교(6절). 사이클 셀 창별 셀 수: " + str(inv["window_cells"]) + ". 그룹 G1/G3/G0는 서로 겹치지 않으며 G2는 캘린더 셀만, G5(WLTP)는 과제 명세 밖이나 참고용으로 포함.\n")
    A("## 4. 동일 SOH에서의 시그니처 (SOC 50%)\n")
    A("각 셀에서 SOH가 목표값을 처음 아래로 가로지르는 두 스펙트럼 사이를 선형 보간(간격 > 0.06이면 제외). 표는 셀 수준 중앙값과 **%p 용량 손실당 신호 변화**(중앙값 / 손실 %p), n = 셀 수 / 조건 수. 검정은 복제 셀의 유사반복을 피하려고 **조건 수준(조건별 중앙값)** 에서 수행.\n")
    rc = R["reach_counts"]
    rows = [[tg] + [f"{rc[str(tg)][g_]['n_cells']} / {rc[str(tg)][g_]['n_reached_checkup']} / {rc[str(tg)][g_]['n_cells_group']}" for g_ in GROUP_ORDER] for tg in SOH_PTS]
    A("매칭 가능 셀 수 / 체크업 용량이 해당 SOH에 도달한 셀 수 / 그룹 전체 (매칭 실패 = 체크업 사이 SOH 하락 간격 > 0.06이라 보간 제외, 또는 EIS 값 결측. 빠르게 열화하는 셀일수록 체크업 간격이 커서 빠지므로 G1 등은 선택 편향 가능):\n")
    A(md_table(["SOH"] + GROUP_ORDER, rows))
    for tg in (0.975, 0.95, 0.90):
        A(f"\n### SOH {tg} (용량 손실 {(1-tg)*100:.1f} %p): 신호/%p, [n 셀/조건]\n")
        rows = []
        for sig in SIG_NAMES:
            row = [SIG_LABEL[sig]]
            for g_ in GROUP_ORDER:
                e = sig_tab[(tg, sig, g_)]
                row.append("-" if e["per_pp"] is None else f"{e['per_pp']:.4g} [{e['n_cells']}/{e['n_cond']}]")
            rows.append(row)
        A(md_table(["신호"] + GROUP_ORDER, rows))
    A("\n(`mean_1mCE`는 손실당 값이 아니라 누적 평균 1-CE의 절대 수준을 손실 %p로 나눈 값이며 사이클 셀만 존재. 펄스 ΔR은 체크업 시점 RT 펄스 블록의 최근접 값.)\n")
    A("## 5. 어느 신호가 먼저 분리되는가\n")
    A("그룹 vs 나머지 전체, 조건 수준 Cliff's delta |d|>=0.5 및 Mann-Whitney p<0.01, 그리고 더 낮은 SOH 점에서도 유지될 때 '지속 분리'로 인정하고 그 첫 SOH를 표기 (+ = 그룹이 더 큼, - = 더 작음, . = 분리 없음). 신호 13개 x 그룹 6개 x SOH 5점 = 390 검정이므로 p<0.01에서도 우연 일치가 몇 개 있을 수 있음.\n")
    rows = []
    for sig in SIG_NAMES:
        row = [SIG_LABEL[sig]]
        for g_ in GROUP_ORDER:
            e = earl[(sig, g_)]
            row.append(".") if e["earliest_persistent"] is None else row.append(f"{e['earliest_persistent']:g}{'+' if e['direction']=='up' else '-'}")
        rows.append(row)
    A(md_table(["신호"] + GROUP_ORDER, rows))
    A("\nKruskal-Wallis (조건 수준, 그룹당 조건 >=3인 그룹만): p 값\n")
    rows = []
    for sig in SIG_NAMES:
        row = [SIG_LABEL[sig]]
        for tg in SOH_PTS:
            k = sep[(sig, tg)]["kw"]
            row.append("-" if k is None else f"{k['p']:.3g} (k={len(k['groups'])})")
        rows.append(row)
    A(md_table(["신호"] + [f"SOH {t}" for t in SOH_PTS], rows))
    A("\n### 조기 특징량(SOH>=0.975)만으로 그룹 분류, 조건 단위 leave-one-condition-out\n")
    et = R["early_classifier"]
    A(f"셀당 SOH>=0.975 구간의 마지막 2개 유효값 평균, 셀 {et['n_cells_early']}개 (구간 내 스펙트럼 >=2개), 마지막 SOH 중앙값 {et['median_soh_last']}, 손실 중앙값 {et['median_caploss_pp']} %p, 그룹별 손실 중앙값 {et['caploss_pp_by_group']}. 그룹 간 도달 손실이 다르므로 용량 손실만 쓴 대조군(C, H)을 같이 표시.\n")
    rows = []
    for c in clf:
        if "balanced_accuracy" not in c:
            rows.append([c.get("label"), c.get("n"), "-", "-", "-", "-", "-", "-"]); continue
        rows.append([c["label"], f"{c['n_cells']}/{c['n_conditions']}", c["model"], c["accuracy"], c["balanced_accuracy"], c["chance_balanced"], c.get("perm_balanced_p95", "-"), c.get("perm_p_value", "-")])
    A(md_table(["모델/특징량", "셀/조건", "모델", "정확도", "균형 정확도", "우연 수준", "순열 95%", "순열 p"], rows))
    A("\n클래스별 재현율: " + "; ".join(f"{c['label']}: {c['recall_by_class']}" for c in clf if "recall_by_class" in c and c["label"].startswith(("B", "F"))) + "\n")
    A("## 6. G4 깊은 창 vs 얕은/중간 창 (같은 SOH에서의 신호 차이, deep - other)\n")
    rows = []
    for tg in (0.975, 0.95, 0.90):
        for sig in ("dR0", "dRarc_mf", "dRarc_fit", "dlogtau_apex", "dLFang", "dR10", "dR1s", "mean_1mCE"):
            for other in ("shallow", "mid"):
                e = g4.get((tg, sig, other))
                if e and e["n_pairs"]:
                    rows.append([tg, SIG_LABEL[sig], other, e["n_pairs"], e["median_diff"], e["frac_deep_higher"], e["p_sign"]])
    A(md_table(["SOH", "신호", "비교", "쌍 수", "중앙 차이", "deep가 클 비율", "부호검정 p"], rows))
    A("\n교란: 창은 DoD 외에 하한 컷오프(2.5 vs 3.249 V)와 상한(4.2 vs 4.092 V)도 다르므로 DoD만 분리되지 않음.\n")
    A("## 7. 신호 안정성 (SOC 50)\n")
    rows = []
    for sig in SIG_NAMES:
        s = R["stability"][sig]
        rows.append([SIG_LABEL[sig], s["frac_available"], s["noise_mad_early"], s["median_at_0p90_cycle"], s.get("snr_median_over_noiseMAD_0p90_cycle"), s["iqr_at_0p90_cycle"], s["n_cells_at_0p90"]])
    A(md_table(["신호", "유효 비율", "초기 노이즈(MAD, SOH>=0.99)", "SOH 0.90 중앙값(사이클 셀)", "|중앙값|/노이즈", "IQR(0.90, 사이클)", "n 셀(0.90)"], rows))
    A("\n## 8. SOC 의존성 (SOC 50 vs 10/30/70/90)\n")
    rows = [[f"{s}%", v["n"], v["R0"], v["R_arc_mf"], v["Rpol_1Hz"], v["fit_Rarc"], v["tau_apex_ms"], v["frac_apex_interior"], v["lf_angle"], v["fit_ok_frac"]] for s, v in R["soc_sensitivity"]["bol_median_by_soc"].items()]
    A(md_table(["SOC", "n 셀(BOL)", "R0", "R_arc 모델프리", "Rpol(1Hz-R0)", "R_arc ECM", "τ apex [ms]", "정점 내부 비율", "LF 각도", "ECM 채택 비율"], rows))
    A("\nSOH 0.95에서 셀 수준 순위 상관 (SOC 50 대비):\n")
    rows = []
    for s, d in R["soc_sensitivity"]["cell_level_rank_agreement_with_soc50_at_soh0.95"].items():
        rows.append([f"{s}%"] + [f"{d[k]['rho']} (n={d[k]['n']})" if k in d else "-" for k in ("dR0", "dRarc_mf", "dRarc_fit", "dlogtau_apex", "dR10")])
    A(md_table(["SOC", "dR0", "dRarc_mf", "dRarc_fit", "dlogtau_apex", "dR10(참고)"], rows))
    A("\n그룹별 SOH 0.95 중앙값 (SOC별, dR0 / dRarc_mf):\n")
    rows = []
    for g_ in GROUP_ORDER:
        row = [g_]
        for s in SOCS:
            e = R["soc_sensitivity"]["matched_0.95_median_by_soc_group"][s][g_]
            row.append(f"{e['dR0']} / {e['dRarc_mf']} (n={e['n']})")
        rows.append(row)
    A(md_table(["그룹"] + [f"SOC {s}%" for s in SOCS], rows))
    A("\n## 9. 문헌 기대와의 비교\n")
    if R["lit_matrix_file"]["exists"]:
        A(f"`lit_mechanism_matrix.json` 존재함: 최상위 키 {R['lit_matrix_file'].get('top_level')}. 아래 표의 '기대'는 이 스크립트에 하드코딩된 표준 기대이며, 파일과의 항목별 대조는 해석 절 참조.\n")
    else:
        A("`lit_mechanism_matrix.json`이 완료 시점에 없어서 표준 기대를 사용하고 모두 **'to be cited'** 로 표시 (출처 미확정).\n")
    rows = [[r["expectation"], r["source"], r["measured"], r["verdict"]] for r in lit]
    A(md_table(["기대", "출처", "측정", "판정"], rows))
    if R.get("lit_vs_data"):
        A("\n### 문헌 매트릭스(`lit_mechanism_matrix.json`, 71행) 항목과 데이터 대조\n")
        A("비교 대비: 플래팅 = G1 vs G0 (같은 저온, 충전율만 다름), SEI = G2a vs G2b (캘린더 고온/고SOC vs 온건), LAM = deep vs shallow 창 (같은 T/충전율/방전율 쌍, 같은 SOH). 판정은 조건 수준 Cliff d>=0.3 & p<0.05를 '유의'로 봄. 문헌 #은 매트릭스의 0-기반 행 번호, refs는 그 파일의 참고문헌 번호.\n")
        rows = [[r["lit_idx"], r["mechanism"], f"{r['diagnostic']} / {r['feature']}", r["lit_direction"], r["refs"], r["specificity"], f"{r['compare']}: {r['signal']}", r["measured"], r["verdict"]] for r in R["lit_vs_data"]]
        A(md_table(["#", "메커니즘", "진단/특징", "문헌 방향", "refs", "문헌 특이성", "데이터 비교", "측정", "판정"], rows))
    A("\n## 10. 한계\n")
    A("- 그룹은 조건 기반 가설일 뿐 메커니즘 직접 측정이 아님. 같은 조건의 복제 셀은 독립이 아니므로 검정은 조건 수준, 그래도 그룹 간 교란(예: G3는 창 3종 x 충전율 3종 혼합, 도달 SOH 상이)이 큼.\n- EIS는 25 C에서만, SOC 50 고정; 스트레스 온도에서의 임피던스가 아님.\n- 'R_arc'는 이 셀에서 수 mOhm로 작고 ECM 파라미터가 경계에 걸리는 경우가 많아 절대 의미를 부여하지 않고 변화량만 사용. DRT는 정규화 강도 고정(λ=0.03), 선택 항목이며 피크 위치는 격자 해상도(약 10점/decade)에 의존.\n- 매칭 SOH 점은 체크업 간격(수일-수십일)에 의존한 선형 보간.\n")
    interp = OUT / "interpretation_ko.md"
    if interp.exists():
        A("\n" + interp.read_text(encoding="utf-8"))
    (OUT / "summary.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
