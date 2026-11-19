# -*- coding: utf-8 -*-
"""
stepD0_selftest.py - parameter-recovery test of the degradation-mode fit (dataset independent).

Synthesises charge/discharge curve pairs from KNOWN (Cp, Cn, nLi) -> true LLI / LAM_PE / LAM_NE relative to a synthetic BOL,
adds nuisance (resistance, hysteresis, charge/discharge asymmetry, white + slowly varying noise, optional OCP mismatch),
runs the SAME pipeline as stepD1 (stepD_core.fit_cell + stepD_pipeline.assemble, averaged curve, fixed window)
and reports recovery errors in percentage points.

    python -I -B -X utf8 stepD0_selftest.py [--n 40] [--jobs 6] [--ocp mohtat_peyman --ocp-module my_ocp.py]

Writes out\\D\\D0_selftest.csv (per case), D0_selftest.md (summary), D0_selftest.png.
"""
import os
import sys
import time
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import pandas as pd
from scipy import optimize

import stepD_core as dc
import stepD_pipeline as dp

OUT_DEFAULT = os.path.join(HERE, "out", "D")
TRUE_BOL = dict(Cp=5.6, Cn=6.2, nLi=5.26)       # Mohtat-like BOL values (from the BOL fits); any chemistry works the same way

SCENARIOS = {
    # name: (nuisance scale, ocp mismatch, fit mode)
    "clean": (0.0, 0.0, "auto"),
    "nuisance": (1.0, 0.0, "auto"),
    "nuisance_x2": (2.0, 0.0, "auto"),
    "ocp_mismatch": (1.0, 1.0, "auto"),
    "ocp_mismatch_x3": (1.0, 3.0, "auto"),
    "ocp_mismatch_x6": (1.0, 6.0, "auto"),
    "chg_only_fallback": (1.0, 0.0, "chg"),
}


def mismatch_ocp(ocp, k=1.0):
    """Generate with a slightly different OCP than the one used for fitting (1 % stoichiometry stretch of graphite, 5 mV ripple on NMC)."""
    Un0, Up0 = ocp.Un, ocp.Up
    return dc.OCPPair(name=ocp.name + "_true", Up=lambda y: Up0(y) + 0.005 * k * np.sin(8.0 * np.asarray(y)),
                      Un=lambda x: Un0((1.0 + 0.01 * k) * np.asarray(x)), x_range=ocp.x_range, y_range=ocp.y_range)


def synth_checkup(cell, k, thr, Cp, Cn, nLi, ocp_gen, cfg, scale, rng):
    """One check-up in the standard long-table format. Returns (DataFrame, true capacity)."""
    st = dc.limit_states(Cp, Cn, nLi, ocp_gen, cfg.vmin, 4.19)    # top state: OCV 4.19 V (4.2 V loaded charge end, rested)
    if st is None:
        return None, np.nan
    xb, xt = st
    yt = (nLi - xt * Cn) / Cp
    smax = (xt - xb) * Cn * 1.15          # data beyond the 3.0 V limit are cut below
    s = np.linspace(0.0, smax, 7000)
    ocv = ocp_gen.Up(yt + s / Cp) - ocp_gen.Un(xt - s / Cn)
    # nuisance: ohmic + kinetic rise at low SOC (rises towards large s), hysteresis, asymmetry, noise
    R = rng.uniform(0.005, 0.020) * scale                       # ohm, times 0.25 A
    asym = 1.0 + rng.uniform(-0.3, 0.3) * (scale > 0)
    hys = rng.uniform(0.004, 0.012) * scale                      # V, charge-discharge half-gap total
    cap_win = (xt - xb) * Cn
    rise = 1.0 + 2.0 * np.exp(-(cap_win - s) / (0.08 * cap_win))   # resistance rise at the low-SOC end
    ir_c = 0.25 * R * asym * rise
    ir_d = 0.25 * R / asym * rise
    hshape = 1.0 + 0.4 * np.sin(2 * np.pi * s / cap_win)
    dV_c = ir_c + 0.5 * hys * hshape
    dV_d = ir_d + 0.5 * hys * hshape
    def noise():
        slow = sum(rng.uniform(0.2, 1.0) * 1e-3 * scale * np.sin(2 * np.pi * (rng.uniform(0.3, 3.0) * s / cap_win) + rng.uniform(0, 6.28)) for _ in range(3))
        return slow + rng.normal(0, 0.4e-3 * scale, len(s))
    Vc = ocv + dV_c + noise()
    Vd = ocv - dV_d + noise()
    # charge branch: from V=3.0 loaded (or later, after a rest) up to s = 0 ; discharge: from s = 0 down to 3.0 V loaded
    c_start_V = 3.0 + rng.uniform(0, 0.25) * (scale > 0)
    ic = np.where(Vc <= max(c_start_V, 3.0))[0]
    sc_max = s[ic[0]] if len(ic) else s[-1]
    mc = s <= sc_max
    idd = np.where(Vd <= 3.0)[0]
    sd_max = s[idd[0]] if len(idd) else s[-1]
    md = s <= sd_max
    # charge branch in time order: Q from 0 (bottom) to sc_max (top). CV tail omitted (flat), fine for the test
    sc, vc = s[mc][::-1], Vc[mc][::-1]
    Qc = sc[0] - sc
    Qd = s[md]
    base = dict(cell=cell, checkup_idx=k, ah_throughput=thr)
    df = pd.concat([pd.DataFrame(dict(base, direction="charge", Q_Ah=Qc, V=vc, T=25.0)),
                    pd.DataFrame(dict(base, direction="discharge", Q_Ah=Qd, V=Vd[md], T=25.0))], ignore_index=True)
    return df, float(sd_max)


def draw_modes(rng, kind):
    if kind == "random":
        return dict(LLI=rng.uniform(0.0, 0.25), LAM_PE=rng.uniform(0.0, 0.15), LAM_NE=rng.uniform(0.0, 0.25))
    base = dict(LLI=0.0, LAM_PE=0.0, LAM_NE=0.0)
    if kind == "LLI_only":
        base["LLI"] = 0.10
    elif kind == "LAM_PE_only":
        base["LAM_PE"] = 0.10
    elif kind == "LAM_NE_only":
        base["LAM_NE"] = 0.10
    elif kind == "LLI+LAM_NE":
        base["LLI"], base["LAM_NE"] = 0.10, 0.10
    elif kind == "all":
        base["LLI"], base["LAM_PE"], base["LAM_NE"] = 0.10, 0.05, 0.08
    return base


def build_cases(n_random, seed):
    rng = np.random.default_rng(seed)
    cases = [("random", i) for i in range(n_random)]
    for kind in ("LLI_only", "LAM_PE_only", "LAM_NE_only", "LLI+LAM_NE", "all"):
        cases += [(kind, i) for i in range(4)]
    return cases


def run_scenario(name, n_random, ocp_name, ocp_module, jobs, seed=1):
    scale, mism, mode = SCENARIOS[name]
    ocp = dp.ensure_ocp(ocp_name, ocp_module)
    ocp_gen = mismatch_ocp(ocp, mism) if mism else ocp
    cfg = dc.FitConfig(n_starts=12, n_starts_ref=24)
    rng = np.random.default_rng(seed + hash(name) % 1000 if False else seed + sum(map(ord, name)))
    cases = build_cases(n_random, seed)
    tabs, truth = [], {}
    for ci, (kind, i) in enumerate(cases):
        cell = f"S{ci:03d}"
        m = draw_modes(rng, kind)
        bol, cap0 = synth_checkup(cell, 0, 0.0, TRUE_BOL["Cp"], TRUE_BOL["Cn"], TRUE_BOL["nLi"], ocp_gen, cfg, scale, rng)
        Cp, Cn, nLi = TRUE_BOL["Cp"] * (1 - m["LAM_PE"]), TRUE_BOL["Cn"] * (1 - m["LAM_NE"]), TRUE_BOL["nLi"] * (1 - m["LLI"])
        aged, cap1 = synth_checkup(cell, 1, 500.0, Cp, Cn, nLi, ocp_gen, cfg, scale, rng)
        if bol is None or aged is None:
            continue
        tabs += [bol, aged]
        cap_win_true = dc.capacity_from_modes(TRUE_BOL, m["LLI"], m["LAM_PE"], m["LAM_NE"], ocp_gen, cfg.vmin, 4.19) / dc.cell_capacity(
            TRUE_BOL["Cp"], TRUE_BOL["Cn"], TRUE_BOL["nLi"], ocp_gen, cfg.vmin, 4.19)
        truth[cell] = dict(kind=kind, cap_ratio_true=cap_win_true, cap_ratio_data=cap1 / cap0, **{f"true_{k}": v for k, v in m.items()})
    curves = pd.concat(tabs, ignore_index=True)
    rows = dp.fit_all(curves, ocp_name, cfg, mode, jobs, ocp_module, log=lambda *a: None)
    df = dp.assemble(rows, ocp, cfg, None, "first_avg", variant=name)
    df = df[(df["checkup_idx"] == 1) & df["rmse_mV"].notna()].copy()
    for k in ("kind", "cap_ratio_true", "cap_ratio_data", "true_LLI", "true_LAM_PE", "true_LAM_NE"):
        df[k] = df["cell"].map(lambda c: truth[c][k])
    # in the single-branch fallback the BOL curve is also charge-only, so the reference is a charge-only fit as well
    df["scenario"] = name
    for k in ("LLI", "LAM_PE", "LAM_NE"):
        df["err_" + k] = 100 * (df[k] - df["true_" + k])
    return df


def summarize(df, outdir, ocp_name):
    L = []
    P = L.append
    P("# D0 self-test: parameter recovery (synthetic, known truth)\n")
    P(f"OCP pair for fitting: `{ocp_name}`; BOL truth Cp={TRUE_BOL['Cp']} Ah, Cn={TRUE_BOL['Cn']} Ah, n_Li={TRUE_BOL['nLi']} Ah. Each case = BOL + one aged check-up with independent noise; "
      "modes are recovered relative to the fitted BOL (as in D1). Error = recovered - true, in percentage points (pp). "
      "'random' cases: LLI 0-25 %, LAM_PE 0-15 %, LAM_NE 0-25 % (uniform, independent); plus 5 structured patterns x 4 repeats.\n")
    P("Nuisance (scale 1): resistance 5-20 mohm x 0.25 A with low-SOC rise, charge/discharge asymmetry +-30 %, hysteresis 4-12 mV with SOC modulation, white 0.4 mV + three slow sinusoids ~1 mV, charge start voltage 3.0-3.25 V. "
      "'nuisance_x2' doubles all; 'ocp_mismatch' (x1, x3, x6) generates with a graphite OCP stretched by 1 % (3 %, 6 %) in stoichiometry and a 5 mV (15, 30 mV) NMC ripple but fits with the nominal OCP (this is the structural-error test; real D1 fits have 2-7 mV RMSE); "
      "'chg_only_fallback' fits the charge branch only with a fitted offset.\n")
    rows = []
    for sc, g in df.groupby("scenario", sort=False):
        for lab, gg in (("all", g), ("random only", g[g["kind"] == "random"])):
            r = dict(scenario=sc, subset=lab, n=len(gg), fit_rmse_median_mV=gg["rmse_mV"].median())
            for k in ("LLI", "LAM_PE", "LAM_NE"):
                e = gg["err_" + k]
                r[f"{k}_bias"] = e.mean()
                r[f"{k}_rmse"] = float(np.sqrt((e ** 2).mean()))
                r[f"{k}_p95"] = float(e.abs().quantile(0.95))
                t = gg["true_" + k] * 100
                r[f"{k}_slope"] = float(np.polyfit(t, 100 * gg[k], 1)[0]) if t.std() > 0 else np.nan
            rows.append(r)
    S = pd.DataFrame(rows)
    S.round(2).to_csv(os.path.join(outdir, "D0_selftest_summary.csv"), index=False)
    P("## Recovery error (pp), random cases and all cases\n")
    cols = ["scenario", "subset", "n", "fit_rmse_median_mV", "LLI_bias", "LLI_rmse", "LLI_p95", "LAM_PE_bias", "LAM_PE_rmse", "LAM_PE_p95", "LAM_NE_bias", "LAM_NE_rmse", "LAM_NE_p95"]
    P(_md(S[cols].round(2)))
    P("\nSlope of recovered vs true (1 = perfect, <1 = shrinkage / poor identifiability), random cases:\n")
    P(_md(S[S["subset"] == "random only"][["scenario", "LLI_slope", "LAM_PE_slope", "LAM_NE_slope"]].round(2)))
    P("\n## Identifiable? (random cases, rule: RMSE < 1.5 pp)\n")
    for sc, g in df[df["kind"] == "random"].groupby("scenario", sort=False):
        ok = {k: float(np.sqrt((g["err_" + k] ** 2).mean())) < 1.5 for k in ("LLI", "LAM_PE", "LAM_NE")}
        P(f"- {sc}: " + ", ".join(f"{k} {'identifiable' if v else 'POOR'}" for k, v in ok.items()))
    P("\n## Structured patterns, scenario `nuisance` (mean recovered vs truth, pp)\n")
    g = df[(df["scenario"] == "nuisance") & (df["kind"] != "random")]
    t = g.groupby("kind").agg(true_LLI=("true_LLI", "mean"), LLI=("LLI", "mean"), true_PE=("true_LAM_PE", "mean"), LAM_PE=("LAM_PE", "mean"),
                              true_NE=("true_LAM_NE", "mean"), LAM_NE=("LAM_NE", "mean"))
    P(_md((100 * t).round(1)))
    P("\n## Collinearity and the LAM_PE question\n")
    g = df[(df["scenario"] == "nuisance") & (df["kind"] == "random")]
    c = np.corrcoef(g["err_LLI"], g["err_LAM_NE"])[0, 1]
    c2 = np.corrcoef(g["err_LAM_PE"], g["err_LLI"])[0, 1]
    c3 = np.corrcoef(g["err_LAM_PE"], g["err_LAM_NE"])[0, 1]
    P(f"error correlations (nuisance, random): LLI-LAM_NE {c:+.2f}, LAM_PE-LLI {c2:+.2f}, LAM_PE-LAM_NE {c3:+.2f}. A strong negative/positive value means the fit trades one mode against the other.")
    # when is LAM_PE worst? split by which electrode limits capacity: Li-limited (LLI large vs LAM_PE) vs not
    g = g.copy()
    g["li_limited"] = g["true_LLI"] > g["true_LAM_PE"] + 0.05
    for flag, gg in g.groupby("li_limited"):
        P(f"- LLI-limited (true LLI > LAM_PE + 5 pp) = {flag}: n={len(gg)}, LAM_PE rmse {np.sqrt((gg['err_LAM_PE']**2).mean()):.2f} pp, LLI rmse {np.sqrt((gg['err_LLI']**2).mean()):.2f}, LAM_NE rmse {np.sqrt((gg['err_LAM_NE']**2).mean()):.2f}")
    P("\n## Capacity-from-modes identity check\n")
    P(f"capacity ratio implied by TRUE modes vs ratio in the synthetic data (3.0 V loaded discharge end): median |diff| {100*(df['cap_ratio_true']-df['cap_ratio_data']).abs().median():.2f} pp, "
      f"max {100*(df['cap_ratio_true']-df['cap_ratio_data']).abs().max():.2f} pp (differences come from the loaded 3.0 V cut-off vs the OCV limit; this is the floor of completion criterion C2).")
    with open(os.path.join(outdir, "D0_selftest.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    return S


def _md(df):
    cols = [str(c) for c in df.reset_index().columns]
    d = df.reset_index()
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in d.iterrows():
        out.append("| " + " | ".join(str(v) for v in r.values) + " |")
    return "\n".join(out)


def figure(df, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scs = ["nuisance", "ocp_mismatch_x3"]
    fig, ax = plt.subplots(len(scs), 3, figsize=(12, 7.5))
    for i, sc in enumerate(scs):
        g = df[(df["scenario"] == sc)]
        for j, (k, lab) in enumerate((("LLI", "LLI"), ("LAM_PE", "LAM_PE"), ("LAM_NE", "LAM_NE"))):
            a = ax[i, j]
            a.scatter(100 * g["true_" + k], 100 * g[k], s=14, c=np.where(g["kind"] == "random", "tab:blue", "tab:red"))
            lim = [-2, 28]
            a.plot(lim, lim, "k--", lw=0.8)
            a.set_xlim(lim); a.set_ylim(lim)
            e = g["err_" + k]
            a.set_title(f"{sc}: {lab}  rmse {np.sqrt((e**2).mean()):.2f} pp", fontsize=9)
            a.set_xlabel("true [%]"); a.set_ylabel("recovered [%]")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "D0_selftest.png"), dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=40, help="random cases per scenario")
    ap.add_argument("--jobs", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--ocp", default="mohtat_peyman")
    ap.add_argument("--ocp-module", default=None)
    ap.add_argument("--out-dir", default=OUT_DEFAULT)
    ap.add_argument("--scenarios", default=",".join(SCENARIOS))
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    t0 = time.time()
    dfs = []
    for sc in a.scenarios.split(","):
        print(time.strftime("%H:%M:%S"), "scenario", sc, flush=True)
        dfs.append(run_scenario(sc, a.n, a.ocp, a.ocp_module, a.jobs))
    df = pd.concat(dfs, ignore_index=True)
    keep = ["scenario", "cell", "kind", "true_LLI", "true_LAM_PE", "true_LAM_NE", "LLI", "LAM_PE", "LAM_NE", "err_LLI", "err_LAM_PE", "err_LAM_NE",
            "rmse_mV", "n_starts_converged", "spread_LAM_PE_pp", "sd_nLi_pct", "cap_ratio_true", "cap_ratio_data"]
    df[keep].to_csv(os.path.join(a.out_dir, "D0_selftest.csv"), index=False, float_format="%.5g")
    S = summarize(df, a.out_dir, a.ocp)
    figure(df, a.out_dir)
    print(S[["scenario", "subset", "n", "LLI_rmse", "LAM_PE_rmse", "LAM_NE_rmse"]].round(2).to_string())
    print("done in %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
