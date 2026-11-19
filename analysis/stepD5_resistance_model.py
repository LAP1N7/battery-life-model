"""D5: 나이퀴스트 저항 성분(R0, R_arc)의 증가를 '시간항 + 사이클항' 수식으로 피팅하고, 온도별로 어느 항이 우세한지 본다.

수식 (성분 i = R0 또는 R_arc, 단위 mOhm):
    ΔR_i(t, EFC) = a_i(T, SOC) * sqrt(t_day)  +  b_i(T, I) * EFC
    a_i = a_ref * exp(-Ea_a/R (1/T - 1/Tref)) * exp(s * (SOC - 50)/50)              시간항: 보관/고온형 (확산 지배 피막·전해질 반응 가정)
    b_i = b_ref * exp(-Ea_b/R (1/T - 1/Tref)) * I^m      사이클항 (Ea_b < 0 이면 추울수록 빠름)
    (두 갈래 사이클항은 KIT·Mohtat에서 식별 안 됨 -> 한 갈래로 축소. V자는 시간항(고온)+사이클항(저온)이 만든다)
    Tref = 298.15 K. 상수(Ea, s, m)는 메커니즘마다 고정, 온도·SOC·C-rate는 각 항의 비중만 바꾼다.

순서: (1) KIT 달력+사이클 셀로 모든 상수 피팅 (달력 셀이 시간항을 고정해 줌)
      (2) Mohtat에 KIT의 Ea·s·m을 고정하고 크기(a_ref, b_ref)만 피팅 -> 자유 피팅과 비교 (이식 검증)
      (3) 온도별 손실 10 % 시점에서 각 항의 비중
입력: out/D/D4_resistance.csv (Mohtat), --kit <spectrum_features_enriched.csv> (KIT, prototype 결과)
출력: out/D/D5_params.csv, D5_shares.csv, D5_fit.png
"""
import sys, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import OUT

plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
D = OUT / "D"; RG = 8.314e-3; TREF = 298.15
NAMES = ["la_ref", "Ea_a", "s", "lb_ref", "Ea_b", "m"]
LO = [-12, 0, -6, -15, -150, -1]; HI = [5, 150, 6, 2, 150, 3]


def model(p, d):
    la, Ea_a, s, lb, Ea_b, m = p
    x = (1 / d["TK"] - 1 / TREF) / RG
    ta = np.exp(la - Ea_a * x + s * (d["soc"] - 50) / 50) * np.sqrt(d["days"])
    hb = np.exp(lb - Ea_b * x) * np.maximum(d["I"], 1e-3) ** m * d["efc"]
    return ta, hb, 0 * hb


def fit(d, y, fixed=None, p0=None):
    """fixed: {이름: 값} 고정. 셀마다 같은 무게가 되도록 1/sqrt(n_cell) 가중."""
    w = 1 / np.sqrt(d.groupby("cell").cell.transform("size").values)
    free = [n for n in NAMES if not fixed or n not in fixed]
    base = dict(zip(NAMES, p0 if p0 is not None else [-3, 50, 0, -6, 0, 0.5]))
    if fixed: base.update(fixed)

    def full(q):
        b = dict(base); b.update(zip(free, q)); return [b[n] for n in NAMES]

    def res(q):
        ta, hb, cb = model(full(q), d); return w * (ta + hb + cb - y)
    lo = [LO[NAMES.index(n)] for n in free]; hi = [HI[NAMES.index(n)] for n in free]
    best = None
    for Ea0 in (20, 70):
        for Ec0 in (-40, 0, 40):
            st = dict(base); st.update(Ea_a=Ea0 if "Ea_a" in free else base["Ea_a"], Ea_b=Ec0 if "Ea_b" in free else base["Ea_b"])
            q0 = np.clip([st[n] for n in free], np.array(lo) + 1e-6, np.array(hi) - 1e-6)
            r = least_squares(res, q0, bounds=(lo, hi), loss="soft_l1", f_scale=0.3, max_nfev=20000)
            if best is None or r.cost < best.cost: best = r
    p = full(best.x); ta, hb, cb = model(p, d)
    rmse = np.sqrt(np.mean((ta + hb + cb - y) ** 2))
    return p, rmse, best


def kit_table(path):
    k = pd.read_csv(path, low_memory=False); k = k[k.soc == 50].sort_values(["cell", "day"])
    out = []
    for c, g in k.groupby("cell"):
        g = g[g.R0.notna() & g.R_arc_mf.notna() & g.soh.notna()]
        b = g[g.soh >= 0.99].head(2)
        if len(g) < 4 or len(b) == 0: continue
        cyc = g.age_type.iloc[0] == 2
        out.append(pd.DataFrame(dict(cell=c, src="KIT", T=g["T"].iloc[0], TK=g["T"].iloc[0] + 273.15,
                                     I=g.chg_rate.iloc[0] if cyc else 0.0, soc=50.0 if cyc else g.cal_soc.iloc[0],
                                     days=g.day.values - g.day.values[0], efc=(g.efc.values - g.efc.values[0]) if cyc else 0.0,
                                     loss=100 * (1 - g.soh.values), age=g.age_type.iloc[0],
                                     dR0=g.R0.values - b.R0.mean(), dRarc=g.R_arc_mf.values - b.R_arc_mf.mean())))
    K = pd.concat(out); return K[K.age.isin([1, 2])]          # 달력(1)+사이클(2), WLTP(3)는 제외


def mohtat_table():
    R = pd.read_csv(D / "D4_resistance.csv")
    chg = {1: 0.2, 4: 1.5, 7: 2.0, 10: 0.2, 13: 0.2, 16: 0.2, 19: 1.5}
    R["I"] = R.cell.map(lambda c: chg[c - (c - 1) % 3]); R["TK"] = R["T"].map({"R": 298.15, "C": 268.15, "H": 318.15})
    R["soc"] = 50.0; R["efc"] = R.ah / (2 * 5.0); R["src"] = "Mohtat"
    return R[~R.cell.isin([15, 16])]                            # 15: BOL 기준 없음, 16: EIS check-up 2개


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--kit", required=True); a = ap.parse_args()
    K = kit_table(a.kit); M = mohtat_table()
    rows, shares, fits = [], [], {}
    for comp in ("dR0", "dRarc"):
        pk, rk, _ = fit(K, K[comp].values)
        rows.append(dict(data="KIT 자유", comp=comp, rmse=rk, n_cell=K.cell.nunique(), **dict(zip(NAMES, pk))))
        pm, rm, _ = fit(M, M[comp].values, fixed=dict(s=0.0))   # Mohtat은 SOC 조건이 없어 s 고정
        rows.append(dict(data="Mohtat 자유", comp=comp, rmse=rm, n_cell=M.cell.nunique(), **dict(zip(NAMES, pm))))
        fx = {n: pk[NAMES.index(n)] for n in ("Ea_a", "Ea_b", "m")}; fx["s"] = 0.0
        pt, rt, _ = fit(M, M[comp].values, fixed=fx)
        rows.append(dict(data="Mohtat, KIT 상수 고정", comp=comp, rmse=rt, n_cell=M.cell.nunique(), **dict(zip(NAMES, pt))))
        fits[comp] = (pm, pt)
        for lab, p in (("자유", pm), ("KIT 고정", pt)):
            ta, hb, cb = model(p, M); tot = ta + hb + cb
            S = M.assign(time=ta, cyc_hot=hb, cyc_cold=cb)
            for c, g in S.groupby("cell"):
                i = (g.loss - 10).abs().idxmin()
                if abs(g.loss[i] - 10) > 3: continue
                r = g.loc[i]; t = r.time + r.cyc_hot + r.cyc_cold
                shares.append(dict(comp=comp, fit=lab, cell=c, T=r["T"], I=r.I, time=r.time / t, cyc_hot=r.cyc_hot / t, cyc_cold=r.cyc_cold / t))
    P = pd.DataFrame(rows); P.to_csv(D / "D5_params.csv", index=False)
    S = pd.DataFrame(shares); S.to_csv(D / "D5_shares.csv", index=False)
    with pd.option_context("display.width", 250):
        print(P.round(3).to_string(index=False))
        print("\n손실 10 % 시점 항 비중 (중앙값):")
        print(S.groupby(["comp", "fit", "T"])[["time", "cyc_hot", "cyc_cold"]].median().round(2).to_string())

    # 그림: Mohtat 셀별 궤적 vs 자유 피팅 (R0, R_arc)
    col = {"C": "#0072B2", "R": "#009E73", "H": "#D55E00"}
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for a_, comp in zip(ax, ("dR0", "dRarc")):
        pm = fits[comp][0]; ta, hb, cb = model(pm, M); M["_fit"] = ta + hb + cb
        for c, g in M.groupby("cell"):
            T = g["T"].iloc[0]; a_.plot(g.days, g[comp], "o", ms=3, color=col[T], alpha=.6); a_.plot(g.days, g._fit, "-", color=col[T], lw=1)
        a_.set_xlabel("경과 일수"); a_.set_ylabel(f"{'ΔR0' if comp == 'dR0' else 'ΔR_arc'} [mΩ]"); a_.grid(alpha=.3)
        a_.set_title(f"{'R0' if comp == 'dR0' else '아크'}: 점 = 측정, 선 = 수식 (파랑 -5, 초록 25, 주황 45 °C)")
    fig.tight_layout(); fig.savefig(D / "D5_fit.png", dpi=130); plt.close(fig)

