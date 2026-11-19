"""D4: 나이퀴스트에서 R0(고주파 절편)와 아크(R_arc) 증가를 셀별로 뽑아 온도·C-rate별 '저항 지문'을 만든다.

입력: common.load_eis (SOC 50% 스펙트럼), out/c20.csv (check-up 용량), out/checkup_time.csv (check-up 시각)
출력: out/D/D4_resistance.csv (check-up별), D4_cells.csv (셀별 요약), D4_arrhenius.csv, D4_nyquist.png, D4_fingerprint.png
정의:
  R0    = -Im = 0 교차점의 Re (HFR)
  R_arc = R0에서 마지막 골(-Im 극소)까지의 실수축 길이 (모델 프리, prototype eis_features와 같은 규칙)
  ΔR    = 첫 check-up 대비 증가량 [mOhm], 같은 용량 손실(5, 10 %)에서 보간
  rate  = ΔR / 경과 일수,  Ea_hot = 25->45 C 겉보기 활성화 에너지, Ea_cold = 25->-5 C (양수 = 추울수록 빠름)
주의: 2전극 상온 EIS라 아크를 음극/양극으로 나눌 수 없다. 이름은 'R0형/아크형 스트레스'로만 쓴다.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))   # python -I 실행에서도 common 찾기
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import load_eis, OUT, TEMP

plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
D = OUT / "D"; D.mkdir(exist_ok=True)
TK = {"R": 298.15, "C": 268.15, "H": 318.15}; RG = 8.314
COL = {"R": "#009E73", "C": "#0072B2", "H": "#D55E00"}
# README Table 1 프로토콜 (셀 번호 -> 충전 C-rate, 방전, DoD). 셀 1-3이 같은 프로토콜, 3 온도.
PROTO = {1: ("0.2", "C/5", 100), 4: ("1.5", "1.5C", 100), 7: ("2.0", "2C", 100), 10: ("0.2", "1.5C", 100),
         13: ("0.2", "C/5", 50), 16: ("0.2", "1.5C", 50), 19: ("1.5", "drive", 50)}


def proto(c):
    return PROTO[c - (c - 1) % 3]


def eis_feat(f, re, mim):
    o = np.argsort(-f); f, re, mim = f[o], re[o], mim[o]
    cr = np.where((mim[:-1] < 0) & (mim[1:] >= 0))[0]
    if len(cr) == 0: return np.nan, np.nan
    j = cr[0]; w = -mim[j] / (mim[j + 1] - mim[j]); R0 = re[j] + w * (re[j + 1] - re[j])
    vall = []
    for i in range(j + 2, len(f) - 2):
        if f[i] < 0.05: break
        if mim[i] <= mim[i - 1] and mim[i + 1] > mim[i] and mim[i + 2] > mim[i + 1]:
            pk = mim[j:i + 1].max()
            if pk - mim[i] >= 0.08 * max(pk, 1e-9): vall.append(i)
    return R0, (re[vall[-1]] - R0) if vall else np.nan


def chg_cap(c):
    """C/20 충전 용량 (check-up별). 고온 셀은 BOL C/20 방전이 없어 c20.csv에 BOL이 빠지므로 보완용."""
    from common import DATA
    d = pd.read_csv(DATA / f"{c:02d}" / "OCV_wExpansion.csv"); d.columns = ["t", "I", "V", "x", "T", "Q", "cap", "cyc", "ah"]
    ch = d[(d.I > 100) & (d.I < 400)]
    return ch.groupby("ah").cap.max()


cap = pd.read_csv(OUT / "c20.csv"); tim = pd.read_csv(OUT / "checkup_time.csv")
rows, spectra = [], {}
for c in range(1, 22):
    spec = load_eis(c)
    cc = cap[cap.cell == c].sort_values("ah"); tt = tim[tim.cell == c].sort_values("ah")
    if cc.ah.min() > 5:     # BOL 방전 용량 없음 -> BOL 충전 용량 x (첫 공통 check-up의 방전/충전 비)
        qc = chg_cap(c); a1 = cc.ah.iloc[0]; r = cc.cap.iloc[0] / qc.iloc[np.argmin(np.abs(qc.index - a1))]
        est = qc.iloc[0] * r
        if cc.cap.iloc[0] <= est <= 5.15:      # 공칭 5 Ah 셀, 다른 셀 BOL 4.9-5.0 Ah
            cc = pd.concat([pd.DataFrame(dict(cell=[c], ah=[0.0], cap=[est])), cc])
            print(f"셀 {c:02d}: BOL 방전 용량 없음 -> 충전 용량으로 추정 {est:.3f} Ah")
        else:
            print(f"셀 {c:02d}: BOL 추정 {est:.3f} Ah 비정상 -> 첫 공통 check-up({a1:.0f} Ah)을 기준으로 씀 (ΔR 과소 가능)")
    for (ah, soc), a in spec.items():
        if soc != 50: continue
        R0, Ra = eis_feat(a[:, 0], a[:, 3], a[:, 4])
        k = np.argmin(np.abs(cc.ah.values - ah)); kt = np.argmin(np.abs(tt.ah.values - ah))
        if abs(cc.ah.values[k] - ah) > 5: continue          # 같은 check-up의 용량이 없으면 제외
        rows.append(dict(cell=c, T=TEMP[c], ah=ah, cap=cc.cap.values[k], days=tt.hours.values[kt] / 24, R0=R0 * 1e3, Rarc=Ra * 1e3))
        spectra[(c, ah)] = a
R = pd.DataFrame(rows).sort_values(["cell", "ah"]).dropna()
R["loss"] = 100 * (1 - R.cap / R.groupby("cell").cap.transform("first"))
R["dR0"] = R.R0 - R.groupby("cell").R0.transform("first"); R["dRarc"] = R.Rarc - R.groupby("cell").Rarc.transform("first")
R.to_csv(D / "D4_resistance.csv", index=False)

cells = []
for c, g in R.groupby("cell"):
    L = np.maximum.accumulate(g.loss.values); ch, dis, dod = proto(c)
    rec = dict(cell=c, T=TEMP[c], chg=ch, dis=dis, dod=dod, R0_bol=g.R0.iloc[0], Rarc_bol=g.Rarc.iloc[0], loss_end=L.max())
    for lv in (5, 10):
        if L.max() >= lv:
            for nm in ("dR0", "dRarc", "days"):
                rec[f"{nm}_{lv}"] = np.interp(lv, L, g[nm].values)
            rec[f"ratio_{lv}"] = rec[f"dRarc_{lv}"] / rec[f"dR0_{lv}"]
            rec[f"rR0_{lv}"] = rec[f"dR0_{lv}"] / rec[f"days_{lv}"]; rec[f"rRarc_{lv}"] = rec[f"dRarc_{lv}"] / rec[f"days_{lv}"]
    t = g.days.values
    for nm in ("dR0", "dRarc"):                                   # ΔR = k t^z (양수 점만)
        ok = (t > 0) & (g[nm].values > 0.02)
        if ok.sum() >= 3:
            z, lk = np.polyfit(np.log(t[ok]), np.log(g[nm].values[ok]), 1); rec[f"z_{nm}"] = z; rec[f"k_{nm}"] = np.exp(lk)
    cells.append(rec)
C = pd.DataFrame(cells); C["proto"] = C.chg + "|" + C.dis + "|" + C.dod.astype(str)
C.to_csv(D / "D4_cells.csv", index=False)

ea = []
for lv in (5, 10):
    for comp in ("R0", "Rarc"):
        P = C.pivot_table(index="proto", columns="T", values=f"r{comp}_{lv}"); P = P[(P > 0).all(axis=1)]
        hot = np.log(P.H / P.R) * RG / (1 / TK["R"] - 1 / TK["H"]) / 1e3
        cold = np.log(P.C / P.R) * RG / (1 / TK["C"] - 1 / TK["R"]) / 1e3
        ea.append(dict(loss=lv, comp=comp, n_proto=len(P), Ea_hot_med=hot.median(), Ea_hot_min=hot.min(), Ea_hot_max=hot.max(),
                       Ea_cold_med=cold.median(), Ea_cold_min=cold.min(), Ea_cold_max=cold.max(),
                       rate_C=P.C.median(), rate_R=P.R.median(), rate_H=P.H.median()))
    Q = C.pivot_table(index="proto", columns="T", values=f"ratio_{lv}").dropna()
    ea.append(dict(loss=lv, comp="ratio C>R / R>H (paired)", n_proto=len(Q), Ea_hot_med=int((Q.R > Q.H).sum()), Ea_cold_med=int((Q.C > Q.R).sum())))
E = pd.DataFrame(ea); E.to_csv(D / "D4_arrhenius.csv", index=False)
print(C[["cell", "T", "proto", "dR0_5", "dRarc_5", "ratio_5", "dR0_10", "dRarc_10", "ratio_10", "z_dR0", "z_dRarc"]].round(3).sort_values(["proto", "T"]).to_string(index=False))
print(E.round(2).to_string(index=False))

# ---- 그림 1: 같은 프로토콜(0.2C 충전, C/5 방전, 100% DoD), 세 온도의 나이퀴스트 BOL vs ~10 % 손실
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
for a_, c in zip(ax, (2, 1, 3)):
    g = R[R.cell == c]; i10 = (g.loss - 10).abs().idxmin()
    for ah, lab, ls in ((g.ah.iloc[0], "BOL", "--"), (g.loc[i10, "ah"], f"손실 {g.loc[i10, 'loss']:.1f}%", "-")):
        s = spectra[(c, ah)]; m = s[:, 0] >= 0.05
        a_.plot(s[m, 3] * 1e3, s[m, 4] * 1e3, ls, color=COL[TEMP[c]], lw=1.6, label=lab)
    a_.set_title(f"셀 {c:02d} ({ {'C': '-5', 'R': '25', 'H': '45'}[TEMP[c]] } °C)"); a_.set_xlabel("Re Z [mΩ]"); a_.legend(); a_.set_aspect("equal", adjustable="datalim"); a_.grid(alpha=.3)
ax[0].set_ylabel("-Im Z [mΩ]")
fig.suptitle("같은 사이클 조건, 세 온도: SOC 50 % 상온 EIS (f ≥ 0.05 Hz)")
fig.tight_layout(); fig.savefig(D / "D4_nyquist.png", dpi=130); plt.close(fig)

# ---- 그림 2: 지문 = 같은 용량 손실에서 (ΔR0, ΔR_arc)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
for lv, a_ in zip((5, 10), ax):
    for T in "CRH":
        s = C[C["T"] == T]
        a_.scatter(s[f"dR0_{lv}"], s[f"dRarc_{lv}"], color=COL[T], s=40, label={"C": "-5 °C", "R": "25 °C", "H": "45 °C"}[T])
        for _, r in s.iterrows():
            a_.annotate(f"{int(r.cell):02d}", (r[f"dR0_{lv}"], r[f"dRarc_{lv}"]), fontsize=7, xytext=(3, 2), textcoords="offset points")
    lim = np.nanmax(C[[f"dR0_{lv}", f"dRarc_{lv}"]].values) * 1.1
    a_.plot([0, lim], [0, lim], "k:", lw=1); a_.set_xlim(0, lim); a_.set_ylim(0, lim)
    a_.set_xlabel("ΔR0 [mΩ]"); a_.set_ylabel("ΔR_arc [mΩ]"); a_.set_title(f"용량 손실 {lv} %에서"); a_.grid(alpha=.3); a_.legend()
fig.suptitle("점선 위 = 아크형, 아래 = R0형")
fig.tight_layout(); fig.savefig(D / "D4_fingerprint.png", dpi=130); plt.close(fig)
