"""D6: 충전 DVA 모드(LLI, LAM_PE, LAM_NE)와 나이퀴스트 저항(ΔR0, ΔR_arc)을 같은 check-up에서 묶어 온도별로 비교.

입력: out/D/D1_dva.csv (DVA 모드, 기준 = ref check-up), out/D/D4_resistance.csv (R0, R_arc)
기준 통일: 고온 셀은 BOL C/20 방전이 없어 DVA 기준이 첫 공통 check-up이다. 저항 증가도 같은 ref check-up 기준으로 다시 잰다.
          (그래서 고온 셀의 시험 직후 R0 급상승은 여기서 빠진다.)
출력: out/D/D6_equal_loss.csv, D6_group.csv, D6_corr.csv, D6_modes_resistance.png
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import OUT

plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
D = OUT / "D"
V = pd.read_csv(D / "D1_dva.csv", dtype={"cell": str}); V["cell"] = V.cell.astype(int)
V = V[(V.variant == "main") & V.pass_rmse.astype(bool)]
R = pd.read_csv(D / "D4_resistance.csv")
COL = {"C": "#0072B2", "R": "#009E73", "H": "#D55E00"}; LAB = {"C": "−5 °C", "R": "25 °C", "H": "45 °C"}

rows = []
for c, v in V.groupby("cell"):
    r = R[R.cell == c].sort_values("ah")
    if len(r) < 3: continue
    v = v.sort_values("ah")
    ref_ah = v.loc[v.is_ref.astype(bool), "ah"]
    ref_ah = ref_ah.iloc[0] if len(ref_ah) else v.ah.iloc[0]
    k0 = np.argmin(np.abs(r.ah.values - ref_ah))
    if abs(r.ah.values[k0] - ref_ah) > 5: continue
    R0r, Rar = r.R0.values[k0], r.Rarc.values[k0]
    for _, x in v.iterrows():
        k = np.argmin(np.abs(r.ah.values - x.ah))
        if abs(r.ah.values[k] - x.ah) > 5: continue
        rows.append(dict(cell=c, T=x.T_group, ah=x.ah, loss=100 * x.loss_meas, LLI=100 * x.LLI, LAM_PE=100 * x.LAM_PE, LAM_NE=100 * x.LAM_NE,
                         dR0=r.R0.values[k] - R0r, dRarc=r.Rarc.values[k] - Rar, ref_not_bol=bool(x.flag_ref_not_bol)))
J = pd.DataFrame(rows).sort_values(["cell", "ah"])
J = J[J.ah >= J.groupby("cell").ah.transform(lambda s: s[J.loc[s.index, "loss"].abs() < 1e-9].min() if (J.loc[s.index, "loss"].abs() < 1e-9).any() else s.min())]

eq = []
for c, g in J.groupby("cell"):
    L = np.maximum.accumulate(g.loss.values)
    for lv in (5, 10):
        if L.max() < lv: continue
        rec = dict(cell=c, T=g["T"].iloc[0], loss=lv, ref_not_bol=g.ref_not_bol.iloc[0])
        for m in ("LLI", "LAM_PE", "LAM_NE", "dR0", "dRarc"):
            rec[m] = np.interp(lv, L, g[m].values)
        eq.append(rec)
E = pd.DataFrame(eq)
E["ratio"] = E.dRarc / E.dR0
E.to_csv(D / "D6_equal_loss.csv", index=False)
G = E.groupby(["loss", "T"])[["LLI", "LAM_PE", "LAM_NE", "dR0", "dRarc", "ratio"]].median().round(3)
G["n"] = E.groupby(["loss", "T"]).size()
G.to_csv(D / "D6_group.csv")
pd.set_option("display.width", 200)
print("같은 손실에서 온도별 중앙값 (모드 %, 저항 mΩ, ref check-up 기준):"); print(G.to_string())

cr = []
for lv in (5, 10):
    s = E[E.loss == lv]
    for res in ("dR0", "dRarc", "ratio"):
        for m in ("LLI", "LAM_PE", "LAM_NE"):
            r_all, p_all = spearmanr(s[res], s[m], nan_policy="omit")
            rec = dict(loss=lv, resistance=res, mode=m, rho_all=r_all, p_all=p_all, n_all=len(s))
            # 온도 효과를 뺀 값: 온도 그룹 안에서 순위를 매겨 합친 상관 (within-T)
            rw = s.groupby("T")[[res, m]].rank()
            r_w, p_w = spearmanr(rw[res], rw[m], nan_policy="omit")
            rec.update(rho_withinT=r_w, p_withinT=p_w)
            cr.append(rec)
C = pd.DataFrame(cr); C.to_csv(D / "D6_corr.csv", index=False)
print("\n같은 손실에서 저항 ↔ 모드 Spearman (all = 온도 포함, withinT = 온도 그룹 안 순위):")
print(C.round(3).to_string(index=False))

fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
s = E[E.loss == 10]
for a_, m in zip(ax, ("LLI", "LAM_NE", "LAM_PE")):
    for T in "CRH":
        q = s[s["T"] == T]
        a_.scatter(q[m], q.dR0, color=COL[T], marker="o", s=40, label=f"{LAB[T]} ΔR0")
        a_.scatter(q[m], q.dRarc, color=COL[T], marker="^", s=40, facecolors="none", label=f"{LAB[T]} ΔR_arc")
    a_.set_xlabel(f"{m} [%] (손실 10 % 시점)"); a_.set_ylabel("저항 증가 [mΩ]"); a_.grid(alpha=.3)
ax[0].legend(fontsize=7, ncol=2)
fig.suptitle("Mohtat: 같은 용량 손실 10 %에서 DVA 모드와 저항 증가 (ref check-up 기준)")
fig.tight_layout(); fig.savefig(D / "D6_modes_resistance.png", dpi=130); plt.close(fig)
