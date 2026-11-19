"""2단계: H-alpha(R_int-SOC 모양 vs dV/dQ 피크 이동)와 HA(같은 SOH, 다른 경로) 검정 -> out/checkups.csv"""
import pickle, numpy as np, pandas as pd
from scipy.signal import find_peaks
from scipy.stats import spearmanr
from common import *

E = pd.read_csv(OUT/"eis_fits.csv"); C = pd.read_csv(OUT/"c20.csv"); CV = pickle.load(open(OUT/"curves.pkl", "rb"))
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
rows = []
for (cell, ah), g in E.groupby(["cell", "ah"]):
    c = C[C.cell == cell]; i = np.argmin(abs(c.ah.values - ah))
    if abs(c.ah.values[i] - ah) > 1: continue
    m = c.iloc[i]; g = g.set_index("soc")
    grid, Vg, dv = CV[cell][m.ah]
    w = (grid > .4*grid.max()) & (grid < .8*grid.max())          # 흑연 stage 피크 탐색 구간
    p, pr = find_peaks(dv[w], prominence=0.001)
    xpk = grid[w][p[np.argmax(pr["prominences"])]] / grid.max() if len(p) else np.nan
    rows.append(dict(cell=cell, ah=ah, temp=TEMP[cell], soh=m.soh, xpk=xpk,
                     Rs50=g.Rs[50], Rint50=g.Rint[50], rlow=g.Rint[20]/g.Rint[50], rhigh=g.Rint[90]/g.Rint[50],
                     soc_Rmin=g.Rint.idxmin()))
F = pd.DataFrame(rows)
F = F[(F.groupby("cell").xpk.transform("first") - .595).abs() < .03]   # BOL 피크 오검출 셀 제외
for col in ["rlow", "rhigh", "xpk", "Rint50"]:
    F["d_"+col] = F[col] - F.groupby("cell")[col].transform("first")
F["d_soh"] = F.soh - 1
F.to_csv(OUT/"checkups.csv", index=False)
A = F[F.ah > 0]
def partial(x, y, z):   # z(SOH) 효과를 2차식으로 제거한 뒤 Spearman
    rx = x - np.polyval(np.polyfit(z, x, 2), z); ry = y - np.polyval(np.polyfit(z, y, 2), z); return spearmanr(rx, ry)
print(f"check-ups {len(A)}, cells {A.cell.nunique()}")
print("R_int 최솟점이 90%(측정 끝)인 비율:", round((F.soc_Rmin == 90).mean(), 2))
for a, b in [("d_rlow", "d_xpk"), ("d_rhigh", "d_xpk"), ("d_Rint50", "d_xpk")]:
    r = spearmanr(A[a], A[b]); p = partial(A[a].values, A[b].values, A.d_soh.values)
    print(f"{a} vs {b}: rho={r.statistic:+.2f} | SOH 보정 후 rho={p.statistic:+.2f} (p={p.pvalue:.1e})")
for t, g in A.groupby("temp"):
    p = partial(g.d_Rint50.values, g.d_xpk.values, g.d_soh.values)
    print(f"  온도 그룹 {t} (n={len(g)}): SOH 보정 후 rho={p.statistic:+.2f}")
