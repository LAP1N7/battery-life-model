import sys, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np, pandas as pd
from scipy.optimize import least_squares
import stepD5_resistance_model as m5
K = m5.kit_table(os.path.join(_HERE, "..", "..", "prototype", "results", "kit_eis", "spectrum_features_enriched.csv"))
M = m5.mohtat_table()
RG = 8.314e-3; TREF = 298.15
tr = M[M["T"].isin(["H", "R"])]; te = M[M["T"] == "C"]
def rm(y, p): return np.sqrt(np.mean((y - p) ** 2))
for comp in ("dR0", "dRarc"):
    # (a) classic accelerated-test model: one process, A exp(-Ea x) t^z, fitted on H+R, extrapolated to C
    def cl(p, d):
        x = (1 / d.TK.values - 1 / TREF) / RG; return np.exp(p[0] - p[1] * x) * d.days.values ** p[2]
    w = 1 / np.sqrt(tr.groupby("cell").cell.transform("size").values)
    r = least_squares(lambda p: w * (cl(p, tr) - tr[comp].values), [-3, 40, 0.7], bounds=([-15, -150, 0.1], [5, 150, 2]), loss="soft_l1", f_scale=0.3)
    pa = r.x
    # (b) mechanism model: Ea_a, Ea_b, m from KIT (incl. 0 C and calendar cells), amplitudes from Mohtat H+R
    pk, _, _ = m5.fit(K, K[comp].values)
    fx = {n: pk[m5.NAMES.index(n)] for n in ("Ea_a", "Ea_b", "m")}; fx["s"] = 0.0
    pb, _, _ = m5.fit(tr, tr[comp].values, fixed=fx)
    ta, hb, cb = m5.model(pb, te); predb = ta + hb + cb
    # (c) mechanism model with all constants fitted on Mohtat H+R only
    pc, _, _ = m5.fit(tr, tr[comp].values, fixed=dict(s=0.0))
    ta, hb, cb = m5.model(pc, te); predc = ta + hb + cb
    predb = pd.Series(np.asarray(predb), index=te.index); y = te[comp].values
    print(f"{comp}: C셀 예측 RMSE [mOhm]  기존(단일 아레니우스 t^z, H+R 학습)={rm(y, cl(pa, te)):.3f} (Ea={pa[1]:.0f}, z={pa[2]:.2f})"
          f" | 메커니즘(KIT 상수)={rm(y, predb):.3f} | 메커니즘(H+R만)={rm(y, predc):.3f} | 측정 평균 크기={np.mean(np.abs(y)):.3f}")
    # end-of-test bias per C cell
    last = te.groupby("cell").tail(1).index
    print("   C셀 마지막 check-up 측정/예측(기존)/예측(KIT 상수):",
          [(int(te.cell[i]), round(te[comp][i], 2), round(float(cl(pa, te.loc[[i]])[0]), 2), round(float(predb.loc[i]), 2)) for i in last])
