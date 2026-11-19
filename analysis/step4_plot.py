"""4단계: 결과 그림 -> out/results.png"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from common import *
E = pd.read_csv(OUT/"eis_fits.csv"); F = pd.read_csv(OUT/"checkups.csv"); R = pd.read_csv(OUT/"ml1.csv", index_col=0)
A = F[F.ah > 0].copy()
for col in ["d_Rint50", "d_xpk"]:
    A["r_"+col] = A[col] - np.polyval(np.polyfit(A.d_soh, A[col], 2), A.d_soh)
col = {"R": "#1d5fd1", "C": "#2a9d8f", "H": "#c0392b"}; nm = {"R": "room", "C": "cold", "H": "hot"}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
for c in [1, 2, 3]:
    ahs = sorted(E[E.cell == c].ah.unique())
    for a, ls, lab in [(ahs[0], "--o", "BOL"), (ahs[-1], "-o", "end")]:
        x = E[(E.cell == c) & (E.ah == a)].sort_values("soc")
        ax[0].plot(x.soc, x.Rint*1e3, ls, color=col[TEMP[c]], ms=3, label=f"cell {c:02d} ({nm[TEMP[c]]}) {lab}")
ax[0].set(xlabel="SOC [%]", ylabel="R_int = R1+R2 [mOhm]", title="R_int(SOC): BOL vs end"); ax[0].legend(fontsize=7)
for t, g in A.groupby("temp"):
    ax[1].scatter(g.r_d_xpk, g.r_d_Rint50*1e3, s=14, color=col[t], label=f"{nm[t]} (n={len(g)})")
ax[1].axhline(0, color=".7", lw=.8); ax[1].axvline(0, color=".7", lw=.8)
ax[1].set(xlabel="dV/dQ peak shift, SOH-adjusted [Q/Qtot]", ylabel="R_int@50% growth, SOH-adjusted [mOhm]", title="Same SOH, different path"); ax[1].legend(fontsize=8)
y = R["GPR"]; ax[2].barh(range(len(y)), y.values, color="0.6"); ax[2].set_yticks(range(len(y))); ax[2].set_yticklabels(y.index)
for i, v in enumerate(y.values): ax[2].text(v+.1, i, f"{v:.2f}", va="center")
ax[2].set(xlabel="SOH RMSE [%p], leave-one-cell-out (GPR)", title="ML-1: which features add information?"); ax[2].invert_yaxis()
plt.tight_layout(); plt.savefig(OUT/"results.png", dpi=140); print("saved", OUT/"results.png")
