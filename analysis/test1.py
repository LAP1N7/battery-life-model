plt.rcParams["font.family"] = "Malgun Gothic"     # 윈도우 기본 한글 글꼴(맑은 고딕)
plt.rcParams["axes.unicode_minus"] = False        # 마이너스 기호 깨짐 방지

# 7개 묶음(온도 R/C/H 3셀씩)의 R_int(SOC): BOL(점선) vs end(실선)
import pandas as pd, matplotlib.pyplot as plt
from common import OUT, TEMP

E = pd.read_csv(OUT / "eis_fits.csv")
col = {"R": "#1d5fd1", "C": "#2a9d8f", "H": "#c0392b"}
cond = ["C/5, 0-100%", "1.5C, 0-100%", "2C, 0-100%", "C/5 충/1.5C 방, 0-100%",
        "C/5, 0-50%", "C/5 충/1.5C 방, 0-50%", "1.5C 충/주행, 0-50%"]
fig, ax = plt.subplots(2, 4, figsize=(18, 8)); ax = ax.ravel()
for g in range(7):
    for c in range(3*g + 1, 3*g + 4):
        x = E[E.cell == c]
        if x.empty: continue
        ahs = sorted(x.ah.unique())
        for a, ls in [(ahs[0], "--"), (ahs[-1], "-")]:
            y = x[x.ah == a].sort_values("soc")
            ax[g].plot(y.soc, y.Rint*1e3, ls, marker="o", ms=3, color=col[TEMP[c]],
                       label=f"{c:02d} {TEMP[c]}" if ls == "-" else None)
    ax[g].set(title=f"cells {3*g+1:02d}-{3*g+3:02d}: {cond[g]}", xlabel="SOC [%]", ylabel="R_int [mOhm]")
    ax[g].legend(fontsize=7)
ax[7].axis("off"); plt.tight_layout(); plt.show()