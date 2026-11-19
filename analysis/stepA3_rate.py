"""트랙 A 3단계: 셀마다 '늙는 속도' 재기 -> out/rates.csv, out/A3_soh_curves.png

모델 (셀마다 따로 맞춤):  SOH = 1 - k * EFC^z
  - EFC = 누적 충방전량(Ah) / (2 x 공칭 5 Ah)  = 완전 충방전 횟수로 환산 (DOD 50% 셀도 같은 잣대)
  - k  = 늙는 속도 (클수록 빨리 늙음)  -> 4단계에서 조건(온도, C-rate)으로 설명할 값
  - z  = 곡선 모양. 0.5면 SEI처럼 점점 느려짐(√), 1이면 일정한 속도
  - z는 두 가지로 맞춰 봄: 0.5 고정 / 1 고정 / 자유. 어느 쪽이 잘 맞는지는 직접 판단
입력: out/c20.csv (step1에서 생성: cell, ah, cap)
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from common import OUT, TEMP
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

Q_NOM = 5.0  # Ah, README: DOD 50% = 2.5 Ah
C = pd.read_csv(OUT / "c20.csv").sort_values(["cell", "ah"])
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
C["efc"] = C.ah / (2 * Q_NOM)

rows = []
for cell, g in C.groupby("cell"):
    x, y = g.efc.values, g.soh.values
    r = dict(cell=cell, temp=TEMP[cell], n=len(g), efc_end=x.max(), soh_end=y.min())
    if len(g) < 3:
        rows.append(r); continue
    for name, z in [("z05", 0.5), ("z1", 1.0)]:
        k, _ = curve_fit(lambda e, k: 1 - k * e**z, x, y, p0=[1e-3])
        r[f"k_{name}"] = k[0]
        r[f"rmse_{name}"] = np.sqrt(np.mean((1 - k[0] * x**z - y)**2))
    try:
        (k, z), _ = curve_fit(lambda e, k, z: 1 - k * e**z, x, y, p0=[1e-3, 0.7], bounds=([0, 0.2], [1, 2]))
        r.update(k_free=k, z_free=z, rmse_free=np.sqrt(np.mean((1 - k * x**z - y)**2)))
    except RuntimeError:
        pass
    rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(OUT / "rates.csv", index=False)
print(R.round(4).to_string(index=False))

# 그림: 7묶음별 SOH vs EFC (점 = 측정, 선 = z 자유 fit)
col = {"R": "#1d5fd1", "C": "#2a9d8f", "H": "#c0392b"}
cond = ["C/5, 0-100%", "1.5C, 0-100%", "2C, 0-100%", "C/5 충전/1.5C 방전, 0-100%",
        "C/5, 0-50%", "C/5 충전/1.5C 방전, 0-50%", "1.5C 충전/주행, 0-50%"]
fig, ax = plt.subplots(2, 4, figsize=(18, 8), sharey=True); ax = ax.ravel()
for gi in range(7):
    for c in range(3*gi + 1, 3*gi + 4):
        g = C[C.cell == c]
        if g.empty: continue
        ax[gi].plot(g.efc, g.soh, "o", ms=3, color=col[TEMP[c]], label=f"{c:02d} {TEMP[c]}")
        r = R[R.cell == c].iloc[0]
        if pd.notna(r.get("k_free")):
            e = np.linspace(0, g.efc.max(), 50)
            ax[gi].plot(e, 1 - r.k_free * e**r.z_free, "-", color=col[TEMP[c]], lw=1)
    ax[gi].set(title=f"cells {3*gi+1:02d}-{3*gi+3:02d}: {cond[gi]}", xlabel="EFC (완전 충방전 환산 횟수)", ylabel="SOH")
    ax[gi].legend(fontsize=7)
ax[7].axis("off"); plt.tight_layout()
plt.savefig(OUT / "A3_soh_curves.png", dpi=130); plt.show()
