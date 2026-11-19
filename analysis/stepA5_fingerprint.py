"""트랙 A 5단계: 초반 EIS 지문이 기준 셀(01)과 얼마나 다른가 -> out/A5_distance.csv, out/A5_fingerprint.png

지문 정의 (셀마다 자기 새것 상태와 비교 -> 셀 간 제조 편차와 측정 온도 차이를 줄이기 위함):
  f(SOC) = ln( X(SOC, 지금) / X(SOC, BOL) )      SOC = 10, 20, ..., 90 %  (9개 점)
  X = |Z| at 1 Hz  (기본, 회로 fitting 없이 측정값 그대로 -> fitting 불안정의 영향 없음)
      R_int = R1 + R2 (확인용, 회로 fitting 값)
  -> "늙으면서 SOC별 저항이 몇 % 변했나"의 곡선

비교 시점: SOH = 0.97, 0.95, 0.92, 0.90  (같은 SOH에서 비교해야 '늙은 정도' 차이가 섞이지 않음)
  - EIS를 잰 check-up들 사이를 SOH 기준으로 직선 보간
  - 그 SOH까지 EIS가 없으면 비움(NaN)

거리 (기준 셀 01과 비교, SOC 9개 점의 RMS):
  d_level : 지문 곡선 전체 차이 (저항이 얼마나 '많이' 변했나 + 모양)
  d_shape : 곡선 평균을 뺀 뒤의 차이 (저항이 'SOC별로 어떻게' 변했나 = 모양만, H-a의 바닥 이동·곡률 포함)
입력: out/eis_fits.csv (step1), out/c20.csv (step1)
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from common import OUT, TEMP

plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

REF = 1                               # 기준 셀 (상온, C/5, 0-100%)
LEVELS = [0.97, 0.95, 0.92, 0.90]     # 비교할 SOH
SOCS = [10, 20, 30, 40, 50, 60, 70, 80, 90]

E = pd.read_csv(OUT / "eis_fits.csv")
C = pd.read_csv(OUT / "c20.csv").sort_values(["cell", "ah"])
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])

def fingerprints(cell, col):
    """check-up별 지문 -> (soh 배열, 지문 행렬 [check-up x SOC])"""
    e = E[E.cell == cell].pivot_table(index="ah", columns="soc", values=col).reindex(columns=SOCS).sort_index()
    e = e.dropna()
    c = C[C.cell == cell]
    soh = np.interp(e.index.values, c.ah.values, c.soh.values)       # EIS check-up 시점의 SOH (C/20 사이 보간)
    f = np.log(e.values / e.values[0])                                # 첫 EIS check-up(BOL) 대비
    return soh, f, e.index.values[0]

def at_level(soh, f, s):
    """SOH = s 에서의 지문 (앞뒤 check-up 직선 보간). 도달 못 했으면 None"""
    if soh.min() > s: return None
    j = np.argmax(soh <= s)
    if j == 0: return f[0]
    w = (soh[j-1] - s) / (soh[j-1] - soh[j])
    return f[j-1] + w * (f[j] - f[j-1])

rows, FP = [], {}
for col, name in [("mag_1Hz", "Z1Hz"), ("Rint", "Rint")]:
    ref_soh, ref_f, _ = fingerprints(REF, col)
    for cell in sorted(E.cell.unique()):
        soh, f, ah0 = fingerprints(cell, col)
        for s in LEVELS:
            a, r = at_level(soh, f, s), at_level(ref_soh, ref_f, s)
            if name == "Z1Hz" and a is not None: FP[(cell, s)] = a
            if a is None or r is None:
                rows.append(dict(x=name, cell=cell, temp=TEMP[cell], soh=s, d_level=np.nan, d_shape=np.nan, n_eis=len(soh), bol_ah=ah0)); continue
            diff = a - r
            rows.append(dict(x=name, cell=cell, temp=TEMP[cell], soh=s, n_eis=len(soh), bol_ah=ah0,
                             d_level=np.sqrt(np.mean(diff**2)),
                             d_shape=np.sqrt(np.mean((diff - diff.mean())**2))))
D = pd.DataFrame(rows)
D.to_csv(OUT / "A5_distance.csv", index=False)

# 표: |Z|1Hz 기준, 모양 거리(d_shape)를 SOH별로
print("기준 셀:", REF, "| 지문 = ln(|Z|1Hz / BOL), 거리 = 기준 셀과의 RMS (SOC 9점)")
print("EIS check-up 수가 적은 셀(3개 이하):", D[(D.x == "Z1Hz") & (D.n_eis <= 3)].cell.unique().tolist())
for x in ["Z1Hz", "Rint"]:
    for k in ["d_shape", "d_level"]:
        t = D[D.x == x].pivot_table(index=["cell", "temp"], columns="soh", values=k).round(3)
        print(f"\n[{x}] {k}  (열 = SOH, 빈칸 = 그 SOH까지 EIS 없음)")
        print(t.to_string())

# 그림: SOH 0.95에서 묶음별 지문 곡선 (기준 셀은 검은 점선)
col = {"R": "#1d5fd1", "C": "#2a9d8f", "H": "#c0392b"}
fig, ax = plt.subplots(2, 4, figsize=(18, 8), sharey=True); ax = ax.ravel()
S = 0.95
for gi in range(7):
    if (REF, S) in FP: ax[gi].plot(SOCS, FP[(REF, S)], "k--", lw=1, label=f"기준 {REF:02d}")
    for c in range(3*gi + 1, 3*gi + 4):
        if (c, S) in FP:
            ax[gi].plot(SOCS, FP[(c, S)], "o-", ms=3, color=col[TEMP[c]], label=f"{c:02d} {TEMP[c]}")
    ax[gi].axhline(0, color=".7", lw=.8)
    ax[gi].set(title=f"cells {3*gi+1:02d}-{3*gi+3:02d} @ SOH {S}", xlabel="SOC [%]", ylabel="ln(|Z|1Hz / BOL)")
    ax[gi].legend(fontsize=7)
ax[7].axis("off"); plt.tight_layout()
plt.savefig(OUT / "A5_fingerprint.png", dpi=130); plt.show()
