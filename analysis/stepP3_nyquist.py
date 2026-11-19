"""P3 보조: check-up별 Nyquist 겹쳐 그리기 — 노화에 따라 '시작점(R0)'과 '반원'이 각각 어떻게 변하는지 눈으로 확인

[무엇을 그리나]
  - 셀마다 한 칸, 같은 SOC에서 측정한 check-up별 스펙트럼을 겹침 (색 = SOH, 노랑 1.0 -> 보라 SOH_MIN)
  - 세로 막대 = R0 (허수부가 0이 되는 지점),  점 = 반원 끝 (-Im의 마지막 국소 최소 = 확산 꼬리 직전)
  - R0와 반원 끝 정의는 stepP3_arc.py 와 동일
[설정] 아래 CELLS, SOC, SOH_MIN 만 바꾸면 다른 셀·SOC도 볼 수 있음
[가정과 한계]
  - 고온 셀은 0 Ah(EIS)와 첫 C/20(약 47 Ah) check-up이 모두 SOH 1.0으로 표시됨 (기준점 artifact)
  - 유도성 구간(-Im < 0)은 대부분 잘라서 그림
입력: ../data/XX/Resistance.csv (common.load_eis), out/c20.csv, out/checkup_time.csv
출력: out/P3_nyquist_SOC{SOC}.png
변경 이력: 2026-11-19 최초 작성 (Claude가 그린 P3_nyquist_overlay.png 를 재현 가능한 스크립트로 정리)
           2026-11-19 반원 끝 = '마지막 국소 최소'로 수정 (stepP3_arc.py 와 동일, 셀 09·15 오검출 확인 후 합의)
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from common import OUT, TEMP, load_eis

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False

# ---------- 설정 ----------
CELLS = [9, 15]
SOC = 10        # 한 번 실행 후 20, 30으로도
#CELLS = [2, 6, 7, 1]          # 저온 02, 고온 06, 상온 급속충전 07, 상온 저속 01
#SOC = 20                      # [%] 10~90 중 하나
SOH_MIN = 0.80                # 이보다 낮은 check-up은 그리지 않음
NAME = {"H": "고온 45 °C", "R": "상온 25 °C", "C": "저온 -5 °C"}

# ---------- check-up 표 (stepP3_arc.py 와 동일) ----------
C = pd.read_csv(OUT / "c20.csv").merge(pd.read_csv(OUT / "checkup_time.csv")[["cell", "ah", "hours"]],
                                         on=["cell", "ah"], how="left")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])

def soh_of(cell, ah):
    if ah < 1: return 1.0
    k = C[(C.cell == cell) & ((C.ah - ah).abs() < 1)]
    return k.soh.iloc[0] if len(k) else np.nan

def r0_and_arc_end(f, re, mim):
    """R0 = 허수부 0 지점의 실수부, 반원 끝 = R0 이후 -Im 국소 최소 중 마지막의 인덱스 (stepP3_arc.py 와 동일)"""
    z = np.where((mim[:-1] < 0) & (mim[1:] >= 0))[0]
    if len(z) == 0: return np.nan, None
    i = z[0]; r0 = np.interp(0, [mim[i], mim[i + 1]], [re[i], re[i + 1]])
    mins = [j for j in range(i + 2, len(mim) - 1)
            if mim[j] <= mim[j - 1] and mim[j] < mim[j + 1] and mim[i + 1:j].max() > mim[j]]
    return r0, (mins[-1] if mins else None)

norm = colors.Normalize(SOH_MIN, 1.0); cmap = cm.viridis
n = len(CELLS); ncol = 2; nrow = int(np.ceil(n / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 5 * nrow), squeeze=False)
for ax, cell in zip(axs.flat, CELLS):
    spec = load_eis(cell); xmax = ymax = 0
    for (ah, soc) in sorted(k for k in spec if k[1] == SOC):
        soh = soh_of(cell, ah)
        if not np.isfinite(soh) or soh < SOH_MIN: continue
        a = spec[(ah, soc)]; f, re, mim = a[:, 0], a[:, 3], a[:, 4]
        col = cmap(norm(soh)); keep = mim >= -0.0002
        ax.plot(re[keep] * 1e3, mim[keep] * 1e3, "-", color=col, lw=1.3)
        r0, j = r0_and_arc_end(f, re, mim)
        if np.isfinite(r0): ax.plot(r0 * 1e3, 0, "|", color=col, ms=14, mew=2)
        if j is not None:
            ax.plot(re[j] * 1e3, mim[j] * 1e3, "o", color=col, ms=5)
            xmax = max(xmax, re[j] * 1e3); ymax = max(ymax, mim[:j + 1].max() * 1e3)
    if xmax > 0: ax.set_xlim(None, xmax * 1.25); ax.set_ylim(-0.2, ymax * 2.2)
    ax.set(title=f"셀 {cell:02d} · {NAME[TEMP[cell]]}", xlabel="Re Z [mΩ]", ylabel="-Im Z [mΩ]"); ax.grid(alpha=.3)
for ax in axs.flat[n:]: ax.axis("off")
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, label=f"SOH (노랑 1.0 -> 보라 {SOH_MIN})", shrink=0.6)
fig.suptitle(f"Nyquist 겹쳐 그리기 (SOC {SOC}%, 25 °C 측정, SOH >= {SOH_MIN})  |  막대 = R0   점 = 반원 끝", fontsize=13)
out = OUT / f"P3_nyquist_SOC{SOC}.png"
fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
print("저장:", out)
