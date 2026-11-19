"""5단계 (직접 해보기): 셀 팽창 데이터로 '저온 셀 = LLI형' 신호를 독립 검증
질문: 같은 SOH에서 저온(C) 셀이 비가역 팽창도 더 큰가? (plating 흔적 후보)
기각 조건(미리 정함): SOH를 맞췄을 때 저온 셀의 비가역 팽창이 다른 그룹과 차이가 없다.
입력: ../data/XX/OCV_wExpansion.csv (Expansion [mu m]), out/c20.csv (step1 결과)
출력: out/expansion.csv, out/expansion.png

변경 이력: 2026-11-19 (C) 45 °C 셀 x_dis 버그 수정. 100–400 mA 샘플 전체 -> 가장 긴 연속 구간(C/20 충전)만 사용.
           1C 충전 CV 꼬리가 bottom 자리에 top 두께로 들어가던 문제. x_irr 기준점을 첫 행 -> BOL(ah = 0)로 명시.
           BOL은 c20.csv와 병합하기 전에 잡는다(45 °C 셀은 c20.csv에 ah = 0 행이 없어 병합 후에는 BOL이 사라짐).
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from common import *


def longest_run(mask):
    """True가 연속된 가장 긴 구간의 위치 인덱스 (CV 꼬리 같은 짧은 구간 제외)"""
    idx = np.flatnonzero(mask)
    if not len(idx): return idx
    runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    return max(runs, key=len)


C = pd.read_csv(OUT/"c20.csv")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
rows = []
for cell in range(1, 22):
    d = pd.read_csv(DATA / f"{cell:02d}" / "OCV_wExpansion.csv")
    d.columns = ["t", "I", "V", "x", "T", "Q", "cap", "cyc", "ah"]
    for ah, v in d.groupby("ah"):
        v = v.sort_values("t")
        run = longest_run(((v.I > 100) & (v.I < 400)).to_numpy())   # 가장 긴 연속 구간 = C/20 충전
        if len(run) < 500: continue
        ch = v.iloc[run]
        x0 = ch.x.iloc[:20].median()    # bottom: C/20 충전 시작 (완전 방전 상태)
        x1 = ch.x.iloc[-20:].median()   # top: C/20 충전 끝 (완전 충전 상태)
        rows.append(dict(cell=cell, ah=ah, temp=TEMP[cell], x_dis=x0, x_rev=x1 - x0))
R = pd.DataFrame(rows)
bol = R[R.ah == 0].set_index("cell").x_dis                     # 기준점을 BOL(ah = 0)로 명시, 병합 전에 잡음
X = R.merge(C[["cell", "ah", "soh"]], on=["cell", "ah"]).sort_values(["cell", "ah"])
X["x_irr"] = X.x_dis - X.cell.map(bol)                         # BOL 대비 비가역 팽창 [um], BOL이 없으면 NaN
if X.x_irr.isna().any(): print("[경고] BOL check-up 없는 셀:", sorted(X[X.x_irr.isna()].cell.unique()))
X.to_csv(OUT/"expansion.csv", index=False)

# ---- TODO 1: 그림 확인 ----------------------------------------------------
# SOH(가로축) 대비 비가역 팽창(세로축)을 온도 그룹별 색으로 그립니다.
col = {"R": "#1d5fd1", "C": "#2a9d8f", "H": "#c0392b"}
for t, g in X.groupby("temp"):
    plt.scatter(g.soh*100, g.x_irr, s=12, color=col[t], label=t)
plt.gca().invert_xaxis(); plt.xlabel("SOH [%]"); plt.ylabel("irreversible expansion vs BOL [um]"); plt.legend()
plt.savefig(OUT/"expansion.png", dpi=130); print("saved out/expansion.png")

# ---- TODO 2: SOH를 맞춘 비교 ----------------------------------------------
# 힌트: SOH 0.80~0.90 구간만 골라 그룹별 x_irr 평균을 비교해 보세요.
#   sub = X[(X.soh > 0.80) & (X.soh < 0.90)]
#   print(sub.groupby("temp").x_irr.describe())

# ---- TODO 3: 판정 ---------------------------------------------------------
# 저온(C) 셀의 x_irr이 다른 그룹보다 뚜렷하게 크면 -> 'LLI형(plating 후보)' 해석과 일치
# 차이가 없으면 -> 기각 조건 충족. 결과를 README의 '1차 결과'에 한 줄로 적어 두세요.
