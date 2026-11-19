"""5단계 (직접 해보기): 셀 팽창 데이터로 '저온 셀 = LLI형' 신호를 독립 검증
질문: 같은 SOH에서 저온(C) 셀이 비가역 팽창도 더 큰가? (plating 흔적 후보)
기각 조건(미리 정함): SOH를 맞췄을 때 저온 셀의 비가역 팽창이 다른 그룹과 차이가 없다.
입력: ../data/XX/OCV_wExpansion.csv (Expansion [mu m]), out/c20.csv (step1 결과)
출력: out/expansion.csv, out/expansion.png
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from common import *

C = pd.read_csv(OUT/"c20.csv")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
rows = []
for cell in range(1, 22):
    d = pd.read_csv(DATA / f"{cell:02d}" / "OCV_wExpansion.csv")
    d.columns = ["t", "I", "V", "x", "T", "Q", "cap", "cyc", "ah"]
    for ah, v in d.groupby("ah"):
        ch = v[(v.I > 100) & (v.I < 400)]
        if len(ch) < 500: continue
        # 완전 방전 상태(C/20 충전 시작 직전)의 두께 = 비가역 팽창의 기준점
        x0 = ch.x.iloc[:20].median()
        # 완전 충전 상태의 두께 - 방전 상태 두께 = 가역 팽창
        x1 = ch.x.iloc[-20:].median()
        rows.append(dict(cell=cell, ah=ah, temp=TEMP[cell], x_dis=x0, x_rev=x1 - x0))
X = pd.DataFrame(rows).merge(C[["cell", "ah", "soh"]], on=["cell", "ah"])
X["x_irr"] = X.x_dis - X.groupby("cell").x_dis.transform("first")   # BOL 대비 비가역 팽창 [um]
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
