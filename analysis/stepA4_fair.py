"""트랙 A 4단계 (공정 비교): 시간 기준 수명식이 '시간 환산 산수' 이상을 설명하는지 확인

문제: 시간 기준 수명은 느린 셀일수록 1회당 시간이 길어서 저절로 길게 나온다.
      평균 기준선은 이걸 모르니 크게 틀리고, 선형식은 이 산수만 배워도 좋아 보인다.
해결: 기준선에 '시간 환산'을 공짜로 준다.
  - h_per_efc = 그 셀의 '첫 번째 check-up까지의 시간 / EFC'  (시험 프로토콜로 정해지는 값, 수명 결과는 안 씀)
  - 공정 기준선 = (나머지 셀들의 횟수 기준 평균 수명) x (그 셀의 h_per_efc)
비교 (모두 leave-one-cell-out, ln 기준 RMSE):
  ① 평균 기준선        : 나머지 셀 평균 수명 [h]
  ② 공정 기준선        : 위 설명 (시간 환산만 반영)
  ③ 선형식 (시간 기준)  : stepA4_time.py와 같은 Ridge 모델
  -> ③이 ②보다 작아야 '조건이 시간 환산 이상으로 수명을 설명한다'고 말할 수 있음
입력: out/c20.csv, out/checkup_time.csv (stepA2_time.py)
"""
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from common import OUT, TEMP

TH = 0.85
Q_NOM = 5.0
GROUP = {0: (0.2, 0.2, 100), 1: (1.5, 1.5, 100), 2: (2.0, 2.0, 100), 3: (0.2, 1.5, 100),
         4: (0.2, 0.2, 50), 5: (0.2, 1.5, 50), 6: (1.5, np.nan, 50)}

C = pd.read_csv(OUT / "c20.csv").sort_values(["cell", "ah"])
Tm = pd.read_csv(OUT / "checkup_time.csv")
C = C.merge(Tm[["cell", "ah", "hours"]], on=["cell", "ah"], how="left")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
C["efc"] = C.ah / (2 * Q_NOM)

def cross(x, s):   # SOH가 처음으로 TH 아래로 내려가는 지점을 직선 보간
    i = np.argmax(s < TH) if (s < TH).any() else 0
    return np.interp(TH, [s[i], s[i-1]], [x[i], x[i-1]]) if i else np.nan

rows = []
for cell, g in C.groupby("cell"):
    g = g.dropna(subset=["hours"])
    first = g[g.efc > 0].iloc[0]                         # 첫 번째 check-up (BOL 다음)
    rows.append(dict(cell=cell, temp=TEMP[cell],
                     efc_life=cross(g.efc.values, g.soh.values),
                     h_life=cross(g.hours.values, g.soh.values),
                     h_per_efc=first.hours / first.efc))
R = pd.DataFrame(rows)
R[["c_ch", "c_dis", "dod"]] = [GROUP[(c - 1) // 3] for c in R.cell]
R = R.dropna(subset=["efc_life", "h_life", "c_dis"]).copy()   # 주행 사이클 셀 제외
R["cold"] = (R.temp == "C").astype(float); R["hot"] = (R.temp == "H").astype(float)
R["ln_ch"] = np.log(R.c_ch); R["ln_dis"] = np.log(R.c_dis)
R["dod50"] = (R.dod == 50).astype(float); R["cold_x_ln_ch"] = R.cold * R.ln_ch
R["y"] = np.log(R.h_life); R["ln_efc"] = np.log(R.efc_life)
FEATS = ["cold", "hot", "ln_ch", "ln_dis", "dod50", "cold_x_ln_ch"]

p1, p2, p3 = [], [], []
for c in R.cell:
    tr, te = R[R.cell != c], R[R.cell == c]
    p1.append(tr.y.mean())                                              # ① 평균
    p2.append(tr.ln_efc.mean() + np.log(te.h_per_efc.iloc[0]))          # ② 공정 기준선
    p3.append(RidgeCV(alphas=np.logspace(-3, 2, 11)).fit(tr[FEATS], tr.y).predict(te[FEATS])[0])  # ③ 선형식
rmse = lambda p: np.sqrt(np.mean((np.array(p) - R.y.values)**2))
print(f"수명 기준 SOH = {TH} | 셀 {len(R)}개 | LOCO 오차 (ln 기준 RMSE, 0.1 ~ 약 10%)")
print(f"  ① 평균 기준선       {rmse(p1):.3f}")
print(f"  ② 공정 기준선       {rmse(p2):.3f}   <- 시간 환산만 반영")
print(f"  ③ 선형식(시간 기준) {rmse(p3):.3f}")
R["pred_fair_h"], R["pred_model_h"] = np.exp(p2), np.exp(p3)
out = R[["cell", "temp", "c_ch", "c_dis", "dod", "h_per_efc", "efc_life", "h_life", "pred_fair_h", "pred_model_h"]].round(1)
print(out.to_string(index=False))
out.to_csv(OUT / "A4_fair.csv", index=False)
