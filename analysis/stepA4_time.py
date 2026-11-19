"""트랙 A 4단계 (기준선 모델, 횟수/시간 두 기준): 조건(온도, C-rate, DOD)으로 '수명'을 설명하는 선형식
LIFE = "efc"   -> 수명 = SOH가 TH까지 떨어질 때의 충방전 환산 횟수
LIFE = "hours" -> 수명 = SOH가 TH까지 떨어질 때의 누적 사이클링 시간 [h] (stepA2_time.py 먼저 실행)

수명 지표: EFC_th = SOH가 TH(기본 0.85)까지 떨어질 때의 완전 충방전 환산 횟수
  - 측정점 사이를 직선 보간해서 구함 (곡선 fit이나 외삽 없음)
  - TH에 도달하지 못한 셀은 자동 제외 (0.85면 21셀 모두 사용, 0.80이면 5셀 제외)
  - 본 분석: TH = 0.85 / 확인용: TH = 0.80 으로 다시 돌려 결론이 같은지 비교

선형식 (ln을 씌우면 곱셈 효과가 덧셈이 됨):
  ln(EFC_th) = b0 + b1*[저온] + b2*[고온] + b3*ln(충전 C-rate) + b4*ln(방전 C-rate) + b5*[DOD 50%] + b6*[저온 x ln(충전 C-rate)]
  - 계수가 음수(-)면 그 조건이 수명을 줄인다는 뜻
  - b6 (저온 x 빠른 충전): H-T1이 맞으면 음수(-)로 커야 함 (저온+빠른 충전이 서로 수명을 깎아 먹음)
  - 주행 사이클 셀(19-21)은 방전 C-rate를 정할 수 없어 기본으로 제외 (USE_DRIVE = False)
검증: leave-one-cell-out (셀 하나씩 빼고 나머지로 식을 만든 뒤 뺀 셀의 수명을 예측)
기준선: 조건을 전혀 안 쓰고 '나머지 셀 평균 수명'으로 예측
입력: out/c20.csv (step1_extract.py에서 생성)
출력: out/A4_life_TH85.csv (또는 TH80), out/A4_pred_TH85.png
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from sklearn.linear_model import RidgeCV
from common import OUT, TEMP

plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

TH = 0.85          # 수명 기준 SOH. 확인용으로 0.80도 돌려 볼 것
LIFE = "hours"     # "efc"(횟수) 또는 "hours"(시간)
USE_DRIVE = False  # 주행 사이클 셀 포함 여부
Q_NOM = 5.0        # Ah, README: DOD 50% = 2.5 Ah

# README Table 1 조건 (묶음 번호 -> 충전 C-rate, 방전 C-rate, DOD)
GROUP = {0: (0.2, 0.2, 100), 1: (1.5, 1.5, 100), 2: (2.0, 2.0, 100), 3: (0.2, 1.5, 100),
         4: (0.2, 0.2, 50), 5: (0.2, 1.5, 50), 6: (1.5, np.nan, 50)}

# 1) 셀마다 SOH = TH가 되는 EFC (처음으로 TH 아래로 내려가는 구간에서 직선 보간)
C = pd.read_csv(OUT / "c20.csv").sort_values(["cell", "ah"])
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
C["efc"] = C.ah / (2 * Q_NOM)
if LIFE == "hours":
    Tm = pd.read_csv(OUT / "checkup_time.csv")
    C = C.merge(Tm[["cell", "ah", "hours"]], on=["cell", "ah"], how="left")
    C["efc"] = C["hours"]          # 아래 계산은 같은 열 이름을 쓰므로 시간으로 바꿔 넣음
    C = C.dropna(subset=["efc"])
rows = []
for cell, g in C.groupby("cell"):
    s, e = g.soh.values, g.efc.values
    i = np.argmax(s < TH) if (s < TH).any() else None
    efc_th = np.interp(TH, [s[i], s[i-1]], [e[i], e[i-1]]) if i else np.nan
    rows.append(dict(cell=cell, temp=TEMP[cell], efc_th=efc_th))
R = pd.DataFrame(rows)
skipped = R[R.efc_th.isna()].cell.tolist()
R = R.dropna(subset=["efc_th"]).copy()
R[["c_ch", "c_dis", "dod"]] = [GROUP[(c - 1) // 3] for c in R.cell]
if not USE_DRIVE:
    R = R[R.c_dis.notna()]
R["cold"] = (R.temp == "C").astype(float)
R["hot"] = (R.temp == "H").astype(float)
R["ln_ch"] = np.log(R.c_ch)
R["ln_dis"] = np.log(R.c_dis)
R["dod50"] = (R.dod == 50).astype(float)
R["cold_x_ln_ch"] = R.cold * R.ln_ch
R["y"] = np.log(R.efc_th)
FEATS = ["cold", "hot", "ln_ch", "ln_dis", "dod50", "cold_x_ln_ch"]

def model():
    return RidgeCV(alphas=np.logspace(-3, 2, 11))   # 셀 수가 적어 계수가 튀지 않게 약하게 묶음

tag = f"TH{int(round(TH*100))}_{LIFE}"
print(f"수명 기준 SOH = {TH} | 사용 셀 {len(R)}개 | {TH} 미도달로 제외: {skipped if skipped else '없음'}")

# 2) 전체 셀로 만든 식 (해석용)
m = model().fit(R[FEATS], R.y)
print(f"Ridge alpha = {m.alpha_:.3g}")
print(f"ln(수명[{LIFE}]) = {m.intercept_:+.3f}" + "".join(f"  {b:+.3f}*{f}" for f, b in zip(FEATS, m.coef_)))

# 3) leave-one-cell-out
pm, pb = [], []
for c in R.cell:
    tr, te = R[R.cell != c], R[R.cell == c]
    pm.append(model().fit(tr[FEATS], tr.y).predict(te[FEATS])[0]); pb.append(tr.y.mean())
R["pred_model"], R["pred_base"] = np.exp(pm), np.exp(pb)
rm = np.sqrt(np.mean((np.array(pm) - R.y)**2)); rb = np.sqrt(np.mean((np.array(pb) - R.y)**2))
print(f"LOCO 오차 (ln 기준 RMSE): 선형식 {rm:.3f} | 기준선(평균) {rb:.3f}   * 0.1 ~ 수명 약 10% 차이")
out = R[["cell", "temp", "c_ch", "c_dis", "dod", "efc_th", "pred_model", "pred_base"]].round(1)
print(out.to_string(index=False))
out.to_csv(OUT / f"A4_life_{tag}.csv", index=False)

# 4) 그림: 실제 수명 vs 예측 수명 (대각선에 가까울수록 잘 맞음)
col = {"R": "#1d5fd1", "C": "#2a9d8f", "H": "#c0392b"}
fig, ax = plt.subplots(figsize=(6, 5.5))
for _, r in R.iterrows():
    ax.plot(r.efc_th, r.pred_model, "o", color=col[r.temp])
    ax.annotate(f"{int(r.cell):02d}", (r.efc_th, r.pred_model), fontsize=8, xytext=(3, 3), textcoords="offset points")
lim = [0, max(R.efc_th.max(), R.pred_model.max()) * 1.1]
ax.plot(lim, lim, "k--", lw=1); ax.set(xlim=lim, ylim=lim)
for t, n in [("R", "room"), ("C", "cold"), ("H", "hot")]:
    ax.plot([], [], "o", color=col[t], label=n)
ax.set(xlabel=f"실제 수명: SOH {TH}까지 {'시간 [h]' if LIFE == 'hours' else 'EFC'}", ylabel="예측 수명 (그 셀을 빼고 만든 식)", title=f"A4 기준선 모델 ({tag})")
ax.legend(); plt.tight_layout()
plt.savefig(OUT / f"A4_pred_{tag}.png", dpi=130); plt.show()
