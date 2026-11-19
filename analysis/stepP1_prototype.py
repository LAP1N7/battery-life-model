"""프로토타입 P1: 메커니즘 2항 수식 + 초반 EIS로 수명(SOH 0.85 도달 시간) 예측

[수식]  용량 손실 Q(t) = A*sqrt(t) + B*t           (t: 노화 사이클링 누적 시간 [h])
        - A*sqrt(t): 시간형 열화 (SEI 성장처럼 시간이 지날수록 느려짐)  -> 가정
        - B*t      : 횟수형/선형 열화 (사이클마다 비슷하게 쌓임)       -> 가정
        수명 = Q(t) = 0.15 가 되는 t  (2차식의 해)

[예측 방법]  셀마다 '초반 시점'(SOH가 처음 0.975 이하가 된 check-up) 1개만 사용
        - 그 시점의 손실 Qe, 시간 te 를 안다.  하지만 점 1개로는 A와 B를 나눌 수 없다 (식별 불가)
        - 그래서 '초반 손실 중 시간형 비율 f'를 정해야 한다:  A = f*Qe/sqrt(te),  B = (1-f)*Qe/te
        - f를 무엇으로 정하느냐가 모델 차이:
            ① 용량만     : f = 모든 셀 공통 상수 (학습 셀에서 최적값)            -> 기준선
            ② M2 (사전)  : f = w * x / (Qe/sqrt(te))   "저항이 빨리 늘면 시간형(SEI)이 크다" 가정
            ③ M2b (사후) : f = clip(a + b*x)              방향을 데이터가 정함 (가중치 2개)
            ④ 대조       : f = clip(a + b*hot)            EIS 대신 '고온 여부' 조건만
        - x = 초반 저항 증가율 = (R_e/R_BOL - 1)/sqrt(te),   R = Rs+R1+R2 (SOC 평균, MATLAB ECM)
          (R1, R2를 따로 쓰지 않고 합만 쓰는 이유: R1이 Rs를 흡수하는 문제에 영향받지 않음)

[평가]  LOCO (셀 하나 빼고 학습 -> 뺀 셀 예측), ln(수명) RMSE.  0.1 ~ 약 10% 오차
        조건 선형식(③ stepA4_fair 결과, A4_fair.csv)과 같은 셀끼리 비교
        섞은 값 시험: x를 셀끼리 무작위로 섞어도 비슷하게 맞으면 우연

입력: out/c20.csv, out/checkup_time.csv, out/A4_fair.csv, ../matlab/out/ecm_matlab.csv
출력: out/P1_pred.csv, out/P1_pred.png
소요: 1~3분 (섞은 값 시험 포함)
변경 이력: 2026-11-18 최초 작성
           2026-11-19 P1_pred.csv 저장을 유효숫자 6자리로 변경 (round(4)가 A·B 값을 뭉개 P3 ② 판정을 훼손함)
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from common import OUT, TEMP

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
TH = 0.85; QL = 1 - TH
ECM = OUT.parent.parent / "matlab" / "out" / "ecm_matlab.csv"
N_PERM = 20

# ---------- 1. 데이터 ----------
C = pd.read_csv(OUT / "c20.csv").merge(pd.read_csv(OUT / "checkup_time.csv")[["cell", "ah", "hours"]],
                                         on=["cell", "ah"], how="left")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
E = pd.read_csv(ECM); E["R"] = E.Rs + E.R1 + E.R2
Rm = E.groupby(["cell", "ah"]).R.mean().reset_index()          # SOC 10~90 평균

def cross(x, s):   # SOH가 처음 TH 아래로 내려가는 지점 직선 보간
    i = np.argmax(s < TH) if (s < TH).any() else 0
    return np.interp(TH, [s[i], s[i - 1]], [x[i], x[i - 1]]) if i else np.nan

rows, skip = [], []
for c, g in C.groupby("cell"):
    g = g.dropna(subset=["hours"]).reset_index(drop=True)
    life = cross(g.hours.values, g.soh.values)
    r = Rm[Rm.cell == c]; rb = r[r.ah < 1]
    k = g[g.soh <= 0.975].head(1)
    if np.isnan(life) or rb.empty or k.empty: skip.append(c); continue
    k = k.iloc[0]; m = r[(r.ah - k.ah).abs() < 1]
    if m.empty: skip.append(c); continue                          # 초반 시점에 EIS 없으면 제외 (미래 값 안 씀)
    gg = g[g.soh >= 0.80]; X = np.c_[np.sqrt(gg.hours), gg.hours]
    ab, *_ = np.linalg.lstsq(X, 1 - gg.soh, rcond=None)           # 전체 궤적 피팅 (형태 점검용, 예측에는 안 씀)
    yh = X @ ab; y_ = 1 - gg.soh
    rows.append(dict(cell=c, temp=TEMP[c], life=life, te=k.hours, Qe=1 - k.soh, soh_e=k.soh,
                     x=(m.R.iloc[0] / rb.R.iloc[0] - 1) / np.sqrt(k.hours),
                     A_full=ab[0], B_full=ab[1], r2=1 - ((y_ - yh) ** 2).sum() / ((y_ - y_.mean()) ** 2).sum()))
D = pd.DataFrame(rows).reset_index(drop=True)
D["hot"] = (D.temp == "H").astype(float)
D["frac"] = D.te / D.life
print("제외된 셀 (수명 미도달 / 초반 EIS 없음):", skip)

# ---------- 2. 수식 ----------
def t85(f, te, Qe):           # 벡터 계산: A sqrt(t) + B t = QL 의 해
    A = f * Qe / np.sqrt(te); B = (1 - f) * Qe / te
    Bs = np.maximum(B, 1e-15)
    s = np.where(B > 1e-12, (-A + np.sqrt(A * A + 4 * Bs * QL)) / (2 * Bs), QL / np.maximum(A, 1e-15))
    return s * s

def lerr(p, life): return np.sqrt(np.mean((np.log(p) - np.log(life)) ** 2))

def loco(fmaker, grid, Dd):
    """fmaker(params, df) -> f 벡터.  grid = 파라미터 후보 목록.  LOCO 예측 반환"""
    P = []
    for c in Dd.cell:
        tr, te = Dd[Dd.cell != c], Dd[Dd.cell == c]
        best = min(grid, key=lambda p: lerr(t85(fmaker(p, tr), tr.te.values, tr.Qe.values), tr.life.values))
        P.append(t85(fmaker(best, te), te.te.values, te.Qe.values)[0])
    return np.array(P)

g1 = [(f,) for f in np.linspace(0, 1, 101)]
g2 = [(w,) for w in np.logspace(-2, 2, 100)]
def g_ab(col, Dd):
    s = 1 / max(Dd[col].abs().max(), 1e-12)
    return [(a, b * s) for a in np.linspace(0, 1, 41) for b in np.linspace(-1, 1, 61)]

f_cap = lambda p, d: np.full(len(d), p[0])
f_M2 = lambda p, d: np.clip(p[0] * d.x.values / (d.Qe.values / np.sqrt(d.te.values)), 0, 1)
f_ab = lambda col: (lambda p, d: np.clip(p[0] + p[1] * d[col].values, 0, 1))

# ---------- 3. 비교 ----------
res = {"① 용량만 (f 공통)": loco(f_cap, g1, D),
       "② M2 사전 가정 (저항↑ → 시간형↑)": loco(f_M2, g2, D),
       "③ M2b 사후 (a + b·x, EIS)": loco(f_ab("x"), g_ab("x", D), D),
       "④ 대조 (a + b·hot, 조건)": loco(f_ab("hot"), g_ab("hot", D), D)}
F = pd.read_csv(OUT / "A4_fair.csv").set_index("cell")
com = D.cell.isin(F.index).values
print(f"\n셀 {len(D)}개 | LOCO 오차 (ln RMSE) | 괄호: 조건 선형식과 같은 {com.sum()}개 셀")
for k, p in res.items():
    print(f"  {k:30s} {lerr(p, D.life.values):.3f}  ({lerr(p[com], D.life.values[com]):.3f})")
print(f"  {'M0 조건 선형식 (A4_fair)':30s}   -    ({lerr(F.loc[D.cell[com]].pred_model_h.values, D.life.values[com]):.3f})")

# ---------- 4. 점검 ----------
print("\n[형태 점검] 전체 궤적을 A√t + Bt 로 피팅한 R²: 중앙값 %.3f, 최소 %.3f" % (D.r2.median(), D.r2.min()))
print("[가정 점검] 초반 저항증가율 x 와 전체 궤적 A(시간형 크기)의 순위상관: %.2f  (사전 가정이 맞으면 +)"
      % D[["x", "A_full"]].corr("spearman").iloc[0, 1])
print("[조기성] 측정 시점 / 수명: 중앙값 %.2f, 범위 %.2f ~ %.2f" % (D.frac.median(), D.frac.min(), D.frac.max()))
nh = D[D.hot == 0].reset_index(drop=True)
print("[고온 제외 %d셀] 용량만 %.3f | M2b %.3f" % (len(nh), lerr(loco(f_cap, g1, nh), nh.life.values),
                                              lerr(loco(f_ab("x"), g_ab("x", nh), nh), nh.life.values)))
rng = np.random.default_rng(0); e_real = lerr(res["③ M2b 사후 (a + b·x, EIS)"], D.life.values); perm = []
for i in range(N_PERM):
    Dp = D.copy(); Dp["x"] = rng.permutation(D.x.values)
    perm.append(lerr(loco(f_ab("x"), g_ab("x", Dp), Dp), Dp.life.values))
    print(f"\r[섞은 값 시험] {i + 1}/{N_PERM}", end="", flush=True)
perm = np.array(perm)
print("\r[섞은 값 시험] x를 섞은 %d회: 중앙값 %.3f, 최소 %.3f | 실제 %.3f 이하인 비율 %.2f"
      % (N_PERM, np.median(perm), perm.min(), e_real, np.mean(perm <= e_real)))

# ---------- 5. 저장 ----------
for k, p in res.items(): D["pred " + k[:2]] = p
D.to_csv(OUT / "P1_pred.csv", index=False, float_format="%.6g")   # 반올림 금지: A·B는 1e-4 크기
fig, ax = plt.subplots(figsize=(6, 6))
for k, mk in [("① 용량만 (f 공통)", "s"), ("③ M2b 사후 (a + b·x, EIS)", "o")]:
    ax.scatter(D.life, res[k], marker=mk, label=k, alpha=.8)
for _, r in D.iterrows(): ax.annotate(f"{r.cell:02d}{r.temp}", (r.life, res["③ M2b 사후 (a + b·x, EIS)"][_]), fontsize=8)
lim = [80, 4000]; ax.plot(lim, lim, "k--", lw=1); ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim,
      xlabel="실제 수명 [h] (SOH 0.85)", ylabel="예측 수명 [h]", title="P1 프로토타입: 초반 1점 + 2항 수식")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "P1_pred.png", dpi=150)
print("\n저장:", OUT / "P1_pred.csv", ",", OUT / "P1_pred.png")
