"""P5: 비이상성(n2) 검증 — 위상 피크 변화를 '반원 크기 비율 r' 말고 '반원 모양 n2'가 추가로 설명하는가

[합의된 사양 2026-11-19 확정]  (H-φ2: 열화 정도·비이상성에 따라 위상 범위와 변화 속도가 달라진다)
 대상: 저 SOC(10~30%), MATLAB ECM separable 스펙트럼 중 n2 흔들림(spreadN2) <= 0.2.  n2가 경계(0.5 또는 1)에 붙은 수는 보고
 (a) 위상 피크 예측: 모델 Z = R0 + R_arc / (1 + (j*w*tau)^n) 에서 |phi|max 를 수치로 계산 (tau는 피크 깊이에 영향 없음)
       n = 1 (이상적 반원)  vs  n = n2 (ECM 값)  -> 실측 피크 깊이(1 Hz~1 kHz)와의 RMSE 비교
       성립 기준: n2 를 넣었을 때 RMSE 가 20% 이상 줄어듦
 (b) 셀별 변화 속도(시간 기울기, SOH >= 0.80, check-up 3개 이상):  d|phi|/dt, dn2/dt, dr/dt
       r 의 영향을 뺀 편순위상관 rho(d|phi|/dt, dn2/dt | dr/dt)
       성립 기준: rho >= 0.5 AND p < 0.05   (부호: n2 감소 -> 반원 납작 -> |phi| 감소 이므로 양의 상관이 기대 방향)
 (c) 탐색(판정 없음): 초기 위상 변화 속도가 P1의 선형 비중(frac)과 관련 있는가 — P1의 x(초기 저항 증가율)와 나란히 순위상관만 보고
 판정 등급(11-19 합의): (a)(b) 모두 충족 = 탐색적 지지 / 하나만 = 부분 / 둘 다 미충족 = 기각.  '확증'은 KIT에서

[예상 문제 — 미리 적어 둠]
 - 고온 셀은 tau1이 10 kHz 경계에 붙어 separable 스펙트럼이 적음 (셀 06·09: 3개, 15: 0개) -> (b)에서 고온 셀이 많이 빠질 수 있음
   -> 빠진 셀 수를 출력. plan B: (b)를 separable 조건 없이 n2 흔들림 조건만으로 다시 보는 '보조 결과'를 함께 출력 (판정 미사용)
 - R_arc는 반원 2개를 합친 지름인데 n2는 두 번째 반원의 값 -> (a)는 근사. 개선이 작으면 '모델 단순화 한계'일 수 있음
 - 셀 수가 적어(<= 14) 편상관 p값이 크게 나오기 쉬움
[그래도 남는 것] 기각돼도 n2의 온도별 변화 방향(고온에서 줄어드는가)은 기술 통계로 남김

입력: out/P3_ratio.csv (stepP3_arc.py, 반원 끝 = 마지막 국소 최소), ../matlab/out/ecm_matlab.csv (stepM1_ecm.m), out/P1_pred.csv
출력: out/P5_spectra.csv, out/P5_cells.csv (판정용), out/P5_cells_aux.csv (보조), out/P5_ideality.png
변경 이력: 2026-11-19 최초 작성 (Claude). 실행 결과는 Claude 환경에서 같은 입력으로 확인함
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata, t as tdist
from common import OUT, TEMP

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
LOW = [10, 20, 30]; SPREAD_MAX = 0.2
HOT = [c for c in range(1, 22) if TEMP[c] == "H"]; COLD = [c for c in range(1, 22) if TEMP[c] == "C"]

# ---------- 입력 합치기 ----------
R = pd.read_csv(OUT / "P3_ratio.csv")
E = pd.read_csv(OUT.parent.parent / "matlab" / "out" / "ecm_matlab.csv")
for d in (R, E): d["ahk"] = d.ah.round(1)
R = R[R.soc.isin(LOW) & (R.found == 1) & (R.r > 0)]
M = R.merge(E[["cell", "ahk", "soc", "n2", "spreadN2", "n2AtBound", "separable"]], on=["cell", "ahk", "soc"], how="inner")
print(f"[데이터] 저 SOC 반원 검출 스펙트럼 {len(R)}개 -> ECM과 짝지어짐 {len(M)}개")

# ---------- (a) 위상 피크 예측: n=1 vs n=n2 ----------
X = np.logspace(-4, 4, 3000)                       # w*tau (무차원)
def phimax(r, n):
    """Z = 1 + r/(1+(jX)^n)  (R0로 정규화) 의 최대 위상 [deg, 음수]"""
    z = 1 + r / (1 + (1j * X) ** n)
    return -np.degrees(np.max(np.arctan2(-z.imag, z.real)))
M["phi_n1"] = [phimax(r, 1.0) for r in M.r]
M["phi_n2"] = [phimax(r, n) for r, n in zip(M.r, M.n2)]
A = M[(M.separable == 1) & (M.spreadN2 <= SPREAD_MAX)].copy()
rm1 = np.sqrt(np.mean((A.phi_n1 - A.pk_depth) ** 2)); rm2 = np.sqrt(np.mean((A.phi_n2 - A.pk_depth) ** 2))
imp = 1 - rm2 / rm1
pass_a = imp >= 0.20
print(f"\n[(a)] 대상 {len(A)}개 (separable & n2 흔들림 <= {SPREAD_MAX}), n2 경계에 붙음 {int(A.n2AtBound.sum())}개")
print(f"  RMSE  n=1: {rm1:.2f}°   n=n2: {rm2:.2f}°   개선 {imp:.0%}  (기준 20% 이상) -> {'충족' if pass_a else '미충족'}")
for t, nm in [("H", "고온"), ("R", "상온"), ("C", "저온")]:
    s = A[A.temp == t]
    if len(s): print(f"    {nm}: {len(s)}개  RMSE n=1 {np.sqrt(np.mean((s.phi_n1-s.pk_depth)**2)):.2f}° -> n2 {np.sqrt(np.mean((s.phi_n2-s.pk_depth)**2)):.2f}°")

# ---------- (b) 셀별 변화 속도와 편순위상관 ----------
def slope(t, y):
    m = np.isfinite(t) & np.isfinite(y)
    return np.polyfit(t[m], y[m], 1)[0] * 100 if m.sum() >= 3 else np.nan     # [단위/100 h]

def cell_rates(D):
    rows = []
    for cell, g in D[D.soh >= 0.80].groupby("cell"):
        k = g.groupby("ahk").agg(hours=("hours", "first"), phi=("pk_depth", "mean"), n2=("n2", "mean"), r=("r", "mean")).dropna()
        rows.append(dict(cell=cell, temp=TEMP[cell], n_ck=len(k), dphi=slope(k.hours.values, -k.phi.values),   # |phi| 기울기
                         dn2=slope(k.hours.values, k.n2.values), dr=slope(k.hours.values, k.r.values)))
    return pd.DataFrame(rows).dropna(subset=["dphi", "dn2", "dr"])

def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    rho = np.corrcoef(ex, ey)[0, 1]; df = len(x) - 3
    p = 2 * tdist.sf(abs(rho) * np.sqrt(df / max(1 - rho ** 2, 1e-12)), df) if df > 0 else np.nan
    return rho, p

def report_b(P, label):
    hot = int(P.cell.isin(HOT).sum()); cold = int(P.cell.isin(COLD).sum())
    rho, p = partial_spearman(P.dphi.values, P.dn2.values, P.dr.values) if len(P) >= 5 else (np.nan, np.nan)
    rho0 = spearmanr(P.dphi, P.dn2)[0] if len(P) >= 5 else np.nan
    print(f"  {label}: 셀 {len(P)}개 (고온 {hot}/7, 저온 {cold}/7)  단순 rho(d|phi|, dn2) = {rho0:.2f}  "
          f"r 제외 편상관 rho = {rho:.2f} (p={p:.3f})")
    return rho, p

print("\n[(b)] 셀별 변화 속도 (SOH >= 0.80, check-up 3개 이상)")
Pb = cell_rates(A)
rho_b, p_b = report_b(Pb, "판정용 (separable & 흔들림 <= 0.2)")
pass_b = bool(rho_b >= 0.5 and p_b < 0.05)
print(f"  -> 기준 rho >= 0.5 AND p < 0.05 : {'충족' if pass_b else '미충족'}")
Pb2 = cell_rates(M[M.spreadN2 <= SPREAD_MAX])
report_b(Pb2, "보조 (separable 조건 없이, 판정 미사용)")

# ---------- 기술 통계: n2 변화 방향 ----------
print("\n[기술 통계] 셀별 n2 기울기 부호 (보조 데이터 기준)  감소 / 증가")
for t, nm in [("H", "고온"), ("R", "상온"), ("C", "저온")]:
    s = Pb2[Pb2.temp == t]; print(f"  {nm}: {int((s.dn2 < 0).sum())} / {int((s.dn2 > 0).sum())}  (셀 {len(s)}개)")

# ---------- (c) 탐색: 초기 위상 변화 속도 vs P1 선형 비중 ----------
P1 = pd.read_csv(OUT / "P1_pred.csv")[["cell", "te", "x", "frac"]]
Rall = pd.read_csv(OUT / "P3_ratio.csv"); Rall = Rall[Rall.soc.isin(LOW)]
early = []
for _, row in P1.iterrows():
    g = Rall[(Rall.cell == row.cell) & (Rall.hours <= row.te + 1e-6)].groupby("ah").agg(h=("hours", "first"), phi=("pk_depth", "mean")).dropna()
    early.append(dict(cell=row.cell, dphi_early=slope(g.h.values, -g.phi.values) if len(g) >= 3 else
                      ((-g.phi.iloc[-1] + g.phi.iloc[0]) / (g.h.iloc[-1] - g.h.iloc[0]) * 100 if len(g) == 2 and g.h.iloc[-1] > g.h.iloc[0] else np.nan)))
C = P1.merge(pd.DataFrame(early), on="cell").dropna()
bias1 = (A.phi_n1 - A.pk_depth).mean(); bias2 = (A.phi_n2 - A.pk_depth).mean()
print(f"  참고: 평균 치우침(예측-실측) n=1 {bias1:.2f}°, n=n2 {bias2:.2f}°  (n2는 같은 스펙트럼에서 피팅한 값이라 (a)는 일부 순환적)")
print(f"\n[(c) 탐색, 판정 없음] 셀 {len(C)}개 — P1 선형 비중(frac)과의 순위상관")
print(f"  초기 저항 증가율 x : rho = {spearmanr(C.x, C.frac)[0]:.2f}   초기 |phi| 변화율 : rho = {spearmanr(C.dphi_early, C.frac)[0]:.2f}")

# ---------- 종합 ----------
grade = "탐색적 지지" if (pass_a and pass_b) else ("부분" if (pass_a or pass_b) else "기각")
print(f"\n[P5 종합] (a) {'충족' if pass_a else '미충족'} / (b) {'충족' if pass_b else '미충족'} -> {grade}  ('확증'은 KIT에서)")

# ---------- 저장·그림 ----------
M.to_csv(OUT / "P5_spectra.csv", index=False, float_format="%.6g")
Pb.to_csv(OUT / "P5_cells.csv", index=False, float_format="%.6g"); Pb2.to_csv(OUT / "P5_cells_aux.csv", index=False, float_format="%.6g")
col = {"H": "tab:red", "R": "tab:gray", "C": "tab:blue"}; name = {"H": "고온 45 °C", "R": "상온 25 °C", "C": "저온 -5 °C"}
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
lim = [A[["pk_depth", "phi_n1", "phi_n2"]].min().min() - 1, 0]
for i, (c, ttl) in enumerate([("phi_n1", f"n = 1 (이상적)  RMSE {rm1:.2f}°"), ("phi_n2", f"n = n2  RMSE {rm2:.2f}°")]):
    for t in "HRC":
        s = A[A.temp == t]; ax[i].scatter(s[c], s.pk_depth, s=8, alpha=.5, color=col[t], label=name[t])
    ax[i].plot(lim, lim, "k--", lw=1); ax[i].set(xlim=lim, ylim=lim, xlabel="예측 피크 [deg]", ylabel="실측 피크 [deg]", title=f"(a) {ttl}")
ax[0].legend(fontsize=8)
for t in "HRC":
    s = Pb2[Pb2.temp == t]; ax[2].scatter(s.dn2, s.dphi, color=col[t], label=name[t])
for _, rr in Pb2.iterrows(): ax[2].annotate(f"{int(rr.cell):02d}", (rr.dn2, rr.dphi), fontsize=7)
ax[2].axhline(0, color="k", lw=.5); ax[2].axvline(0, color="k", lw=.5)
ax[2].set(xlabel="dn2/dt [1/100 h]", ylabel="d|phi|/dt [deg/100 h]", title=f"(b) 셀별 변화 속도 (보조 데이터)  편상관(판정용) {rho_b:.2f}")
ax[2].legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "P5_ideality.png", dpi=150); plt.close(fig)
print("\n저장:", OUT / "P5_spectra.csv", OUT / "P5_cells.csv", OUT / "P5_cells_aux.csv", OUT / "P5_ideality.png")
