"""P7 ①: SOC 의존성으로 DRT 봉우리 정체 좁히기 — 고주파(측정 범위 밖) 저항은 접촉 저항인가, DRT_shift 봉우리는 양극 전하전달인가

[합의된 사양 2026-11-19 확정]
 대상: 21셀, 첫·마지막 check-up (SOH >= 0.80, SOC 10~90% 9개가 모두 있는 check-up), DRT 설정은 P6와 동일 (drt_core.py)
 지표 (스펙트럼마다): R_HF  = log10 tau < -4.8 구간 저항 (10 kHz보다 빠름, 외삽값)
                      R_inf = DRT의 순수 저항,  R_HF + R_inf 도 함께 보고 (외삽 구간과 R_inf가 값을 주고받을 수 있어서)
                      R_main = log10 tau -3 ~ -0.7 구간 저항 (DRT_shift 봉우리),  tau_main = 그 구간의 면적 가중 중심
 SOC 의존 정도: 셀 x check-up 마다 SOC 10~90% 의 최대/최소 비
 판정 단위: 셀 x check-up (최대 42개).  lambda = 1e-2 로 판정, 3e-3 / 3e-2 는 강건성 보고
 예측 A (R_HF = 접촉 저항): 최대/최소 비 <= 1.5 인 단위가 70% 이상 -> 접촉 저항 또는 표면막 해석 지지
 예측 B (DRT_shift 봉우리 = 양극 전하전달): R_main 최대/최소 비 >= 2 이고 최대값이 SOC 10 또는 20% 인 단위가 70% 이상 -> 지지
   근거: 전체 셀의 전하전달은 낮은 SOC에서 양극이 지배 (Wildfeuer et al., Iurilli 2022에서 인용)
 반증: A - R_HF가 SOC에 따라 2배 이상 변하는 단위가 다수 / B - R_main이 SOC에 둔감
 경계값(1.5배, 2배)은 Claude 제안값이며 문헌 근거 없음. 결과를 본 뒤 바꾸지 않음
 보고 전용(판정 미사용): 최대/최소 비의 분포, 여러 경계값에서의 충족 비율 (KIT 경계값 결정용), 첫 -> 마지막 SOC 프로필 변화
 온도: 한 check-up의 SOC 9개는 같은 시기에 측정 -> SOC 비교는 check-up 간 온도 잡음 영향이 작음. 셀 13의 36.6 °C check-up은 SOC 9개가 없어 원래 제외됨

[예상 문제] R_HF는 외삽값 / 높은 SOC에서 음극 전하전달이 주 봉우리 구간에 겹칠 수 있음 / SOC 10%는 확산 오염 가능
[plan B] A 기각 시 인덕턴스 처리 비교로 artifact 점검.  [그래도 남는 것] 봉우리별 SOC 프로필 표 (KIT 사전 등록 정의)

입력: ../data/XX/Resistance.csv, out/P3_ratio.csv (check-up별 SOH)
출력: out/P7_soc_spectra.csv, out/P7_soc_units.csv, out/P7_soc.png
변경 이력: 2026-11-19 최초 작성 (Claude). Claude 환경에서 같은 입력으로 실행 확인
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from common import OUT, TEMP, load_eis, SOCS
from drt_core import drt, TAU

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
LAMS = [3e-3, 1e-2, 3e-2]; LAM_J = 1e-2; A_MAX = 1.5; B_MIN = 2.0; NEED = 0.70
LT = np.log10(TAU); HF = LT < -4.8; MAIN = (LT >= -3) & (LT <= -0.7)
NAME = {"H": "45 °C 사이클", "R": "25 °C 사이클", "C": "-5 °C 사이클"}

R3 = pd.read_csv(OUT / "P3_ratio.csv"); R3["ahk"] = R3.ah.round(1)
rows = []
for cell in range(1, 22):
    spec = load_eis(cell); soh = R3[R3.cell == cell].groupby("ahk").soh.first()
    aks = sorted({round(a, 1) for (a, s) in spec}); aks = [a for a in aks if a in soh.index and soh[a] >= 0.80]
    if len(aks) < 2: print(f"cell {cell:02d}: 대상 check-up 부족 ({len(aks)})"); continue
    for tag, ak in [("first", aks[0]), ("last", aks[-1])]:
        a_true = [a for (a, s) in spec if round(a, 1) == ak][0]
        for soc in SOCS:
            d = spec[(a_true, soc)]
            for lam in LAMS:
                g, Ri, L, res = drt(d[:, 0], d[:, 3], d[:, 4], lam)
                Rm = g[MAIN].sum()
                rows.append(dict(cell=cell, temp=TEMP[cell], tag=tag, ahk=ak, soh=soh[ak], soc=soc, lam=lam, res=res,
                                 R_HF=g[HF].sum(), R_inf=Ri, R_main=Rm, logtau_main=np.sum(LT[MAIN] * g[MAIN]) / Rm if Rm > 0 else np.nan))
S = pd.DataFrame(rows); S["R_HF_inf"] = S.R_HF + S.R_inf
S.to_csv(OUT / "P7_soc_spectra.csv", index=False, float_format="%.6g")

def units(lam):
    out = []
    for (cell, tag), g in S[S.lam == lam].groupby(["cell", "tag"]):
        g = g.set_index("soc").sort_index()
        r = lambda c: g[c].max() / max(g[c].min(), 1e-12)
        out.append(dict(cell=cell, temp=TEMP[cell], tag=tag, ratio_HF=r("R_HF"), ratio_HF_inf=r("R_HF_inf"), ratio_main=r("R_main"),
                        argmax_main=int(g.R_main.idxmax()), argmax_HF=int(g.R_HF.idxmax())))
    U = pd.DataFrame(out)
    U["A_ok"] = U.ratio_HF <= A_MAX
    U["B_ok"] = (U.ratio_main >= B_MIN) & U.argmax_main.isin([10, 20])
    return U

U = units(LAM_J); n = len(U)
fa, fb = U.A_ok.mean(), U.B_ok.mean()
print(f"[데이터] 셀 {U.cell.nunique()}개, 판정 단위(셀 x check-up) {n}개, 스펙트럼 {len(S[S.lam == LAM_J])}개, 상대 잔차 중앙값 {S.res.median():.4f}")
print(f"\n[예측 A: R_HF = 접촉 저항] SOC 최대/최소 비 <= {A_MAX}: {int(U.A_ok.sum())} / {n} ({fa:.0%}, 기준 70%) -> {'지지' if fa >= NEED else '미충족'}")
print(f"   R_HF 비 중앙값 {U.ratio_HF.median():.2f} (범위 {U.ratio_HF.min():.2f}~{U.ratio_HF.max():.2f}),  R_HF+R_inf 비 중앙값 {U.ratio_HF_inf.median():.2f}")
print(f"[예측 B: DRT_shift 봉우리 = 양극 전하전달] 비 >= {B_MIN} 이고 최대가 SOC 10/20%: {int(U.B_ok.sum())} / {n} ({fb:.0%}, 기준 70%) -> {'지지' if fb >= NEED else '미충족'}")
print(f"   R_main 비 중앙값 {U.ratio_main.median():.2f} (범위 {U.ratio_main.min():.2f}~{U.ratio_main.max():.2f}),  최대 SOC 분포 {U.argmax_main.value_counts().sort_index().to_dict()}")
for t in "HRC":
    s = U[U.temp == t]; print(f"   {NAME[t]}: A {int(s.A_ok.sum())}/{len(s)}, B {int(s.B_ok.sum())}/{len(s)}")
for lam in LAMS:
    if lam == LAM_J: continue
    V = units(lam); print(f"  [강건성] lambda={lam:g}: A {V.A_ok.mean():.0%}, B {V.B_ok.mean():.0%}")

print("\n[보고 전용, 판정 미사용] 경계값별 충족 비율 (KIT 경계값 결정용)")
print("  A (R_HF 비 <= x):  " + "  ".join(f"x={x}: {(U.ratio_HF <= x).mean():.0%}" for x in [1.2, 1.5, 2, 3]))
print("  B (R_main 비 >= x & 최대 SOC 10/20): " + "  ".join(f"x={x}: {((U.ratio_main >= x) & U.argmax_main.isin([10, 20])).mean():.0%}" for x in [1.5, 2, 3, 5]))
U.to_csv(OUT / "P7_soc_units.csv", index=False, float_format="%.6g")

# ---------- 그림: SOC 프로필 (lambda 1e-2), 사이클 온도별 중앙값 ----------
P = S[S.lam == LAM_J]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
col = {"H": "tab:red", "R": "tab:gray", "C": "tab:blue"}
for k, (c, ttl) in enumerate([("R_HF", "R_HF (log τ < -4.8, 측정 범위 밖)"), ("R_main", "R_main (DRT_shift 봉우리)"), ("logtau_main", "DRT_shift 봉우리 위치 log10 τ")]):
    for t in "HRC":
        for tag, ls in [("first", "--"), ("last", "-")]:
            m = P[(P.temp == t) & (P.tag == tag)].groupby("soc")[c].median()
            ax[k].plot(m.index, m.values * (1e3 if c != "logtau_main" else 1), ls, color=col[t], label=f"{NAME[t]} {'첫' if tag == 'first' else '마지막'}")
    ax[k].set(xlabel="SOC [%]", ylabel="[mΩ]" if c != "logtau_main" else "log10 τ [s]", title=ttl)
ax[0].legend(fontsize=7); fig.suptitle("SOC 프로필 (셀 중앙값, 점선 = 첫 check-up, 실선 = 마지막)"); fig.tight_layout()
fig.savefig(OUT / "P7_soc.png", dpi=150); plt.close(fig)
print("\n저장:", OUT / "P7_soc_spectra.csv", OUT / "P7_soc_units.csv", OUT / "P7_soc.png")
