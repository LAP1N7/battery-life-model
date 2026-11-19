"""P8: DRT_shift 봉우리의 SOC 구간별 증가 <-> DVA 열화 모드(LLI, LAM_PE, LAM_NE) 대조 (트랙 C 결과 사용)

[사전 등록 2026-11-19 14:55, 업무 전달서 기록 — DVA 결과를 보기 전에 고정]
 입력: 트랙 C D1_dva.csv (셀별 LLI·LAM_PE·LAM_NE, C의 기준 check-up 대비),  P7 ①의 R_main (lambda = 1e-2)
 시점: 셀마다 P7의 '마지막' check-up (SOH >= 0.80). D1에서 같은 check-up(+-1 Ah)의 모드를 씀
 DRT 쪽: SOC 구간별 R_main 증가율 = (마지막 / 첫 check-up) - 1,  구간 = 낮음 10-20 % / 중간 40-60 % / 높음 80-90 % (구간 안 SOC 평균)
 혼동 통제: 용량 손실 Q = 1 - SOH(마지막 check-up)를 통제한 편순위상관 (순위를 Q 순위에 회귀한 잔차끼리의 상관)
 예측: (1) 높음 구간 <-> LAM_NE  rho >= 0.5, p < 0.05
       (2) 낮음 구간 <-> LAM_PE  rho >= 0.5, p < 0.05
       (3) 높음 구간은 LAM_PE보다 LAM_NE와 더 강하게 (보고)
 판정: (1)(2) 모두 = U자의 양극·음극 분해 탐색적 지지 / 하나 = 부분 / 없음 = 기각
 보고 전용: LLI와의 상관, 고주파 총저항(R_HF + R_inf) 증가율과의 상관, 통제 없는 단순 상관

[예상 문제 — 실행 전에 적어 둠]
 - 트랙 C 결과: 같은 온도 안에서 세 모드가 거의 같이 움직임(r 0.85~0.995) -> 모드별 연결을 구분하지 못할 수 있음. 모드끼리의 상관을 함께 출력
 - LAM_PE는 식별성이 낮음 (C의 이전 탐색 대비 순위상관 0.18) -> (2)는 약한 시험
 - 45 °C 셀은 C의 기준 check-up이 BOL이 아님 (47~294 Ah). P7의 첫 check-up과 기준이 다를 수 있음
 - 셀 수 약 21, 경계값 0.5는 제안값
[그래도 남는 것] 모드끼리의 공선성 수치와 구간별 증가율 표 -> KIT 시험 설계 근거

입력: out/P7_soc_spectra.csv (stepP7_soc.py), out/D/D1_dva.csv (트랙 C, 브랜치 track-c-mechanism 8bdf52c에서 가져옴)
출력: out/P8_drt_dva.csv, out/P8_drt_dva.png
변경 이력: 2026-11-19 최초 작성 (Claude). Claude 환경에서 같은 입력으로 실행 확인
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata, t as tdist
from common import OUT, TEMP

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
LAM = 1e-2; BANDS = {"낮음": [10, 20], "중간": [40, 50, 60], "높음": [80, 90]}; RHO = 0.5; PMAX = 0.05

S = pd.read_csv(OUT / "P7_soc_spectra.csv"); S = S[S.lam == LAM]
M = pd.read_csv(OUT / "D" / "D1_dva.csv")
rows = []
for cell, g in S.groupby("cell"):
    rec = dict(cell=cell, temp=TEMP[cell])
    for band, socs in BANDS.items():
        f = g[(g.tag == "first") & g.soc.isin(socs)].R_main.mean(); l = g[(g.tag == "last") & g.soc.isin(socs)].R_main.mean()
        rec[f"dR_{band}"] = l / f - 1
    f = g[g.tag == "first"].R_HF_inf.mean(); l = g[g.tag == "last"].R_HF_inf.mean(); rec["dR_HF_tot"] = l / f - 1
    ak = g[g.tag == "last"].ahk.iloc[0]; rec["ah_last"] = ak; rec["Q"] = 1 - g[g.tag == "last"].soh.iloc[0]
    m = M[(M.cell == cell) & ((M.ah - ak).abs() < 1)]
    if len(m):
        m = m.iloc[0]
        rec.update(LLI=m.LLI, LAM_PE=m.LAM_PE, LAM_NE=m.LAM_NE, ref_not_bol=bool(m.flag_ref_not_bol), dva_rmse_mV=m.rmse_mV)
    rows.append(rec)
D = pd.DataFrame(rows).dropna(subset=["LLI"])
print(f"[데이터] 셀 {len(D)}개 (DVA 짝 없음 {len(rows) - len(D)}개), 기준이 BOL 아닌 셀 {int(D.ref_not_bol.sum())}개")

def pcorr(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r = np.corrcoef(ex, ey)[0, 1]; df = len(x) - 3
    return r, 2 * tdist.sf(abs(r) * np.sqrt(df / max(1 - r * r, 1e-12)), df)

print("\n[모드끼리 상관 (예상 문제 점검)]  순위상관:  " + ",  ".join(
    f"{a}-{b} {spearmanr(D[a], D[b])[0]:.2f}" for a, b in [("LLI", "LAM_NE"), ("LLI", "LAM_PE"), ("LAM_PE", "LAM_NE")]))
print("   용량 손실 Q 통제 후:  " + ",  ".join(
    f"{a}-{b} {pcorr(D[a], D[b], D.Q)[0]:.2f}" for a, b in [("LLI", "LAM_NE"), ("LLI", "LAM_PE"), ("LAM_PE", "LAM_NE")]))

tab = []
for band in BANDS:
    for mode in ["LAM_NE", "LAM_PE", "LLI"]:
        r, p = pcorr(D[f"dR_{band}"], D[mode], D.Q); r0, p0 = spearmanr(D[f"dR_{band}"], D[mode])
        tab.append(dict(band=band, mode=mode, rho_partial=r, p_partial=p, rho_simple=r0, p_simple=p0))
for mode in ["LAM_NE", "LAM_PE", "LLI"]:
    r, p = pcorr(D.dR_HF_tot, D[mode], D.Q); r0, p0 = spearmanr(D.dR_HF_tot, D[mode])
    tab.append(dict(band="고주파 총저항", mode=mode, rho_partial=r, p_partial=p, rho_simple=r0, p_simple=p0))
T = pd.DataFrame(tab)
print("\n[편순위상관 (Q 통제)]  행 = R_main 증가율 구간, 괄호 = 통제 없는 단순 상관")
for _, r in T.iterrows():
    print(f"  {r.band:6s} <-> {r['mode']:6s}: rho = {r.rho_partial:+.2f} (p={r.p_partial:.3f})   [단순 {r.rho_simple:+.2f}, p={r.p_simple:.3f}]")

g = lambda b, m: T[(T.band == b) & (T["mode"] == m)].iloc[0]
c1 = g("높음", "LAM_NE"); c2 = g("낮음", "LAM_PE")
p1 = bool(c1.rho_partial >= RHO and c1.p_partial < PMAX); p2 = bool(c2.rho_partial >= RHO and c2.p_partial < PMAX)
p3 = bool(c1.rho_partial > g("높음", "LAM_PE").rho_partial)
grade = "탐색적 지지" if (p1 and p2) else ("부분" if (p1 or p2) else "기각")
print(f"\n[판정] (1) 높음<->LAM_NE {'충족' if p1 else '미충족'} / (2) 낮음<->LAM_PE {'충족' if p2 else '미충족'} -> {grade}")
print(f"       (3) 높음 구간: LAM_NE {c1.rho_partial:+.2f} vs LAM_PE {g('높음', 'LAM_PE').rho_partial:+.2f} -> {'LAM_NE 쪽이 강함' if p3 else 'LAM_NE 쪽이 강하지 않음'} (보고)")

D.to_csv(OUT / "P8_drt_dva.csv", index=False, float_format="%.6g")
col = {"H": "tab:red", "R": "tab:gray", "C": "tab:blue"}; NAME = {"H": "45 °C 사이클", "R": "25 °C 사이클", "C": "-5 °C 사이클"}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
for k, (band, mode) in enumerate([("높음", "LAM_NE"), ("낮음", "LAM_PE"), ("중간", "LLI")]):
    for t in "HRC":
        s = D[D.temp == t]; ax[k].scatter(s[mode] * 100, s[f"dR_{band}"] * 100, color=col[t], label=NAME[t])
    r = g(band, mode)
    ax[k].set(xlabel=f"{mode} [%]", ylabel=f"R_main 증가율, {band} SOC [%]", title=f"{band} SOC <-> {mode}: Q 통제 rho {r.rho_partial:+.2f}")
ax[0].legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "P8_drt_dva.png", dpi=150); plt.close(fig)
print("\n저장:", OUT / "P8_drt_dva.csv", OUT / "P8_drt_dva.png")
