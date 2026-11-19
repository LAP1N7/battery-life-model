"""R1: 충전 속도 특성 열화 <-> 음극 활물질 손실(LAM_NE) 사전 등록 시험 (H-R)

[사전 등록 2026-11-19 오후 확정, 업무 전달서 기록 — 속도 시험 결과를 보기 전에 고정]
 데이터: Mohtat data/XX/Crate_wExpansion.csv. check-up마다 C/10(0.5 A)·C/5(1 A)·C/2(2.5 A)·1C(5 A) CC 충전, 각 충전 뒤 C/2 방전
 지표: r = Q충전(1C) / Q충전(C/10),  dr = r(마지막 check-up) - r(첫 check-up)
 H-R1: LAM_NE가 큰 셀일수록 dr이 더 떨어진다. 용량 손실 Q = 1 - SOH 통제 편순위상관 rho <= -0.4, 단측 p < 0.05
 H-R2 (핵심): R0 증가에 따른 옴 강하를 보정해도 H-R1이 유지된다.
       보정 = 첫 check-up 곡선에서 충전 상한 전압을 I x dR0만큼 내려 다시 읽은 용량으로 r(첫)을 계산 (마지막 check-up은 그대로)
 등급: H-R1·H-R2 모두 = 탐색 이상(확증 미만) / H-R1만 = 저항 효과로 설명 가능 / 둘 다 아님 = 기각
 제외: 측정 온도(중앙값) 30 °C 초과 check-up.  시점 맞춤: 누적 Ah 최근접, 안 되면 보간 (판정 기준은 그대로)
 C가 이 파일을 본 적이 있으면 비블라인드로 표기 (C 답변 대기)

[구현 세부 — 결과 보기 전에 고정, 사양표에 없던 부분]
 - 충전 구간: 전류가 명목값(500/1000/2500/5000 mA)의 ±8 % 안에 연속으로 있는 가장 긴 구간, 0.1 h 이상. Q = 전류 적분
 - 충전 상한 전압 V_cut: 첫 check-up의 같은 속도 구간 최대 전압
 - R0: P5_spectra.csv의 R0 (SOC 10~90 % 중앙값). Crate check-up의 누적 Ah에 가장 가까운 EIS check-up(차이 <= 5 % 또는 5 Ah), 없으면 누적 Ah 선형 보간
 - dR0 < 0 이면 0으로 둠 (보정 안 함, 개수 보고)
 - LAM_NE와 Q(=1-SOH): 트랙 C D1_dva.csv에서 Crate 마지막 check-up과 누적 Ah가 가장 가까운 행 (같은 조건)
 - 마지막 check-up = 네 속도 충전이 모두 있고 30 °C 이하인 가장 늦은 check-up

[예상 문제]
 - 저항이 늘면 CC 충전량이 주는 것은 당연 -> H-R2로 분리
 - C의 LAM_NE는 비블라인드, 용량 정합 통과 24 % -> 순위 잡음
 - LAM_NE와 LLI가 같이 움직임(Q 통제 후 0.02로 분리됨, P8) -> Q 통제로 대응
[그래도 남는 것] 온도별 속도 특성 열화 표, 수명식 후보 지표

입력: data/XX/Crate_wExpansion.csv, out/P5_spectra.csv (stepP5_ideality.py), out/D/D1_dva.csv (트랙 C)
출력: out/R1_rate_checkups.csv, out/R1_rate_cells.csv, out/R1_rate.png
변경 이력: 2026-11-19 최초 작성 (Claude). 합성 데이터로만 시험함, 실제 데이터 미실행
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata, t as tdist
from common import DATA, OUT, TEMP

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
RATES = {"C10": 500, "C5": 1000, "C2": 2500, "1C": 5000}   # mA
TOL = 0.08; MIN_H = 0.1; TMAX = 30.0; RHO = -0.4; PMAX = 0.05


def charges(ck):
    """한 check-up에서 속도별 CC 충전 구간 -> {rate: (Q[Ah] 누적 배열, V 배열, Vmax, Q총)}"""
    t = ck["Time [s]"].to_numpy(float); I = ck["Current [mA]"].to_numpy(float); V = ck["Voltage [V]"].to_numpy(float)
    out = {}
    for name, nom in RATES.items():
        m = np.abs(I - nom) <= TOL * nom
        if not m.any():
            continue
        idx = np.flatnonzero(m); br = np.flatnonzero(np.diff(idx) > 1)
        runs = np.split(idx, br + 1); run = max(runs, key=lambda r: t[r[-1]] - t[r[0]])
        if (t[run[-1]] - t[run[0]]) / 3600 < MIN_H:
            continue
        tt, ii, vv = t[run], I[run], V[run]
        q = np.concatenate([[0], np.cumsum(0.5 * (ii[1:] + ii[:-1]) * np.diff(tt))]) / 3.6e6   # mA*s -> Ah
        out[name] = (q, vv, vv.max(), q[-1])
    return out


def q_at(q, v, vlim):
    """전압이 처음 vlim에 닿을 때까지의 충전량"""
    k = np.flatnonzero(v >= vlim)
    return q[-1] if not len(k) else q[k[0]]


def pcorr(x, y, z):   # stepP8과 같은 정의
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r = np.corrcoef(ex, ey)[0, 1]; df = len(x) - 3
    return r, 2 * tdist.sf(abs(r) * np.sqrt(df / max(1 - r * r, 1e-12)), df)


def r0_at(P, cell, ah):
    g = P[P.cell == cell].groupby("ah").R0.median().sort_index()
    if g.empty:
        return np.nan, "없음"
    d = (g.index.to_numpy() - ah); k = np.argmin(np.abs(d))
    if abs(d[k]) <= max(5.0, 0.05 * max(ah, 1)):
        return g.iloc[k], "최근접"
    if g.index.min() <= ah <= g.index.max():
        return float(np.interp(ah, g.index, g.values)), "보간"
    return np.nan, "범위 밖"


P = pd.read_csv(OUT / "P5_spectra.csv"); M = pd.read_csv(OUT / "D" / "D1_dva.csv")
P["cell"] = P.cell.astype(int); M["cell"] = M.cell.astype(int)
ck_rows, cell_rows = [], []
for cell in range(1, 22):
    f = DATA / f"{cell:02d}" / "Crate_wExpansion.csv"
    if not f.exists():
        print(f"[경고] 셀 {cell:02d}: Crate 파일 없음"); continue
    d = pd.read_csv(f)
    valid = []
    for cyc, ck in d.groupby("Cycle number"):
        ah = ck["Total Ah throughput [Ah]"].iloc[0]; T = ck["Temperature [C]"].median(); ch = charges(ck)
        ok = len(ch) == 4 and T <= TMAX
        ck_rows.append(dict(cell=cell, cyc=cyc, ah=ah, T_med=T, n_rates=len(ch), used=ok,
                            **{f"Q_{k}": v[3] for k, v in ch.items()}))
        if ok:
            valid.append((ah, ch))
    if len(valid) < 2:
        print(f"[경고] 셀 {cell:02d}: 쓸 수 있는 check-up {len(valid)}개 -> 제외"); continue
    (ah0, c0), (ah1, c1) = valid[0], valid[-1]
    r_first = c0["1C"][3] / c0["C10"][3]; r_last = c1["1C"][3] / c1["C10"][3]
    R0a, s0 = r0_at(P, cell, ah0); R0b, s1 = r0_at(P, cell, ah1); dR0 = R0b - R0a
    neg = bool(dR0 < 0); dR0c = max(dR0, 0.0) if np.isfinite(dR0) else np.nan
    if np.isfinite(dR0c):
        q1 = q_at(c0["1C"][0], c0["1C"][1], c0["1C"][2] - 5.0 * dR0c)
        q10 = q_at(c0["C10"][0], c0["C10"][1], c0["C10"][2] - 0.5 * dR0c)
        r_first_corr = q1 / q10
    else:
        r_first_corr = np.nan
    m = M[M.cell == cell]
    if m.empty:
        lam = Q = np.nan; dah = np.nan
    else:
        k = (m.ah - ah1).abs().idxmin(); lam, Q, dah = m.loc[k, "LAM_NE"], 1 - m.loc[k, "soh"], m.loc[k, "ah"] - ah1
    cell_rows.append(dict(cell=cell, temp=TEMP[cell], ah_first=ah0, ah_last=ah1, r_first=r_first, r_last=r_last,
                          dr=r_last - r_first, dR0=dR0, dR0_neg=neg, R0_match=f"{s0}/{s1}",
                          r_first_corr=r_first_corr, dr_corr=r_last - r_first_corr,
                          LAM_NE=lam, Q=Q, dva_ah_gap=dah))

pd.DataFrame(ck_rows).to_csv(OUT / "R1_rate_checkups.csv", index=False)
D = pd.DataFrame(cell_rows); D.to_csv(OUT / "R1_rate_cells.csv", index=False)
D = D.dropna(subset=["LAM_NE", "Q", "dr"])
print(f"[데이터] 셀 {len(D)}개 사용, 30 °C 초과로 뺀 check-up {sum((not r['used']) and r['n_rates'] == 4 for r in ck_rows)}개, "
      f"dR0<0 (보정 0) {int(D.dR0_neg.sum())}셀, R0 시점 맞춤 {D.R0_match.value_counts().to_dict()}, "
      f"DVA 누적 Ah 차이 최대 {D.dva_ah_gap.abs().max():.1f} Ah")


def judge(col):
    x = D.dropna(subset=[col]); r, p2 = pcorr(x.LAM_NE, x[col], x.Q); p1 = p2 / 2 if r < 0 else 1 - p2 / 2
    return r, p1, len(x), (r <= RHO) and (p1 < PMAX)


r1, p1, n1, ok1 = judge("dr"); r2, p2, n2, ok2 = judge("dr_corr")
print(f"\n[H-R1] LAM_NE vs dr (Q 통제):          rho {r1:+.2f}, 단측 p {p1:.3f}, n {n1} -> {'충족' if ok1 else '미충족'}")
print(f"[H-R2] LAM_NE vs dr_보정 (Q 통제):     rho {r2:+.2f}, 단측 p {p2:.3f}, n {n2} -> {'충족' if ok2 else '미충족'}")
grade = "탐색 이상 지지 (확증 미만)" if ok1 and ok2 else ("저항 효과로 설명 가능" if ok1 else ("기각" if not ok2 else "H-R2만 충족 (사전 등록에 없는 경우 -> 논의)"))
print(f"[판정] {grade}")
print(f"[보고 전용] 통제 없는 순위상관: dr {spearmanr(D.LAM_NE, D.dr)[0]:+.2f}, dr_보정 {spearmanr(D.LAM_NE, D.dr_corr)[0]:+.2f}")
print("[보고 전용] 온도별 dr 중앙값: " + ", ".join(f"{t} {g.dr.median():+.4f}" for t, g in D.groupby("temp")))

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
cols = {"R": "tab:green", "C": "tab:blue", "H": "tab:red"}
for a, col, ttl in [(ax[0], "dr", f"H-R1: rho {r1:+.2f}"), (ax[1], "dr_corr", f"H-R2 (옴 보정): rho {r2:+.2f}")]:
    for t, g in D.groupby("temp"):
        a.scatter(g.LAM_NE, g[col], c=cols.get(t, "k"), label={"R": "25 °C", "C": "−5 °C", "H": "45 °C"}.get(t, t))
    a.set_xlabel("LAM_NE (트랙 C DVA)"); a.set_ylabel("Δr = r(마지막) − r(첫)"); a.set_title(ttl); a.axhline(0, c="gray", lw=0.5)
ax[0].legend(title="사이클 온도"); fig.tight_layout(); fig.savefig(OUT / "R1_rate.png", dpi=130)
print(f"\n저장: {OUT / 'R1_rate_cells.csv'}, {OUT / 'R1_rate.png'}")
