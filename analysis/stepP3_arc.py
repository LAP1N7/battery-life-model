"""P3: 반원 검증 — ① 위상 피크 변화가 '반원 저항 / 고주파 저항' 비율 때문인가 (모델 없이)
                  ② √t 성분(A)·선형 성분(B)이 반원별 성장(R1, R2)과 대응하는가

[합의된 사양 2026-11-18 확정]
 ① 비율 (모델 없이, 원본 스펙트럼에서 직접)
    R0    = 허수부가 0이 되는 점(유도성 -> 용량성 전환)의 실수부         [Ohm]
    R_arc = 반원 끝(확산 꼬리 직전 = R0 이후 -Im의 국소 최소 중 마지막)의 실수부 - R0   [Ohm]
            (11-19 수정: 원래 '첫 국소 최소'였으나 반원이 2개일 때 두 반원 사이 골짜기(약 200 Hz)를 잡는 오검출 확인
             -> 꼬리 직전 골짜기(0.1~5 Hz)를 잡도록 '마지막 국소 최소'로 수정. 판정 기준은 그대로)
    r     = R_arc / R0                                                    [-]
    위상 예측: tan|phi|max = r / (2*sqrt(1+r))  (Rs + R||C 반원 하나 가정)
    실측 피크 깊이 = 1 Hz~1 kHz 위상 최솟값 [deg]
    판정 (저 SOC 10~30%, 셀별 r·피크의 시간 기울기):
      고온(45 °C) 7셀 중 5개 이상 r 감소  AND  저온(-5 °C) 7셀 중 5개 이상 r 증가
      AND  r 변화와 위상 변화 방향 일치(r↑ -> 피크 더 음수) 셀 비율 80% 이상
 ② A·B <-> R1·R2 (MATLAB ECM, separable 스펙트럼만, SOC 10~30%)
    성장 속도 = (R/R_ref - 1)의 시간 기울기 x 100  [%/100 h],  R_ref = 첫 separable check-up
    A, B = P1_pred.csv 의 A_full, B_full (A<0 셀 07·12·21 제외)
    판정: rho(A,R1) >= 0.5 AND rho(B,R2) >= 0.5 AND rho(A,R1) > rho(A,R2) AND rho(B,R2) > rho(B,R1)
 공통: 기울기는 SOH >= 0.80 구간, check-up 3개 이상인 셀만. 제외 수는 모두 출력.

[가정과 한계 — 미리 적어 둠]
 - R0에는 10 kHz 밖의 빠른 과정과 인덕턴스 영향이 섞임 -> 결론은 "고주파 쪽 저항"까지
 - 위상 예측식은 반원 하나 가정 (실제는 반원 2개 + 확산) -> 예측·실측은 경향 비교용
 - 고온 셀 SOH 기준점은 첫 C/20(약 47 Ah) -> A·B에 artifact 가능
 - ② 는 약 14셀: rho 약 0.53 이 유의성 경계 -> 경계 근처는 '판정 보류'
 - 사양과 다른 점: 위상 피크 깊이는 P2_window.csv 대신 같은 스펙트럼에서 직접 계산 (정의 동일, 짝 맞추기 단순화)

입력: ../data/XX/Resistance.csv (common.load_eis), out/c20.csv, out/checkup_time.csv,
      out/P1_pred.csv, ../matlab/out/ecm_matlab.csv
출력: out/P3_ratio.csv (스펙트럼별), out/P3_cell.csv (셀별 기울기), out/P3_AB.csv, out/P3_ratio.png, out/P3_AB.png
변경 이력: 2026-11-19 최초 작성
           2026-11-19 (1) 범례 마이너스 기호 수정(사용자)  (2) CSV 저장을 유효숫자 6자리로 (반올림이 B 값을 뭉갬)
                      (3) 반원 끝 오검출 의심 스펙트럼 집계 추가 — 보고만 하고 판정 기준·대상은 바꾸지 않음
           2026-11-19 (4) 반원 끝 = '마지막 국소 최소'로 수정 (사용자 합의). 근거: 저 SOC 600개 중 국소 최소가 2개인
                      16개는 첫 최소가 130~425 Hz(두 반원 사이), 마지막이 0.16~0.35 Hz. 1개인 584개는 모두 0.1~5 Hz.
                      수정 전 결과는 P3_ratio_v1.csv 등으로 보관하지 않음 (변경 이력에 수치 기록: 고온 5/7, 저온 7/7, 89%)
                      판정 문구에 3단계 등급(확증/탐색적 지지/기각) 안내 추가 — 이 스크립트는 반원끝 방식 하나라 '확증'은 못 냄
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from common import OUT, TEMP, load_eis

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
CELLS = range(1, 22)
LOW = [10, 20, 30]
HOT = [c for c in CELLS if TEMP[c] == "H"]; COLD = [c for c in CELLS if TEMP[c] == "C"]
ECM = OUT.parent.parent / "matlab" / "out" / "ecm_matlab.csv"

# ---------- check-up 표: 누적 Ah -> 시간 [h], SOH ----------
C = pd.read_csv(OUT / "c20.csv").merge(pd.read_csv(OUT / "checkup_time.csv")[["cell", "ah", "hours"]],
                                         on=["cell", "ah"], how="left")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])

def ck(cell, ah):
    """EIS check-up의 (시간, SOH). BOL(0 Ah)은 (0, 1.0)."""
    if ah < 1: return 0.0, 1.0
    k = C[(C.cell == cell) & ((C.ah - ah).abs() < 1)]
    return (k.hours.iloc[0], k.soh.iloc[0]) if len(k) else (np.nan, np.nan)

# ---------- ① 스펙트럼별 R0, R_arc, r, 위상 ----------
def arc_metrics(a):
    """a = [f, |Z|, phase(deg), Re, -Im], 주파수 내림차순"""
    f, ph, re, mim = a[:, 0], a[:, 2], a[:, 3], a[:, 4]
    band = (f >= 1) & (f <= 1000)
    pk = ph[band].min() if band.any() else np.nan
    out = dict(pk_depth=pk, R0=np.nan, R_arc=np.nan, found=0)
    z = np.where((mim[:-1] < 0) & (mim[1:] >= 0))[0]           # 유도성 -> 용량성
    if len(z) == 0: return out
    i = z[0]; out["R0"] = np.interp(0, [mim[i], mim[i + 1]], [re[i], re[i + 1]])
    m = mim[i + 1:]; r_ = re[i + 1:]
    mins = [j for j in range(1, len(m) - 1)                    # 국소 최소 후보
            if m[j] <= m[j - 1] and m[j] < m[j + 1] and m[:j].max() > m[j]]
    if mins:                                                    # 마지막 국소 최소 = 확산 꼬리 직전 = 반원 끝
        j = mins[-1]; out.update(R_arc=r_[j] - out["R0"], found=1, n_min=len(mins))
    return out

rows = []
for cell in CELLS:
    try: spec = load_eis(cell)
    except Exception as ex: print(f"cell {cell:02d}: 읽기 실패 - {ex}"); continue
    for (ah, soc), a in spec.items():
        h, soh = ck(cell, ah)
        rows.append(dict(cell=cell, temp=TEMP[cell], ah=ah, hours=h, soh=soh, soc=soc, **arc_metrics(a)))
R = pd.DataFrame(rows)
R["r"] = R.R_arc / R.R0
R["phi_pred"] = -np.degrees(np.arctan(R.r / (2 * np.sqrt(1 + R.r))))
# 반원 끝 오검출 의심 (보고 전용, 판정에 사용하지 않음): 같은 셀·SOC에서 직전 check-up 대비 r이 40% 넘게 급락
R = R.sort_values(["cell", "soc", "ah"])
R["r_prev"] = R.groupby(["cell", "soc"]).r.shift()
R["suspect"] = ((R.found == 1) & (R.r < 0.6 * R.r_prev)).astype(int)
R.to_csv(OUT / "P3_ratio.csv", index=False, float_format="%.6g")

low = R[R.soc.isin(LOW)]
print(f"[① 데이터] 스펙트럼 {len(R)}개 | 저 SOC {len(low)}개 중 반원 끝 검출 {int(low.found.sum())}개 "
      f"(미검출 {int((low.found == 0).sum())}개 제외)")
ok = R[(R.found == 1) & (R.r > 0)]
if len(ok) > 2:
    rho, p = spearmanr(ok.phi_pred, ok.pk_depth)
    print(f"[① 점검] 예측 위상 깊이 vs 실측 피크 깊이 (전체 SOC, {len(ok)}개): 순위상관 {rho:.2f} (p={p:.1e})")

def slope(t, y):
    m = np.isfinite(t) & np.isfinite(y)
    return np.polyfit(t[m], y[m], 1)[0] if m.sum() >= 3 else np.nan

cells1 = []
for cell, g in low[(low.found == 1) & (low.soh >= 0.80)].groupby("cell"):
    k = g.groupby("ah").agg(hours=("hours", "first"), r=("r", "mean"), pk=("pk_depth", "mean"), n=("soc", "count"))
    k = k[k.n >= 2].dropna(subset=["hours"])                   # 저 SOC 3개 중 2개 이상 검출된 check-up
    cells1.append(dict(cell=cell, temp=TEMP[cell], n_ck=len(k),
                       r_slope=slope(k.hours.values, k.r.values) * 100,      # [1/100 h]
                       pk_slope=slope(k.hours.values, k.pk.values) * 100))   # [deg/100 h]
P = pd.DataFrame(cells1).dropna(subset=["r_slope", "pk_slope"])
sus = low[(low.suspect == 1) & (low.soh >= 0.80)]
print(f"[오검출 의심, 보고 전용] 판정 구간(SOH>=0.80, 저 SOC)에서 r 급락(직전 대비 -40% 초과) 스펙트럼 {len(sus)}개, "
      f"해당 셀 {sorted(sus.cell.unique().tolist())}  -> 판정은 사전 기준 그대로 유지")
P["agree"] = np.sign(P.r_slope) == -np.sign(P.pk_slope)         # r↑ -> 피크 더 음수
hot_dn = int((P[P.cell.isin(HOT)].r_slope < 0).sum()); n_hot = int(P.cell.isin(HOT).sum())
cold_up = int((P[P.cell.isin(COLD)].r_slope > 0).sum()); n_cold = int(P.cell.isin(COLD).sum())
agree = P.agree.mean() if len(P) else np.nan
pass1 = hot_dn >= 5 and cold_up >= 5 and agree >= 0.8
print(f"\n[① 판정] 분석 셀 {len(P)}개 (check-up 3개 미만 등으로 제외 {21 - len(P)}개)")
print(f"  고온 r 감소: {hot_dn} / {n_hot}  (기준 5/7 이상)")
print(f"  저온 r 증가: {cold_up} / {n_cold}  (기준 5/7 이상)")
print(f"  r 변화와 위상 변화 방향 일치: {agree:.0%}  (기준 80% 이상)")
print(f"  -> 반원끝 방식 기준 {'충족' if pass1 else '미충족'}  (11-19 합의: 이 결과 하나로는 '탐색적 지지'까지. '확증'은 KIT에서 두 정의 일치 시)")
print(f"  참고: 국소 최소가 2개 이상이라 정의 수정의 영향을 받은 저 SOC 스펙트럼 {int((low.get('n_min', 0) >= 2).sum())}개")
print(P.round(4).to_string(index=False))

# ---------- ② A·B <-> R1·R2 ----------
E = pd.read_csv(ECM)
E = E[(E.separable == 1) & E.soc.isin(LOW)]
g2 = []
for cell, g in E.groupby("cell"):
    k = g.groupby("ah").agg(R1=("R1", "mean"), R2=("R2", "mean")).reset_index()
    k[["hours", "soh"]] = [ck(cell, a) for a in k.ah]
    k = k[(k.soh >= 0.80)].dropna().sort_values("hours")
    if len(k) < 3: continue
    g2.append(dict(cell=cell, n_ck=len(k),
                   R1_rate=slope(k.hours.values, k.R1.values / k.R1.iloc[0] - 1) * 100 * 100,   # [%/100 h]
                   R2_rate=slope(k.hours.values, k.R2.values / k.R2.iloc[0] - 1) * 100 * 100))
G = pd.DataFrame(g2)
Pp = pd.read_csv(OUT / "P1_pred.csv")[["cell", "temp", "A_full", "B_full"]]
AB = Pp.merge(G, on="cell", how="inner")
AB = AB[AB.A_full > 0]
print(f"\n[② 데이터] P1 셀 {len(Pp)}개 -> separable check-up 3개 이상 + A>0: {len(AB)}개")
res = {}
for a_, r_ in [("A_full", "R1_rate"), ("B_full", "R2_rate"), ("A_full", "R2_rate"), ("B_full", "R1_rate")]:
    rho, p = spearmanr(AB[a_], AB[r_]) if len(AB) > 3 else (np.nan, np.nan)
    res[(a_, r_)] = rho
    print(f"  rho({a_[0]}, {r_[:2]}) = {rho:.2f}  (p={p:.3f})")
rA1, rB2, rA2, rB1 = res[("A_full", "R1_rate")], res[("B_full", "R2_rate")], res[("A_full", "R2_rate")], res[("B_full", "R1_rate")]
pass2 = rA1 >= 0.5 and rB2 >= 0.5 and rA1 > rA2 and rB2 > rB1
print(f"[② 판정] {'성립' if pass2 else '기각'}  (기준: 주 비교 둘 다 >= 0.5 이고 교차 비교보다 강함)")
AB.to_csv(OUT / "P3_AB.csv", index=False, float_format="%.6g")
P.to_csv(OUT / "P3_cell.csv", index=False, float_format="%.6g")

# ---------- 그림 ----------
col = {"H": "tab:red", "R": "tab:gray", "C": "tab:blue"}; name = {"H": "고온 45 °C", "R": "상온 25 °C", "C": "저온 -5 °C"}
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for cell, g in low[(low.found == 1)].groupby("cell"):
    k = g.groupby("ah").agg(soh=("soh", "first"), r=("r", "mean")).dropna().sort_values("soh", ascending=False)
    t = TEMP[cell]; ax[0].plot(k.soh, k.r, "-o", ms=3, color=col[t], alpha=.7)
for t in "HRC": ax[0].plot([], [], color=col[t], label=name[t])
ax[0].invert_xaxis(); ax[0].set(xlabel="SOH", ylabel="r = R_arc / R0 (저 SOC 평균)", title="① 반원/고주파 저항 비율의 변화"); ax[0].legend()
for t in "HRC":
    s = ok[ok.temp == t]; ax[1].scatter(s.phi_pred, s.pk_depth, s=6, color=col[t], alpha=.5, label=name[t])
ax[1].set(xlabel="예측 위상 깊이 [deg] (r로 계산)", ylabel="실측 피크 깊이 [deg] (1 Hz–1 kHz)", title="① 예측 vs 실측 위상"); ax[1].legend()
fig.tight_layout(); fig.savefig(OUT / "P3_ratio.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(2, 2, figsize=(10, 8))
for (i, j), (a_, r_) in zip([(0, 0), (0, 1), (1, 0), (1, 1)], [("A_full", "R1_rate"), ("A_full", "R2_rate"), ("B_full", "R1_rate"), ("B_full", "R2_rate")]):
    for t in "HRC":
        s = AB[AB.temp == t]; ax[i, j].scatter(s[r_], s[a_], color=col[t], label=name[t])
    for _, rr in AB.iterrows(): ax[i, j].annotate(f"{int(rr.cell):02d}", (rr[r_], rr[a_]), fontsize=7)
    ax[i, j].set(xlabel=f"{r_[:2]} 성장 [%/100 h]", ylabel=f"{a_[0]} ({'√t' if a_[0] == 'A' else '선형'})",
                 title=f"rho = {res[(a_, r_)]:.2f}" + ("  [주 비교]" if (a_, r_) in [("A_full", "R1_rate"), ("B_full", "R2_rate")] else ""))
ax[0, 0].legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "P3_AB.png", dpi=150); plt.close(fig)
print("\n저장:", OUT / "P3_ratio.csv", OUT / "P3_cell.csv", OUT / "P3_AB.csv", OUT / "P3_ratio.png", OUT / "P3_AB.png")
