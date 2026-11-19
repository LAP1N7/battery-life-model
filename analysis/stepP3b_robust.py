"""P3b: P3 ① 판정의 강건성 — 반원 끝 검출 방법에 좌우되는지 두 방법 + 신뢰구간으로 재판정

[합의된 사양 2026-11-19 확정]
  방법 A : P3와 같은 모델 없는 비율  r_A = R_arc / R0        (P3_ratio.csv 사용, 반원 끝 = -Im 첫 국소 최소)
  방법 B : ECM 기반 비율             r_B = R2 / (Rs + R1)    (ecm_matlab.csv, separable 스펙트럼만)
           고주파 저항은 Rs+R1로 묶음 (11-18 합의: Rs와 R1은 원리적으로 분리 불가) -> 반원 끝 검출을 쓰지 않음
  공통   : 저 SOC(10~30%) 평균 (check-up마다 2개 이상 SOC 필요), SOH >= 0.80, check-up 3개 이상 셀만
  신뢰구간: 셀별 r의 시간 기울기에 잔차 부트스트랩 95% 구간 (2000회, 난수 고정)
           구간 상한 < 0 -> '감소', 하한 > 0 -> '증가', 0 포함 -> '변화 없음'(감소·증가로 세지 않음)
  판정(방법마다): 고온(45 °C) 7셀 중 '감소' 5개 이상 AND 저온(-5 °C) 7셀 중 '증가' 5개 이상
                 AND 방향 일치(r 증가 <-> 위상 피크 더 음수) 80% 이상 [변화 없음 셀 제외하고 계산]
  종합: 강건 성립 = A와 B 모두 판정 충족 / 부분 성립 = 저온 기준만 A·B 모두 충족 / 그 외 기각
  보조(판정 미사용): 중 SOC(40~60), 고 SOC(70~90)에서 방법 A 기울기 부호
  plan B: 고온 부분이 무너지면 '저온 = 반원 성장형'만 결론, 고온은 KIT(B6)에서 재시험

[예상 문제 — 미리 적어 둠]
  - 신뢰구간 적용으로 기준이 엄격해져 셀 수가 줄 수 있음
  - 방법 B는 separable(약 60%)만 써서 고온 셀 일부가 빠질 수 있음 -> 빠진 셀 출력
  - check-up이 적은 셀(예: 09)은 구간이 넓어 '변화 없음'이 되기 쉬움

입력: out/P3_ratio.csv (stepP3_arc.py 결과), out/c20.csv, out/checkup_time.csv, ../matlab/out/ecm_matlab.csv
출력: out/P3b_cells.csv, out/P3b_compare.png
변경 이력: 2026-11-19 최초 작성
           2026-11-19 전달 전 시험 실행에서 찾은 수정: (1) 위상 기울기를 전체 SOC가 아닌 저 SOC로 계산 (버그)
                      (2) 방법 B의 check-up 짝 맞추기 (Ah 반올림 차이로 방향 일치가 계산되지 않던 버그)
                      (3) 보조 B2 = (R1+R2)/Rs 와 ECM 저항별 진단 추가 — 사후 추가라 종합 판정에 쓰지 않음
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from common import OUT, TEMP

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
LOW, MID, HIGH = [10, 20, 30], [40, 50, 60], [70, 80, 90]
HOT = [c for c in range(1, 22) if TEMP[c] == "H"]; COLD = [c for c in range(1, 22) if TEMP[c] == "C"]
NB = 2000; rng = np.random.default_rng(0)
ECM = OUT.parent.parent / "matlab" / "out" / "ecm_matlab.csv"

# ---------- check-up 표 (P3와 동일) ----------
C = pd.read_csv(OUT / "c20.csv").merge(pd.read_csv(OUT / "checkup_time.csv")[["cell", "ah", "hours"]],
                                         on=["cell", "ah"], how="left")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
def ck(cell, ah):
    if ah < 1: return 0.0, 1.0
    k = C[(C.cell == cell) & ((C.ah - ah).abs() < 1)]
    return (k.hours.iloc[0], k.soh.iloc[0]) if len(k) else (np.nan, np.nan)

# ---------- 두 방법의 check-up별 r ----------
A = pd.read_csv(OUT / "P3_ratio.csv")
A = A[(A.found == 1) & (A.r > 0)][["cell", "ah", "hours", "soh", "soc", "r", "pk_depth"]]
E = pd.read_csv(ECM)
E = E[E.separable == 1].copy()
E["r"] = E.R2 / (E.Rs + E.R1)
E[["hours", "soh"]] = [ck(c, a) for c, a in zip(E.cell, E.ah)]
# 보조 B2 (사후 추가, 판정 미사용): 보이는 반원 전체(R1+R2) / Rs  — 사전 정의 B는 R1을 분모(고주파)에 넣어 A와 정의가 어긋남
E2 = E.copy(); E2["r"] = (E2.R1 + E2.R2) / E2.Rs
A["ahk"] = A.ah.round(1)
pk = A[A.soc.isin(LOW)].groupby(["cell", "ahk"]).pk_depth.mean()      # 위상 피크 (저 SOC, 측정값, 두 방법 공통)

def per_checkup(D, socs):
    g = D[D.soc.isin(socs) & (D.soh >= 0.80)].dropna(subset=["hours"])
    k = g.groupby(["cell", "ah"]).agg(hours=("hours", "first"), r=("r", "mean"), n=("soc", "count")).reset_index()
    return k[k.n >= 2]

def boot_slope(t, y):
    """잔차 부트스트랩: 기울기 [1/100 h]와 95% 구간"""
    b, a = np.polyfit(t, y, 1); res = y - (a + b * t)
    bs = [np.polyfit(t, a + b * t + rng.choice(res, len(res)), 1)[0] for _ in range(NB)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return b * 100, lo * 100, hi * 100

rows = []
for name, D in [("A", A), ("B", E), ("B2", E2)]:
    K = per_checkup(D, LOW)
    for cell in range(1, 22):
        k = K[K.cell == cell].sort_values("hours")
        if len(k) < 3: rows.append(dict(method=name, cell=cell, temp=TEMP[cell], n_ck=len(k), cls="데이터 부족")); continue
        s, lo, hi = boot_slope(k.hours.values, k.r.values)
        cls = "감소" if hi < 0 else ("증가" if lo > 0 else "변화 없음")
        p = pk.reindex(list(zip(k.cell, k.ah.round(1)))).values
        m = np.isfinite(p)
        pks = np.polyfit(k.hours.values[m], p[m], 1)[0] * 100 if m.sum() >= 3 else np.nan
        rows.append(dict(method=name, cell=cell, temp=TEMP[cell], n_ck=len(k), slope=s, lo=lo, hi=hi, cls=cls,
                         pk_slope=pks, agree=(np.sign(s) == -np.sign(pks)) if cls != "변화 없음" and np.isfinite(pks) else np.nan))
R = pd.DataFrame(rows)

# ---------- 판정 ----------
res = {}
for name in ["A", "B", "B2"]:
    r = R[R.method == name]
    hot_dn = int(((r.cell.isin(HOT)) & (r.cls == "감소")).sum()); cold_up = int(((r.cell.isin(COLD)) & (r.cls == "증가")).sum())
    ag = r.agree.dropna().astype(bool); agree = ag.mean() if len(ag) else np.nan
    full = hot_dn >= 5 and cold_up >= 5 and agree >= 0.8; cold_ok = cold_up >= 5
    res[name] = (full, cold_ok)
    lack = r[r.cls == "데이터 부족"].cell.tolist()
    print(f"[방법 {name}] 고온 감소 {hot_dn}/7 | 저온 증가 {cold_up}/7 | 방향 일치 {agree:.0%} ({len(ag)}셀) "
          f"-> {'충족' if full else '미충족'}   (데이터 부족 셀: {lack})")
print("  (B2는 사후 추가 보조 정의라 종합 판정에 쓰지 않음)")
verdict = "강건 성립" if res["A"][0] and res["B"][0] else ("부분 성립 (저온만)" if res["A"][1] and res["B"][1] else "기각")
print(f"\n[P3b 종합 판정] {verdict}")
if verdict != "강건 성립": print("  -> 합의된 plan B: '저온 = 반원 성장형'만 결론, 고온은 KIT(B6)에서 재시험" if verdict.startswith("부분") else "  -> ① 결론을 보류하고 KIT에서 재시험")

piv = R.pivot_table(index=["cell", "temp"], columns="method", values="cls", aggfunc="first")
print("\n셀별 분류 (A: 반원 끝 기반, B: ECM R2/(Rs+R1), B2: ECM (R1+R2)/Rs 보조)"); print(piv.to_string())

# ---------- 진단: ECM 저항별 변화 (저 SOC, BOL -> 마지막 check-up, 판정 미사용) ----------
print("\n[진단] ECM 저항별 변화율과 tau1 중앙값 (tau1 약 1.6e-5 s = 10 kHz 측정 한계)")
Ed = pd.read_csv(ECM); Ed = Ed[Ed.soc.isin(LOW)].copy()
Ed[["hours", "soh"]] = [ck(c, a) for c, a in zip(Ed.cell, Ed.ah)]
for cell, g in Ed[Ed.soh >= 0.80].groupby("cell"):
    g = g.groupby("ah").agg(Rs=("Rs", "mean"), R1=("R1", "mean"), R2=("R2", "mean"), tau1=("tau1", "median")).reset_index()
    if len(g) < 2: continue
    b, e = g.iloc[0], g.iloc[-1]
    print(f"  {cell:02d} {TEMP[cell]}  Rs {e.Rs/b.Rs-1:+.2f}  R1 {e.R1/b.R1-1:+.2f}  R2 {e.R2/b.R2-1:+.2f}  tau1 {g.tau1.median():.1e} s")

# ---------- 보조: 다른 SOC 대역 (방법 A, 판정 미사용) ----------
print("\n[보조, 판정 미사용] 방법 A 기울기 부호 — 셀 수 (감소 / 증가)")
for bname, socs in [("중 SOC 40~60", MID), ("고 SOC 70~90", HIGH)]:
    K = per_checkup(A, socs); out = {"H": [0, 0], "C": [0, 0]}
    for cell, k in K.groupby("cell"):
        if len(k) < 3 or TEMP[cell] == "R": continue
        s = np.polyfit(k.hours, k.r, 1)[0]; out[TEMP[cell]][0 if s < 0 else 1] += 1
    print(f"  {bname}: 고온 {out['H'][0]} / {out['H'][1]}   저온 {out['C'][0]} / {out['C'][1]}")

R.to_csv(OUT / "P3b_cells.csv", index=False, float_format="%.6g")

# ---------- 그림: 셀별 기울기와 95% 구간 ----------
col = {"H": "tab:red", "R": "tab:gray", "C": "tab:blue"}
fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
order = sorted(range(1, 22), key=lambda c: ("CRH".index(TEMP[c]), c))
for ax, name in zip(axs, ["A", "B", "B2"]):
    r = R[R.method == name].set_index("cell")
    for i, c in enumerate(order):
        if c not in r.index or r.loc[c, "cls"] == "데이터 부족": ax.text(0, i, "데이터 부족", fontsize=7, va="center", color="0.5"); continue
        x = r.loc[c]; ax.errorbar(x.slope, i, xerr=[[x.slope - x.lo], [x.hi - x.slope]], fmt="o", color=col[TEMP[c]], capsize=3)
    ax.axvline(0, color="k", lw=.8); ax.set_yticks(range(len(order))); ax.set_yticklabels([f"{c:02d} {TEMP[c]}" for c in order])
    ax.set(title=f"방법 {name}: " + {"A": "R_arc / R0 (반원 끝 기반)", "B": "R2 / (Rs+R1) (ECM, 사전 정의)", "B2": "(R1+R2) / Rs (ECM, 사후 보조)"}[name],
           xlabel="r 기울기 [1/100 h] (점 = 추정값, 막대 = 95% 구간)")
fig.suptitle(f"P3b: 저 SOC 비율 r 의 변화 — 종합 판정: {verdict}"); fig.tight_layout()
fig.savefig(OUT / "P3b_compare.png", dpi=150); plt.close(fig)
print("\n저장:", OUT / "P3b_cells.csv", ",", OUT / "P3b_compare.png")
