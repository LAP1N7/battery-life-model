"""C6: DVA pass-rate notation, offset-corrected gap, 45 degC reference-point comparison.
Change history:
  2026-11-19  v1  first version (team C6). Reads out/D/D1_dva.csv (variant=='main') and
                  data/NN/OCV_wExpansion.csv (H cells, charge capacity only). Does not modify
                  common.py, the DVA reference, or any D1 output. No Crate_wExpansion / new datasets.
Outputs: out/D/C6_dva_reference.csv, out/D/C6_dva_reference.md
"""
import numpy as np, pandas as pd
from pathlib import Path
from common import DATA, OUT

THR_PP = 1.0            # stepD_core.pass_criteria: cap_thr=0.01 -> |loss_implied - loss_meas| < 1 pp (strict)
K_HOT = 0.9745          # team definition: BOL discharge = BOL C/20 charge x 0.9745
HOT = [3, 6, 9, 12, 15, 18, 21]
D = OUT / "D"

d = pd.read_csv(D / "D1_dva.csv"); d = d[d.variant == "main"].copy()
for c in ["is_ref", "pass_cap", "pass_rmse", "cap_evaluable"]:
    d[c] = d[c].astype(bool)

# --- 1. recount (reference check-ups excluded: their gap is 0 by construction)
nr = d[~d.is_ref]; ev = nr[nr.cap_evaluable]
tot = dict(rows=len(d), ref=int(d.is_ref.sum()), nonref=len(nr), evaluable=len(ev),
           cap=int(ev.pass_cap.sum()), rmse=int(nr.pass_rmse.sum()), both=int((ev.pass_cap & ev.pass_rmse).sum()))

# --- 2. offset correction: subtract the reference check-up's loss_gap_pp from later check-ups
ref_gap = d[d.is_ref].set_index("cell").loss_gap_pp
d["ref_gap_pp"] = d.cell.map(ref_gap)
d["gap_corr_pp"] = d.loss_gap_pp - d.ref_gap_pp
d["pass_cap_corr"] = d.gap_corr_pp.abs() < THR_PP
# supplementary sensitivity (NOT the requested definition): subtract the first non-ref evaluable gap of the cell
first = d[(~d.is_ref) & d.cap_evaluable].sort_values(["cell", "ah"]).groupby("cell").loss_gap_pp.first()
d["first_gap_pp"] = d.cell.map(first)
d["gap_corr_first_pp"] = d.loss_gap_pp - d.first_gap_pp
d["pass_cap_corr_first"] = d.gap_corr_first_pp.abs() < THR_PP

# --- 3. 45 degC: SOH on c20_hotbol basis (recomputed: no c20_hotbol.csv in repo)
rows = []
for c in HOT:
    r = pd.read_csv(DATA / f"{c:02d}" / "OCV_wExpansion.csv"); r.columns = ["t", "I", "V", "x", "T", "Q", "cap", "cyc", "ah"]
    ah0 = r.ah.min(); ch = r[(r.ah == ah0) & (r.I > 100) & (r.I < 400)]
    rows.append(dict(cell=c, bol_chg=ch.cap.max(), cap_hotbol=ch.cap.max() * K_HOT))
hot = pd.DataFrame(rows)
d = d.merge(hot, on="cell", how="left")
d["soh_hotbol"] = d.cap_dis_meas / d.cap_hotbol
d["diff_pp"] = 100 * (d.soh_hotbol - d.soh_vs_ref)
refcap = d[d.is_ref].set_index("cell")[["cap_dis_meas", "checkup_idx"]]
hot["ref_idx"] = hot.cell.map(refcap.checkup_idx); hot["ref_cap"] = hot.cell.map(refcap.cap_dis_meas)
hot["ref_soh_hotbol_pct"] = 100 * hot.ref_cap / hot.cap_hotbol
d.to_csv(D / "C6_dva_reference.csv", index=False, float_format="%.6g")

# --- md
def tbl(df): return df.to_markdown(index=False) if hasattr(df, "to_markdown") else df.to_string(index=False)
try:
    import tabulate  # noqa
except ImportError:
    def tbl(df):
        h = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns) + "\n"
        return h + "\n".join("| " + " | ".join(str(v) for v in r) + " |" for r in df.astype(object).itertuples(index=False))
def grp(col):
    g = ev.assign(p=d.loc[ev.index, col]).groupby("T_group")
    return g.p.agg(["sum", "size"])
g = []
for k, name in [("T", None)]: pass
e = d[(~d.is_ref) & d.cap_evaluable]
ts = []
for tg, x in list(e.groupby("T_group")) + [("전체", e)]:
    n = len(x); a = int(x.pass_cap.sum()); b = int(x.pass_cap_corr.sum()); f = int(x.pass_cap_corr_first.sum())
    ts.append([tg, n, f"{a}/{n} ({100*a/n:.1f}%)", f"{b}/{n} ({100*b/n:.1f}%)", f"{f}/{n} ({100*f/n:.1f}%)",
               f"{x.loss_gap_pp.median():+.2f}"])
t2 = pd.DataFrame(ts, columns=["T_group", "평가 가능 n", "원래 pass_cap", "보정(기준 gap 차감)", "참고: 첫 비기준 gap 차감", "gap 중앙값 pp"])
h = d[d.cell.isin(HOT) & ~d.is_ref & d.soh_vs_ref.notna()]
t3 = h.pivot_table(index="checkup_idx", columns="cell", values="diff_pp").round(2).reset_index()
t3["checkup_idx"] = t3.checkup_idx.astype(int)
t3.columns = ["check-up"] + [f"셀{int(c)}" for c in t3.columns[1:]]
t3 = t3.fillna("-")
sm = h.groupby("cell").diff_pp.agg(["min", "median", "max"]).round(2).reset_index(); sm["cell"] = sm.cell.astype(int)
t4 = hot.assign(bol_chg=hot.bol_chg.round(4), cap_hotbol=hot.cap_hotbol.round(4), ref_cap=hot.ref_cap.round(4),
                ref_soh_hotbol_pct=hot.ref_soh_hotbol_pct.round(2))
t4["cell"] = t4.cell.astype(int); t4["ref_idx"] = t4.ref_idx.astype(int)
t4.columns = ["셀", "BOL 충전용량 Ah", "c20_hotbol BOL 방전 Ah", "DVA 기준 check-up idx", "DVA 기준 용량 Ah", "기준점의 hotbol SOH %"]
md = f"""# C6 DVA 기준점: 통과율 표기, 오프셋 보정, 45 °C 기준점

## 변경 이력
- 2026-11-19 v1: 최초 작성 (`stepD_c6_reference.py`). D1 출력, common.py, DVA 기준은 바꾸지 않았다. Crate_wExpansion과 새 데이터는 열지 않았다.

## 정의
- 대상: `out/D/D1_dva.csv`, variant==main, 전체 {tot['rows']}행 (21셀). is_ref 행 {tot['ref']}개.
- loss_meas = 1 - cap_dis_meas / (기준 check-up의 cap_dis_meas); loss_implied = 모드(LLI, LAM)로 모델이 낸 용량 손실; loss_gap_pp = 100 x (loss_implied - loss_meas).
- pass_cap = |loss_implied - loss_meas| < 1 pp (엄격 부등호, `stepD_core.pass_criteria`, cap_thr=0.01). pass_rmse = rmse_mV < 10.
- 분모: 기준 check-up 제외(기준은 정의상 gap 0이라 자동 통과) 후 cap_evaluable(C/20 방전 용량 있음)인 행. 이 정의로 분모 {tot['evaluable']}.
- c20_hotbol.csv는 저장소에 없다. 그래서 `OCV_wExpansion.csv`에서 BOL(ah 최소) C/20 충전용량(전류 100~400 mA 구간 cap 최대, load_c20과 같은 방식)을 구해 x {K_HOT}로 **재계산**했다.

## 1. 통과 수 재집계
- pass_cap: **{tot['cap']}/{tot['evaluable']}** ({100*tot['cap']/tot['evaluable']:.1f} %). 53/220 표기 맞음.
- pass_rmse: {tot['rmse']}/{tot['nonref']} (분모 = 기준 제외 전체 {tot['nonref']}행, 평가 가능 여부와 무관).
- 둘 다 통과: {tot['both']}/{tot['evaluable']} (분모 = 평가 가능 {tot['evaluable']}).
- 기준 행 포함 시 pass_cap {tot['cap']+tot['ref']}/{tot['evaluable']+tot['ref']}로 나오므로, 보고 때 분모에 기준 제외를 명시할 것.

## 2. 오프셋 보정 gap
셀마다 기준 check-up의 loss_gap_pp를 이후 check-up에서 뺐다. **전 셀에서 기준 행 gap이 정확히 0**(loss_meas = 0이고 loss_implied도 기준 모드 대비 0이 되도록 계산되기 때문)이라, 이 보정은 항등이고 통과 수는 원래와 같다. 즉 "기준 행의 gap을 오프셋으로 본다"는 정의로는 보정 효과가 없다.
참고로 오프셋을 다르게 잡은 민감도(요청 정의 아님): 각 셀의 첫 비기준 평가 가능 check-up의 gap을 뺐다. 첫 점은 정의상 항상 통과가 되므로 상한 쪽으로 편향된다.

{tbl(t2)}

## 3. 45 °C 기준점 (c20_hotbol 대비 DVA)
{tbl(t4)}

SOH 차이 = 100 x (SOH_hotbol - soh_vs_ref), pp. SOH_hotbol = cap_dis_meas / (BOL 충전 x {K_HOT}), soh_vs_ref = cap_dis_meas / DVA 기준 용량. check-up별(행 = check-up idx, 기준 check-up 제외):

{tbl(t3)}

셀별 요약 (pp):

{tbl(sm)}

DVA 기준이 BOL(0)이 아니라 check-up 1 또는 2인 이유는 H 셀 BOL에 C/20 방전이 없기 때문이다(D1 flag_ref_not_bol). 그래서 기준점 자체가 이미 hotbol SOH로 100 %가 아니고(위 표 마지막 열), 이후 모든 check-up에서 차이는 (기준점 hotbol SOH - 100) x soh_vs_ref 꼴로 열화와 함께 줄어든다.

## 4. 의견 (기준은 바꾸지 않음)
- 맞추는 쪽 장점: SOH 표기가 21셀 전부 "BOL = 100 %"로 통일된다. 위 표처럼 45 °C 셀은 기준점이 이미 열화된 상태(hotbol SOH {t4['기준점의 hotbol SOH %'].min():.1f} ~ {t4['기준점의 hotbol SOH %'].max():.1f} %)라 DVA SOH를 그대로 쓰면 열화를 과소평가한다. 최대 격차 {h.diff_pp.abs().max():.1f} pp.
- 단점: 0.9745는 충전/방전 효율 가정이라 BOL 방전을 실측하지 못한 값이고, 셀 간 불확실성이 그대로 SOH 오프셋이 된다. 또 DVA 모드(LLI, LAM)의 기준 상태는 모델 피팅이 가능한 check-up이어야 하므로 BOL로 옮기면 모드 기준과 용량 기준이 갈라진다.
- 제안: 모드 분해의 기준(is_ref)은 그대로 두고, 보고용 SOH 열만 c20_hotbol 기준으로 따로 병기한다. 통과율(pass_cap)은 기준 행을 빼도 보정 효과가 없으므로 오프셋 보정 대신 현재 정의 유지가 낫다.
"""
(D / "C6_dva_reference.md").write_text(md, encoding="utf-8")
print(md)
