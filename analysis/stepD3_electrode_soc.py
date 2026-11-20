# -*- coding: utf-8 -*-
"""
stepD3_electrode_soc.py - 같은 셀 SOC에서의 전극 리튬화(x_NE, y_PE) 정렬표 (과업 C4).

변경 이력
  2026-11-19  최초 작성. D1_dva.csv(variant=='main')에서 셀 SOC 10..90 %별 x_NE, y_PE 계산, 자기 검증, 요약 작성.
              common.py / D1 입력 / 판정 기준은 수정하지 않음 (읽기 전용).

실행:  python -I -X utf8 stepD3_electrode_soc.py
정의 (stepD_core.py 머리말 L12-L14, L182-L188 기준):
  x_NE = 음극(흑연) 리튬 분율, y_PE = 양극(NMC) 리튬 분율. 완전 방전(OCV=vmin)에서 x 작음 / y 큼, 완전 충전에서 x 큼 / y 작음.
  Q는 완전 방전 상태부터 센 용량:  x = x0 + Q/Cn,  y = y0 - Q/Cp   (Li 보존).
  셀 SOC = Q / Q_cell,  Q_cell = (xT - x0) * Cn  (D1 모델의 완전 방전 -> 완전 충전 구간).
  => x_NE(soc) = x0 + soc*(xT - x0),  y_PE(soc) = y0 - soc*(xT - x0)*Cn/Cp.
"""
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

CONFIG = dict(
    in_csv=HERE / "out" / "D" / "D1_dva.csv",
    out_csv=HERE / "out" / "D" / "D3_electrode_soc.csv",
    out_md=HERE / "out" / "D" / "D3_electrode_soc.md",
    variant="main",
    soc_list=[10, 20, 30, 40, 50, 60, 70, 80, 90],   # [%]
    soc_primary=50,                                   # 혼동 크기 및 1차 검증 SOC [%]
    pass_ratio=0.90,                                  # 검증 ② 통과 기준 (바꾸지 말 것)
    lo=0.0, hi=1.0,                                   # 검증 ① 구간
)


def build_table(m, cfg):
    rows = []
    for r in m.itertuples(index=False):
        span = r.xT - r.x0
        for s in cfg["soc_list"]:
            f = s / 100.0
            rows.append(dict(cell=r.cell, ah=r.ah, checkup_idx=r.checkup_idx, T_group=r.T_group, soc_cell=s,
                             x_NE=r.x0 + f * span, y_PE=r.y0 - f * span * r.Cn / r.Cp,
                             LLI=r.LLI, pass_rmse=r.pass_rmse, pass_cap=r.pass_cap, is_ref=r.is_ref))
    return pd.DataFrame(rows)


def check_direction(t, cfg, subset_mask=None, soc_only=None):
    """기준 check-up 대비 이후(ΔLLI>0) check-up에서 Δx_NE<0 비율. 단위 = (check-up, SOC) 쌍. subset_mask: 사용할 행."""
    n_tot = n_ok = 0
    fails = []
    for cell, g in t.groupby("cell"):
        ref = g[g["is_ref"]]
        if ref.empty:
            continue
        ref_ah = ref["ah"].iloc[0]
        ref = ref.set_index("soc_cell")
        rlli = ref["LLI"].iloc[0]
        gg = g[~g["is_ref"] & (g["ah"] > ref_ah)]
        if subset_mask is not None:
            gg = gg[subset_mask.loc[gg.index]]
        if soc_only is not None:
            gg = gg[gg["soc_cell"].isin(soc_only)]
        for r in gg.itertuples():
            dl = r.LLI - rlli
            if not dl > 0:
                continue
            dx = r.x_NE - ref.loc[r.soc_cell, "x_NE"]
            n_tot += 1
            if dx < 0:
                n_ok += 1
            else:
                fails.append(dict(cell=cell, ah=r.ah, checkup_idx=r.checkup_idx, soc_cell=r.soc_cell, dLLI_pp=100 * dl, dx_NE=dx))
    return n_tot, n_ok, pd.DataFrame(fails)


def md_table(df):
    cols = [str(c) for c in df.columns]
    L = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        L.append("| " + " | ".join(f"{v:.4g}" if isinstance(v, float) else str(v) for v in r.values) + " |")
    return "\n".join(L)


def main():
    cfg = CONFIG
    d = pd.read_csv(cfg["in_csv"], dtype={"cell": str})
    m = d[d["variant"] == cfg["variant"]].copy()
    m["T_group"] = m["T_group"].astype(str)
    t = build_table(m, cfg)
    # 같은 (cell, ah) 행 정렬 보조: itertuples 행 순서 유지
    t.to_csv(cfg["out_csv"], index=False, float_format="%.6g")

    # ---- 정합 점검: soc=100 %에서 y가 yT와 일치하는가
    span = m["xT"] - m["x0"]
    y_full = m["y0"] - span * m["Cn"] / m["Cp"]
    max_y_mis = float((y_full - m["yT"]).abs().max())
    qcell = span * m["Cn"]
    cap_gap = (qcell - m["cap_model"]) / m["cap_model"] * 100   # SOC 분모(fit 완전 충전 기준) vs OCV-한계 기준 cap_model

    # ---- ① BOL / 기준 check-up 범위
    ref_t = t[t["is_ref"]]
    c0_t = t[t["checkup_idx"] == 0]
    def rng_ok(x):
        return int(((x["x_NE"] >= cfg["lo"]) & (x["x_NE"] <= cfg["hi"]) & (x["y_PE"] >= cfg["lo"]) & (x["y_PE"] <= cfg["hi"])).sum()), len(x)
    ok_ref, n_ref = rng_ok(ref_t)
    ok_c0, n_c0 = rng_ok(c0_t)
    ok_all, n_all = rng_ok(t)
    bad_ref = ref_t[~((ref_t["x_NE"].between(0, 1)) & (ref_t["y_PE"].between(0, 1)))]
    bad_all = t[~((t["x_NE"].between(0, 1)) & (t["y_PE"].between(0, 1)))]

    # ---- ② 방향성
    S = cfg["soc_primary"]
    res = {}
    for name, mask in (("전체 check-up", None), ("pass_rmse & pass_cap 통과분", (t["pass_rmse"] & t["pass_cap"]))):
        n50, k50, f50 = check_direction(t, cfg, mask, [S])
        na, ka, fa = check_direction(t, cfg, mask, None)
        res[name] = (n50, k50, f50, na, ka, fa)
    # 기준 이전(ΔLLI<0 → Δx>0 기대) 참고: H 셀
    pre_n = pre_ok = 0
    for cell, g in t.groupby("cell"):
        ref = g[g["is_ref"]]
        if ref.empty:
            continue
        rah = ref["ah"].iloc[0]
        rlli = ref["LLI"].iloc[0]
        for r in g[(g["ah"] < rah) & (g["soc_cell"] == S)].itertuples():
            if r.LLI - rlli < 0:
                pre_n += 1
                pre_ok += int(r.x_NE - ref[ref["soc_cell"] == S]["x_NE"].iloc[0] > 0)

    # ---- ④ 혼동 크기: SOC 50 %에서 기준 대비 x_NE 이동량 |Δx| (기준 이후 check-up)
    mv = []
    for cell, g in t[t["soc_cell"] == S].groupby("cell"):
        ref = g[g["is_ref"]]
        if ref.empty:
            continue
        g2 = g[g["ah"] > ref["ah"].iloc[0]].copy()
        if g2.empty:
            continue
        g2["dx"] = g2["x_NE"] - ref["x_NE"].iloc[0]
        g2["dy"] = g2["y_PE"] - ref["y_PE"].iloc[0]
        g2["dLLI_pp"] = 100 * (g2["LLI"] - ref["LLI"].iloc[0])
        i = g2["dx"].abs().idxmax()
        mv.append(dict(cell=cell, T_group=g2["T_group"].iloc[0], max_abs_dx=abs(g2.loc[i, "dx"]), dx_at_max=g2.loc[i, "dx"],
                       dy_at_max=g2.loc[i, "dy"], ah_at_max=g2.loc[i, "ah"], dLLI_pp_at_max=g2.loc[i, "dLLI_pp"],
                       max_abs_dy=g2["dy"].abs().max(), n_checkups=len(g2)))
    mv = pd.DataFrame(mv)
    med = mv.groupby("T_group")[["max_abs_dx", "max_abs_dy", "dLLI_pp_at_max"]].median()
    med["n_cells"] = mv.groupby("T_group").size()
    # x_NE의 가능한 전체 SOC 범위 대비 (xT - x0 중앙값) 참고
    span_ref = float((m[m["is_ref"]]["xT"] - m[m["is_ref"]]["x0"]).median())

    # ---- 통과 현황 (한계)
    nr = m[~m["is_ref"]]
    n_nr = len(nr)
    n_ev = int(nr["cap_evaluable"].sum())
    n_both = int((nr["pass_rmse"] & nr["pass_cap"]).sum())

    # ---- 보고서
    L = []
    P = L.append
    P("# D3 전극 리튬화 정렬표 (stepD3_electrode_soc.py 자동 생성)\n")
    P("## 변경 이력\n- 2026-11-19 최초 작성 (입력 D1_dva.csv은 읽기 전용, common.py 및 기준 미변경)\n")
    P("## 정의 (코드 근거)\n"
      "- x_NE = 음극(흑연) Li 분율, y_PE = 양극(NMC) Li 분율. 모델 V(Q)=Up(y0-Q/Cp)-Un(x0+Q/Cn), Q는 빈 상태부터 (stepD_core.py L11-L12).\n"
      "- 완전 충전 기준 파라미터화: V(s)=Up(yT+s/Cp)-Un(xT-s/Cn), s=완전 충전 아래 용량, n_Li=yT*Cp+xT*Cn=y0*Cp+x0*Cn (L13-L14, L182-L188).\n"
      "- 방향: 충전 시 x 증가, y 감소. Up은 y에 대해 감소(Up_nmc 4.3452-1.6518y+..., stepD_ocp_mohtat.py L30-L33), Un도 x에 대해 감소(L19-L27, exp(-75x) 항) -> 만충 x 크고 y 작음. D1 값: x0 약 0.006, xT 약 0.79, y0 약 0.85, yT 약 0.03.\n"
      "- (y0, x0) = 모델 OCV가 vmin(3.0 V)이 되는 빈 상태 (empty_state, stepD_core.py L301-L316; 호출 stepD_pipeline.py L112). (yT, xT) = 곡선 s=0인 충전 끝 (curve top, L104-L115, L175).\n"
      "- LLI = 1 - nLi/nLi_ref, LAM_PE = 1-Cp/Cp_ref, LAM_NE = 1-Cn/Cn_ref (modes_vs_ref, L319-L320). 기준(is_ref) = 셀의 첫 averaged 정상 check-up (choose_reference 'first_avg', stepD_pipeline.py L60-L69, 호출 L273 in stepD1_dva.py).\n"
      "- 용량 축: Q_cell = (xT - x0)*Cn. 셀 SOC = Q/Q_cell (0 %=완전 방전, 100 %=완전 충전). x_NE=x0+soc*(xT-x0), y_PE=y0-soc*(xT-x0)*Cn/Cp.\n"
      f"- 일관성 점검: soc=100 %에서 y가 yT와 어긋나는 최대값 {max_y_mis:.2e} (Li 보존 확인). Q_cell과 OCV 한계(3.0/4.2 V) 기준 cap_model 차이(%): 중앙값 {cap_gap.median():.2f}, 최대 절대 {cap_gap.abs().max():.2f}.\n")
    P(f"입력: variant=='{cfg['variant']}' {len(m)}행, 셀 {m['cell'].nunique()}개 (T_group: "
      + ", ".join(f"{k} {v}" for k, v in m.groupby('T_group')['cell'].nunique().items()) + f"). 출력 {len(t)}행 -> {cfg['out_csv'].name}\n")
    P("## 검증 ① 기준 check-up의 x_NE, y_PE가 [0,1] 안인가\n")
    P(f"- is_ref 행: {ok_ref}/{n_ref} 통과 (셀 {m[m['is_ref']]['cell'].nunique()}개 x SOC 9점)")
    P(f"- checkup_idx==0 행: {ok_c0}/{n_c0} 통과")
    P(f"- 참고, 모든 check-up: {ok_all}/{n_all} 통과")
    if len(bad_ref):
        P("- is_ref 실패 행:\n" + md_table(bad_ref[["cell", "ah", "soc_cell", "x_NE", "y_PE"]]))
    if len(bad_all):
        P("- 전체 중 범위 밖 행(상위 20):\n" + md_table(bad_all[["cell", "ah", "checkup_idx", "soc_cell", "x_NE", "y_PE"]].head(20)))
    P("")
    P("## 검증 ② LLI가 클수록 x_NE가 낮아지는가 (기준 check-up 대비 이후 check-up, ΔLLI>0만)\n")
    P(f"기준 {cfg['pass_ratio']*100:.0f} %. 단위 = (check-up, SOC) 쌍; ΔLLI = LLI - LLI_기준, Δx_NE = x_NE - x_NE_기준 (같은 SOC).\n")
    ver = {}
    for name, (n50, k50, f50, na, ka, fa) in res.items():
        r50 = k50 / n50 if n50 else np.nan
        ra = ka / na if na else np.nan
        ver[name] = (r50, ra)
        P(f"- {name}: SOC {S} % {k50}/{n50} = {100*r50:.1f} % ({'통과' if r50 >= cfg['pass_ratio'] else '실패'}); "
          f"SOC 10-90 % 합계 {ka}/{na} = {100*ra:.1f} % ({'통과' if ra >= cfg['pass_ratio'] else '실패'})")
    P(f"- 참고: 기준 이전 check-up(H 셀, ΔLLI<0)에서 SOC {S} %의 Δx_NE>0 일치 {pre_ok}/{pre_n}.\n")
    n50, k50, f50, na, ka, fa = res["전체 check-up"]
    if len(fa):
        P("### 실패 목록 (전체 check-up, Δx_NE>=0인 쌍)\n")
        P(f"실패 {len(fa)}쌍, check-up {fa.groupby(['cell','checkup_idx']).ngroups}개, 셀 {fa['cell'].nunique()}개. SOC 50 % 실패 {len(f50)}쌍.\n")
        fa2 = fa.merge(m[["cell", "checkup_idx", "pass_rmse", "pass_cap", "T_group", "LAM_NE", "LAM_PE"]], on=["cell", "checkup_idx"], how="left")
        fa2["LAM_NE"] = 100 * fa2["LAM_NE"]
        fa2["LAM_PE"] = 100 * fa2["LAM_PE"]
        fa2["LAMNE_minus_LLI_pp"] = fa2["LAM_NE"] - fa2["dLLI_pp"]
        g = fa2.groupby(["cell", "T_group", "ah", "checkup_idx"]).agg(
            n_soc_fail=("soc_cell", "size"), soc_failed=("soc_cell", lambda v: ",".join(str(int(i)) for i in v)),
            dLLI_pp=("dLLI_pp", "first"), LAM_NE=("LAM_NE", "first"), LAM_PE=("LAM_PE", "first"),
            LAMNE_minus_LLI_pp=("LAMNE_minus_LLI_pp", "first"), pass_rmse=("pass_rmse", "first"), pass_cap=("pass_cap", "first")).reset_index()
        P("check-up 단위 요약 (실패한 SOC 목록 포함; 전체 값은 D3_electrode_soc.csv):\n")
        P(md_table(g.round(2)) + "\n")
        P(f"실패 check-up 중 LAM_NE > ΔLLI: {int((g['LAMNE_minus_LLI_pp'] > 0).sum())}/{len(g)}.\n")
    else:
        P("실패 없음.\n")
    n50p, k50p, f50p, nap, kap, fap = res["pass_rmse & pass_cap 통과분"]
    P("### 통과분만으로 재계산 (기준은 그대로)\n")
    P(f"- SOC {S} %: {k50p}/{n50p} = {100*ver['pass_rmse & pass_cap 통과분'][0]:.1f} %; SOC 10-90 % 합계 {kap}/{nap} = {100*ver['pass_rmse & pass_cap 통과분'][1]:.1f} %. 실패 {len(fap)}쌍.\n")
    if len(fap):
        gp = fap.groupby(["cell", "ah", "checkup_idx"]).agg(n_soc_fail=("soc_cell", "size"), soc_failed=("soc_cell", lambda v: ",".join(str(int(i)) for i in v)), dLLI_pp=("dLLI_pp", "first")).reset_index()
        P(md_table(gp.round(2)) + "\n")
    P("## LLI 혼동 크기 (SOC 50 %, 기준 대비 x_NE 이동량)\n")
    P(f"기준 check-up 이후 중 |Δx_NE|가 최대인 check-up (셀별). 참고: 기준 check-up의 셀 SOC 0-100 % 사이 x 폭(xT-x0) 중앙값 {span_ref:.3f}.\n")
    P(md_table(mv.round(4).sort_values("cell")) + "\n")
    P("T_group별 중앙값 (셀별 최대값들의 중앙값; dx는 |Δx_NE|):\n")
    P(md_table(med.reset_index().round(4)) + "\n")
    P(f"전체 셀 중앙값 {mv['max_abs_dx'].median():.4f}, 최대 {mv['max_abs_dx'].max():.4f}.\n")
    P("## 한계\n"
      f"- DVA 용량 정합: 기준 제외 check-up {n_nr}개 중 용량 평가 가능 {n_ev}, RMSE<10 mV와 용량 1 pp 이내를 모두 통과한 것 {n_both}개 (D1 요약의 53/220 계열; 여기서는 variant=main 기준 재계산 값). 통과 못한 check-up의 LLI는 신뢰도가 낮다.\n"
      "- 45 degC(H) 셀은 기준 check-up이 BOL이 아니다(check-up 1 또는 2, flag_ref_not_bol). 따라서 H 셀의 LLI와 x_NE 이동은 '기준 이후' 값이며 R/C 셀의 BOL 기준값과 직접 비교할 수 없다.\n"
      "- 셀 SOC는 DVA 모델의 빈~만충 구간(Q_cell=(xT-x0)*Cn)이며 측정 용량이나 OCV 한계 기준 cap_model과 약간 다르다(위 차이 참조). EIS 측정 SOC와 맞추려면 별도 보정이 필요.\n"
      "- 절대 stoichiometry는 OCP 모델(Mohtat/PyBaMM)의 파라미터화 산물이며 y의 절대 오프셋은 아티팩트다 (stepD_ocp_mohtat.py L43 notes). 이동량(Δ)만 의미를 두는 것이 안전하다.\n"
      "- 'x_NE 이동량'은 LLI 외에 LAM_NE, LAM_PE의 영향도 함께 포함한다(LLI만 분리한 값이 아님).\n")
    cfg["out_md"].write_text("\n".join(L), encoding="utf-8")
    print("wrote", cfg["out_csv"], cfg["out_md"])


if __name__ == "__main__":
    main()
