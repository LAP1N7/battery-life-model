"""P6: DRT(시간상수 분포) — n2 상승이 '겹쳐 있던 반원들의 분리' 때문인가, 피팅 artifact인가 (모델 가정 없이 확인)

[합의된 사양 2026-11-19]  (사용자 가설: 초기에는 여러 반원이 겹쳐 n이 낮았고, 열화되며 분리되거나 하나가 커져 n이 올랐다)
 방법: Z = R_inf + jwL + sum_k g_k/(1+jw*tau_k),  tau 격자 1e-7~1e3 s (121점), 2차 미분 Tikhonov 정규화 + 비음수 최소제곱(NNLS)
       |Z| 가중, 피팅은 0.2 Hz 이상 (확산 꼬리 영향 줄임)
 봉우리: 반원 영역 tau = 1/(2*pi*f), f 0.8 Hz~10 kHz 안의 내부 극대, prominence >= 창 안 최대값의 10%
 대상: 저 SOC(10·20·30%), 셀별 첫 check-up과 마지막 check-up (SOH >= 0.80, 3개 SOC가 모두 있는 check-up)
 셀 판정(분리 증가): 3개 SOC 중 2개 이상에서 [봉우리 수 증가] 또는 [가장 큰 두 봉우리 간격이 0.2 decade 이상 증가]
 강건성: lambda 3개(3e-3, 1e-2, 3e-2). 판정은 1e-2, 나머지는 같은 판정이 나오는지 보고
 성립 기준: 분리 증가 셀이 대상 셀의 70% 이상 (20셀이면 14셀) -> 탐색적 지지
 반증 조건: 봉우리 수·간격이 그대로인데 n2만 오름 -> n2 상승은 피팅 artifact
 보조(판정 미사용): g 질량의 log(tau) 사분위 범위(IQR, 분포 폭) 변화, ECM n2 변화와의 대조
 사후 추가(그림을 본 뒤, 판정 미사용): 저주파 주 봉우리의 폭(FWHM)·비중 — 가설의 다른 갈래 '하나가 비대해짐' 탐색

[예상 문제 — 합성 데이터로 미리 확인함]
 - 잡음 0.2%에서 lambda가 작으면 가짜 봉우리 생김 -> 판정 lambda 1e-2, prominence 10%로 고정
 - 넓은 반원 하나(n=0.7)와 가까운 두 반원(간격 0.6 decade)은 DRT로도 구분이 흔들림 (합성 시험에서 seed에 따라 봉우리 1~3개)
   -> 봉우리 '수'보다 '간격'을 함께 봄. 그래도 구분이 안 되면 '판정 불가'로 기록
 - 10 kHz보다 빠른 과정(고온 tau1 경계)은 보이지 않음
[그래도 남는 것] 기각이면 "n2 상승 = 피팅 artifact 가능성 큼 -> n을 수명 모델에 쓰지 않음" 결정이 남음

입력: ../data/XX/Resistance.csv (common.load_eis), out/P3_ratio.csv (check-up별 SOH), ../matlab/out/ecm_matlab.csv
출력: out/P6_drt_spectra.csv, out/P6_drt_cells.csv, out/P6_drt.png
변경 이력: 2026-11-19 최초 작성 (Claude). 합성 데이터 시험 후 파라미터 고정, 실데이터는 Claude 환경에서 같은 입력으로 실행 확인
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from common import OUT, TEMP, load_eis
from drt_core import drt, peaks, TAU

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
LOW = [10, 20, 30]; LAMS = [3e-3, 1e-2, 3e-2]; LAM_J = 1e-2; SEP_MIN = 0.2

R3 = pd.read_csv(OUT / "P3_ratio.csv")
E = pd.read_csv(OUT.parent.parent / "matlab" / "out" / "ecm_matlab.csv")
for d in (R3, E): d["ahk"] = d.ah.round(1)

def sep2(lt, h):
    if len(lt) < 2: return 0.0
    i = np.argsort(h)[-2:]; return float(abs(lt[i[0]] - lt[i[1]]))

def iqr_logtau(g):
    w = (TAU >= 1 / (2 * np.pi * 1e4)) & (TAU <= 1 / (2 * np.pi * 0.8)); gw = g[w]; lt = np.log10(TAU[w])
    c = np.cumsum(gw) / gw.sum(); return float(np.interp(0.75, c, lt) - np.interp(0.25, c, lt))

rows = []
for cell in range(1, 22):
    try: spec = load_eis(cell)
    except Exception as ex: print(f"cell {cell:02d}: 읽기 실패 - {ex}"); continue
    soh = R3[R3.cell == cell].groupby("ahk").soh.first()
    aks = sorted({round(a, 1) for (a, s) in spec if s in LOW and all((a, x) in spec for x in LOW)})
    aks = [a for a in aks if a in soh.index and soh[a] >= 0.80]
    if len(aks) < 2: print(f"cell {cell:02d}: 대상 check-up 부족 ({len(aks)})"); continue
    for tag, ak in [("first", aks[0]), ("last", aks[-1])]:
        a_true = [a for (a, s) in spec if round(a, 1) == ak][0]
        for soc in LOW:
            d = spec[(a_true, soc)]
            for lam in LAMS:
                g, Ri, L, res = drt(d[:, 0], d[:, 3], d[:, 4], lam); lt, h = peaks(g)
                rows.append(dict(cell=cell, temp=TEMP[cell], tag=tag, ahk=ak, soh=soh[ak], soc=soc, lam=lam, res=res,
                                 n_pk=len(lt), sep=sep2(lt, h), iqr=iqr_logtau(g),
                                 pk_logtau=";".join(f"{x:.2f}" for x in lt)))
S = pd.DataFrame(rows)
S = S.merge(E[["cell", "ahk", "soc", "n2"]], on=["cell", "ahk", "soc"], how="left")
S.to_csv(OUT / "P6_drt_spectra.csv", index=False, float_format="%.6g")
print(f"[데이터] 스펙트럼 {S[S.lam == LAM_J].shape[0]}개 (셀 {S.cell.nunique()}개 x 첫/마지막 x SOC 3), 상대 잔차 중앙값 {S.res.median():.4f}")

def cell_verdict(lam):
    out = []
    for cell, g in S[S.lam == lam].groupby("cell"):
        f = g[g.tag == "first"].set_index("soc"); l = g[g.tag == "last"].set_index("soc")
        up = [(l.n_pk[s] > f.n_pk[s]) or (l.sep[s] - f.sep[s] >= SEP_MIN) for s in LOW]
        out.append(dict(cell=cell, temp=TEMP[cell], split=int(sum(up) >= 2), n_soc_up=int(sum(up)),
                        npk_f=f.n_pk.mean(), npk_l=l.n_pk.mean(), sep_f=f.sep.mean(), sep_l=l.sep.mean(),
                        iqr_f=f.iqr.mean(), iqr_l=l.iqr.mean(), n2_f=f.n2.mean(), n2_l=l.n2.mean()))
    return pd.DataFrame(out)

V = cell_verdict(LAM_J)
need = int(np.ceil(0.7 * len(V))); k = int(V.split.sum())
print(f"\n[P6 판정, lambda={LAM_J:g}] 분리 증가 셀 {k} / {len(V)}  (기준 {need} 이상) -> {'탐색적 지지' if k >= need else '기각'}")
for t, nm in [("H", "고온"), ("R", "상온"), ("C", "저온")]:
    s = V[V.temp == t]; print(f"  {nm}: {int(s.split.sum())} / {len(s)}")
for lam in LAMS:
    if lam == LAM_J: continue
    Vx = cell_verdict(lam); same = (Vx.set_index("cell").split == V.set_index("cell").split).mean()
    print(f"  [강건성] lambda={lam:g}: 분리 증가 {int(Vx.split.sum())} / {len(Vx)}, 셀 판정 일치 {same:.0%}")
print(f"\n[보조] 평균 봉우리 수 {V.npk_f.mean():.2f} -> {V.npk_l.mean():.2f},  두 봉우리 간격 {V.sep_f.mean():.2f} -> {V.sep_l.mean():.2f} decade,"
      f"  분포 폭(IQR) {V.iqr_f.mean():.2f} -> {V.iqr_l.mean():.2f} decade")
print(f"       ECM n2 {V.n2_f.mean():.3f} -> {V.n2_l.mean():.3f}  (n2 오른 셀 {int((V.n2_l > V.n2_f).sum())} / {len(V)})")
art = V[(V.split == 0) & (V.n2_l > V.n2_f)]
print(f"       반증 조건 해당(분리 증가 없음 + n2 상승) 셀: {len(art)}개 {art.cell.tolist()}")
print(V[['cell','temp','split','n_soc_up','npk_f','npk_l','sep_f','sep_l','n2_f','n2_l']].round(3).to_string(index=False))
# ---------- 사후 추가 (11-19, 그림을 본 뒤 추가 — 판정 미사용, '탐색'으로만 보고) ----------
# 사용자 가설의 다른 갈래 "하나가 비대해져서 n이 올랐다": 저주파 주 봉우리(tau 1e-3~0.2 s, 두 번째 반원)의 폭(FWHM)과 비중
LT = np.log10(TAU); WM = (LT >= -3) & (LT <= -0.7); WA = (LT >= np.log10(1 / (2 * np.pi * 1e4))) & (LT <= -0.7)
def main_peak(g):
    gw = g[WM]; x = LT[WM]; i = int(np.argmax(gw)); h = gw[i] / 2; l = r = i
    while l > 0 and gw[l] > h: l -= 1
    while r < len(gw) - 1 and gw[r] > h: r += 1
    return x[r] - x[l], gw[i] / g[WA].sum(), x[i]
mp = []
for cell in V.cell:
    spec = load_eis(cell); rec = dict(cell=cell, temp=TEMP[cell])
    for tag in ["first", "last"]:
        ak = S[(S.cell == cell) & (S.tag == tag)].ahk.iloc[0]; a_true = [a for (a, s) in spec if round(a, 1) == ak][0]
        for lam in LAMS:
            v = np.mean([main_peak(drt(spec[(a_true, s)][:, 0], spec[(a_true, s)][:, 3], spec[(a_true, s)][:, 4], lam)[0]) for s in LOW], axis=0)
            rec[f"fw_{tag}_{lam:g}"], rec[f"share_{tag}_{lam:g}"], rec[f"pos_{tag}_{lam:g}"] = v
    mp.append(rec)
MP = pd.DataFrame(mp); V = V.merge(MP, on=["cell", "temp"])
print(f"\n[사후 탐색, 판정 미사용] 저주파 주 봉우리(두 번째 반원)가 좁아진(FWHM 감소) 셀:")
for lam in LAMS:
    k2 = int((V[f"fw_last_{lam:g}"] < V[f"fw_first_{lam:g}"]).sum())
    print(f"  lambda={lam:g}: {k2} / {len(V)}   (비중 증가 {int((V[f'share_last_{lam:g}'] > V[f'share_first_{lam:g}']).sum())} / {len(V)})")
print(f"  lambda={LAM_J:g} 평균 FWHM {V[f'fw_first_{LAM_J:g}'].mean():.2f} -> {V[f'fw_last_{LAM_J:g}'].mean():.2f} decade, "
      f"위치 log10 tau {V[f'pos_first_{LAM_J:g}'].mean():.2f} -> {V[f'pos_last_{LAM_J:g}'].mean():.2f}")
V.to_csv(OUT / "P6_drt_cells.csv", index=False, float_format="%.6g")

# ---------- 그림: 대표 셀 첫/마지막 DRT (SOC 20) ----------
pick = [c for c in [2, 6, 7, 1] if c in V.cell.values]
fig, axs = plt.subplots(1, len(pick), figsize=(4.2 * len(pick), 3.8), squeeze=False)
NAME = {"H": "고온", "R": "상온", "C": "저온"}
for ax, cell in zip(axs[0], pick):
    spec = load_eis(cell)
    for tag, c in [("first", "gold"), ("last", "purple")]:
        ak = S[(S.cell == cell) & (S.tag == tag)].ahk.iloc[0]
        a_true = [a for (a, s) in spec if round(a, 1) == ak][0]; d = spec[(a_true, 20)]
        g, *_ = drt(d[:, 0], d[:, 3], d[:, 4], LAM_J)
        ax.plot(np.log10(TAU), g * 1e3, color=c, label=f"{'첫' if tag == 'first' else '마지막'} ({ak:.0f} Ah)")
    ax.axvspan(np.log10(1 / (2 * np.pi * 1e4)), np.log10(1 / (2 * np.pi * 0.8)), color="0.92", zorder=0)
    ax.set(xlim=(-6, 1), xlabel="log10 tau [s]", ylabel="g [mΩ/격자]", title=f"셀 {cell:02d} {NAME[TEMP[cell]]} SOC 20%")
    ax.legend(fontsize=8)
fig.suptitle("DRT 첫 vs 마지막 check-up (회색 = 봉우리 판정 구간 0.8 Hz~10 kHz)"); fig.tight_layout()
fig.savefig(OUT / "P6_drt.png", dpi=150); plt.close(fig)
print("\n저장:", OUT / "P6_drt_spectra.csv", OUT / "P6_drt_cells.csv", OUT / "P6_drt.png")
