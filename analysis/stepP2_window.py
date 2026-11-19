"""P2: 열화 구간(SOH 2% 창)별 저항 증가율과 EIS 위상 지표 -> 열화 모드 전환 지점 탐색

[합의된 사양 2026-11-18]
  창     : SOH 2% 폭 (1.00-0.98, 0.98-0.96, ...)
  SOC    : 9개 모두 (10~90%). 그림은 히트맵(창 x SOC) + 3구간(저 10-30 / 중 40-60 / 고 70-90)
  지표   : (1) dR/dQ  = 창 안에서 [R/R_BOL 증가 %] / [용량 손실 %]   R = Rs+R1+R2 (MATLAB ECM, ecm_matlab.csv)
           (2) 위상 a = 반원 피크 깊이: 1 Hz-1 kHz 대역 위상 최솟값 [deg] 과 그 주파수 [Hz]
           (3) 위상 b = 저주파 위상: sweep 최저 주파수(약 0.01 Hz)의 위상 [deg]
           (참고) 전체 sweep 위상 max-min [deg]  (고주파 끝은 인덕턴스, 저주파 끝은 확산이 좌우함)
  판정   : dR/dQ 가 창마다 중앙값 대비 +-30% 이내 -> 'SEI 지배 가설 성립(일정)'
           dR/dQ 가 직전 창들 중앙값 대비 30% 넘게 벗어난 첫 창 -> '꺾임 후보'
           위상의 '같이 변함' 기준은 미합의 -> 표와 그림으로만 보고 (자동 판정 안 함)

[방법과 가정]
  - check-up(용량 C/20 + EIS)마다 SOH, R(SOC), 위상 지표를 구함
  - 창 경계(SOH=1.00, 0.98, ...)에서의 값은 이웃 check-up 사이 직선 보간 (손실 Q 기준)
  - 창 안에 실제 check-up이 없으면 interp=1 표시 (그 창의 기울기는 이웃 창과 같아질 수 있음)
  - EIS BOL(누적 Ah=0)의 SOH는 1.0으로 둠. 고온 셀은 C/20 첫 측정이 약 47 Ah 지점이라
    그 사이 손실이 빠져 있음 (알려진 한계)
  - EIS와 C/20 check-up은 누적 Ah 차이 1 Ah 이내면 같은 check-up으로 봄

입력: out/c20.csv, ../matlab/out/ecm_matlab.csv, ../data/XX/Resistance.csv (common.load_eis)
출력: out/P2_window.csv, out/P2/cellXX.png
변경 이력: 2026-11-18 최초 작성
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from common import OUT, TEMP, load_eis, SOCS

plt.rcParams["font.family"] = "Malgun Gothic"; plt.rcParams["axes.unicode_minus"] = False
CELLS = range(1, 22)
STEP = 0.02; TOL = 0.30; BAND = (1.0, 1000.0)
BANDS = {"저 SOC(10-30)": [10, 20, 30], "중 SOC(40-60)": [40, 50, 60], "고 SOC(70-90)": [70, 80, 90]}
ECM = OUT.parent.parent / "matlab" / "out" / "ecm_matlab.csv"
FIG = OUT / "P2"; FIG.mkdir(exist_ok=True)

C = pd.read_csv(OUT / "c20.csv")
C["soh"] = C.groupby("cell").cap.transform(lambda s: s / s.iloc[0])
E = pd.read_csv(ECM); E["R"] = E.Rs + E.R1 + E.R2

def phase_metrics(a):
    """a = [f, |Z|, phase(deg), Re, -Im] (주파수 내림차순). 위상 음수 = 용량성"""
    f, ph = a[:, 0], a[:, 2]
    m = (f >= BAND[0]) & (f <= BAND[1])
    i = np.argmin(ph[m])
    return dict(pk_depth=ph[m][i], pk_f=f[m][i], pk_edge=int(i in (0, m.sum() - 1)),
                ph_low=ph[np.argmin(f)], f_low=f.min(), rng_all=ph.max() - ph.min())

def checkup_table(cell):
    """check-up x SOC 표: soh, Qloss, R, 위상 지표"""
    spec = load_eis(cell); c = C[C.cell == cell]; e = E[E.cell == cell]; rows = []
    for (ah, soc), a in spec.items():
        if ah < 1: soh = 1.0
        else:
            k = c[(c.ah - ah).abs() < 1]
            if k.empty: continue
            soh = k.soh.iloc[0]
        r = e[((e.ah - ah).abs() < 1) & (e.soc == soc)]
        rows.append(dict(cell=cell, ah=ah, soc=soc, soh=soh, R=r.R.iloc[0] if len(r) else np.nan, **phase_metrics(a)))
    T = pd.DataFrame(rows)
    if T.empty: return T
    T = T.sort_values(["soc", "ah"])
    T["R_rel"] = T.R / T.groupby("soc").R.transform("first") - 1      # BOL 대비 증가율
    T["Q"] = 1 - T.soh
    return T

def windows(T):
    """SOC별로 창 경계에서 보간 -> 창별 지표"""
    out = []
    for soc, g in T.groupby("soc"):
        g = g.sort_values("Q").drop_duplicates("Q")
        if len(g) < 2: continue
        Qmax = g.Q.max()
        for lo_soh in np.arange(1.0, 1 - Qmax - 1e-9, -STEP):
            q0, q1 = 1 - lo_soh, 1 - lo_soh + STEP
            if q1 > Qmax + 1e-9: break
            iv = lambda col, q: np.interp(q, g.Q, g[col])
            n_in = int(((g.Q > q0 + 1e-9) & (g.Q < q1 - 1e-9)).sum())
            out.append(dict(soc=soc, win=f"{lo_soh:.2f}-{lo_soh - STEP:.2f}", soh_hi=lo_soh, n_ck=n_in,
                            interp=int(n_in == 0 and not (np.isclose(g.Q, q0).any() and np.isclose(g.Q, q1).any())),
                            dRdQ=(iv("R_rel", q1) - iv("R_rel", q0)) / STEP,      # (R 증가 %)/(손실 %)
                            pk_depth=iv("pk_depth", (q0 + q1) / 2), pk_f=iv("pk_f", (q0 + q1) / 2),
                            ph_low=iv("ph_low", (q0 + q1) / 2), rng_all=iv("rng_all", (q0 + q1) / 2)))
    return pd.DataFrame(out)

def judge(s):
    """창 순서대로 놓인 dR/dQ 시리즈 -> (일정 여부, 꺾임 후보 창 번호)"""
    s = s.dropna().values
    if len(s) < 3: return np.nan, np.nan
    const = bool(np.all(np.abs(s / np.median(s) - 1) <= TOL)) if np.median(s) > 0 else False
    brk = next((i for i in range(2, len(s)) if np.median(s[:i]) > 0 and abs(s[i] / np.median(s[:i]) - 1) > TOL), np.nan)
    return const, brk

allW, summ = [], []
for cell in CELLS:
    try: T = checkup_table(cell)
    except Exception as ex: print(f"cell {cell:02d}: 읽기 실패 - {ex}"); continue
    if T.empty: print(f"cell {cell:02d}: 맞는 check-up 없음"); continue
    W = windows(T)
    if W.empty: print(f"cell {cell:02d}: 창 없음 (EIS check-up {T.ah.nunique()}개)"); continue
    W.insert(0, "cell", cell); W.insert(1, "temp", TEMP[cell]); allW.append(W)
    nck = T.ah.nunique()
    for bname, socs in BANDS.items():
        b = W[W.soc.isin(socs)].groupby("soh_hi", sort=False).agg(dRdQ=("dRdQ", "mean"), interp=("interp", "max"))
        const, brk = judge(b.dRdQ)
        summ.append(dict(cell=cell, temp=TEMP[cell], band=bname, n_eis_ck=nck, n_win=len(b),
                         n_interp=int(b.interp.sum()), dRdQ_med=np.median(b.dRdQ) if len(b) else np.nan,
                         const=const, brk_win=b.index[int(brk)] if brk == brk else None))
    # ---- 그림: 히트맵 3개 + 3구간 선 ----
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for j, (col, lab) in enumerate([("dRdQ", "dR/dQ [%R / %Q]"), ("pk_depth", "반원 피크 위상 [deg]"), ("ph_low", "저주파 위상 [deg]")]):
        H = W.pivot_table(index="soc", columns="soh_hi", values=col).sort_index(axis=1, ascending=False)
        im = ax[0, j].imshow(H.values, aspect="auto", origin="lower", cmap="viridis")
        ax[0, j].set_xticks(range(H.shape[1])); ax[0, j].set_xticklabels([f"{v:.2f}" for v in H.columns], rotation=90, fontsize=7)
        ax[0, j].set_yticks(range(H.shape[0])); ax[0, j].set_yticklabels(H.index)
        ax[0, j].set(title=lab, xlabel="창 시작 SOH", ylabel="SOC [%]"); fig.colorbar(im, ax=ax[0, j])
        for bname, socs in BANDS.items():
            b = W[W.soc.isin(socs)].groupby("soh_hi").agg(v=(col, "mean"), it=("interp", "max")).sort_index(ascending=False)
            ln, = ax[1, j].plot(b.index, b.v, "-o", ms=3, label=bname)
            ax[1, j].plot(b.index[b.it == 1], b.v[b.it == 1], "o", mfc="white", color=ln.get_color())
        ax[1, j].invert_xaxis(); ax[1, j].set(xlabel="창 시작 SOH", ylabel=lab); ax[1, j].legend(fontsize=7)
    fig.suptitle(f"cell {cell:02d} ({TEMP[cell]}) | EIS check-up {nck}개 | 빈 원 = 창 안에 실측 check-up 없음(보간)")
    fig.tight_layout(); fig.savefig(FIG / f"cell{cell:02d}.png", dpi=120); plt.close(fig)
    print(f"cell {cell:02d}: EIS check-up {nck}개, 창 {W.soh_hi.nunique()}개", flush=True)

pd.concat(allW).round(5).to_csv(OUT / "P2_window.csv", index=False)
S = pd.DataFrame(summ)
print("\n[판정 요약] dR/dQ 기준 (위상 동반 여부는 그림으로 확인, 자동 판정 없음)")
print("  const=True: 모든 창이 중앙값 +-30% 이내 | brk_win: 직전 창들 대비 30% 넘게 벗어난 첫 창 (창 시작 SOH)")
print("  n_interp: 실측 없이 보간된 창 수 (많으면 판정 신뢰도 낮음)")
print(S.round(3).to_string(index=False))
S.round(4).to_csv(OUT / "P2_summary.csv", index=False)
print("\n저장:", OUT / "P2_window.csv", ",", OUT / "P2_summary.csv", ",", FIG)
