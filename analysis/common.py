"""공통: 경로, 데이터 로딩, ECM fitting, C/20 dV/dQ.
데이터 위치: 이 폴더(analysis)와 같은 위치의 ../data/01 ~ 21
단위: 저항 Ohm, 주파수 Hz, 용량 Ah, SOC %
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.signal import savgol_filter

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE / "out"; OUT.mkdir(exist_ok=True)
SOCS = [90, 80, 70, 60, 50, 40, 30, 20, 10]  # check-up 안에서 EIS 측정 순서 (100%에서 10%씩 방전)

def load_eis(cell):
    """Resistance.csv -> {(ah, soc): array[f, |Z|, phase, Re, -Im]}. ah = 누적 Ah (check-up 구분용)"""
    d = pd.read_csv(DATA / f"{cell:02d}" / "Resistance.csv")
    d.columns = ["t", "I", "V", "f", "mag", "ph", "re", "mim", "T", "Q", "cyc", "ah"]
    e = d[d.f > 0].copy()
    e["sid"] = ((e.f.diff() > 0) | (e.t.diff() > 600)).cumsum()   # 주파수가 다시 커지면 새 스펙트럼
    spec = {}
    for ah, g in e.groupby("ah"):
        sids = sorted(g.sid.unique())
        if len(sids) != 9:          # SOC 9개가 다 있는 check-up만 사용
            continue
        for soc, s in zip(SOCS, sids):
            x = g[g.sid == s].sort_values("f", ascending=False)
            spec[(ah, soc)] = x[["f", "mag", "ph", "re", "mim"]].values
    return spec

# ECM: Rs + L + R1||CPE1 + R2||CPE2 + 확산 CPE(Qw, nw)
LO = [0, 0, 0, -3, 0.5, 0, -3, 0.5, -2, 0.3]
HI = [0.05, 1e-5, 0.05, 4, 1, 0.05, 4, 1, 5, 1]
def Zm(p, w):
    Rs, L, R1, lQ1, n1, R2, lQ2, n2, lQw, nw = p
    jw = 1j * w
    return Rs + jw*L + R1/(1 + R1*10**lQ1*jw**n1) + R2/(1 + R2*10**lQ2*jw**n2) + 1/(10**lQw*jw**nw)

def fit_ecm(f, re, mim, p0=None):
    """|Z| 가중 최소제곱. 반환: 파라미터, 상대 RMS 오차.
    주의: R1/R2 분리는 식별성이 낮음 -> R_int = R1+R2 만 신뢰."""
    w = 2*np.pi*f; Z = re - 1j*mim
    def r(p):
        d = (Zm(p, w) - Z) / np.abs(Z); return np.r_[d.real, d.imag]
    inits = [[re.min(), 1e-7, .001, 0, .8, .002, 1, .8, 2, .6], [re.min(), 1e-7, .002, 1, .8, .001, 0, .8, 2, .6]]
    if p0 is not None: inits.insert(0, p0)
    best = None
    for i0 in inits:
        i0 = np.clip(i0, np.array(LO)+1e-9, np.array(HI)-1e-9)
        s = least_squares(r, i0, bounds=(LO, HI), max_nfev=3000)
        if best is None or s.cost < best.cost: best = s
    return best.x, np.sqrt(2*best.cost/len(f))

def load_c20(cell):
    """OCV_wExpansion.csv -> check-up별 C/20 방전용량, 충전곡선의 dV/dQ"""
    d = pd.read_csv(DATA / f"{cell:02d}" / "OCV_wExpansion.csv")
    d.columns = ["t", "I", "V", "x", "T", "Q", "cap", "cyc", "ah"]
    rows, curves = [], {}
    for ah, v in d.groupby("ah"):
        ch = v[(v.I > 100) & (v.I < 400)]; dc = v[(v.I < -100) & (v.I > -400)]   # C/20 = 약 250 mA
        if len(ch) < 500 or len(dc) < 500: continue
        q, V = ch.cap.values, ch.V.values
        grid = np.arange(0.02, q.max()-0.02, 0.01)
        Vg = np.interp(grid, q, V)
        curves[ah] = (grid, Vg, savgol_filter(Vg, 31, 3, deriv=1, delta=0.01))
        rows.append(dict(cell=cell, ah=ah, cap=dc.cap.max()))
    return pd.DataFrame(rows), curves

TEMP = {c: "RCH"[(c-1) % 3] for c in range(1, 22)}   # README Table 1: R=상온, C=저온, H=고온
