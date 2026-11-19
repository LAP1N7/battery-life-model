import numpy as np
from scipy.optimize import nnls
from scipy.signal import find_peaks

TAU = np.logspace(-7, 3, 121)          # [s] 시간상수 격자 (측정 창보다 넓게)

def drt(f, re, mim, lam, fmin=0.2):
    """Z = R_inf + j*w*L + sum g_k/(1+j*w*tau_k). 반환: g [Ohm/격자], R_inf, L, 상대 잔차"""
    m = f >= fmin; f, re, mim = f[m], re[m], mim[m]
    w = 2 * np.pi * f; Z = re - 1j * mim; s = np.abs(Z)
    K = 1 / (1 + 1j * np.outer(w, TAU))
    A = np.c_[K, np.ones(len(w)), 1j * w * 1e-6]                 # L 단위 uH
    A = A / s[:, None]; b = Z / s
    Ar = np.r_[A.real, A.imag]; br = np.r_[b.real, b.imag]
    sc = np.median(re)                                             # 저항 규모로 정규화
    Ar[:, :len(TAU)] *= sc
    n = len(TAU); D = np.zeros((n - 2, n + 2))                     # 2차 미분 정규화 (g에만)
    for i in range(n - 2): D[i, i:i + 3] = [1, -2, 1]
    x, _ = nnls(np.r_[Ar, np.sqrt(lam) * D], np.r_[br, np.zeros(n - 2)], maxiter=5000)
    g = x[:n] * sc; Rinf, L = x[n], x[n + 1] * 1e-6
    res = np.sqrt(np.mean((Ar @ x - br) ** 2))
    return g, Rinf, L, res

def peaks(g, fmax=1e4, fmin=0.8, prom=0.10):
    """반원 영역(0.5 Hz~10 kHz, 확산 꼬리 제외)에 해당하는 tau 안의 봉우리: 위치 log10(tau), 높이"""
    lo, hi = 1 / (2 * np.pi * fmax), 1 / (2 * np.pi * fmin)
    w = (TAU >= lo) & (TAU <= hi); gw = g[w]; lt = np.log10(TAU[w])
    if gw.max() <= 0: return np.array([]), np.array([])
    pk, _ = find_peaks(gw, prominence=prom * gw.max())         # 창 가장자리는 봉우리로 세지 않음 (내부 극대만)
    return lt[pk], gw[pk]
