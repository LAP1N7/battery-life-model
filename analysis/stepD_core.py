# -*- coding: utf-8 -*-
"""
stepD_core.py - dataset-agnostic degradation-mode (DVA half-cell) fitting.

No dataset-specific code lives here. A dataset provides
  * a standard long table (see stepD_loaders.py): cell, checkup_idx, ah_throughput, direction, Q_Ah, V, T[, thickness_um]
  * an OCP pair registered with register_ocp() (chemistry matched to the cell)
Everything else (curve averaging, bounded multistart fit, modes, capacity-from-modes check, identifiability)
is generic.

Model (per check-up):   V(Q) = Up(y0 - Q/Cp) - Un(x0 + Q/Cn)
  Q counted from the empty state. Internally the top-of-charge parametrisation is used
  (s = capacity below the full state):  V(s) = Up(yT + s/Cp) - Un(xT - s/Cn),  (yT, xT) = stoichiometry at the full state.
  n_Li = yT*Cp + xT*Cn = y0*Cp + x0*Cn  (Li conservation, independent of where along the charge it is evaluated).

Only numpy / pandas / scipy are used. Python 3.13/3.14 compatible.
"""
from __future__ import annotations
import warnings
from dataclasses import dataclass, field
import numpy as np
from scipy import optimize

# ----------------------------------------------------------------------------------------------
# OCP registry
# ----------------------------------------------------------------------------------------------
@dataclass
class OCPPair:
    name: str
    Up: callable            # positive electrode OCP vs Li fraction y (vectorised, V vs Li/Li+)
    Un: callable            # negative electrode OCP vs Li fraction x
    x_range: tuple = (-0.05, 1.10)   # physically meaningful graphite stoichiometry range for this OCP
    y_range: tuple = (-0.05, 1.20)   # positive electrode stoichiometry range
    source: str = ""
    chemistry: str = ""
    notes: str = ""


OCP_REGISTRY: dict[str, OCPPair] = {}


def register_ocp(pair: OCPPair) -> OCPPair:
    OCP_REGISTRY[pair.name] = pair
    return pair


def get_ocp(name: str) -> OCPPair:
    if name not in OCP_REGISTRY:
        raise KeyError(f"OCP pair '{name}' not registered. Registered: {sorted(OCP_REGISTRY)}. "
                       "Import the module that registers it (e.g. stepD_ocp_mohtat) before calling get_ocp.")
    return OCP_REGISTRY[name]


# ----------------------------------------------------------------------------------------------
# configuration
# ----------------------------------------------------------------------------------------------
@dataclass
class FitConfig:
    vmin: float = 3.0              # cycler / OCV window limits [V]
    vmax: float = 4.2
    fit_vlo: float = 3.30          # fixed voltage window used for the fit, on the averaged curve [V]
    fit_vhi: float = 4.18
    n_grid: int = 300              # points of the resampled curve inside the fit window
    top_tol: float = 0.05          # a branch must reach within this of vmax at the top to be 'complete'
    bot_tol: float = 0.25          # ...and within this of vmin at the bottom (loose: charge may start after a rest)
    min_points: int = 200          # minimum raw samples of a branch
    n_starts: int = 16             # multistart count (random starts; warm starts are added)
    n_starts_ref: int = 40         # more starts for the first (reference) fit of a cell
    conv_tol_mV: float = 0.5       # a start is 'converged' if its RMSE is within this of the best RMSE
    cap_bounds: tuple = (0.6, 2.2)  # Cp, Cn bounds as multiples of the measured capacity
    voff_max: float = 0.08         # |constant polarisation offset| bound for single-branch fits [V]
    seed: int = 0
    avg_voltage_offset: bool = False   # if True also fit a constant offset in averaged mode (not recommended)


# ----------------------------------------------------------------------------------------------
# curve preparation
# ----------------------------------------------------------------------------------------------
def _clean_branch(Q, V):
    Q = np.asarray(Q, float)
    V = np.asarray(V, float)
    m = np.isfinite(Q) & np.isfinite(V)
    return Q[m], V[m]


def branch_status(Q, V, direction, cfg: FitConfig):
    """Return (complete: bool, reason)."""
    Q, V = _clean_branch(Q, V)
    if len(Q) < cfg.min_points:
        return False, "too_few_points"
    if direction == "charge":
        if V.max() < cfg.vmax - cfg.top_tol:
            return False, "charge_not_reaching_top"
        if V[:20].min() > cfg.vmin + cfg.bot_tol + 0.1 and V.min() > cfg.vmin + 0.5:
            return False, "charge_start_too_high"
    else:
        if V[:20].max() < cfg.vmax - cfg.top_tol - 0.05:
            return False, "discharge_not_from_top"
        if V.min() > cfg.vmin + cfg.bot_tol:
            return False, "discharge_not_reaching_bottom"
    return True, "ok"


def _s_curve(Q, V, direction):
    """Capacity below the top state, ascending, unique."""
    Q, V = _clean_branch(Q, V)
    if direction == "charge":
        s = Q[-1] - Q
    else:
        s = Q - Q[0]
    o = np.argsort(s, kind="mergesort")
    s, V = s[o], V[o]
    s, ui = np.unique(s, return_index=True)
    return s, V[ui]


def build_curve(ck, cfg: FitConfig, mode="auto"):
    """
    ck: DataFrame of one check-up (standard long table rows).
    mode: 'auto' (avg if both branches complete, else chg, else dis), or 'avg' | 'chg' | 'dis'.
    Returns dict(kind, s, V, Qcap, flags...) or dict(kind=None, reason=...).
    """
    br = {}
    for d in ("charge", "discharge"):
        g = ck[ck["direction"] == d]
        if len(g):
            br[d] = (g["Q_Ah"].values, g["V"].values)
    st = {d: branch_status(*br[d], d, cfg) if d in br else (False, "missing") for d in ("charge", "discharge")}
    cap_d = float(br["discharge"][0][-1]) if "discharge" in br else np.nan
    # optional: capacity as reported by the cycler for this discharge (column cap_cycler_Ah); preferred for soh so that it
    # matches the team's c20.csv; the integrated value is kept as cap_dis_int
    g_d = ck[ck["direction"] == "discharge"]
    cap_rep = float(g_d["cap_cycler_Ah"].iloc[0]) if ("cap_cycler_Ah" in g_d.columns and len(g_d)) else np.nan
    cap_d_int = cap_d
    if np.isfinite(cap_rep):
        cap_d = cap_rep
    cap_c = float(br["charge"][0][-1]) if "charge" in br else np.nan
    info = dict(charge_status=st["charge"][1], discharge_status=st["discharge"][1], cap_dis=cap_d, cap_dis_int=cap_d_int, cap_chg=cap_c)
    if mode == "auto":
        if st["charge"][0] and st["discharge"][0]:
            kind = "avg"
        elif st["charge"][0]:
            kind = "chg"
        elif st["discharge"][0]:
            kind = "dis"
        else:
            return dict(kind=None, reason=f"charge:{st['charge'][1]};discharge:{st['discharge'][1]}", **info)
    else:
        kind = mode
        need = {"avg": ("charge", "discharge"), "chg": ("charge",), "dis": ("discharge",)}[kind]
        if not all(st[d][0] for d in need):
            return dict(kind=None, reason=f"{kind}_unavailable", **info)
    if kind == "avg":
        sc, vc = _s_curve(*br["charge"], "charge")
        sd, vd = _s_curve(*br["discharge"], "discharge")
        smax = min(sc[-1], sd[-1])
        sg = np.linspace(0.0, smax, 1500)
        v = 0.5 * (np.interp(sg, sc, vc) + np.interp(sg, sd, vd))
        info["v_gap_mV_median"] = float(np.median(np.interp(sg, sc, vc) - np.interp(sg, sd, vd)) * 1e3)
        s_all, V_all = sg, v
    elif kind == "chg":
        s_all, V_all = _s_curve(*br["charge"], "charge")
    else:
        s_all, V_all = _s_curve(*br["discharge"], "discharge")
    # fixed voltage window on the (averaged or single) curve; the window is defined from the top, V decreases with s
    m = (V_all <= cfg.fit_vhi) & (V_all >= cfg.fit_vlo)
    if m.sum() < 30:
        return dict(kind=None, reason="window_too_short", **info)
    idx = np.where(m)[0]
    # contiguous from the first index inside the window to the last (V is monotone except noise)
    i0, i1 = idx[0], idx[-1]
    sg2 = np.linspace(s_all[i0], s_all[i1], cfg.n_grid)
    vg2 = np.interp(sg2, s_all, V_all)
    Qcap = np.nanmax([cap_d, cap_c]) if np.isfinite([cap_d, cap_c]).any() else float(s_all[-1])
    return dict(kind=kind, s=sg2, V=vg2, Qcap=float(Qcap), s_lo=float(s_all[i0]), s_hi=float(s_all[i1]),
                v_lo=float(V_all[i1]), v_hi=float(V_all[i0]), s_end=float(s_all[-1]), **info)


# ----------------------------------------------------------------------------------------------
# model and fit
# ----------------------------------------------------------------------------------------------
def model_V(p, s, ocp: OCPPair):
    """p = (Cp, Cn, yT, xT[, Voff]);  V(s) with s = capacity below the full state."""
    Cp, Cn, yT, xT = p[:4]
    v = ocp.Up(yT + s / Cp) - ocp.Un(xT - s / Cn)
    if len(p) > 4:
        v = v + p[4]
    return v


def _bounds(curve, ocp, cfg, kind):
    Q = curve["Qcap"]
    lo = [cfg.cap_bounds[0] * Q, cfg.cap_bounds[0] * Q, ocp.y_range[0], 0.3]
    hi = [cfg.cap_bounds[1] * Q, cfg.cap_bounds[1] * Q, 0.6, ocp.x_range[1]]
    if kind == "chg":       # charging: terminal voltage above OCV -> Voff >= 0
        lo.append(0.0); hi.append(cfg.voff_max)
    elif kind == "dis":
        lo.append(-cfg.voff_max); hi.append(0.0)
    elif kind == "avg" and cfg.avg_voltage_offset:
        lo.append(-0.02); hi.append(0.02)
    return np.array(lo), np.array(hi)


def _residual(p, s, V, ocp, s_hi):
    r = (model_V(p, s, ocp) - V) * 1e3
    # soft domain penalty: stoichiometry at the bottom of the data must stay inside the OCP's valid range
    Cp, Cn, yT, xT = p[:4]
    xb = xT - s_hi / Cn
    yb = yT + s_hi / Cp
    pen = max(0.0, ocp.x_range[0] - xb) + max(0.0, yb - ocp.y_range[1]) + max(0.0, ocp.y_range[0] - yT) + max(0.0, xT - ocp.x_range[1])
    return np.append(r, 2e4 * pen)


@dataclass
class FitResult:
    kind: str
    p: np.ndarray                  # (Cp, Cn, yT, xT[, Voff])
    rmse_mV: float
    n_converged: int
    n_starts: int
    solutions: list = field(default_factory=list)   # parameter vectors of the converged starts
    at_bound: list = field(default_factory=list)
    sd: np.ndarray | None = None   # Jacobian-based 1-sigma for each parameter (white-noise assumption)
    cond: float = np.nan
    corr: np.ndarray | None = None
    success: bool = True

    @property
    def Cp(self): return float(self.p[0])
    @property
    def Cn(self): return float(self.p[1])
    @property
    def yT(self): return float(self.p[2])
    @property
    def xT(self): return float(self.p[3])
    @property
    def nLi(self): return float(self.p[2] * self.p[0] + self.p[3] * self.p[1])


def fit_halfcell(curve, ocp: OCPPair, cfg: FitConfig, warm_starts=(), rng=None, n_starts=None):
    """Bounded multistart least squares. Returns FitResult."""
    kind = curve["kind"]
    s, V = curve["s"], curve["V"]
    s_hi = float(s[-1])
    lo, hi = _bounds(curve, ocp, cfg, kind)
    npar = len(lo)
    rng = rng or np.random.default_rng(cfg.seed)
    n_starts = cfg.n_starts if n_starts is None else n_starts
    Q = curve["Qcap"]
    starts = []
    for w in warm_starts:
        w = np.asarray(w, float)
        if len(w) < npar:
            w = np.r_[w, np.zeros(npar - len(w))]
        starts.append(np.clip(w[:npar], lo + 1e-9, hi - 1e-9))
    for _ in range(n_starts):
        st = [rng.uniform(0.8, 1.6) * Q, rng.uniform(0.8, 1.6) * Q, rng.uniform(0.0, 0.3), rng.uniform(0.6, 1.0)]
        if npar > 4:
            st.append(rng.uniform(lo[4], hi[4]))
        starts.append(np.clip(np.array(st), lo + 1e-9, hi - 1e-9))
    xs = np.r_[Q, Q, 0.1, 0.1, 0.02][:npar]
    sols = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for st in starts:
            try:
                r = optimize.least_squares(_residual, st, args=(s, V, ocp, s_hi), bounds=(lo, hi), x_scale=xs,
                                           max_nfev=300, xtol=1e-10, ftol=1e-10)
            except Exception:
                continue
            rm = float(np.sqrt(np.mean(r.fun[:-1] ** 2)))
            sols.append((rm, r))
    if not sols:
        return FitResult(kind, np.full(npar, np.nan), np.nan, 0, len(starts), success=False)
    sols.sort(key=lambda t: t[0])
    best_rm, best = sols[0]
    conv = [r.x for rm, r in sols if rm <= best_rm + cfg.conv_tol_mV]
    ab = [i for i in range(npar) if (best.x[i] - lo[i]) < 1e-3 * (hi[i] - lo[i]) or (hi[i] - best.x[i]) < 1e-3 * (hi[i] - lo[i])]
    # Jacobian-based uncertainty (white noise, sigma from residual)
    sd = np.full(npar, np.nan)
    cond = np.nan
    corr = None
    try:
        J = best.jac[:-1] / xs[None, :]
        sv = np.linalg.svd(J, compute_uv=False)
        cond = float(sv[0] / max(sv[-1], 1e-300))
        sig2 = best_rm ** 2
        cov = sig2 * np.linalg.pinv(J.T @ J)
        sd = np.sqrt(np.clip(np.diag(cov), 0, None)) * xs
        d = np.sqrt(np.clip(np.diag(cov), 1e-300, None))
        corr = cov / np.outer(d, d)
    except Exception:
        pass
    return FitResult(kind, best.x.copy(), best_rm, len(conv), len(starts), solutions=conv, at_bound=ab, sd=sd,
                     cond=cond, corr=corr, success=True)


# ----------------------------------------------------------------------------------------------
# states, modes, capacity
# ----------------------------------------------------------------------------------------------
def empty_state(p, ocp: OCPPair, cfg: FitConfig, s_hi_data):
    """(y0, x0, s_empty) at the empty state (model OCV = vmin); NaN if not found."""
    Cp, Cn, yT, xT = p[:4]
    f = lambda s: ocp.Up(yT + s / Cp) - ocp.Un(xT - s / Cn) - cfg.vmin
    smax = min(Cn * (xT - ocp.x_range[0]), Cp * (ocp.y_range[1] - yT))
    a = max(s_hi_data * 0.5, 0.0)
    try:
        grid = np.linspace(a, smax, 200)
        fv = f(grid)
        j = np.where(np.sign(fv[:-1]) != np.sign(fv[1:]))[0]
        if len(j) == 0:
            return np.nan, np.nan, np.nan
        se = optimize.brentq(f, grid[j[0]], grid[j[0] + 1])
        return yT + se / Cp, xT - se / Cn, se
    except Exception:
        return np.nan, np.nan, np.nan


def modes_vs_ref(Cp, Cn, nLi, ref):
    return dict(LLI=1.0 - nLi / ref["nLi"], LAM_PE=1.0 - Cp / ref["Cp"], LAM_NE=1.0 - Cn / ref["Cn"])


def cell_capacity(Cp, Cn, nLi, ocp: OCPPair, vmin=3.0, vmax=4.2):
    """
    Capacity [Ah] between the OCV limits determined by (Cp, Cn, nLi) alone. Li is conserved: y = (nLi - x*Cn)/Cp.
    V(x) is monotonic increasing in x, so each limit has one root. Returns NaN if a limit lies outside the OCP domain.
    """
    st = limit_states(Cp, Cn, nLi, ocp, vmin, vmax)
    return np.nan if st is None else float((st[1] - st[0]) * Cn)


def limit_states(Cp, Cn, nLi, ocp: OCPPair, vmin=3.0, vmax=4.2):
    """(x_bottom, x_top) graphite stoichiometry at the OCV limits for given (Cp, Cn, nLi); None if outside the OCP domain."""
    xl = max(ocp.x_range[0], (nLi - ocp.y_range[1] * Cp) / Cn)
    xh = min(ocp.x_range[1], (nLi - ocp.y_range[0] * Cp) / Cn)
    if not xh > xl:
        return None
    Vx = lambda x: ocp.Up((nLi - x * Cn) / Cp) - ocp.Un(x)
    try:
        if Vx(xl) > vmin or Vx(xh) < vmax:
            return None
        xb = optimize.brentq(lambda x: Vx(x) - vmin, xl, xh)
        xt = optimize.brentq(lambda x: Vx(x) - vmax, xl, xh)
        return float(xb), float(xt)
    except Exception:
        return None


def capacity_from_modes(ref, LLI, LAM_PE, LAM_NE, ocp: OCPPair, vmin=3.0, vmax=4.2):
    """Cell capacity implied by degradation modes applied to a reference parameter set {'Cp','Cn','nLi'}."""
    return cell_capacity(ref["Cp"] * (1 - LAM_PE), ref["Cn"] * (1 - LAM_NE), ref["nLi"] * (1 - LLI), ocp, vmin, vmax)


# ----------------------------------------------------------------------------------------------
# identifiability report
# ----------------------------------------------------------------------------------------------
def identifiability(fit: FitResult):
    """Dict of relative Jacobian uncertainties (percent) and condition number for one fit."""
    if fit.sd is None or not np.all(np.isfinite(fit.sd)):
        return dict(sd_Cp_pct=np.nan, sd_Cn_pct=np.nan, sd_nLi_pct=np.nan, sd_yT=np.nan, sd_xT=np.nan, cond=np.nan,
                    corr_Cp_Cn=np.nan)
    Cp, Cn, yT, xT = fit.p[:4]
    sd = fit.sd
    # delta method for nLi = yT*Cp + xT*Cn using the covariance implied by corr and sd
    if fit.corr is not None:
        cov = fit.corr[:4, :4] * np.outer(sd[:4], sd[:4])
        grad = np.array([yT, xT, Cp, Cn])      # d nLi / d (Cp, Cn, yT, xT)
        # order of p: Cp, Cn, yT, xT -> d/dCp = yT, d/dCn = xT, d/dyT = Cp, d/dxT = Cn
        var_n = float(grad @ cov @ grad)
        sd_n = np.sqrt(max(var_n, 0.0))
        cpcn = float(fit.corr[0, 1])
    else:
        sd_n, cpcn = np.nan, np.nan
    nLi = yT * Cp + xT * Cn
    return dict(sd_Cp_pct=100 * sd[0] / Cp, sd_Cn_pct=100 * sd[1] / Cn, sd_nLi_pct=100 * sd_n / nLi, sd_yT=sd[2],
                sd_xT=sd[3], cond=fit.cond, corr_Cp_Cn=cpcn)


# ----------------------------------------------------------------------------------------------
# per-cell pipeline
# ----------------------------------------------------------------------------------------------
def fit_cell(tbl, ocp: OCPPair, cfg: FitConfig, mode="auto", seed=None):
    """
    tbl: standard long table of ONE cell. Fits every check-up. Returns list of dict rows (one per check-up, in order)
    with the fit and its auxiliary results; modes are NOT computed here (reference chosen by assign_modes).
    """
    rng = np.random.default_rng(cfg.seed if seed is None else seed)
    rows = []
    prev = {}        # kind -> last params
    last_any = None
    for k, ck in tbl.groupby("checkup_idx", sort=True):
        thr = float(ck["ah_throughput"].max())
        cv = build_curve(ck, cfg, mode)
        row = dict(checkup_idx=int(k), ah=thr, T_mean=float(np.nanmean(ck["T"].values)) if "T" in ck else np.nan,
                   kind=cv.get("kind"), cap_dis=cv.get("cap_dis"), cap_chg=cv.get("cap_chg"),
                   reason=cv.get("reason", ""), fit=None, curve=cv)
        if cv.get("kind") is None:
            rows.append(row)
            continue
        ws = []
        if cv["kind"] in prev:
            ws.append(prev[cv["kind"]])
        if last_any is not None:
            ws.append(last_any)
        first_of_kind = cv["kind"] not in prev
        ns = cfg.n_starts_ref if first_of_kind else cfg.n_starts
        fit = fit_halfcell(cv, ocp, cfg, warm_starts=ws, rng=rng, n_starts=ns)
        row["fit"] = fit
        if fit.success:
            prev[cv["kind"]] = fit.p
            last_any = fit.p[:4]
        rows.append(row)
    return rows


def pass_criteria(rms_mV, cap_loss_implied, cap_loss_meas, rms_thr=10.0, cap_thr=0.01):
    c1 = bool(np.isfinite(rms_mV) and rms_mV < rms_thr)
    c2 = bool(np.isfinite(cap_loss_implied) and np.isfinite(cap_loss_meas) and abs(cap_loss_implied - cap_loss_meas) < cap_thr)
    return c1, c2
