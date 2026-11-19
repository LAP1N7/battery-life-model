# -*- coding: utf-8 -*-
"""
stepD_pipeline.py - dataset-agnostic orchestration shared by stepD0 (selftest) and stepD1: fit all cells (optionally in
parallel), choose the reference check-up, compute modes, capacity-from-modes check, pass flags.
Still no dataset-specific code; the OCP pair is looked up by registry name.
"""
import os
import importlib
import importlib.util
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

import stepD_core as dc

OCP_MODULES = {"mohtat_peyman": "stepD_ocp_mohtat"}


def ensure_ocp(name, ocp_module=None):
    """Make sure the OCP pair is registered (import its module). ocp_module: path to a .py file that registers pairs (YOUR code)."""
    if name in dc.OCP_REGISTRY:
        return dc.get_ocp(name)
    if ocp_module:
        spec = importlib.util.spec_from_file_location("stepD_user_ocp", os.path.abspath(ocp_module))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    elif name in OCP_MODULES:
        importlib.import_module(OCP_MODULES[name])
    return dc.get_ocp(name)


def _worker(args):
    cell, tbl, ocp_name, ocp_module, cfg, mode = args
    ocp = ensure_ocp(ocp_name, ocp_module)
    rows = dc.fit_cell(tbl, ocp, cfg, mode=mode)
    return cell, rows


def fit_all(curves, ocp_name, cfg, mode="auto", jobs=1, ocp_module=None, log=print):
    cells = list(dict.fromkeys(curves["cell"].tolist()))
    args = [(c, curves[curves["cell"] == c], ocp_name, ocp_module, cfg, mode) for c in cells]
    out = {}
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for cell, rows in ex.map(_worker, args):
                out[cell] = rows
                log(f"  fitted {cell} ({len(rows)} check-ups)")
    else:
        for a in args:
            cell, rows = _worker(a)
            out[cell] = rows
            log(f"  fitted {cell} ({len(rows)} check-ups)")
    return out


def _ref_dict(fit):
    return dict(Cp=fit.Cp, Cn=fit.Cn, nLi=fit.nLi, yT=fit.yT, xT=fit.xT)


def choose_reference(rows, policy="first_avg"):
    """Return index into rows of the reference check-up (or None)."""
    ok = [i for i, r in enumerate(rows) if r["fit"] is not None and r["fit"].success]
    if not ok:
        return None
    if policy == "first_avg":
        for i in ok:
            if rows[i]["kind"] == "avg":
                return i
        return ok[0]
    if policy == "second_avg":      # first averaged check-up AFTER the first check-up (avoids the fresh-cell first charge)
        for i in ok:
            if rows[i]["kind"] == "avg" and rows[i]["checkup_idx"] > rows[0]["checkup_idx"]:
                return i
        return None
    if policy == "bol":
        return ok[0] if rows[ok[0]]["checkup_idx"] == rows[0]["checkup_idx"] else None
    raise ValueError(policy)


def assemble(rows_by_cell, ocp, cfg, meta=None, ref_policy="first_avg", ref_override=None, variant="main"):
    """
    rows_by_cell: output of fit_all. ref_override: {cell: dict(Cp,Cn,nLi,cap)} absolute reference (e.g. same-batch BOL of another cell).
    Returns DataFrame (one row per fitted check-up).
    """
    meta_idx = meta.set_index("cell") if meta is not None else None
    out = []
    for cell, rows in rows_by_cell.items():
        ri = choose_reference(rows, ref_policy)
        if ref_override and cell in ref_override:
            ref = ref_override[cell]
            ref_idx_label = ref.get("label", -1)
            ref_cap_meas = ref.get("cap_meas", np.nan)
            ref_not_bol = True
        elif ri is None:
            continue
        else:
            ref = _ref_dict(rows[ri]["fit"])
            ref_idx_label = rows[ri]["checkup_idx"]
            ref_cap_meas = rows[ri]["cap_dis"]
            ref_not_bol = rows[ri]["checkup_idx"] != rows[0]["checkup_idx"]
        cap_first = next((r['cap_dis'] for r in rows if r['cap_dis'] is not None and np.isfinite(r['cap_dis']) and r['cap_chg'] is not None and np.isfinite(r['cap_chg'])), np.nan)
        cap_ref_model = dc.cell_capacity(ref["Cp"], ref["Cn"], ref["nLi"], ocp, cfg.vmin, cfg.vmax)
        for i, r in enumerate(rows):
            f = r["fit"]
            rec = dict(cell=cell, checkup_idx=r["checkup_idx"], ah=r["ah"], T_mean=r["T_mean"], variant=variant,
                       mode_used=r["kind"], reason=r["reason"])
            if f is None or not f.success:
                rec.update(is_ref=False, flag_no_fit=True)
                out.append(rec)
                continue
            cv = r["curve"]
            y0, x0, s_e = dc.empty_state(f.p, ocp, cfg, cv["s_hi"])
            mo = dc.modes_vs_ref(f.Cp, f.Cn, f.nLi, ref)
            # multistart spread of the modes (pp) among converged starts
            sp = {"LLI": [], "LAM_PE": [], "LAM_NE": []}
            for sol in f.solutions:
                m = dc.modes_vs_ref(sol[0], sol[1], sol[2] * sol[0] + sol[3] * sol[1], ref)
                for k in sp:
                    sp[k].append(m[k])
            sd = {k: (100 * float(np.std(v)) if len(v) > 1 else np.nan) for k, v in sp.items()}
            idf = dc.identifiability(f)
            cap_model = dc.cell_capacity(f.Cp, f.Cn, f.nLi, ocp, cfg.vmin, cfg.vmax)
            cap_impl = dc.capacity_from_modes(ref, mo["LLI"], mo["LAM_PE"], mo["LAM_NE"], ocp, cfg.vmin, cfg.vmax)
            loss_impl = 1 - cap_impl / cap_ref_model if np.isfinite(cap_impl) and np.isfinite(cap_ref_model) else np.nan
            cap_meas = r["cap_dis"]
            soh_ref = cap_meas / ref_cap_meas if np.isfinite(cap_meas) and np.isfinite(ref_cap_meas) else np.nan   # relative to the mode reference
            soh = cap_meas / cap_first if np.isfinite(cap_meas) and np.isfinite(cap_first) else np.nan           # team convention (c20.csv): first check-up with both C/20 branches
            loss_meas = 1 - soh_ref if np.isfinite(soh_ref) else np.nan
            c1, c2 = dc.pass_criteria(f.rmse_mV, loss_impl, loss_meas)
            is_ref = (ref_idx_label == r["checkup_idx"]) and not (ref_override and cell in ref_override)
            rec.update(
                soh=soh, soh_vs_ref=soh_ref, Cp=f.Cp, Cn=f.Cn, y0=y0, x0=x0, LLI=mo["LLI"], LAM_PE=mo["LAM_PE"], LAM_NE=mo["LAM_NE"],
                rmse_mV=f.rmse_mV, n_starts_converged=f.n_converged, n_starts=f.n_starts,
                spread_LLI_pp=sd["LLI"], spread_LAM_PE_pp=sd["LAM_PE"], spread_LAM_NE_pp=sd["LAM_NE"],
                sd_Cp_pct=idf["sd_Cp_pct"], sd_Cn_pct=idf["sd_Cn_pct"], sd_nLi_pct=idf["sd_nLi_pct"], cond=idf["cond"],
                corr_Cp_Cn=idf["corr_Cp_Cn"],
                nLi=f.nLi, yT=f.yT, xT=f.xT, Voff_mV=(1e3 * f.p[4] if len(f.p) > 4 else np.nan),
                cap_dis_meas=cap_meas, cap_model=cap_model, loss_meas=loss_meas, loss_implied=loss_impl,
                loss_gap_pp=(100 * (loss_impl - loss_meas) if np.isfinite(loss_impl) and np.isfinite(loss_meas) else np.nan),
                pass_rmse=c1, pass_cap=c2, cap_evaluable=bool(np.isfinite(loss_impl) and np.isfinite(loss_meas)),
                is_ref=bool(is_ref), ref_checkup_idx=ref_idx_label,
                flag_ref_not_bol=bool(ref_not_bol), flag_curve_not_avg=bool(r["kind"] != "avg"),
                flag_at_bound=bool(len(f.at_bound) > 0),
                flag_multimodal=bool(np.nanmax(list(sd.values())) > 1.0) if np.isfinite(list(sd.values())).any() else False,
                flag_no_empty_state=bool(not np.isfinite(y0)),
                v_fit_lo=cv.get("v_lo"), v_fit_hi=cv.get("v_hi"), s_fit_hi=cv.get("s_hi"),
            )
            out.append(rec)
    df = pd.DataFrame(out)
    return df


def matched_value(g, col, level=0.10, xcol="loss_meas"):
    """
    Value of `col` at the point where the measured capacity loss first reaches `level` (linear interpolation between check-ups).
    g: DataFrame of one cell sorted by ah. Returns (value, how) with how in {'interp', 'last_not_reached', 'na'}.
    """
    h = g[[xcol, col, "ah"]].dropna().sort_values("ah")
    if len(h) < 2:
        return np.nan, "na"
    x, y = h[xcol].values, h[col].values
    if x.max() < level:
        return float(y[-1]), "last_not_reached"
    j = int(np.argmax(x >= level))
    if j == 0:
        return float(y[0]), "interp"
    w = (level - x[j - 1]) / (x[j] - x[j - 1])
    return float(y[j - 1] + w * (y[j] - y[j - 1])), "interp"
