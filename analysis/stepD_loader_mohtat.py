# -*- coding: utf-8 -*-
"""
Mohtat adapter: reads <data_dir>/<NN>/OCV_wExpansion.csv (READ ONLY; nothing in the data folder is executed) and returns
the standard tables used by stepD_core / stepD*.

Check-up structure in OCV_wExpansion.csv (one 'Cycle number' per check-up): optional 1C discharge, optional rest,
C/20 charge (CC, then CV tail and rest at the top), C/20 discharge. Current in mA, +charge.
"""
import os
import numpy as np
import pandas as pd

C_NOM_AH = 5.0
I_C20_MA = 0.05 * C_NOM_AH * 1000.0      # 250 mA
TEMP_GROUP = {0: "R", 1: "C", 2: "H"}      # README Table 1 ordering within each condition triple
TEMP_NOMINAL = {"R": 25.0, "C": -5.0, "H": 45.0}
COND = {1: (100, 0.2, "C/5", 0.2), 2: (100, 1.5, "1.5C", 1.5), 3: (100, 2.0, "2C", 2.0), 4: (100, 0.2, "1.5C", 1.5),
        5: (50, 0.2, "C/5", 0.2), 6: (50, 0.2, "1.5C", 1.5), 7: (50, 1.5, "drive", None)}
CELLS = [f"{i:02d}" for i in range(1, 22)]
REQ_COLS = ["Capacity [Ah]", "Time [s]", "Current [mA]", "Voltage [V]", "Expansion [mu m]", "Temperature [C]",
            "Cycle number", "Total Ah throughput [Ah]"]

REST_I_MA = 20.0


def cell_meta(cell):
    n = int(cell)
    g = (n - 1) // 3 + 1
    t = TEMP_GROUP[(n - 1) % 3]
    dod, chg, dlab, dis = COND[g]
    # same-batch room-temperature mate: R cell of the same condition group
    mate = f"{3 * (g - 1) + 1:02d}"
    return dict(cell=cell, group=g, T_group=t, T_nom=TEMP_NOMINAL[t], dod=dod, chg_rate=chg, dis_label=dlab,
                room_mate=mate, c_nom_ah=C_NOM_AH)


def _runs(mask):
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    b = np.r_[0, np.where(np.diff(idx) > 1)[0] + 1, len(idx)]
    return [(int(idx[a]), int(idx[c - 1])) for a, c in zip(b[:-1], b[1:])]


def _cum_ah(t, I_mA):
    return np.concatenate([[0.0], np.cumsum(0.5 * (I_mA[1:] + I_mA[:-1]) * np.diff(t))]) / 3.6e6


def load_cell(cell, data_dir):
    """Returns (curves_df, rest_df) for one cell."""
    path = os.path.join(data_dir, cell, "OCV_wExpansion.csv")
    o = pd.read_csv(path, usecols=REQ_COLS)
    curves, rests = [], []
    for k, (cy, s) in enumerate(o.groupby("Cycle number", sort=True)):
        s = s.reset_index(drop=True)
        t = s["Time [s]"].values
        I = s["Current [mA]"].values
        V = s["Voltage [V]"].values
        X = s["Expansion [mu m]"].values
        T = s["Temperature [C]"].values
        thr = float(s["Total Ah throughput [Ah]"].max())
        base = dict(cell=cell, checkup_idx=k, cycle_number=int(cy), ah_throughput=thr)
        # same current windows as common.load_c20: charge 100<I<400 mA, discharge -400<I<-100 mA, >500 samples
        chg_runs = [g for g in _runs((I > 100) & (I < 400)) if g[1] - g[0] > 500]
        dis_runs = [g for g in _runs((I < -100) & (I > -400)) if g[1] - g[0] > 500]
        e_chg = None
        if chg_runs:
            a, b = max(chg_runs, key=lambda g: g[1] - g[0])
            # extend through the CV tail (current decaying, still charging, V ~ top)
            e = b
            while e + 1 < len(I) and I[e + 1] > REST_I_MA and V[e + 1] > V[b] - 0.01:
                e += 1
            e_chg = e
            q = _cum_ah(t[a:e + 1], I[a:e + 1])
            curves.append(pd.DataFrame(dict(base, direction="charge", Q_Ah=q, V=V[a:e + 1], T=T[a:e + 1],
                                            thickness_um=X[a:e + 1], t_s=t[a:e + 1] - t[a])))
        if dis_runs:
            a, b = max(dis_runs, key=lambda g: g[1] - g[0])
            q = -_cum_ah(t[a:b + 1], I[a:b + 1])
            cap_rep = float(s["Capacity [Ah]"].values[a:b + 1].max())
            curves.append(pd.DataFrame(dict(base, direction="discharge", Q_Ah=q, cap_cycler_Ah=cap_rep, V=V[a:b + 1], T=T[a:b + 1],
                                            thickness_um=X[a:b + 1], t_s=t[a:b + 1] - t[a])))
        if e_chg is not None and V[e_chg] > 4.15:      # only a charge that reached the top yields a top-of-charge rest
            r0 = e_chg + 1
            while r0 < len(I) and abs(I[r0]) >= REST_I_MA and I[r0] > 0:     # remaining CV samples
                r0 += 1
            r1 = r0
            while r1 < len(I) and abs(I[r1]) < REST_I_MA:
                r1 += 1
            if r1 - r0 >= 5:
                rests.append(pd.DataFrame(dict(base, phase="top_rest", t_rest_s=t[r0:r1] - t[r0], V=V[r0:r1],
                                               thickness_um=X[r0:r1], T=T[r0:r1])))
    curves = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    rests = pd.concat(rests, ignore_index=True) if rests else pd.DataFrame()
    return curves, rests


def load(data_dir, cells=None):
    """Standard dataset object: dict(curves, rest, meta)."""
    cells = cells or CELLS
    cs, rs, ms = [], [], []
    for c in cells:
        cu, re_ = load_cell(c, data_dir)
        cs.append(cu)
        rs.append(re_)
        ms.append(cell_meta(c))
    return dict(curves=pd.concat(cs, ignore_index=True), rest=pd.concat(rs, ignore_index=True), meta=pd.DataFrame(ms),
                info=dict(dataset="mohtat", vmin=3.0, vmax=4.2, c_nom_ah=C_NOM_AH, rate="C/20 (250 mA)", ocp="mohtat_peyman",
                          t_col="Temperature [C] (cell surface)"))
