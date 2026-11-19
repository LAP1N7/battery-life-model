"""Tidy check-up table for the KIT LG INR18650HG2 dataset (Luh & Blank, Sci. Data 2024, RADAR v2 result data).

One row per check-up capacity measurement:
  cell, P (parameter id), age_type (1 calendar, 2 cyclic, 3 profile), T, cal_soc, v_min, v_max,
  chg_rate, dis_rate, profile, day, efc (discharged Ah / 3 Ah), cap_Ah, q (cap / first cap)
"""
import pathlib
import numpy as np, pandas as pd

E = pathlib.Path(__file__).resolve().parent.parent / "data" / "kit_radar" / "result_v2" / "Unnamed entity"
OUT = pathlib.Path(__file__).resolve().parent / "results" / "kit_checkups.csv"
QNOM = 3.0


def build():
    rows = []
    for cf in sorted((E / "cfg").glob("cell_cfg_P*.csv")):
        cfg = pd.read_csv(cf, sep=";").iloc[0]
        if int(cfg["cell_used"]) != 1 or int(cfg["parameter_id"]) == 0:
            continue
        ef = E / "cell_eocv2" / cf.name.replace("cell_cfg_", "cell_eocv2_")
        if not ef.exists():
            continue
        d = pd.read_csv(ef, sep=";")
        if d.empty:
            continue
        cu = d[(d["cyc_condition"] == 2) & (d["cyc_charged"] == 0)].copy()     # check-up discharges
        cu["cap"] = cu["delta_q_dischg_Ah"].abs()
        cu = cu[cu["cap"] > 0.5]
        if cu.empty:
            continue
        # one value per check-up block: keep the largest discharge (full 1/3 C capacity discharge)
        cu["blk"] = (cu["timestamp_s"].diff().fillna(1e9) > 6 * 3600).cumsum()
        cu = cu.loc[cu.groupby("blk")["cap"].idxmax()]
        t0 = d["timestamp_s"].iloc[0]
        for _, r in cu.iterrows():
            rows.append(dict(cell=cf.stem.replace("cell_cfg_", ""), P=int(cfg["parameter_id"]), age_type=int(cfg["age_type"]),
                             T=int(cfg["age_temp"]), cal_soc=int(cfg["age_soc"]), v_min=float(cfg["V_min_cyc_V"]), v_max=float(cfg["V_max_cyc_V"]),
                             chg_rate=float(cfg["age_chg_rate"]), dis_rate=float(cfg["age_dischg_rate"]), profile=int(cfg["age_profile"]),
                             day=(r["timestamp_s"] - t0) / 86400.0, efc=abs(r["total_q_dischg_sum_Ah"]) / QNOM, cap_Ah=r["cap"]))
    df = pd.DataFrame(rows).sort_values(["cell", "day"])
    df["q"] = df["cap_Ah"] / df.groupby("cell")["cap_Ah"].transform("first")
    df.to_csv(OUT, index=False)
    return df


if __name__ == "__main__":
    df = build()
    print(len(df), "check-ups,", df["cell"].nunique(), "cells")
    s = df.groupby(["age_type", "P", "T", "cal_soc", "v_min", "v_max", "chg_rate", "dis_rate", "profile"]).agg(
        cells=("cell", "nunique"), days=("day", "max"), efc=("efc", "max"), q_end=("q", "last")).reset_index()
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 200)
    print(s.to_string(index=False))
