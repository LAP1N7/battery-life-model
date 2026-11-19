"""Stage 1 (quick): refit a BLAST-Lite-style calendar + cycle model to KIT LG INR18650HG2 data.

Model form (same structure as BLAST-Lite NCA/NMC models, parameters refitted here):
  q_loss_cal = exp(lnA - Ea/R (1/T - 1/Tref)) * exp(c * SOC) * (t_days)^z
  q_loss_cyc = exp(lnB) * (exp(-Eh/R (1/T - 1/Tref)) + exp(lnW + El/R (1/T - 1/Tref)))   # high-T and low-T branch
               * (1 + kd*DoD) * (1 + kc*C_chg + kdis*C_dis) * EFC^zc
Calendar cells are fitted first; cycle cells are fitted on (measured loss - calendar loss during the test).
Known simplifications: check-up cycling of calendar cells is ignored; post-knee points (q < 0.7) are dropped.
"""
import json, pathlib, numpy as np, pandas as pd
from scipy.optimize import least_squares

HERE = pathlib.Path(__file__).resolve().parent
R, TREF = 8.314, 298.15
df = pd.read_csv(HERE / "results" / "kit_checkups.csv")
df["loss"] = 1 - df["q"]
df = df[df["q"] >= 0.7]

WIN = {(2.5, 4.2): (1.0, 0.50), (3.249, 4.2): (0.9, 0.55), (3.249, 4.092): (0.8, 0.50)}   # DoD, mean SOC
cal = df[df.age_type == 1].copy(); cal["soc"] = cal.cal_soc / 100
cyc = df[df.age_type == 2].copy()
cyc["dod"] = [WIN[(round(a, 3), round(b, 3))][0] for a, b in zip(cyc.v_min, cyc.v_max)]
cyc["soc"] = [WIN[(round(a, 3), round(b, 3))][1] for a, b in zip(cyc.v_min, cyc.v_max)]


def arr(T, Ea):
    return np.exp(-Ea / R * (1 / (T + 273.15) - 1 / TREF))


def f_cal(p, t, T, soc):
    return np.exp(p[0]) * arr(T, p[1]) * np.exp(p[2] * soc) * np.power(np.maximum(t, 0), p[3])


def f_cyc(p, efc, T, dod, cc, cd):
    temp = arr(T, p[1]) + np.exp(p[2]) * arr(T, -p[3])
    return np.exp(p[0]) * temp * (1 + p[4] * dod) * (1 + p[5] * cc + p[6] * cd) * np.power(np.maximum(efc, 0), p[7])


rc = least_squares(lambda p: f_cal(p, cal.day.values, cal["T"].values, cal.soc.values) - cal.loss.values,
                   [np.log(0.003), 30000, 1.0, 0.5], bounds=([-20, 0, -5, 0.2], [2, 150000, 10, 1.2]))
pc = rc.x
cyc["loss_cyc"] = cyc.loss - f_cal(pc, cyc.day.values, cyc["T"].values, cyc.soc.values)
ry = least_squares(lambda p: f_cyc(p, cyc.efc.values, cyc["T"].values, cyc.dod.values, cyc.chg_rate.values, cyc.dis_rate.values) - cyc.loss_cyc.values,
                   [np.log(0.002), 20000, np.log(0.5), 40000, 0.5, 0.5, 0.5, 0.6],
                   bounds=([-20, 0, -10, 0, -0.9, -0.5, -0.5, 0.2], [5, 80000, 10, 120000, 3, 10, 10, 1.5]))
py = ry.x
rm_c = np.sqrt(np.mean(rc.fun ** 2)); rm_y = np.sqrt(np.mean(ry.fun ** 2))
cal_share = (f_cal(pc, cyc.day.values, cyc["T"].values, cyc.soc.values) / cyc.loss.clip(lower=1e-3)).clip(0, 1)

# held-out check: WLTP profile cells (never fitted). mean SOC/DoD/C-rates approximated from the profile definition.
prof = df[df.age_type == 3].copy()
PROF = {2: (0.9, 0.55, 0.33, 0.35), 3: (0.8, 0.50, 0.33, 0.35), 4: (0.8, 0.50, 1.67, 1.0)}   # dod, soc, C_chg, C_dis(approx)
pr = []
for (P, T, pf), g in prof.groupby(["P", "T", "profile"]):
    dod, soc, cc, cd = PROF[pf]
    pred = f_cal(pc, g.day.values, g["T"].values, np.full(len(g), soc)) + f_cyc(py, g.efc.values, g["T"].values, dod, cc, cd)
    pr.append(dict(P=int(P), T=int(T), profile=int(pf), n=len(g), rmse=float(100 * np.sqrt(np.mean((pred - g.loss.values) ** 2)))))

out = dict(cal=dict(lnA=pc[0], Ea=pc[1], c=pc[2], z=pc[3], rmse_pp=100 * rm_c, n_points=len(cal)),
           cyc=dict(lnB=py[0], Eh=py[1], lnW=py[2], El=py[3], kd=py[4], kc=py[5], kdis=py[6], zc=py[7], rmse_pp=100 * rm_y, n_points=len(cyc)),
           cal_share_in_cycle_tests=dict(median=float(np.median(cal_share)), p10=float(np.percentile(cal_share, 10)), p90=float(np.percentile(cal_share, 90))),
           wltp_check=pr)
(HERE / "results" / "stage1_kit_params.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(f"calendar: Ea={pc[1]/1000:.1f} kJ/mol, SOC coef c={pc[2]:.2f} (SOC 0->1 x{np.exp(pc[2]):.2f}), z={pc[3]:.2f}, RMSE {100*rm_c:.2f} %p ({len(cal)} pts)")
print(f"cycle: Eh={py[1]/1000:.1f}, El={py[3]/1000:.1f} kJ/mol, kd={py[4]:.2f}, kc={py[5]:.2f}, kdis={py[6]:.2f}, zc={py[7]:.2f}, RMSE {100*rm_y:.2f} %p ({len(cyc)} pts)")
print("calendar share inside cycle tests: median {:.0%} (p10 {:.0%}, p90 {:.0%})".format(*out["cal_share_in_cycle_tests"].values()))
print("WLTP cells (held out), RMSE %p:", [(d["T"], d["profile"], round(d["rmse"], 1)) for d in pr])
