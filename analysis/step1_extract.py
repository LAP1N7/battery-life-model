"""1단계: 21개 셀의 EIS fitting + C/20 특징 추출 -> out/eis_fits.csv, out/c20.csv, out/curves.pkl
소요: 수 분"""
import pickle, numpy as np, pandas as pd
from common import *

E, C, CV = [], [], {}
for cell in range(1, 22):
    spec = load_eis(cell); p0 = None
    for (ah, soc), x in spec.items():
        f, mag, ph, re, mim = x.T
        p, err = fit_ecm(f, re, mim, p0)
        lf = np.log10(f[::-1])
        E.append(dict(cell=cell, ah=ah, soc=soc, Rs=p[0], Rint=p[2]+p[5], err=err,
                      ph_10Hz=np.interp(1, lf, ph[::-1]), ph_1Hz=np.interp(0, lf, ph[::-1]),
                      ph_01Hz=np.interp(-1, lf, ph[::-1]), ph_001=ph[-1],
                      mag_1Hz=np.interp(0, lf, mag[::-1]), mag_01Hz=np.interp(-1, lf, mag[::-1])))
    c, cv = load_c20(cell); C.append(c); CV[cell] = cv
    print(f"cell {cell:02d}: EIS {len(spec)} spectra, C/20 {len(c)} check-ups", flush=True)
E = pd.DataFrame(E); C = pd.concat(C)
E.to_csv(OUT/"eis_fits.csv", index=False); C.to_csv(OUT/"c20.csv", index=False)
pickle.dump(CV, open(OUT/"curves.pkl", "wb"))
print("max fit error:", round(E.err.max(), 4))
