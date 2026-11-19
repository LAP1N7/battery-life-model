"""3단계: ML-1 정보량 판정. 셀 하나씩 빼는 교차검증(leave-one-cell-out)으로 SOH 예측 오차 비교 -> out/ml1.csv"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as Ck
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from common import *

E = pd.read_csv(OUT/"eis_fits.csv"); F = pd.read_csv(OUT/"checkups.csv")
e50 = E[E.soc == 50][["cell", "ah", "mag_1Hz", "mag_01Hz", "ph_10Hz", "ph_1Hz", "ph_01Hz", "ph_001"]]
D = F.merge(e50, on=["cell", "ah"]).dropna(subset=["soh"])
SETS = {"Ah throughput (baseline)": ["ah"],
        "|Z| (1, 0.1 Hz)": ["mag_1Hz", "mag_01Hz"],
        "phase (10..0.01 Hz)": ["ph_10Hz", "ph_1Hz", "ph_01Hz", "ph_001"],
        "physics (Rs, Rint, shape)": ["Rs50", "Rint50", "rlow", "rhigh"],
        "physics + phase": ["Rs50", "Rint50", "rlow", "rhigh", "ph_10Hz", "ph_1Hz", "ph_01Hz", "ph_001"]}
def gpr(k): return make_pipeline(StandardScaler(), GaussianProcessRegressor(Ck()*RBF(np.ones(k)) + WhiteKernel(), normalize_y=True))
def ridge(k): return make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 13)))
res = []
y = D.soh.values * 100
for name, cols in SETS.items():
    X = D[cols].values
    for mname, mk in [("GPR", gpr), ("Ridge", ridge)]:
        pred = np.full(len(y), np.nan)
        for c in D.cell.unique():
            tr = (D.cell != c).values
            m = mk(len(cols)).fit(X[tr], y[tr]); pred[~tr] = m.predict(X[~tr])
        res.append(dict(features=name, model=mname, RMSE=np.sqrt(np.mean((pred - y)**2))))
R = pd.DataFrame(res).pivot(index="features", columns="model", values="RMSE").reindex(SETS.keys()).round(2)
print(R.to_string()); print(f"n={len(D)} check-ups, {D.cell.nunique()} cells")
R.to_csv(OUT/"ml1.csv")
