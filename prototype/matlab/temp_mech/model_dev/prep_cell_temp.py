import numpy as np, pandas as pd, os
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "data")   # Mothat/data
J = pd.read_csv(r"..\mohtat_modes\mohtat_modes_joined.csv")
amb = J.groupby("cell")["T_C"].first()
rows = []
for c in range(1, 22):
    f = os.path.join(D, "%02d" % c, "cycling_wExpansion.csv")
    sw = sc = nw = nc = 0.0; n = 0; prev_t = None
    for ch in pd.read_csv(f, usecols=["Time [s]", "Current [mA]", "Temperature [C]"], chunksize=500000):
        t = ch["Time [s]"].values.astype(float); i = ch["Current [mA]"].values.astype(float); T = ch["Temperature [C]"].values.astype(float)
        n += len(t)
        # 각 샘플의 가중치 = 다음 샘플까지 시간 간격 (청크 경계 처리 위해 이전 마지막 시간 이어붙임)
        tt = np.r_[t, t[-1]] if False else t
        if prev_t is not None:
            tt = np.r_[prev_t[0], t]; dt = np.diff(tt)  # 이전 청크 마지막 샘플의 간격은 별도 처리 불가 -> 후행 간격 사용
            # 단순화: 샘플 j의 가중치 = t[j]-t[j-1] (후행 차분)
        else:
            dt = np.r_[0.0, np.diff(t)]
        prev_t = [t[-1]]
        dt = np.clip(dt, 0, 600)  # 큰 결측 구간 제한
        ok = np.isfinite(T) & (T > -40) & (T < 100) & np.isfinite(dt)
        sw += (T[ok] * dt[ok]).sum(); nw += dt[ok].sum()
        okc = ok & (i > 0)
        sc += (T[okc] * dt[okc]).sum(); nc += dt[okc].sum()
    rows.append(dict(cell=c, T_eff_mean=sw / nw, T_eff_chg_mean=sc / nc, T_amb=amb[c], n_rows=n))
df = pd.DataFrame(rows); df["dev"] = df.T_eff_mean - df.T_amb
print(df.round(2).to_string(index=False))
df.drop(columns="dev").to_csv("cell_temp.csv", index=False)
print("anomalies:", df[df.dev.abs() > 5].cell.tolist())
