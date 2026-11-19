import pathlib, pandas as pd
E = pathlib.Path(__file__).resolve().parent.parent / "data" / "kit_radar" / "result_v2" / "Unnamed entity"
x = pd.read_excel(E / "Data Structure v05.3.xlsx", sheet_name=None, header=None)
keys = ["cyc_condition", "cyc_charged", "age_type", "age_profile", "cap_aged_est", "soh_cap", "timestamp_s", "age_soc"]
for s, df in x.items():
    for _, r in df.iterrows():
        vals = [str(v) for v in r.tolist() if str(v) != "nan"]
        line = " | ".join(vals)
        if any(k in line for k in keys):
            print(f"[{s}] {line[:500]}")
