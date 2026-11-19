"""트랙 A 2단계 보강: check-up마다 '실제로 노화 사이클을 돈 시간(시간 단위)' 구하기 -> out/checkup_time.csv

왜 필요한가: 수명을 '충방전 횟수(EFC)'로만 재면, 느리게 충방전한 셀은 1회당 시간이 길어서
             그동안 시간에 비례해 진행되는 열화(SEI 성장 등, 특히 고온)가 더 쌓인다. 그래서 '시간' 기준도 본다.
방법:
  - OCV_wExpansion.csv(check-up)의 'Cycle number' = 그 check-up 직전까지 돈 노화 사이클 번호
  - cycling_wExpansion.csv(노화 사이클)의 'Time [s]'는 사이클링 내내 누적되는 시간 (셀 01에서 확인: 되감김 없음)
  - 두 파일을 사이클 번호로 이어서, 각 check-up 시점의 누적 사이클링 시간을 구함
주의: check-up 자체에 걸린 시간은 포함하지 않음 (모든 셀이 비슷한 check-up을 받으므로 비교에는 큰 영향 없음, 확인 필요)
소요: 셀마다 100 MB 안팎 파일을 읽어서 몇 분 걸릴 수 있음
"""
import pandas as pd
from common import DATA, OUT

rows = []
for cell in range(1, 22):
    try:
        o = pd.read_csv(DATA / f"{cell:02d}" / "OCV_wExpansion.csv", usecols=[7, 8])
        o.columns = ["cyc", "ah"]
        ck = o.groupby("ah").cyc.max().reset_index()                      # check-up별 (누적 Ah, 사이클 번호)
        c = pd.read_csv(DATA / f"{cell:02d}" / "cycling_wExpansion.csv", usecols=[0, 7])
        c.columns = ["t", "cyc"]
        tmax = c.groupby("cyc").t.max() / 3600                             # 사이클 번호별 누적 시간 [h]
        for _, r in ck.iterrows():
            n = int(r.cyc)
            hours = 0.0 if n == 0 else float(tmax[tmax.index <= n].iloc[-1]) if (tmax.index <= n).any() else float("nan")
            rows.append(dict(cell=cell, ah=r.ah, cyc=n, hours=hours))
        print(f"cell {cell:02d}: check-up {len(ck)}개, 마지막 {rows[-1]['hours']:.0f} h", flush=True)
    except Exception as e:
        print(f"cell {cell:02d}: 실패 - {e}")
pd.DataFrame(rows).to_csv(OUT / "checkup_time.csv", index=False)
print("저장:", OUT / "checkup_time.csv")
