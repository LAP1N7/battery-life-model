"""load_c20 수정 제안 (common.py는 수정하지 않음; 이 파일은 common.py를 import하지 않는 독립 사본).

변경 이력
  2026-11-19 (C5) 제안 작성. common.load_c20(54-67줄)과 동일한 함수에 '가장 긴 연속 구간(longest run)'만 적용.
    - 문제: 100-400 mA / -400..-100 mA 필터에 1C 충전 CV 꼬리(전류 100-400 mA로 감쇠, V = 4.2 V)가 섞임.
      45 C 셀(3,6,9,12,15,18,21)의 ah > 0 check-up 중 42개에서 이 꼬리(약 60-70 샘플)가 C/20 충전 '앞'에 위치하고
      cap 열이 3.4 Ah 부근 -> 0 으로 리셋되어 충전 q가 단조 증가가 아님 -> np.interp(grid, q, V)가 오염.
      (|dV| 최대 977 mV, dV/dQ 차이 최대 6.3 V/Ah, 흑연 피크 위치 xpk가 15개 check-up에서 변함)
    - 방전 쪽은 run 밖 샘플(전환 구간 1-81개)이 있어도 dc.cap.max() 차이가 3.3e-4 Ah(0.007 %) 이하라 사실상 무영향이지만,
      충전과 동일 규칙을 적용해 일관성 확보.
    - 파일 순서 = 시간 순서임을 확인함(모든 check-up에서 t 단조 증가). 그래도 run 판정이 순서에 의존하므로 t 기준 안정 정렬 추가.
    - 그 외(열 이름, 500 샘플 기준, grid, savgol 31/3, 반환 형식)는 동일. 500 샘플 기준은 run 길이에 적용;
      현 데이터에서 통과/탈락 check-up 집합은 바뀌지 않음(496개 중 차이 0).
    - 주의: 셀 9, ah = 1430.9 check-up은 C/20 충전이 3.71 V(0.93 Ah)에서 끊긴 불완전 곡선. longest run을 써도 곡선이 부분 구간뿐이므로
      step2처럼 grid.max() 비율에 의존하는 분석에서는 이 check-up을 제외하는 것을 권장(여기서는 필터를 추가하지 않음).
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"


def longest_run(mask):
    """True가 연속된 가장 긴 구간의 위치 인덱스 (CV 꼬리/전환 샘플 같은 짧은 구간 제외). step5_expansion_TODO.longest_run과 동일."""
    idx = np.flatnonzero(mask)
    if not len(idx): return idx
    runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    return max(runs, key=len)


def load_c20(cell):
    """OCV_wExpansion.csv -> check-up별 C/20 방전용량, 충전곡선의 dV/dQ (C/20 충/방전 = 가장 긴 연속 구간만)"""
    d = pd.read_csv(DATA / f"{cell:02d}" / "OCV_wExpansion.csv")
    d.columns = ["t", "I", "V", "x", "T", "Q", "cap", "cyc", "ah"]
    rows, curves = [], {}
    for ah, v in d.groupby("ah"):
        v = v.sort_values("t", kind="stable")
        ch = v.iloc[longest_run(((v.I > 100) & (v.I < 400)).to_numpy())]      # C/20 = 약 250 mA
        dc = v.iloc[longest_run(((v.I < -100) & (v.I > -400)).to_numpy())]
        if len(ch) < 500 or len(dc) < 500: continue
        q, V = ch.cap.values, ch.V.values
        grid = np.arange(0.02, q.max()-0.02, 0.01)
        Vg = np.interp(grid, q, V)
        curves[ah] = (grid, Vg, savgol_filter(Vg, 31, 3, deriv=1, delta=0.01))
        rows.append(dict(cell=cell, ah=ah, cap=dc.cap.max()))
    return pd.DataFrame(rows), curves
