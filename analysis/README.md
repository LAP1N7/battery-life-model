# Mohtat 데이터 분석 (2026-11-13 1차)

## 실행 방법
1. Python 3.10 이상에서 패키지 설치
   pip install numpy pandas scipy scikit-learn matplotlib
2. 이 analysis 폴더에서 순서대로 실행 (data 폴더는 한 단계 위에 있어야 함)
   python step1_extract.py   # EIS fitting + C/20 특징 (수 분)
   python step2_paths.py     # H-alpha / 같은 SOH 다른 경로 검정
   python step3_ml1.py       # ML-1 정보량 판정
   python step4_plot.py      # 그림 out/results.png
결과 파일은 out 폴더에 저장됩니다.

## 1차 결과 (Claude 실행 기준, 다시 돌리면 수치가 약간 다를 수 있음)
- H-alpha: R_int(SOC)는 저SOC에서 커지는 모양, 최솟점이 90%(측정 끝)인 경우 25%.
  모양 지표(R20/R50, R90/R50)와 dV/dQ 피크 이동의 상관은 SOH 보정 후 유의하지 않음 -> 지지 안 됨.
- 같은 SOH, 다른 경로: SOH 보정 후 R_int 증가와 dV/dQ 피크 이동이 반비례 (전체 rho -0.63).
  단 온도 그룹별로는 저온 -0.74, 고온 -0.01, 상온 +0.16 -> 주로 저온 셀과 그룹 간 차이가 만든 관계.
- ML-1 (GPR, SOH RMSE %p): Ah 4.25 / |Z| 4.90 / phase 3.78 / physics 7.57 / physics+phase 3.91

## 한계
- LLI는 정밀 DVA가 아니라 dV/dQ 흑연 피크 위치로 대신한 간이 지표
- R1/R2 분리는 식별성 낮음 -> R_int = R1+R2 사용
- 일부 셀은 EIS check-up 수가 적음, BOL 피크 오검출 셀은 제외
