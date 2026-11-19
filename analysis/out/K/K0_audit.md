# K0_audit — KIT 데이터 점검 결과 (트랙 B, 과업 B2)

> `stepK0_audit.py`가 자동으로 작성. **"확인 필요"는 사람이 확인해야 함.** 열 이름 전체 목록: `out/K/K0_columns.txt`

## 0. 내려받은 파일
- 폴더: `C:\Users\user\Documents\GitHub\battery-life-model\analysis\out\K`
- EOC 240개, EIS 468개, CFG 240개, 남은 압축 파일 6개

## ① SOH 0.85 도달 셀 수
- SOH = 각 check-up 용량 / 첫 check-up 용량 (사전 등록과 동일). 용량 열: `cap_aged_est_Ah`
- **도달 셀 수: 187 / 240**  (기준 10개 이상)
- check-up 구분: 조건 열 `cyc_condition`의 값 중 'cu'가 들어간 행. 조건 열 값 상위: `{}`

## ② check-up별 시간·누적 충전량
- 시간 열: `timestamp_s` (있음)
- 누적 충전량 열: `total_q_chg_CU_RT_Ah` (있음 — EFC = 누적 Ah / (2 x 3 Ah)로 계산 가능)

## ③ 첫 check-up 상온 EIS
- **첫 check-up에 상온 EIS가 있는 셀: 456 / 468**
- 자동으로 찾은 EIS 열: `{'시간': 'timestamp_s', '주파수': 'freq_Hz', '온도': 't_avg_degC', 'SoC': 'soc_nom', '실수부': 'z_re_comp_mOhm', '허수부': 'z_im_comp_mOhm', '조건': None}`
- 가정: 상온 = 평균 온도 20.0~30.0 °C, 첫 check-up = 첫 EIS 측정 후 3일 이내

## ④ 셀별 조건표
- `out/K/K0_conditions.csv` (셀, 파라미터 ID, cfg 열, check-up 수, 마지막 SOH, 0.85 도달, 첫 상온 EIS)

## 발견한 문제 / 확인 필요
- 압축 파일 6개가 남아 있음 (풀었는지 확인)
- EIS '조건' 열을 자동으로 못 찾음 -> 확인 필요
- 온도 이상 구간(2022-11-27~29, 2023-04 25 °C 풀) 표시: 아직 안 함 — 로더(B3)에서 처리

## 결론 (정오 sync용)
- KIT 사용: **진행** — SOH 0.85 도달 187셀
