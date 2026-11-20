# 트랙 C: 메커니즘 (DVA · 팽창 · 저항)

Mohtat 21셀(UMBL 5 Ah NMC/흑연 파우치, −5 / 25 / 45 °C)을 대상으로, check-up마다 열화 모드(LLI, LAM_PE, LAM_NE)를 분리하고 저항(R0, 아크)과 연결한다.
작업은 두 갈래로 나뉜다.

| 갈래 | 질문 | 스크립트 | 산출물 |
|---|---|---|---|
| ① 원래 트랙 C | A의 √t·선형 성분이 각각 어느 모드와 맞나 (H-D1, H-D2) | stepD0 ~ stepD3 | D0_*, D1_*, D2_*, D3_* |
| ② 온도 메커니즘 | 같은 용량 손실에서 온도별로 저항과 모드가 어떻게 다른가, 식으로 쓰면 무엇인가 | stepD4 ~ stepD6 | D4_*, D5_*, D6_* |

공용 모듈: `stepD_core.py`(DVA 피팅 핵심), `stepD_loaders.py`, `stepD_loader_mohtat.py`, `stepD_ocp_mohtat.py`(전극 OCP), `stepD_pipeline.py`.
원본 데이터(`data/`)는 올리지 않았다. 각자 내려받아 같은 위치에 둔다.

## 스크립트와 산출물

| 스크립트 | 하는 일 | 주요 산출물 |
|---|---|---|
| stepD0_selftest.py | 합성 곡선으로 DVA 피팅 자기검증 (잡음, OCP 불일치) | D0_selftest.csv, D0_selftest_summary.csv, D0_selftest.png |
| stepD1_dva.py | check-up별 DVA 피팅: Cp, Cn, y0, x0 → LLI, LAM_PE, LAM_NE | D1_dva.csv, D1_sensitivity.csv, D1_modes_at_10pct.csv, D1_fit_cell01~03.png |
| stepD2_expansion.py | 같은 SOC·휴지 조건에서 비가역 두께 변화 | D2_expansion.csv |
| stepD3_compare_prep.py | A의 P1_pred(A_full, B_full)와 D1을 셀 단위로 묶고 H-D1·H-D2 순위상관 | D3_join.csv, D3_correlations.csv (+ _standin) |
| stepD4_resistance.py | EIS에서 R0(−Im = 0 교차)와 아크 저항 추출, 온도별 지문 | D4_resistance.csv, D4_cells.csv, D4_arrhenius.csv, D4_nyquist.png, D4_fingerprint.png |
| stepD5_resistance_model.py, stepD5b_holdout.py | 저항 증가 모델 피팅, 셀 holdout 검증 | D5_params.csv, D5_shares.csv, D5_fit.png |
| stepD6_mode_resistance.py | 같은 손실(5/10/20 %)에서 모드와 저항 비교 | D6_equal_loss.csv, D6_group.csv, D6_corr.csv, D6_modes_resistance.png |

전극 OCP 출처: PyBaMM 파라미터 세트 Mohtat2020 (`graphite_ocp_PeymanMPM`, `NMC_ocp_PeymanMPM`). 이 파일은 셀을 NMC532로 적고 있어, 데이터 논문 쪽 NMC111 표기와 다르다.

## 지금까지의 결과

### ① 원래 트랙 C
- **DVA 피팅:** 248 check-up 중 247개가 RMSE < 10 mV (중앙값 3.0 mV). 예외는 cell 07 check-up 6 (10.2 mV).
- **용량 손실 일치 (1 %p 이내):** 기준 check-up 제외 53/220 (24 %)만 통과 (차이 중앙값 +1.24 %p, 항상 양수). 기준 check-up을 포함한 74/241은 기준 행이 정의상 0이라 부풀려진 값이다. 제안 완료 기준을 충족하지 못했다.
  - 주원인: BOL에서 모델 용량이 실측보다 약 1 % 크고 그 오프셋이 이후로 전파된다. 기준을 BOL 다음 check-up으로 옮기면 193/220 통과. OCP 차이(NMC532 vs 111)와 전압창은 부분 원인.
  - LAM_NE 셀 간 순위는 설정 변형에 안정적(ρ 0.97–0.99)이나, 45 °C 셀과 나머지의 기준 처리 차이로 그룹 간 약 3 %p 편향이 있다.
- **고온 셀 64행:** 첫 C/20이 BOL이 아니다 (`flag_ref_not_bol`). 첫 측정을 기준으로 삼았다.
- **H-D1 (√t 비중 ~ LLI 비중) / H-D2a (선형 비중 ~ LAM 비중):**
  - 전체 ρ = 0.69, p = 0.002 (n 17)로 제안 기준(ρ ≥ 0.5, p < 0.05)을 통과한다.
  - 두 검정은 서로 1 − x 관계라 독립 근거는 하나다.
  - 온도 안에서 보면 상온 ρ = 0.86 (p 0.01). 저온은 −0.56 (n 5)으로 판단할 수 없다.
- **H-D2b (비가역 팽창 ~ 선형 성분):** ρ = −0.23, p = 0.37. 지지되지 않는다.

### ② 온도 메커니즘 (탐색적)
1. 같은 용량 손실에서 LLI는 세 온도가 거의 같다 (용량 손실 ≈ LLI). 온도가 바꾸는 것은 함께 생기는 LAM과 저항이다.
2. 손실당 저항 증가는 고온이 가장 크다. KIT 데이터에서도 같은 방향이다.
3. ΔR = s_cell · L^p_T에서 온도에 따라 가장 크게 바뀌는 상수는 p다. 저온 약 1.7 (가속), 상온 0.9–1.2, 고온 0.75–0.88 (둔화). Mohtat에서만 확인했다.
4. "저온 = 아크형" 패턴은 Mohtat에서만 보이고, KIT에서는 반대 방향이다.
5. 기울기 s(T)를 아레니우스 식 하나로 쓰는 것은 두 데이터셋에 공통으로 맞지 않는다.
6. 세 모드가 같이 움직여 (온도 안 상관 r 0.85–0.995) 모드별 저항 계수는 분리할 수 없다.

### ③ 모델 개발 D1–D3 (사전 등록: `prototype/matlab/temp_mech/model_dev/prereg.md`, Mohtat 21셀)
- **D1 기준점 민감도:** 저온·상온 셀의 기준을 같은 프로토콜 고온 셀의 기준 Ah 이후로 옮겨 다시 맞췄다.
  - p 저온 − 고온: R0 0.91 [0.15, 1.52], 아크 0.71 [0.29, 1.02]로 기준을 맞춰도 남는다 → **저온 가속은 기준점 탓이 아니다 (지지).**
  - p 상온 − 고온: CI가 0을 포함 → 기준점 효과와 구분 불가.
  - 시간 지수 z: 기준을 맞추면 LLI·LAM_NE의 온도 차이는 사라지고 LAM_PE는 부호가 뒤집힌다 → **z의 온도 차이 주장은 철회.** z는 시간 원점에 따라 크게 바뀌므로 원점이 다른 셀끼리 비교하지 않는다.
- **D2 구조 모델:** 셀별 s 대신 온도·충전 속도로 쓴 식을 셀 하나 빼기(LOCO)로 검증했다.
  - 아크: 충전 속도를 넣은 모델의 LOCO 오차가 온도만 쓴 모델의 0.84–0.85배 → **충전 속도가 셀 간 차이를 설명한다 (지지).** R0는 설명하지 못했다(0.99–1.03배).
  - 실측 셀 온도는 주변 온도보다 낫지 않았다(0.995배). 상온 급속충전 셀 4·7번은 사이클 중 실제로 6–7 °C 높았다.
  - 21셀 공통 식(상수 5개, T °C, chg C-rate, L %, ΔR mΩ):
    ΔR_arc = exp(−4.05 + 0.070·T − 0.148·ln chg) · L^(1.69 − 0.0231·T),
    ΔR0 = exp(−3.97 + 0.067·T − 0.088·ln chg) · L^(1.50 − 0.0192·T).
    셀별 s 모델(R² 0.97)보다 설명력은 낮다(R² 0.90).
- **D3 고온 전이:** 고온 셀로만 맞춘 식으로 상온·저온 셀을 예측했다(예측/실측, 0.8–1.25면 재현).
  - LLI: 0.82–1.13 → **재현됨.**
  - 저항: 1.4–3.5배 → **재현 안 됨 (과대 예측).** 같은 손실 기준에서 고온 시험은 저항을 크게 잡는다. 시간 축 기준(D5b, 저온 0.17–0.22배 과소)과는 축이 다르니 보고할 때 축을 밝힌다.
  - LAM_PE 재현 안 됨, LAM_NE 저온 재현 안 됨, 상온은 대체로 재현.
- **정리:** 주장은 "고온 가속은 용량(LLI)은 재현하지만 저항과 모드 구성은 재현하지 못한다"로 좁힌다.

E1~E5 식 피팅, s(T), KIT 비교, D1–D3에 쓴 MATLAB 코드는 `prototype/matlab/temp_mech/`에 있다.

## 해야 할 것

- [ ] **C2 사양표 A와 합의:** 입력 곡선, 맞출 파라미터, 판정 기준. 합의 기록이 아직 없다.
- [x] **용량 손실 일치 미달 원인 정리:** 위 ①에 정리. 표기는 53/220으로 통일하고 기준 오프셋을 뺀 차이를 함께 보고하기로 A와 합의.
- [ ] **BOL 모델 용량이 실측보다 약 1 % 큰 원인 확인:** 첫 사이클 상태, 전압창 밖 외삽, BOL 첫 낙폭 중 무엇인지.
- [ ] **OCP 불일치 자기검증 보강:** 불일치 상황에서 적합한 모드로 계산한 용량 손실까지 비교. 양극이 NMC111인지 532인지 데이터 논문에서 확인.
- [ ] **H-D1·H-D2 판정을 sync용으로 정리:** 독립 근거가 하나라는 점과 저온 반대 방향을 함께 적는다.
- [x] **step5 팽창 45 °C 버그 수정** (가장 긴 C/20 충전 구간, 기준점 BOL).
- [ ] **독립 데이터로 D3 반복:** 흑연 단독 셀 데이터(Medina 2026, Röder 2025 후보)에서 "LLI 재현, 저항 과대"가 나오는지. 범위는 실행 전에 적어 둔다.
- [ ] **R0 구조 모델:** 아크와 달리 충전 속도로 설명되지 않음. 다른 설명 변수 탐색(탐색적).
- [ ] **Mohtat EIS 측정 온도 확인** (`Resistance.csv`의 T 열).
- [ ] **팀 M2 판정과 온도별 분류표를 sync 공유용으로 정리.**
- 보류: KIT 관련 추가 분석(B의 사전 등록 시험이 끝날 때까지 열지 않음). Mohtat와 KIT는 음극(흑연 vs 흑연+SiO)과 저온 조건(−5 vs 0 °C)이 달라 저온 결과는 직접 비교 대상이 아니다.

## 실행

```
cd analysis
python stepD0_selftest.py
python stepD1_dva.py
python stepD2_expansion.py
python stepD3_compare_prep.py
python stepD4_resistance.py
python stepD5_resistance_model.py
python stepD5b_holdout.py
python stepD6_mode_resistance.py
```
결과는 모두 `analysis/out/D/`에 저장된다. `common.py`는 수정하지 않았다.
