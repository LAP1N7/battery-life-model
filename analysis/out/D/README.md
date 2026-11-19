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
- **용량 손실 일치 (1 %p 이내):** 241개 중 74개만 통과 (차이 중앙값 1.2 %p). 제안 완료 기준을 충족하지 못했다.
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

E1~E5 식 피팅, s(T), KIT 비교에 쓴 MATLAB 코드는 팀 규칙(`matlab/`은 트랙 A)에 따라 레포 밖 개인 `prototype/matlab/temp_mech/`에 있다.

## 해야 할 것

- [ ] **C2 사양표 A와 합의:** 입력 곡선, 맞출 파라미터, 판정 기준. 합의 기록이 아직 없다.
- [ ] **용량 손실 일치 미달 (74/241) 원인 정리:** 기준은 바꾸지 않고 그대로 보고한다. 원인 후보는 OCP 불일치(NMC532 vs 111)와 전압 창.
- [ ] **고온 셀 BOL 부재 처리:** 같은 배치 상온 셀 BOL로 보완할지 sync에서 결정한다.
- [ ] **H-D1·H-D2 판정을 sync용으로 정리:** 독립 근거가 하나라는 점과 저온 반대 방향을 함께 적는다.
- [ ] **KIT에서도 p가 온도에 따라 바뀌는지 확인:** 결론 3을 두 데이터셋 근거로 올릴 수 있는지를 정한다.
- [ ] **Mohtat와 KIT의 저온 결과가 반대인 이유:** 셀 종류, 충전 속도, EIS 측정 조건을 비교한다.
- [ ] **팀 M2 판정과 온도별 분류표를 sync 공유용으로 정리.**
- [ ] **temp_mech MATLAB 재실행 확인:** 폴더 이동 뒤 경로 3줄을 고쳤으나 아직 실행 검증은 하지 않았다.

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
