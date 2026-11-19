# prototype/ (트랙 C, 레포 밖 작업 사본)

트랙 C가 레포 밖 `prototype/`에서 한 분석 코드와 작은 결과 파일을 올린 것이다.
원본 데이터와 1 MB 이상 파일(파생 데이터 표 포함)은 올리지 않았다. 코드는 원래 폴더 구조(`켄/prototype`, `켄/data`, `켄/Mothat/data`) 기준 상대경로라 그대로 실행하려면 같은 구조가 필요하다.

| 경로 | 데이터 | 내용 |
|---|---|---|
| matlab/temp_mech/mohtat_modes | Mohtat | ICA 특징, 모드 시간 거동, 모드–저항 결합 |
| matlab/temp_mech/mohtat_eqs | Mohtat | E1–E5: 모드 대 손실, 시간 지수, ΔR = s·L^p (M1–M4) |
| matlab/temp_mech/model_dev | Mohtat | D1 기준점 민감도, D2 s의 구조 모델(S0–S4), D3 고온→상온·저온 전이. 사전 등록: prereg.md |
| matlab/temp_mech/temp_eqs | KIT + Mohtat | 같은 손실 비교(K), s(T) 아레니우스(S), 모드 부분상관(P) |
| matlab/temp_mech/kit_models | KIT | 아크 분해, 저항 모델 M1–M4 온도 하나 빼기 검증, WLTP 예측 |
| matlab/temp_mech/kit_optemp | KIT | 운전온도 대 25 °C 겉보기 활성화에너지 |
| kit_data.py, kit_inspect.py, stage1_kit_fit.py | KIT | 체크업 용량 추출, BLAST-Lite형 용량손실 모델, WLTP 예측 |
| kit_eis_signatures.py, mech_crosscheck.py | KIT (+M50T) | EIS·펄스·CE 서명 비교 |
| results/ | KIT | 위 스크립트의 요약 결과(1 MB 미만만) |

**주의:** KIT 관련 파일에는 결과가 들어 있다. B의 KIT 사전 등록 시험이 끝나기 전에는 B가 이 결과를 열지 않도록 한다.
MATLAB 코드는 Gemini(neon)가 작성하고 트랙 C가 실행·수정했다. 수정한 줄에는 `% 수정:` 주석이 있다.
