function P = lifeParams()
% lifeParams.m : 수명 수식 모델(stepM2_life.m) 설정값을 한 곳에 모음 (재현성: 값을 바꾸면 변경 이력에 기록)
% 단위: 시간 [h], 용량 손실 Q = 1 - SOH [-], 저항 [Ohm]
% 변경 이력: 2026-11-19 최초 작성 (Claude)
%            2026-11-19 excludeCk 추가: 셀 13의 2682 Ah check-up은 36.6 °C에서 측정된 이상치 -> knee 진단에서 제외 (사용자 합의)
here = fileparts(mfilename('fullpath'));
P.anaOut  = fullfile(here, '..', 'analysis', 'out');   % Python 결과 폴더 (c20.csv, checkup_time.csv, P1_pred.csv)
P.ecmCsv  = fullfile(here, 'out', 'ecm_matlab.csv');   % stepM1_ecm.m 결과
P.outDir  = fullfile(here, 'out');
P.TH      = 0.85;                 % 수명 기준 SOH
P.SOH_E   = 0.975;                % 초기 측정 시점: SOH가 처음 이 값 이하가 된 check-up (P1과 동일)
P.SOH_E_SENS = [0.98 0.975 0.96]; % 민감도 점검 값
P.gridF   = linspace(0, 1, 101);  % ① 용량만: f 후보
P.gridA   = linspace(0, 1, 41);   % ③ f = clip(a + b*x): a 후보
P.gridB   = linspace(-1, 1, 61);  %                     b 후보 (x 최대 절댓값으로 나눠 씀)
P.grid3   = {linspace(0,1,21), linspace(-1,1,31), linspace(-1,1,21)};  % 확장(가중치 3개): a, b, c
P.nPerm   = 99;                   % 섞은 값 시험 횟수 (p = (k+1)/(nPerm+1), 최소 0.01)
P.seed    = 0;
P.accept  = 0.10;                 % 확장 항 채택: LOCO 오차가 10% 이상(상대) 줄고
P.pMax    = 0.05;                 %               섞은 값 시험 p < 0.05
P.spreadN2Max = 0.2;              % E1에서 쓰는 n2: 흔들림 20% 이하만
P.kneeGain = 0.005;
P.excludeCk = [13 2682];         % [셀, 누적 Ah] — 이 check-up(+-1 Ah)은 용량 궤적에서 뺌 (측정 온도 이상치, EIS_temperature.csv)               % E3 knee 진단: R^2가 이만큼 이상 오르면 'knee 신호'
end
