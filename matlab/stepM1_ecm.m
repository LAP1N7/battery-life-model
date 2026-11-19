%% stepM1_ecm.m : 21개 셀 EIS -> ECM fitting -> out\ecm_matlab.csv + Python 결과와 대조
% 위치: ...\Mothat\matlab\  (데이터는 ...\Mothat\data, Python 결과는 ...\Mothat\analysis\out)
% 필요: Optimization Toolbox (lsqnonlin)
% 소요: 10분 안팎 예상 (스펙트럼 약 1,900개 x 시작점 6개, 미검증)
% 출력 열: cell, ah, soc, Rs, R1, R2, Rint(=R1+R2), tau1, tau2, ratio, spreadR1, spreadR2, nGood, atEdge, separable, relErr,
%          n1, n2, nw, spreadN2, n2AtBound   (11-19 추가: P5 비이상성 검증용. 기존 열과 값은 그대로)
% R1/R2 사용 규칙(11-18 합의): separable = 1 인 스펙트럼만 R1, R2를 따로 해석. 나머지는 Rint만.
% 검증: Python(analysis\out\eis_fits.csv)과 같은 스펙트럼끼리 Rs, Rint 상대 차이를 출력

here    = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', 'data');
outDir  = fullfile(here, 'out');
if ~exist(outDir, 'dir'), mkdir(outDir); end

CELLS = 1:21;     % 처음에는 1:2 로 빠르게 시험해 보고, 대조 결과가 괜찮으면 1:21 로
rows = [];
for c = CELLS
    spec = loadEIS(dataDir, c);
    for k = 1:numel(spec)
        s = spec(k);
        o = fitECM_mohtat(s.f, s.re, s.mim);
        rows(end+1, :) = [c, s.ah, s.soc, o.Rs, o.R1, o.R2, o.Rint, o.tau1, o.tau2, o.ratio, ...
                          o.spreadR1, o.spreadR2, o.nGood, o.atEdge, o.separable, o.relErr, ...
                          o.n1, o.n2, o.nw, o.spreadN2, o.n2AtBound]; %#ok<AGROW>
    end
    fprintf('cell %02d: 스펙트럼 %d개\n', c, numel(spec));
end
T = array2table(rows, 'VariableNames', {'cell', 'ah', 'soc', 'Rs', 'R1', 'R2', 'Rint', 'tau1', 'tau2', ...
    'ratio', 'spreadR1', 'spreadR2', 'nGood', 'atEdge', 'separable', 'relErr', 'n1', 'n2', 'nw', 'spreadN2', 'n2AtBound'});
writetable(T, fullfile(outDir, 'ecm_matlab.csv'));
fprintf('저장: %s\n', fullfile(outDir, 'ecm_matlab.csv'));
fprintf('fit 상대 오차: 중앙값 %.4f, 최대 %.4f\n', median(T.relErr), max(T.relErr));
fprintf('R1/R2 분리 가능(separable) 비율: %.1f%% (%d / %d)\n', 100*mean(T.separable), sum(T.separable), height(T));
fprintf('  탈락 사유: tau 비율 < 20 → %d개, 흔들림 > 20%% → %d개, tau가 측정 창 경계 → %d개\n', sum(T.ratio < 20), ...
        sum(T.ratio >= 20 & (T.spreadR1 > 0.2 | T.spreadR2 > 0.2)), sum(T.atEdge == 1));
fprintf('  Rs가 0에 붙은 스펙트럼(Rs < 0.1 mOhm): %d개 (수정 전 417개, 0에 가까워야 정상)\n', sum(T.Rs < 1e-4));
bol = T.ah < 1;
fprintf('  새것(BOL) 스펙트럼 중 분리 가능: %d / %d\n', sum(T.separable & bol), sum(bol));
use = T.separable == 1 & T.spreadN2 <= 0.2;   % P5 사용 조건 (11-19 합의)
fprintf('\n[P5용 n2] 중앙값 %.3f (범위 %.3f~%.3f) | separable 중 흔들림 20%% 이하: %d / %d | 경계(0.5 또는 1)에 붙음: %d개\n', ...
        median(T.n2), min(T.n2), max(T.n2), sum(use), sum(T.separable), sum(T.n2AtBound));

%% Python 결과와 대조 (같은 cell, ah, soc)
pyFile = fullfile(here, '..', 'analysis', 'out', 'eis_fits.csv');
if exist(pyFile, 'file')
    P = readtable(pyFile);
    P = P(:, {'cell', 'ah', 'soc', 'Rs', 'Rint', 'err'});
    P.Properties.VariableNames(4:6) = {'Rs_py', 'Rint_py', 'err_py'};
    T2 = T;  T2.ah = round(T2.ah, 6);  P.ah = round(P.ah, 6);       % 소수점 표기 차이로 짝이 안 맞는 것 방지
    J = innerjoin(T2, P, 'Keys', {'cell', 'ah', 'soc'});
    rsM = J.Rs; riM = J.Rint; rsP = J.Rs_py; riP = J.Rint_py;
    dRs = abs(rsM - rsP) ./ rsP;  dRi = abs(riM - riP) ./ riP;
    fprintf('\n[Python 대조] 짝지어진 스펙트럼 %d개\n', height(J));
    fprintf('  fit 오차: MATLAB이 더 좋음 %d / 같음(±1%%) %d / 더 나쁨 %d\n', sum(J.relErr < 0.99*J.err_py), ...
            sum(abs(J.relErr ./ J.err_py - 1) <= 0.01), sum(J.relErr > 1.01*J.err_py));
    fprintf('  Rs + Rint 합 상대 차이: 중앙값 %.2f%%\n', 100*median(abs((J.Rs + J.Rint) ./ (J.Rs_py + J.Rint_py) - 1)));
    fprintf('  Rs   상대 차이: 중앙값 %.2f%%, 5%% 넘는 경우 %d개\n', 100*median(dRs), sum(dRs > 0.05));
    fprintf('  Rint 상대 차이: 중앙값 %.2f%%, 5%% 넘는 경우 %d개\n', 100*median(dRi), sum(dRi > 0.05));
    ok = J.Rs_py >= 1e-4;   % Python에서도 Rs가 0에 붙지 않은 스펙트럼만 비교 (Python 412개는 같은 문제가 있었음)
    fprintf('  [Python Rs 정상 스펙트럼 %d개만] Rs 상대 차이 중앙값 %.2f%%, Rint 상대 차이 중앙값 %.2f%%\n', sum(ok), ...
            100*median(dRs(ok)), 100*median(dRi(ok)));
    fprintf('  참고: 이번 수정으로 tau를 측정 창 안으로 제한했기 때문에 fit 오차가 Python보다 약간 커질 수 있음 (제약을 둔 대가)\n');
else
    fprintf('Python 결과 파일이 없어 대조를 건너뜀: %s\n', pyFile);
end
