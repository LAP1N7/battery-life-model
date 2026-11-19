function spec = loadEIS(dataDir, cell)
%LOADEIS  Mohtat Resistance.csv -> check-up별 EIS 스펙트럼
%   spec = loadEIS(dataDir, cell)
%   입력  dataDir : ...\Mothat\data 폴더 경로 (char/string)
%         cell    : 셀 번호 (1~21)
%   출력  spec    : 구조체 배열, 원소 하나 = 스펙트럼 하나
%                   .ah   누적 충방전량 [Ah] (check-up 구분용)
%                   .soc  SOC [%] (90, 80, ..., 10)
%                   .f    주파수 [Hz] (높은 주파수 -> 낮은 주파수)
%                   .re   Re(Z) [Ohm]
%                   .mim  -Im(Z) [Ohm]
%   규칙 (Python common.load_eis와 동일):
%     - f > 0 인 행만 사용
%     - 주파수가 다시 커지거나 시간 간격이 600 s를 넘으면 새 스펙트럼
%     - check-up 하나에 스펙트럼 9개(SOC 90 -> 10%)가 다 있는 경우만 사용
%   열 순서 (Resistance.csv): t, I, V, f, |Z|, phase, Re, -Im, T, Q, cycle, Ah

SOCS = [90 80 70 60 50 40 30 20 10];
M = readmatrix(fullfile(dataDir, sprintf('%02d', cell), 'Resistance.csv'));
t = M(:,1); f = M(:,4); re = M(:,7); mim = M(:,8); ah = M(:,12);

keep = f > 0;
t = t(keep); f = f(keep); re = re(keep); mim = mim(keep); ah = ah(keep);

newSpec = [false; diff(f) > 0 | diff(t) > 600];   % 새 스펙트럼 시작 표시
sid = cumsum(newSpec);

spec = struct('ah', {}, 'soc', {}, 'f', {}, 're', {}, 'mim', {});
ahList = unique(ah);
for a = ahList'
    inA = (ah == a);
    sids = unique(sid(inA));
    if numel(sids) ~= 9, continue; end              % SOC 9개가 다 있는 check-up만
    for k = 1:9
        r = inA & (sid == sids(k));
        [fs, order] = sort(f(r), 'descend');
        reK = re(r); mimK = mim(r);
        spec(end+1) = struct('ah', a, 'soc', SOCS(k), 'f', fs, ...
                             're', reK(order), 'mim', mimK(order)); %#ok<AGROW>
    end
end
end
