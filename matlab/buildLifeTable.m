function D = buildLifeTable(P, sohE)
% buildLifeTable.m : 셀별 [수명, 초기 시점 te, 초기 손실 Qe, 초기 저항 증가율 x, 초기 n2 변화율 z, 온도] 표 만들기
%   입력  P: lifeParams(),  sohE: 초기 측정 시점 기준 SOH (예: 0.975)
%   출력  D: 구조체 (필드마다 셀 수 길이의 열 벡터)
%     cell, temp('R'/'C'/'H'), life [h] (SOH 0.85 도달, 직선 보간), te [h], Qe [-],
%     x = (R_e/R_BOL - 1)/sqrt(te) [1/sqrt(h)],  R = Rs+R1+R2 의 SOC 평균 (MATLAB ECM)
%     z = (n2_e - n2_BOL)/sqrt(te)  (n2 흔들림 <= spreadN2Max 인 스펙트럼의 SOC 평균, 없으면 NaN)
%   규칙 (Python stepP1_prototype.py 와 동일)
%     - SOH = check-up 용량 / 그 셀의 첫 check-up 용량,  시간 = checkup_time.csv 의 hours
%     - 초기 시점 = SOH가 처음 sohE 이하인 check-up,  그 시점 +-1 Ah 안에 EIS가 없으면 셀 제외 (미래 값 사용 금지)
%     - 수명이 0.85에 도달하지 않은 셀 제외
% 변경 이력: 2026-11-19 최초 작성
[Hc, Mc] = readNumCsv(fullfile(P.anaOut, 'c20.csv'));
[Ht, Mt] = readNumCsv(fullfile(P.anaOut, 'checkup_time.csv'));
[He, Me] = readNumCsv(P.ecmCsv);
cC = getCol(Hc,Mc,'cell'); aC = getCol(Hc,Mc,'ah'); cap = getCol(Hc,Mc,'cap');
cT = getCol(Ht,Mt,'cell'); aT = getCol(Ht,Mt,'ah'); hT = getCol(Ht,Mt,'hours');
cE = getCol(He,Me,'cell'); aE = getCol(He,Me,'ah');
RE = getCol(He,Me,'Rs') + getCol(He,Me,'R1') + getCol(He,Me,'R2');
n2E = getCol(He,Me,'n2'); okN = getCol(He,Me,'spreadN2') <= P.spreadN2Max;
[tf, loc] = ismember([cC aC], [cT aT], 'rows');
hours = nan(size(cC)); hours(tf) = hT(loc(tf));
TEMPS = 'RCH'; QL = 1 - P.TH;
D = struct('cell',[],'temp',[],'life',[],'te',[],'Qe',[],'x',[],'z',[]);
for c = unique(cC)'
    i = find(cC == c); soh = cap(i) / cap(i(1)); h = hours(i);
    keep = ~isnan(h); soh = soh(keep); h = h(keep); a = aC(i(keep));
    j = find(soh < P.TH, 1);                                  % 수명: 처음 0.85 아래로 내려간 점과 그 전 점 사이 보간
    if isempty(j) || j == 1, continue; end
    life = h(j-1) + (P.TH - soh(j-1)) * (h(j) - h(j-1)) / (soh(j) - soh(j-1));
    k = find(soh <= sohE, 1); if isempty(k), continue; end
    ie = find(cE == c); [ag, ~, gi] = unique(aE(ie));         % ECM: (cell, ah)별 SOC 평균
    Rg = accumarray(gi, RE(ie), [], @mean);
    nOK = accumarray(gi, okN(ie)); nS = accumarray(gi, n2E(ie) .* okN(ie)); n2g = nS ./ nOK;   % 흔들림 조건 만족분 평균
    b = find(ag < 1, 1); m = find(abs(ag - a(k)) < 1, 1);
    if isempty(b) || isempty(m), continue; end
    te = h(k);
    D.cell(end+1,1) = c; D.temp(end+1,1) = TEMPS(mod(c-1,3)+1); D.life(end+1,1) = life;
    D.te(end+1,1) = te; D.Qe(end+1,1) = 1 - soh(k);
    D.x(end+1,1) = (Rg(m) / Rg(b) - 1) / sqrt(te);
    D.z(end+1,1) = (n2g(m) - n2g(b)) / sqrt(te);
end
D.hot = double(D.temp == 'H');
end
