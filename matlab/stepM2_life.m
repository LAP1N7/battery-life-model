%% stepM2_life.m : 수명 수식 모델 (MATLAB) — ① P1 재현(교차검증) ② 기준점 민감도 ③ 확장 항 시험 ④ knee 진단
% 위치: ...\Mothat\matlab\   실행: MATLAB 명령창에서  stepM2_life
% 필요: MATLAB 기본 기능만 (툴박스 없음). 함수 파일: lifeParams, lifeModel, buildLifeTable, locoGrid, fGrid, readNumCsv, getCol, ternary, subD
% 입력: ..\analysis\out\c20.csv, checkup_time.csv, P1_pred.csv(대조용),  out\ecm_matlab.csv (stepM1_ecm.m)
% 출력: out\M2_life.csv (셀별 예측), out\M2_summary.txt (콘솔과 같은 요약)
%
% [합의된 사양 2026-11-19]
%  수식: Q(t) = A*sqrt(t) + B*t,  A = f*Qe/sqrt(te), B = (1-f)*Qe/te,  수명 = Q = 0.15 되는 t
%  ① 재현: f = 상수(용량만), f = clip(a + b*x)(EIS) 를 Python P1과 같은 격자·LOCO로 -> ln 예측 차이 0.01 이하면 통과
%  ② 민감도: 초기 시점 기준 SOH 0.98 / 0.975 / 0.96
%  ③ 확장 (결과 보기 전 등록): E1 f에 초기 n2 변화율 z 추가,  E2 f에 고온 항 추가 (25/45 °C 셀만, -5 °C 제외)
%     채택 기준: LOCO 오차가 기준 모델보다 10% 이상(상대) 작고 AND 섞은 값 시험 p < 0.05
%     (E2 참고: 온도가 2개뿐이면 Arrhenius 계수 exp(-Ea/R*(1/T - 1/Tref)) 는 '고온 여부' 하나의 배수와 같아짐)
%  ④ knee 진단(예측에 안 씀): 전체 궤적에 Q = A*sqrt(t) + B*t + K*max(0, t - tk) 를 A,B,K >= 0 로 피팅 (lsqnonneg, tk는 격자)
% [사양과 다른 점] 사양표의 lsqnonlin 대신: ①~③은 Python과 같은 격자 탐색(재현이 목적), ④는 선형 부분이라 lsqnonneg로 정확히 풂
%                  -> Optimization Toolbox가 필요 없음
% [가정·한계] 21셀 중 수명 도달·초기 EIS가 있는 셀만 (P1과 같은 17셀 예상), 고온 셀 SOH 기준점 artifact 남음
% 검증: Claude가 같은 입력으로 Octave에서 실행해 Python 결과와 대조함 (MATLAB R2026b 실행은 사용자 확인 필요)
% 변경 이력: 2026-11-19 최초 작성 (Claude)
%            2026-11-19 ④ knee 진단에서 P.excludeCk(셀 13, 2682 Ah, 36.6 °C 측정) 제외

P = lifeParams(); QL = 1 - P.TH;
if ~exist(P.outDir, 'dir'), mkdir(P.outDir); end
lerr = @(p, life) sqrt(mean((log(p(:)) - log(life(:))).^2));
rng(P.seed);
out = {};   % 요약 문장 모음
say = @(varargin) fprintf([varargin{1} '\n'], varargin{2:end});

%% ① P1 재현
D = buildLifeTable(P, P.SOH_E); n = numel(D.life);
p1 = locoGrid(fGrid(P, D, 'cap'), D, QL);
p3 = locoGrid(fGrid(P, D, 'ab', D.x), D, QL);
e1 = lerr(p1, D.life); e3 = lerr(p3, D.life);
say('[① 재현] 셀 %d개 | LOCO 오차(ln RMSE): 용량만 %.3f, EIS(a + b*x) %.3f   (Python P1: 0.388, 0.168)', n, e1, e3);
[Hp, Mp] = readNumCsv(fullfile(P.anaOut, 'P1_pred.csv'));
cp = getCol(Hp, Mp, 'cell'); [~, ip] = ismember(D.cell, cp);
d1 = max(abs(log(p1) - log(Mp(ip, find(strcmp(Hp, 'pred ①'), 1)))));
d3 = max(abs(log(p3) - log(Mp(ip, find(strcmp(Hp, 'pred ③'), 1)))));
ok1 = d1 <= 0.01 && d3 <= 0.01;
say('  Python 셀별 예측과의 최대 차이 (ln): 용량만 %.2g, EIS %.2g  -> %s', d1, d3, ternary(ok1, '통과 (기준 0.01 이하)', '불일치 - 원인 확인 필요'));
% 섞은 값 시험 (EIS 모델): x를 셀끼리 섞어 같은 LOCO
k = 0;
for r = 1:P.nPerm
    Dp = D; Dp.x = D.x(randperm(n));
    k = k + (lerr(locoGrid(fGrid(P, Dp, 'ab', Dp.x), Dp, QL), Dp.life) <= e3);
end
say('  섞은 값 시험 %d회: 실제 이하 %d회 -> p = %.2f', P.nPerm, k, (k + 1) / (P.nPerm + 1));

%% ② 기준점 민감도
say('\n[② 민감도] 초기 시점 기준 SOH | 셀 수 | 측정 시점/수명 중앙값 | 용량만 | EIS');
for s = P.SOH_E_SENS
    Ds = buildLifeTable(P, s);
    a1 = lerr(locoGrid(fGrid(P, Ds, 'cap'), Ds, QL), Ds.life);
    a3 = lerr(locoGrid(fGrid(P, Ds, 'ab', Ds.x), Ds, QL), Ds.life);
    say('  %.3f | %2d | %.2f | %.3f | %.3f', s, numel(Ds.life), median(Ds.te ./ Ds.life), a1, a3);
end

%% ③ 확장 항
say('\n[③ 확장 항] 채택 = LOCO 오차 %.0f%% 이상 감소 AND p < %.2f', P.accept * 100, P.pMax);
ext = {'E1 초기 n2 변화율 z', 'z', true(n,1) & ~isnan(D.z); ...
       'E2 고온 항 (25/45 °C 셀만)', 'hot', D.temp ~= 'C'};
for q = 1:size(ext, 1)
    m = ext{q, 3}; Dq = subD(D, m); v2 = Dq.(ext{q, 2});
    eb = lerr(locoGrid(fGrid(P, Dq, 'ab', Dq.x), Dq, QL), Dq.life);
    ee = lerr(locoGrid(fGrid(P, Dq, 'abc', Dq.x, v2), Dq, QL), Dq.life);
    gain = 1 - ee / eb; kk = 0;
    if gain >= P.accept          % 개선이 기준을 넘을 때만 섞은 값 시험 (시간 절약)
        for r = 1:P.nPerm
            vp = v2(randperm(numel(v2)));
            kk = kk + (lerr(locoGrid(fGrid(P, Dq, 'abc', Dq.x, vp), Dq, QL), Dq.life) <= ee);
        end
        pq = (kk + 1) / (P.nPerm + 1);
    else
        pq = NaN;
    end
    acc = gain >= P.accept && pq < P.pMax;
    say('  %s: 셀 %d개 | 기준 %.3f -> 확장 %.3f (변화 %+.0f%%) | p = %s -> %s', ext{q,1}, numel(Dq.life), eb, ee, ...
        -gain * 100, num2str(pq, '%.2f'), ternary(acc, '채택', '탈락 (단순 모델 유지)'));
end

%% ④ knee 진단 (전체 궤적, 예측에 사용 안 함)
[Hc, Mc] = readNumCsv(fullfile(P.anaOut, 'c20.csv')); [Ht, Mt] = readNumCsv(fullfile(P.anaOut, 'checkup_time.csv'));
cC = getCol(Hc,Mc,'cell'); aC = getCol(Hc,Mc,'ah'); cap = getCol(Hc,Mc,'cap');
[tf, loc] = ismember([cC aC], [getCol(Ht,Mt,'cell') getCol(Ht,Mt,'ah')], 'rows');
hrs = nan(size(cC)); hT = getCol(Ht,Mt,'hours'); hrs(tf) = hT(loc(tf));
say('\n[④ knee 진단] 전체 궤적 R^2: 2항 -> 2항 + knee  (R^2 %.3f 이상 오르고 K > 0 이면 knee 신호, 제외 check-up %d개)', P.kneeGain, size(P.excludeCk,1));
K = [];
for c = D.cell'
    bad = false(size(cC));
    for e = 1:size(P.excludeCk, 1), bad = bad | (cC == P.excludeCk(e,1) & abs(aC - P.excludeCk(e,2)) < 1); end
    i = find(cC == c & ~isnan(hrs) & ~bad); t = hrs(i); y = 1 - cap(i) / cap(find(cC == c, 1));
    X2 = [sqrt(t) t]; b2 = lsqnonneg(X2, y); r2a = 1 - sum((y - X2*b2).^2) / sum((y - mean(y)).^2);
    best = [r2a NaN 0];
    for tk = linspace(0.2, 0.9, 15) * max(t)
        X3 = [X2 max(0, t - tk)]; b3 = lsqnonneg(X3, y);
        r2 = 1 - sum((y - X3*b3).^2) / sum((y - mean(y)).^2);
        if r2 > best(1), best = [r2 tk b3(3)]; end
    end
    K(end+1, :) = [c r2a best]; %#ok<AGROW>
end
sig = (K(:,3) - K(:,2) >= P.kneeGain) & K(:,5) > 0;
say('  knee 신호 셀 %d / %d: %s', sum(sig), size(K,1), mat2str(K(sig,1)'));
say('  2항 R^2 중앙값 %.3f -> knee 포함 %.3f', median(K(:,2)), median(K(:,3)));

%% 저장
T = [D.cell, D.te, D.Qe, D.x, D.z, D.life, p1, p3];
fid = fopen(fullfile(P.outDir, 'M2_life.csv'), 'w');
fprintf(fid, 'cell,te_h,Qe,x,z,life_h,pred_cap,pred_eis\n'); fprintf(fid, '%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', T'); fclose(fid);
say('\n저장: %s', fullfile(P.outDir, 'M2_life.csv'));
