function out = fitECM_mohtat(f, re, mim, opt)
%FITECM  등가회로(ECM) fitting + R1/R2 식별성 점검
%   회로: Rs + jwL + R1||CPE1 + R2||CPE2 + 확산 CPE(Qw, nw)
%   R||CPE는 시간상수 tau로 표현:  R / (1 + (jw*tau)^n)   (Q = tau^n / R 와 같은 식)
%   tau는 측정 주파수 창 안으로 제한: 1/(2*pi*f_max) <= tau <= 1/(2*pi*f_min)  (11-18 수정)
%     이유: 제한이 없으면 R1의 tau가 측정 범위(10 kHz)보다 빨라져 R1이 Rs를 흡수함
%           (전체 실행에서 1,836개 중 417개가 Rs = 0, R1 = Rs + R1 이 됨. Python 결과도 412개 같은 문제)
%   out = fitECM(f, re, mim)
%   out = fitECM(f, re, mim, opt)   opt.ratioMin (기본 20), opt.spreadMax (기본 0.20), opt.costTol (기본 0.05)
%
%   입력  f [Hz], re = Re(Z) [Ohm], mim = -Im(Z) [Ohm]  (열 벡터)
%   출력  out 구조체
%     .p        최적 파라미터 [Rs, L, R1, log10tau1, n1, R2, log10tau2, n2, log10Qw, nw] (정렬 후)
%     .Rs, .R1, .R2, .Rint   [Ohm]   (Rint = R1 + R2)
%     .tau1, .tau2           [s]     정렬 규칙: tau1 < tau2 (R1 = 빠른 과정)
%     .atEdge                tau1 또는 tau2가 측정 창 경계에 붙음 (그 RC는 측정 범위 밖 -> 분리 불확실)
%     .ratio                 tau2 / tau1
%     .spreadR1, .spreadR2   거의 같은 품질로 맞는 해들 사이의 R1, R2 흔들림 = (최대 - 최소) / 중앙값
%     .nGood                 그런 해의 개수
%     .separable             ratio >= ratioMin, spreadR1·spreadR2 <= spreadMax, 그리고 atEdge가 아님
%     .relErr                |Z|로 나눈 상대 잔차의 RMS (Python fit_ecm의 err와 같은 정의)
%     .n1, .n2, .nw          CPE 이상성 지수 (1 = 이상적 커패시터, 1보다 작을수록 반원이 눌림)  (11-19 추가, P5용)
%     .spreadN2              거의 같은 품질로 맞는 해들 사이 n2 흔들림 = (최대 - 최소) / 중앙값
%     .n2AtBound             n2가 피팅 범위 끝(0.5 또는 1)에서 0.01 이내 -> 경계에 붙음
%   변경 이력: 11-18 tau 측정 창 제한·파라미터 스케일링 / 11-19 n1·n2·nw·spreadN2·n2AtBound 출력 추가 (피팅 자체는 그대로)
%
%   안전장치 (11-18 합의)
%     ① 순서 고정: 해마다 tau로 정렬해 R1 = 빠른 과정, R2 = 느린 과정으로 이름을 고정 (라벨 바뀜 방지)
%     ② 분리도:    tau2/tau1 < ratioMin(20) 이면 두 반원이 겹쳐 분리 불확실
%     ③ 흔들림:    시작점 6개 중 최저 비용의 (1 + costTol) 이내인 해들에서 R1, R2가 spreadMax(20%) 넘게 달라지면 불확실
%   분리 불확실한 스펙트럼은 해석에 R1, R2 대신 Rint만 쓸 것
%   필요: Optimization Toolbox (lsqnonlin)

if nargin < 4, opt = struct(); end
if ~isfield(opt, 'ratioMin'),  opt.ratioMin  = 20;   end
if ~isfield(opt, 'spreadMax'), opt.spreadMax = 0.20; end
if ~isfield(opt, 'costTol'),   opt.costTol   = 0.05; end

f = f(:); re = re(:); mim = mim(:);
ltLo = log10(1 / (2*pi*max(f)));   ltHi = log10(1 / (2*pi*min(f)));   % 측정 창 [log10 s]
LO = [0,    0,    0,    ltLo, 0.5, 0,    ltLo, 0.5, -2, 0.3];
HI = [0.05, 1e-5, 0.05, ltHi, 1,   0.05, ltHi, 1,    5, 1  ];

w = 2*pi*f;
Z = re - 1i*mim;
% 파라미터 크기 맞추기: 저항은 mOhm, 인덕턴스는 uH 단위로 바꿔서 최적화 (11-18 수정)
%   이유: 저항이 1e-3 Ohm, L이 1e-7 H 수준이라 기본 허용오차로는 최적화가 일찍 멈춤
%         (첫 실행에서 Python보다 fit 오차가 큰 경우 766/1836개 발생)
S = [1e-3, 1e-6, 1e-3, 1, 1, 1e-3, 1, 1, 1, 1];
resid = @(x) residual(x .* S, w, Z);

r0 = min(re);
% 시작점 6개: 앞의 2개는 Python fit_ecm과 동일, 나머지는 R과 Q를 바꿔 다른 해를 찾기 위함
% 시작점 6개 (log10 tau [s]): 빠른/느린 조합을 바꿔 다른 해를 찾기 위함
inits = [r0, 1e-7, 0.001,  -3,   0.8, 0.002,  -1,   0.8, 2, 0.6;
         r0, 1e-7, 0.002,  -1,   0.8, 0.001,  -3,   0.8, 2, 0.6;
         r0, 1e-7, 0.001,  -4,   0.8, 0.001,   0,   0.8, 2, 0.6;
         r0, 1e-7, 0.001,   0,   0.8, 0.001,  -4,   0.8, 2, 0.6;
         r0, 1e-7, 0.003,  -2,   0.9, 0.0005, -0.5, 0.7, 2, 0.6;
         r0, 1e-7, 0.0005, -0.5, 0.7, 0.003,  -2,   0.9, 2, 0.6];

opts = optimoptions('lsqnonlin', 'Display', 'off', 'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000, ...
                    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10);
nS = size(inits, 1);
sols = zeros(nS, 10); costs = inf(nS, 1);
for k = 1:nS
    q0 = min(max(inits(k,:), LO + 1e-9), HI - 1e-9);   % 시작점을 경계 안으로
    [x, resnorm] = lsqnonlin(resid, q0 ./ S, LO ./ S, HI ./ S, opts);
    sols(k,:) = sortByTau(x .* S); costs(k) = resnorm;
end

[best, ib] = min(costs);
p = sols(ib,:);
good = costs <= best * (1 + opt.costTol);               % 거의 같은 품질로 맞는 해들
R1s = sols(good, 3); R2s = sols(good, 6); N2s = sols(good, 8);
spread = @(x) (max(x) - min(x)) / max(median(x), eps);

[t1, t2] = taus(p);
out.p = p;
out.Rs = p(1); out.R1 = p(3); out.R2 = p(6); out.Rint = p(3) + p(6);
out.tau1 = t1; out.tau2 = t2; out.ratio = t2 / max(t1, eps);
out.spreadR1 = spread(R1s); out.spreadR2 = spread(R2s); out.nGood = sum(good);
edgeTol = 0.05;   % log10 기준으로 경계에서 0.05 이내면 '붙었다'고 봄
out.atEdge = any(abs(p([4 7]) - ltLo) < edgeTol) || any(abs(p([4 7]) - ltHi) < edgeTol);
out.separable = out.ratio >= opt.ratioMin && out.spreadR1 <= opt.spreadMax && out.spreadR2 <= opt.spreadMax && ~out.atEdge;
out.relErr = sqrt(best / numel(f));
out.n1 = p(5); out.n2 = p(8); out.nw = p(10);                       % 11-19 추가 (P5)
out.spreadN2 = spread(N2s);
out.n2AtBound = abs(p(8) - LO(8)) < 0.01 || abs(p(8) - HI(8)) < 0.01;
end

function q = sortByTau(q)
% ① 순서 고정: tau가 작은(빠른) RC를 R1 자리에 둠
[t1, t2] = taus(q);
if t1 > t2
    q([3 4 5 6 7 8]) = q([6 7 8 3 4 5]);
end
end

function [t1, t2] = taus(q)
t1 = 10^q(4);   % R1||CPE1 의 시간상수 [s]
t2 = 10^q(7);   % R2||CPE2 의 시간상수 [s]
end

function r = residual(q, w, Z)
d = (ecmZ(q, w) - Z) ./ abs(Z);
r = [real(d); imag(d)];
end

function Zm = ecmZ(q, w)
Rs = q(1); L = q(2); R1 = q(3); t1 = 10^q(4); n1 = q(5);
R2 = q(6); t2 = 10^q(7); n2 = q(8); Qw = 10^q(9); nw = q(10);
jw = 1i*w;
Zm = Rs + jw*L + R1 ./ (1 + (jw*t1).^n1) + R2 ./ (1 + (jw*t2).^n2) + 1 ./ (Qw*jw.^nw);
end
