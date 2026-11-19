function t = lifeModel(f, te, Qe, QL)
% lifeModel.m : 2항 수식 Q(t) = A*sqrt(t) + B*t 에서 Q = QL 이 되는 시간 (수명)
%   입력  f  : 초기 손실 중 sqrt(t)형 비율 [0~1]   (행렬 가능: 후보 x 셀)
%         te : 초기 측정 시점 [h] (행 벡터, 셀별)
%         Qe : 초기 측정 시점 손실 1-SOH [-] (행 벡터)
%         QL : 수명 기준 손실 (0.15 = SOH 0.85)
%   출력  t  : 수명 [h]  (f와 같은 크기)
%   가정  A = f*Qe/sqrt(te), B = (1-f)*Qe/te  (초기 1점 통과 조건),  A, B >= 0
% 변경 이력: 2026-11-19 최초 작성 (Python stepP1_prototype.py 의 t85 와 같은 식)
A  = f .* Qe ./ sqrt(te);
B  = (1 - f) .* Qe ./ te;
Bs = max(B, 1e-15);
s  = (-A + sqrt(A.^2 + 4 .* Bs .* QL)) ./ (2 .* Bs);
lin = B <= 1e-12;                       % B = 0 이면 sqrt(t)항만:  s = QL / A
sA = QL ./ max(A, 1e-15);
s(lin) = sA(lin);
t = s.^2;
end
