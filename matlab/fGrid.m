function F = fGrid(P, D, kind, v1, v2)
% fGrid.m : 가중치 후보 격자 -> f 행렬 [후보 x 셀]
%   kind 'cap' : f = 상수                       (① 용량만)
%        'ab'  : f = clip(a + b*v1)             (③ EIS, v1 = x)
%        'abc' : f = clip(a + b*v1 + c*v2)      (확장: v2 = z 또는 hot)
%   b, c 후보는 각 변수의 최대 절댓값으로 나눠 씀 (Python과 동일한 격자)
switch kind
    case 'cap'
        F = repmat(P.gridF(:), 1, numel(D.life));
    case 'ab'
        s = 1 / max(max(abs(v1)), 1e-12);
        [b, a] = ndgrid(P.gridB * s, P.gridA);           % 순서: a 바깥, b 안쪽 (Python과 같은 동점 처리)
        F = min(max(a(:) + b(:) * v1(:)', 0), 1);
    case 'abc'
        s = 1 / max(max(abs(v1)), 1e-12); s2 = 1 / max(max(abs(v2)), 1e-12);
        [c, b, a] = ndgrid(P.grid3{3} * s2, P.grid3{2} * s, P.grid3{1});
        F = min(max(a(:) + b(:) * v1(:)' + c(:) * v2(:)', 0), 1);
end
end
