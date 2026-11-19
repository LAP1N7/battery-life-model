function pred = locoGrid(F, D, QL)
% locoGrid.m : LOCO(셀 하나 빼고 학습 -> 뺀 셀 예측). 학습 = 후보 격자 중 ln(수명) RMSE 최소
%   입력  F : [후보 수 x 셀 수] 의 f 행렬 (후보마다 셀별 f),  D: buildLifeTable 결과,  QL: 0.15
%   출력  pred : 셀별 LOCO 예측 수명 [h]
n = numel(D.life); pred = nan(n, 1);
te = D.te'; Qe = D.Qe'; L = log(D.life');
T = log(lifeModel(F, te, Qe, QL));                % 후보 x 셀
for i = 1:n
    tr = true(1, n); tr(i) = false;
    e = sqrt(mean((T(:, tr) - L(tr)).^2, 2));
    [~, j] = min(e);
    pred(i) = exp(T(j, i));
end
end
