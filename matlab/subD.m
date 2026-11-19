function Dq = subD(D, m)
% subD.m : buildLifeTable 구조체에서 m(논리 벡터)에 해당하는 셀만 남김
f = fieldnames(D); Dq = D;
for i = 1:numel(f), Dq.(f{i}) = D.(f{i})(m); end
end
