function v = getCol(H, M, name)
% getCol.m : readNumCsv 결과에서 열 이름으로 열 꺼내기
k = find(strcmp(H, name), 1); assert(~isempty(k), ['열 없음: ' name]); v = M(:, k);
end
