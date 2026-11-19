function [H, M] = readNumCsv(path)
% readNumCsv.m : 숫자 CSV 읽기 (첫 줄 = 열 이름). 툴박스·table 없이 동작 (MATLAB/Octave 공통)
%   출력 H: 열 이름 cell, M: 숫자 행렬 (숫자가 아닌 칸은 NaN, 예: temp 열의 R/C/H 문자는 NaN)
% 변경 이력: 2026-11-19 최초 작성
fid = fopen(path, 'r', 'n', 'UTF-8'); assert(fid > 0, ['파일 없음: ' path]);
hdr = fgetl(fid); if double(hdr(1)) == 65279, hdr = hdr(2:end); end   % BOM 제거
H = strtrim(strsplit(hdr, ','));
C = textscan(fid, repmat('%s', 1, numel(H)), 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);
C = C{1}; M = str2double(C);
end
