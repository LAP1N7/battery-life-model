%% 0 설정
clear; clc; close all;

t_start = tic;

% 입력 경로 및 출력 폴더 설정
JOINED  = 'mohtat_modes_joined.csv';
RAW     = '..\..\..\..\Mothat\data';   % temp_mech로 이동하며 경로 수정
OUT_DIR = 'out';
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

% 로그 파일 열기
log_path = fullfile(OUT_DIR, 'log.txt');
log_fid = fopen(log_path, 'w');
if log_fid < 0
    error('로그 파일을 열 수 없습니다: %s', log_path);
end
cleanup_obj = onCleanup(@() safe_fclose(log_fid));

log_msg(log_fid, '=== Mohtat 2020 열화 모드 / ICA / 저항 통합 분석 시작 ===\n');
log_msg(log_fid, '시작 시각: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
log_msg(log_fid, 'JOINED: %s\n', JOINED);
log_msg(log_fid, 'RAW   : %s\n\n', RAW);

% 온도 그룹, 색상, 모드, ICA 창 정의
T_ORDER = ["C", "R", "H"];
T_LABEL = ["C (-5 °C)", "R (25 °C)", "H (45 °C)"];
COL_MAP = struct('C', [0 0.45 0.70], 'R', [0 0.62 0.45], 'H', [0.84 0.37 0]);
MODES   = ["LLI", "LAM_PE", "LAM_NE"];
WINS    = [3.45 3.56; 3.56 3.74; 3.76 3.90];   % 실제 피크 위치(약 3.50, 3.65, 3.82 V)에 맞춤
WIN_NAMES = ["W1", "W2", "W3"];

% ICA 전압 격자 및 가우시안 스무딩 파라미터
V_GRID  = (3.45:0.002:4.18)';
DV_GRID = 0.002;
SIGMA_G = 4;

% JOINED 테이블 읽기
if ~exist(JOINED, 'file')
    log_msg(log_fid, '[오류] JOINED 파일을 찾을 수 없습니다: %s\n', JOINED);
    error('JOINED 파일 없음: %s', JOINED);
end
J = readtable(JOINED, 'TextType', 'string');
% UTF-8 BOM 등으로 첫 컬럼명에 특수문자가 섞인 경우 정리
vnames = J.Properties.VariableNames;
vnames{1} = regexprep(vnames{1}, '^[^A-Za-z0-9]+', '');
J.Properties.VariableNames = vnames;

J.T      = strtrim(string(J.T));
J.cell   = to_double_vec(J.cell);
J.T_C    = to_double_vec(J.T_C);
J.ah     = to_double_vec(J.ah);
J.days   = to_double_vec(J.days);
J.efc    = to_double_vec(J.efc);
J.loss   = to_double_vec(J.loss);
J.LLI    = to_double_vec(J.LLI);
J.LAM_PE = to_double_vec(J.LAM_PE);
J.LAM_NE = to_double_vec(J.LAM_NE);
J.dR0    = to_double_vec(J.dR0);
J.dRarc  = to_double_vec(J.dRarc);

log_msg(log_fid, 'JOINED 로드 완료: %d 행, %d 셀\n\n', height(J), numel(unique(J.cell)));


%% 1 ICA
log_msg(log_fid, '--- [1] ICA 분석 시작 ---\n');

all_cells = sort(unique(J.cell));
ica_cell_data = struct([]);

for ci = 1:numel(all_cells)
    c = all_cells(ci);
    fpath = fullfile(RAW, sprintf('%02d', c), 'OCV_wExpansion.csv');

    entry = struct();
    entry.cell = c;
    entry.ok = false;
    entry.checkups = struct([]);

    try
        if ~exist(fpath, 'file')
            error('파일이 존재하지 않습니다: %s', fpath);
        end
        log_msg(log_fid, 'Cell %02d 읽는 중: %s\n', c, fpath);
        raw = readmatrix(fpath, 'NumHeaderLines', 1);
        if isempty(raw) || size(raw, 2) < 9
            error('컬럼 수 부족 또는 빈 파일 (size = %dx%d)', size(raw, 1), size(raw, 2));
        end
        % 필요한 컬럼만 유지: 2=current_mA, 3=voltage_V, 7=capacity_Ah, 9=total_Ah_throughput
        raw = raw(:, [2, 3, 7, 9]);
        valid_rows = all(isfinite(raw), 2);
        raw = raw(valid_rows, :);
        if isempty(raw)
            error('유효한 수치 행이 없습니다.');
        end

        % total_Ah_throughput 기준 체크업 그룹화 (상대 허용오차 1e-6)
        [uAh, ~, grpIdx] = uniquetol(raw(:, 4), 1e-6);
        n_grp = numel(uAh);

        chks = struct([]);
        k_valid = 0;

        for g = 1:n_grp
            mask_g = (grpIdx == g);
            I_mA = raw(mask_g, 1);
            V_v  = raw(mask_g, 2);
            Q_ah = raw(mask_g, 3);
            ah_g = mean(raw(mask_g, 4));

            % 충전 구간: 100 < current_mA < 400, 최소 500 행 이상
            chg_mask = (I_mA > 100) & (I_mA < 400);
            if sum(chg_mask) < 500
                continue;
            end

            Vc = V_v(chg_mask);
            Qc = Q_ah(chg_mask);

            % Q 기준 정렬
            [Qs, ord] = sort(Qc);
            Vs = Vc(ord);
            Q_total = Qs(end) - Qs(1);

            % 비증가 V 잡음 제거: 단조 포락선(cummax) 유지 후 엄격 증가점 추출
            V_env = cummax(Vs);
            keep_mono = [true; diff(V_env) > 0];
            Vu = V_env(keep_mono);
            Qu = Qs(keep_mono);
            if numel(Vu) < 10
                continue;
            end

            % 균일 V 격자 보간, 중앙 차분 dQ/dV, 가우시안 스무딩 (sigma = 4)
            Q_grid = interp1(Vu, Qu, V_GRID, 'linear', NaN);
            dqdv_smooth = compute_dqdv_smooth(Q_grid, DV_GRID, SIGMA_G);
            if sum(isfinite(dqdv_smooth)) < 20
                continue;
            end

            k_valid = k_valid + 1;
            chks(k_valid).ah = ah_g;
            chks(k_valid).Q_total = Q_total;
            chks(k_valid).dqdv = dqdv_smooth;
        end

        if k_valid == 0
            error('유효한 C/20 충전 체크업이 없습니다.');
        end

        % 기준 체크업: JOINED에서 해당 셀의 최소 ah에 가장 가까운 total_Ah
        j_sub = J(J.cell == c, :);
        min_ah_j = min(j_sub.ah);
        chk_ahs = [chks.ah]';
        [~, ref_idx] = min(abs(chk_ahs - min_ah_j));

        % 기준 곡선에서 prominence >= 5% of max 인 국소 최대값 탐색 (창당 최고 1개)
        dqdv_ref = chks(ref_idx).dqdv;
        max_ref = max(dqdv_ref, [], 'omitnan');
        [p_idx, p_prom] = find_local_maxima_prom(dqdv_ref);
        prom_ok = p_prom >= 0.05 * max_ref;
        p_idx = p_idx(prom_ok);

        ref_h = nan(1, 3);
        ref_pos = nan(1, 3);
        ref_area = nan(1, 3);

        for w = 1:3
            w_mask = (V_GRID >= WINS(w, 1)) & (V_GRID <= WINS(w, 2));
            Vw = V_GRID(w_mask);
            yw = dqdv_ref(w_mask);
            fin_w = isfinite(yw);
            if sum(fin_w) >= 2
                ref_area(w) = trapz(Vw(fin_w), yw(fin_w));
            end

            % 해당 창 내 prominence 통과 피크 중 최고점
            p_in_w = p_idx(V_GRID(p_idx) >= WINS(w, 1) & V_GRID(p_idx) <= WINS(w, 2));
            if ~isempty(p_in_w)
                [pk_val, best_k] = max(dqdv_ref(p_in_w));
                best_idx = p_in_w(best_k);
                ref_h(w) = pk_val;
                ref_pos(w) = V_GRID(best_idx);
            else
                % 예외 처리: 창 내 prominence 기준 피크가 없으면 창 내 최대값 사용
                if any(fin_w)
                    [pk_val, local_k] = max(yw, [], 'omitnan');
                    ref_h(w) = pk_val;
                    ref_pos(w) = Vw(local_k);
                    log_msg(log_fid, '  [주의] Cell %02d 기준곡선 %s 창 내 prominence>=5%% 피크 없음 -> 창 최대값 사용\n', c, WIN_NAMES(w));
                end
            end
        end

        % 모든 체크업에 대해 창별 특성 및 정규화 값 계산
        for k = 1:k_valid
            yk = chks(k).dqdv;
            h_k = nan(1, 3);
            pos_k = nan(1, 3);
            area_k = nan(1, 3);

            for w = 1:3
                w_mask = (V_GRID >= WINS(w, 1)) & (V_GRID <= WINS(w, 2));
                Vw = V_GRID(w_mask);
                yw = yk(w_mask);
                fin_w = isfinite(yw);

                if sum(fin_w) >= 2
                    area_k(w) = trapz(Vw(fin_w), yw(fin_w));
                end
                if k == ref_idx
                    h_k(w) = ref_h(w);
                    pos_k(w) = ref_pos(w);
                elseif any(fin_w)
                    [pk_val, local_k] = max(yw, [], 'omitnan');
                    h_k(w) = pk_val;
                    pos_k(w) = Vw(local_k);
                end
            end

            chks(k).is_ref = (k == ref_idx);
            chks(k).height = h_k;
            chks(k).pos = pos_k;
            chks(k).area = area_k;
            chks(k).height_rel = h_k ./ ref_h;
            chks(k).area_rel = area_k ./ ref_area;
            chks(k).shift_mV = 1000 * (pos_k - ref_pos);
        end

        entry.ok = true;
        entry.ref_idx = ref_idx;
        entry.checkups = chks;
        log_msg(log_fid, '  Cell %02d 완료: 유효 체크업 %d/%d, ref_ah = %.2f (JOINED min = %.2f)\n', ...
            c, k_valid, n_grp, chks(ref_idx).ah, min_ah_j);

    catch ME
        log_msg(log_fid, '  [경고] Cell %02d 처리 실패: %s\n', c, ME.message);
    end

    if isempty(ica_cell_data)
        ica_cell_data = entry;
    else
        ica_cell_data(end + 1) = entry; %#ok<SAGROW>
    end
end

% JOINED 테이블과 ICA 특성 조인 (cell 및 최근접 ah, |dAh| <= 5)
nJ = height(J);
ica_ah        = nan(nJ, 1);
ica_ah_diff   = nan(nJ, 1);
Q_total_col   = nan(nJ, 1);
h_mat         = nan(nJ, 3);
pos_mat       = nan(nJ, 3);
area_mat      = nan(nJ, 3);
hrel_mat      = nan(nJ, 3);
arel_mat      = nan(nJ, 3);
shift_mat     = nan(nJ, 3);
matched_k     = nan(nJ, 1);

cell_list_ica = [ica_cell_data.cell];
n_matched = 0;

for i = 1:nJ
    c = J.cell(i);
    idx_c = find(cell_list_ica == c, 1);
    if isempty(idx_c) || ~ica_cell_data(idx_c).ok
        continue;
    end
    chks = ica_cell_data(idx_c).checkups;
    ahs = [chks.ah]';
    [min_d, best_k] = min(abs(ahs - J.ah(i)));
    if min_d <= 5
        n_matched = n_matched + 1;
        ica_ah(i)      = chks(best_k).ah;
        ica_ah_diff(i) = chks(best_k).ah - J.ah(i);
        Q_total_col(i) = chks(best_k).Q_total;
        h_mat(i, :)     = chks(best_k).height;
        pos_mat(i, :)   = chks(best_k).pos;
        area_mat(i, :)  = chks(best_k).area;
        hrel_mat(i, :)  = chks(best_k).height_rel;
        arel_mat(i, :)  = chks(best_k).area_rel;
        shift_mat(i, :) = chks(best_k).shift_mV;
        matched_k(i)    = best_k;
    end
end

ICA_TBL = J;
ICA_TBL.ica_ah        = ica_ah;
ICA_TBL.ica_ah_diff   = ica_ah_diff;
ICA_TBL.Q_total       = Q_total_col;
ICA_TBL.W1_height     = h_mat(:, 1);
ICA_TBL.W1_pos        = pos_mat(:, 1);
ICA_TBL.W1_area       = area_mat(:, 1);
ICA_TBL.W1_height_rel = hrel_mat(:, 1);
ICA_TBL.W1_area_rel   = arel_mat(:, 1);
ICA_TBL.W1_shift_mV   = shift_mat(:, 1);
ICA_TBL.W2_height     = h_mat(:, 2);
ICA_TBL.W2_pos        = pos_mat(:, 2);
ICA_TBL.W2_area       = area_mat(:, 2);
ICA_TBL.W2_height_rel = hrel_mat(:, 2);
ICA_TBL.W2_area_rel   = arel_mat(:, 2);
ICA_TBL.W2_shift_mV   = shift_mat(:, 2);
ICA_TBL.W3_height     = h_mat(:, 3);
ICA_TBL.W3_pos        = pos_mat(:, 3);
ICA_TBL.W3_area       = area_mat(:, 3);
ICA_TBL.W3_height_rel = hrel_mat(:, 3);
ICA_TBL.W3_area_rel   = arel_mat(:, 3);
ICA_TBL.W3_shift_mV   = shift_mat(:, 3);

writetable(ICA_TBL, fullfile(OUT_DIR, 'ica_features.csv'));
log_msg(log_fid, 'ICA 조인 완료: %d / %d 행 매칭 (|dAh| <= 5). 저장 -> out\\ica_features.csv\n', n_matched, nJ);

% Spearman 및 loss 통제 편상관(partial Spearman) 계산
ICA_FEAT_NAMES = [ ...
    "W1_height_rel", "W1_area_rel", "W1_shift_mV", ...
    "W2_height_rel", "W2_area_rel", "W2_shift_mV", ...
    "W3_height_rel", "W3_area_rel", "W3_shift_mV"];
SCOPES = ["all", "C", "R", "H"];

n_corr_rows = numel(ICA_FEAT_NAMES) * numel(MODES) * numel(SCOPES);
c_feat  = strings(n_corr_rows, 1);
c_mode  = strings(n_corr_rows, 1);
c_scope = strings(n_corr_rows, 1);
c_rho   = nan(n_corr_rows, 1);
c_rhop  = nan(n_corr_rows, 1);
c_n     = zeros(n_corr_rows, 1);

row_i = 0;
for si = 1:numel(SCOPES)
    sc = SCOPES(si);
    if sc == "all"
        mask_sc = ICA_TBL.loss > 0;
    else
        mask_sc = (ICA_TBL.loss > 0) & (ICA_TBL.T == sc);
    end
    loss_sc = ICA_TBL.loss(mask_sc);

    for fi = 1:numel(ICA_FEAT_NAMES)
        fname = ICA_FEAT_NAMES(fi);
        x_sc = ICA_TBL.(fname)(mask_sc);
        for mi = 1:numel(MODES)
            mname = MODES(mi);
            y_sc = ICA_TBL.(mname)(mask_sc);

            valid = isfinite(x_sc) & isfinite(y_sc) & isfinite(loss_sc);
            n_v = sum(valid);

            row_i = row_i + 1;
            c_feat(row_i)  = fname;
            c_mode(row_i)  = mname;
            c_scope(row_i) = sc;
            c_n(row_i)     = n_v;
            if n_v >= 4
                c_rho(row_i)  = spearman_r(x_sc(valid), y_sc(valid));
                c_rhop(row_i) = partial_spearman_r(x_sc(valid), y_sc(valid), loss_sc(valid));
            end
        end
    end
end

ICA_MODES_TBL = table(c_feat, c_mode, c_scope, c_rho, c_rhop, c_n, ...
    'VariableNames', {'feature', 'mode', 'scope', 'rho', 'rho_partial', 'n'});
writetable(ICA_MODES_TBL, fullfile(OUT_DIR, 'ica_vs_modes.csv'));
log_msg(log_fid, 'ICA-모드 상관 분석 저장 -> out\\ica_vs_modes.csv\n');

% 전체(all) 범위 주요 결과 요약 출력
log_msg(log_fid, '\n[요약] ICA 특성 vs DVA 모드 Spearman 상관 (scope = all, loss > 0):\n');
log_msg(log_fid, '%-14s | %-6s | %7s | %11s | %4s\n', 'feature', 'mode', 'rho', 'rho_partial', 'n');
log_msg(log_fid, '%s\n', repmat('-', 1, 54));
sub_all = ICA_MODES_TBL(ICA_MODES_TBL.scope == "all", :);
for i = 1:height(sub_all)
    log_msg(log_fid, '%-14s | %-6s | %7.3f | %11.3f | %4d\n', ...
        sub_all.feature(i), sub_all.mode(i), sub_all.rho(i), sub_all.rho_partial(i), sub_all.n(i));
end
log_msg(log_fid, '\n');


%% 2 모드 진행식
log_msg(log_fid, '--- [2] 모드 진행식(Mode kinetics) 분석 시작 ---\n');

% [중요 주의사항]
% 일정 프로토콜 단일 셀 내부에서는 efc가 days에 비례하므로,
% a (sqrt-shape 항)와 b (linear-shape 항)는 시간 열화 vs 사이클 열화가 아니라
% 오직 곡선 형태(제곱근형 vs 선형)에 의해서만 구분됩니다.
log_msg(log_fid, ['[중요 참고] 단일 일정-프로토콜 셀에서는 efc가 days에 비례하므로, ' ...
    'a(sqrt-shape 항)와 b(linear-shape 항)는 시간 vs 사이클의 물리적 분리가 아니라 ' ...
    '오직 곡선 형태(sqrt vs linear)로만 구분됩니다.\n\n']);

n_kin_rows = numel(all_cells) * numel(MODES);
k_cell   = nan(n_kin_rows, 1);
k_T      = strings(n_kin_rows, 1);
k_mode   = strings(n_kin_rows, 1);
k_npts   = zeros(n_kin_rows, 1);
k_a      = nan(n_kin_rows, 1);
k_b      = nan(n_kin_rows, 1);
k_rmse   = nan(n_kin_rows, 1);
k_r2     = nan(n_kin_rows, 1);
k_share  = nan(n_kin_rows, 1);
k_z      = nan(n_kin_rows, 1);
k_Apow   = nan(n_kin_rows, 1);

% 그림용 피팅 정보 보관
kin_fit_store = struct();
row_i = 0;

for ci = 1:numel(all_cells)
    c = all_cells(ci);
    sub = J(J.cell == c, :);
    [~, ord_d] = sort(sub.days);
    sub = sub(ord_d, :);
    T_c = sub.T(1);

    for mi = 1:numel(MODES)
        mname = MODES(mi);
        y_all = sub.(mname);
        d_all = sub.days;
        e_all = sub.efc;

        valid = isfinite(y_all) & isfinite(d_all) & isfinite(e_all) & (d_all >= 0);
        d = d_all(valid);
        e = e_all(valid);
        y = y_all(valid);

        row_i = row_i + 1;
        k_cell(row_i) = c;
        k_T(row_i)    = T_c;
        k_mode(row_i) = mname;
        k_npts(row_i) = numel(y);

        if numel(y) >= 2
            X = [sqrt(d), e];
            coef = lsqnonneg(X, y);
            a_val = coef(1); % sqrt-shape 항
            b_val = coef(2); % linear-shape 항

            yhat = X * coef;
            rmse_val = sqrt(mean((y - yhat).^2));
            ss_res = sum((y - yhat).^2);
            ss_tot = sum((y - mean(y)).^2);
            if ss_tot > 0
                r2_val = 1 - ss_res / ss_tot;
            else
                r2_val = NaN;
            end

            % 마지막 체크업에서의 sqrt-shape 비중
            d_end = d(end);
            e_end = e(end);
            t_sqrt = a_val * sqrt(d_end);
            t_lin  = b_val * e_end;
            if (t_sqrt + t_lin) > 0
                share_val = t_sqrt / (t_sqrt + t_lin);
            else
                share_val = NaN;
            end

            % 형태 전용 대안 모델: m = A * days^z (m > 0.05 & days > 0 로그-선형 회귀)
            pmask = (y > 0.05) & (d > 0);
            if sum(pmask) >= 2
                lx = log(d(pmask));
                ly = log(y(pmask));
                Zmat = [ones(numel(lx), 1), lx];
                if rank(Zmat) == 2
                    pbeta = Zmat \ ly;
                    Apow_val = exp(pbeta(1));
                    z_val    = pbeta(2);
                else
                    Apow_val = NaN;
                    z_val    = NaN;
                end
            else
                Apow_val = NaN;
                z_val    = NaN;
            end

            k_a(row_i)     = a_val;
            k_b(row_i)     = b_val;
            k_rmse(row_i)  = rmse_val;
            k_r2(row_i)    = r2_val;
            k_share(row_i) = share_val;
            k_z(row_i)     = z_val;
            k_Apow(row_i)  = Apow_val;

            % 궤적 그림용 저장
            fld = sprintf('c%d_%s', c, mname);
            kin_fit_store.(fld) = struct('days', d, 'efc', e, 'y', y, 'a', a_val, 'b', b_val, 'T', T_c);
        end
    end
end

KIN_CELLS_TBL = table(k_cell, k_T, k_mode, k_npts, k_a, k_b, k_rmse, k_r2, k_share, k_z, k_Apow, ...
    'VariableNames', {'cell', 'T', 'mode', 'n_pts', 'a_sqrt_shape', 'b_linear_shape', ...
    'rmse', 'r2', 'sqrt_share', 'z', 'A_pow'});
writetable(KIN_CELLS_TBL, fullfile(OUT_DIR, 'mode_kinetics_cells.csv'));

% 온도 x 모드별 그룹 요약 (sqrt-share 및 z의 중앙값 및 범위, 셀 수)
n_grp_rows = numel(T_ORDER) * numel(MODES);
g_T         = strings(n_grp_rows, 1);
g_mode      = strings(n_grp_rows, 1);
g_ncells    = zeros(n_grp_rows, 1);
g_sh_med    = nan(n_grp_rows, 1);
g_sh_min    = nan(n_grp_rows, 1);
g_sh_max    = nan(n_grp_rows, 1);
g_sh_rng    = nan(n_grp_rows, 1);
g_z_med     = nan(n_grp_rows, 1);
g_z_min     = nan(n_grp_rows, 1);
g_z_max     = nan(n_grp_rows, 1);
g_z_rng     = nan(n_grp_rows, 1);
g_r2_med    = nan(n_grp_rows, 1);

row_i = 0;
for ti = 1:numel(T_ORDER)
    T_cur = T_ORDER(ti);
    for mi = 1:numel(MODES)
        mname = MODES(mi);
        sub = KIN_CELLS_TBL(KIN_CELLS_TBL.T == T_cur & KIN_CELLS_TBL.mode == mname, :);
        row_i = row_i + 1;
        g_T(row_i)      = T_cur;
        g_mode(row_i)   = mname;
        g_ncells(row_i) = height(sub);

        sh_v = sub.sqrt_share(isfinite(sub.sqrt_share));
        if ~isempty(sh_v)
            g_sh_med(row_i) = median(sh_v);
            g_sh_min(row_i) = min(sh_v);
            g_sh_max(row_i) = max(sh_v);
            g_sh_rng(row_i) = max(sh_v) - min(sh_v);
        end

        z_v = sub.z(isfinite(sub.z));
        if ~isempty(z_v)
            g_z_med(row_i) = median(z_v);
            g_z_min(row_i) = min(z_v);
            g_z_max(row_i) = max(z_v);
            g_z_rng(row_i) = max(z_v) - min(z_v);
        end

        r2_v = sub.r2(isfinite(sub.r2));
        if ~isempty(r2_v)
            g_r2_med(row_i) = median(r2_v);
        end
    end
end

KIN_GRP_TBL = table(g_T, g_mode, g_ncells, g_sh_med, g_sh_min, g_sh_max, g_sh_rng, ...
    g_z_med, g_z_min, g_z_max, g_z_rng, g_r2_med, ...
    'VariableNames', {'T', 'mode', 'n_cells', ...
    'sqrt_share_median', 'sqrt_share_min', 'sqrt_share_max', 'sqrt_share_range', ...
    'z_median', 'z_min', 'z_max', 'z_range', 'r2_median'});
writetable(KIN_GRP_TBL, fullfile(OUT_DIR, 'mode_kinetics_groups.csv'));

log_msg(log_fid, '[요약] 온도 x 모드별 진행식 형태 (sqrt-shape share & power-law z):\n');
log_msg(log_fid, '%-3s | %-6s | %2s | %19s | %19s | %6s\n', ...
    'T', 'mode', 'n', 'sqrt_share med[rng]', 'z med[min, max]', 'R2_med');
log_msg(log_fid, '%s\n', repmat('-', 1, 68));
for i = 1:height(KIN_GRP_TBL)
    log_msg(log_fid, '%-3s | %-6s | %2d | %5.2f [%4.2f, %4.2f] | %5.2f [%4.2f, %4.2f] | %6.3f\n', ...
        KIN_GRP_TBL.T(i), KIN_GRP_TBL.mode(i), KIN_GRP_TBL.n_cells(i), ...
        KIN_GRP_TBL.sqrt_share_median(i), KIN_GRP_TBL.sqrt_share_min(i), KIN_GRP_TBL.sqrt_share_max(i), ...
        KIN_GRP_TBL.z_median(i), KIN_GRP_TBL.z_min(i), KIN_GRP_TBL.z_max(i), ...
        KIN_GRP_TBL.r2_median(i));
end
log_msg(log_fid, '저장 -> out\\mode_kinetics_cells.csv, out\\mode_kinetics_groups.csv\n\n');


%% 3 모드-저항 연결
log_msg(log_fid, '--- [3] 모드 -> 저항 변화(dR0, dRarc) 연결 분석 시작 ---\n');

TARGETS = ["dR0", "dRarc"];
N_BOOT  = 1000;

% (a) 온도별 다중선형회귀 (무절편) + 셀 클러스터 부트스트랩
n_link_rows = numel(TARGETS) * numel(T_ORDER) * numel(MODES);
l_target  = strings(n_link_rows, 1);
l_T       = strings(n_link_rows, 1);
l_mode    = strings(n_link_rows, 1);
l_coef    = nan(n_link_rows, 1);
l_cilo    = nan(n_link_rows, 1);
l_cihi    = nan(n_link_rows, 1);
l_clear   = false(n_link_rows, 1);
l_r2      = nan(n_link_rows, 1);
l_nrows   = zeros(n_link_rows, 1);
l_ncells  = zeros(n_link_rows, 1);
l_cond    = nan(n_link_rows, 1);
l_rLP     = nan(n_link_rows, 1);
l_rLN     = nan(n_link_rows, 1);
l_rPN     = nan(n_link_rows, 1);

row_i = 0;
for tgi = 1:numel(TARGETS)
    tname = TARGETS(tgi);
    for ti = 1:numel(T_ORDER)
        T_cur = T_ORDER(ti);

        mask = (J.T == T_cur) & (J.loss > 0) & isfinite(J.(tname)) & ...
            isfinite(J.LLI) & isfinite(J.LAM_PE) & isfinite(J.LAM_NE);
        sub = J(mask, :);
        y = sub.(tname);
        X = [sub.LLI, sub.LAM_PE, sub.LAM_NE];
        cells_g = sub.cell;

        n_rows_g  = numel(y);
        u_cells_g = unique(cells_g);
        n_cells_g = numel(u_cells_g);

        if n_rows_g >= 3 && rank(X) == 3
            beta = X \ y;
            yhat = X * beta;
            ss_res = sum((y - yhat).^2);
            ss_tot = sum((y - mean(y)).^2);
            r2_g = 1 - ss_res / max(ss_tot, eps);
            cond_X = cond(X);
            rLP = pearson_r(X(:, 1), X(:, 2));
            rLN = pearson_r(X(:, 1), X(:, 3));
            rPN = pearson_r(X(:, 2), X(:, 3));

            % 셀 클러스터 부트스트랩 (1000회, rng(1))
            rng(1);
            cell_idx_map = cell(n_cells_g, 1);
            for k = 1:n_cells_g
                cell_idx_map{k} = find(cells_g == u_cells_g(k));
            end

            boot_betas = nan(N_BOOT, 3);
            for b = 1:N_BOOT
                pick = randi(n_cells_g, n_cells_g, 1);
                idx_b = vertcat(cell_idx_map{pick});
                Xb = X(idx_b, :);
                yb = y(idx_b);
                if size(Xb, 1) >= 3 && rank(Xb) == 3 && cond(Xb) < 1e10
                    boot_betas(b, :) = (Xb \ yb)';
                end
            end

            ci_lo_vec = nan(3, 1);
            ci_hi_vec = nan(3, 1);
            for m = 1:3
                ci_lo_vec(m) = my_prctile(boot_betas(:, m), 2.5);
                ci_hi_vec(m) = my_prctile(boot_betas(:, m), 97.5);
            end
        else
            beta = nan(3, 1);
            r2_g = NaN; cond_X = NaN;
            rLP = NaN; rLN = NaN; rPN = NaN;
            ci_lo_vec = nan(3, 1);
            ci_hi_vec = nan(3, 1);
        end

        for mi = 1:3
            row_i = row_i + 1;
            l_target(row_i) = tname;
            l_T(row_i)      = T_cur;
            l_mode(row_i)   = MODES(mi);
            l_coef(row_i)   = beta(mi);
            l_cilo(row_i)   = ci_lo_vec(mi);
            l_cihi(row_i)   = ci_hi_vec(mi);
            l_clear(row_i)  = isfinite(ci_lo_vec(mi)) && isfinite(ci_hi_vec(mi)) && ...
                ((ci_lo_vec(mi) > 0) || (ci_hi_vec(mi) < 0));
            l_r2(row_i)     = r2_g;
            l_nrows(row_i)  = n_rows_g;
            l_ncells(row_i) = n_cells_g;
            l_cond(row_i)   = cond_X;
            l_rLP(row_i)    = rLP;
            l_rLN(row_i)    = rLN;
            l_rPN(row_i)    = rPN;
        end
    end
end

LINK_TBL = table(l_target, l_T, l_mode, l_coef, l_cilo, l_cihi, l_clear, ...
    l_r2, l_nrows, l_ncells, l_cond, l_rLP, l_rLN, l_rPN, ...
    'VariableNames', {'target', 'T', 'mode', 'coef_mOhm_per_pctp', 'ci_lo', 'ci_hi', 'clear', ...
    'r2', 'n_rows', 'n_cells', 'cond_X', 'corr_LLI_PE', 'corr_LLI_NE', 'corr_PE_NE'});
writetable(LINK_TBL, fullfile(OUT_DIR, 'mode_resistance_link.csv'));

log_msg(log_fid, '\n[요약] 온도별 모드 -> 저항 계수 (mOhm per %%p) 및 95%% 셀-부트스트랩 CI:\n');
log_msg(log_fid, '%-5s | %-2s | %-6s | %8s | %17s | %-5s | %6s | %6s | %4s/%2s\n', ...
    'dR', 'T', 'mode', 'coef', '95% CI', 'clear', 'R2', 'condX', 'nrow', 'nc');
log_msg(log_fid, '%s\n', repmat('-', 1, 86));
for i = 1:height(LINK_TBL)
    clr_str = "no";
    if LINK_TBL.clear(i), clr_str = "CLEAR"; end
    log_msg(log_fid, '%-5s | %-2s | %-6s | %8.4f | [%7.4f, %7.4f] | %-5s | %6.3f | %6.1f | %4d/%2d\n', ...
        LINK_TBL.target(i), LINK_TBL.T(i), LINK_TBL.mode(i), ...
        LINK_TBL.coef_mOhm_per_pctp(i), LINK_TBL.ci_lo(i), LINK_TBL.ci_hi(i), ...
        clr_str, LINK_TBL.r2(i), LINK_TBL.cond_X(i), LINK_TBL.n_rows(i), LINK_TBL.n_cells(i));
end

% 다중공선성(모드 간 상관계수) 진단 출력
log_msg(log_fid, '\n[진단] 온도별 모드 간 상관계수 및 조건수 (dR0 기준 유효 행):\n');
log_msg(log_fid, '%-5s | %-2s | %7s | %11s | %11s | %11s\n', ...
    'dR', 'T', 'cond_X', 'r(LLI,PE)', 'r(LLI,NE)', 'r(PE,NE)');
log_msg(log_fid, '%s\n', repmat('-', 1, 60));
for tgi = 1:numel(TARGETS)
    for ti = 1:numel(T_ORDER)
        idx1 = find(LINK_TBL.target == TARGETS(tgi) & LINK_TBL.T == T_ORDER(ti), 1);
        log_msg(log_fid, '%-5s | %-2s | %7.2f | %11.3f | %11.3f | %11.3f\n', ...
            LINK_TBL.target(idx1), LINK_TBL.T(idx1), LINK_TBL.cond_X(idx1), ...
            LINK_TBL.corr_LLI_PE(idx1), LINK_TBL.corr_LLI_NE(idx1), LINK_TBL.corr_PE_NE(idx1));
    end
end

% (b) 통합(pooled) 모델 비교: 공유 계수 (k=3) vs 온도별 계수 (k=9), AIC = n*log(RSS/n) + 2*k
log_msg(log_fid, '\n[비교] 통합 모델 AIC (공유 계수 k=3 vs 온도별 계수 k=9, AIC = n*log(RSS/n) + 2k):\n');
log_msg(log_fid, '%-5s | %4s | %10s | %10s | %10s | %10s | %10s | %-12s\n', ...
    'dR', 'n', 'RSS_shared', 'RSS_Tspec', 'AIC_shared', 'AIC_Tspec', 'dAIC(S-T)', 'preferred');
log_msg(log_fid, '%s\n', repmat('-', 1, 90));

for tgi = 1:numel(TARGETS)
    tname = TARGETS(tgi);
    mask = (J.loss > 0) & isfinite(J.(tname)) & isfinite(J.LLI) & isfinite(J.LAM_PE) & isfinite(J.LAM_NE);
    sub = J(mask, :);
    y = sub.(tname);
    n_p = numel(y);
    X_sh = [sub.LLI, sub.LAM_PE, sub.LAM_NE];

    % 공유 계수 모델 (k = 3)
    b_sh = X_sh \ y;
    rss_sh = sum((y - X_sh * b_sh).^2);
    k_sh = 3;
    aic_sh = n_p * log(rss_sh / n_p) + 2 * k_sh;

    % 온도별 계수 모델 (k = 9)
    iC = double(sub.T == "C");
    iR = double(sub.T == "R");
    iH = double(sub.T == "H");
    X_ts = [X_sh .* iC, X_sh .* iR, X_sh .* iH];
    b_ts = X_ts \ y;
    rss_ts = sum((y - X_ts * b_ts).^2);
    k_ts = 9;
    aic_ts = n_p * log(rss_ts / n_p) + 2 * k_ts;

    daic = aic_sh - aic_ts; % 양수이면 온도별 계수 모델 우세
    if aic_ts < aic_sh
        pref = "T-specific";
    else
        pref = "Shared";
    end
    log_msg(log_fid, '%-5s | %4d | %10.3f | %10.3f | %10.2f | %10.2f | %+10.2f | %-12s\n', ...
        tname, n_p, rss_sh, rss_ts, aic_sh, aic_ts, daic, pref);
end

% (c) 온도별 단순 민감도 d(dR)/d(loss) (원점 통과 기울기, mOhm per % capacity loss) + 부트스트랩 CI
n_sens_rows = numel(TARGETS) * numel(T_ORDER);
s_target = strings(n_sens_rows, 1);
s_T      = strings(n_sens_rows, 1);
s_slope  = nan(n_sens_rows, 1);
s_cilo   = nan(n_sens_rows, 1);
s_cihi   = nan(n_sens_rows, 1);
s_clear  = false(n_sens_rows, 1);
s_r2     = nan(n_sens_rows, 1);
s_nrows  = zeros(n_sens_rows, 1);
s_ncells = zeros(n_sens_rows, 1);

row_i = 0;
for tgi = 1:numel(TARGETS)
    tname = TARGETS(tgi);
    for ti = 1:numel(T_ORDER)
        T_cur = T_ORDER(ti);
        mask = (J.T == T_cur) & (J.loss > 0) & isfinite(J.(tname)) & isfinite(J.loss);
        sub = J(mask, :);
        y = sub.(tname);
        x = sub.loss;
        cells_g = sub.cell;

        n_rows_g  = numel(y);
        u_cells_g = unique(cells_g);
        n_cells_g = numel(u_cells_g);

        row_i = row_i + 1;
        s_target(row_i) = tname;
        s_T(row_i)      = T_cur;
        s_nrows(row_i)  = n_rows_g;
        s_ncells(row_i) = n_cells_g;

        if n_rows_g >= 2 && sum(x.^2) > 0
            sl = (x' * y) / (x' * x);
            yhat = sl * x;
            ss_res = sum((y - yhat).^2);
            ss_tot = sum((y - mean(y)).^2);
            r2_s = 1 - ss_res / max(ss_tot, eps);

            rng(1);
            cell_idx_map = cell(n_cells_g, 1);
            for k = 1:n_cells_g
                cell_idx_map{k} = find(cells_g == u_cells_g(k));
            end

            boot_sl = nan(N_BOOT, 1);
            for b = 1:N_BOOT
                pick = randi(n_cells_g, n_cells_g, 1);
                idx_b = vertcat(cell_idx_map{pick});
                xb = x(idx_b);
                yb = y(idx_b);
                denom_b = xb' * xb;
                if denom_b > 1e-12
                    boot_sl(b) = (xb' * yb) / denom_b;
                end
            end

            cilo = my_prctile(boot_sl, 2.5);
            cihi = my_prctile(boot_sl, 97.5);

            s_slope(row_i) = sl;
            s_cilo(row_i)  = cilo;
            s_cihi(row_i)  = cihi;
            s_clear(row_i) = isfinite(cilo) && isfinite(cihi) && ((cilo > 0) || (cihi < 0));
            s_r2(row_i)    = r2_s;
        end
    end
end

SENS_TBL = table(s_target, s_T, s_slope, s_cilo, s_cihi, s_clear, s_r2, s_nrows, s_ncells, ...
    'VariableNames', {'target', 'T', 'slope_mOhm_per_pct_loss', 'ci_lo', 'ci_hi', 'clear', ...
    'r2', 'n_rows', 'n_cells'});
writetable(SENS_TBL, fullfile(OUT_DIR, 'resistance_per_loss.csv'));

log_msg(log_fid, '\n[요약] 용량 손실당 저항 증가 민감도 d(dR)/d(loss) [mOhm per %% capacity loss]:\n');
log_msg(log_fid, '%-5s | %-2s | %8s | %17s | %-5s | %6s | %4s/%2s\n', ...
    'dR', 'T', 'slope', '95% CI', 'clear', 'R2', 'nrow', 'nc');
log_msg(log_fid, '%s\n', repmat('-', 1, 66));
for i = 1:height(SENS_TBL)
    clr_str = "no";
    if SENS_TBL.clear(i), clr_str = "CLEAR"; end
    log_msg(log_fid, '%-5s | %-2s | %8.4f | [%7.4f, %7.4f] | %-5s | %6.3f | %4d/%2d\n', ...
        SENS_TBL.target(i), SENS_TBL.T(i), SENS_TBL.slope_mOhm_per_pct_loss(i), ...
        SENS_TBL.ci_lo(i), SENS_TBL.ci_hi(i), clr_str, SENS_TBL.r2(i), ...
        SENS_TBL.n_rows(i), SENS_TBL.n_cells(i));
end
log_msg(log_fid, '저장 -> out\\mode_resistance_link.csv, out\\resistance_per_loss.csv\n\n');


%% 4 그림
log_msg(log_fid, '--- [4] 그림 생성 시작 ---\n');

% (1) fig_ica_curves.png: 1x3 타일, 셀 2 (C), 1 (R), 3 (H)
fig1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1350 420]);
tl1 = tiledlayout(fig1, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
rep_cells = [2, 1, 3];
rep_Ts    = ["C", "R", "H"];
win_shades = [0.92 0.94 0.98; 0.96 0.93 0.92; 0.92 0.96 0.93];

for p = 1:3
    ax = nexttile(tl1, p);
    hold(ax, 'on');
    c = rep_cells(p);
    T_c = rep_Ts(p);
    col = COL_MAP.(T_c);

    idx_c = find(cell_list_ica == c, 1);
    if isempty(idx_c) || ~ica_cell_data(idx_c).ok
        title(ax, sprintf('Cell %d (%s) - No ICA data', c, T_c));
        continue;
    end

    chks = ica_cell_data(idx_c).checkups;
    ref_k = ica_cell_data(idx_c).ref_idx;
    y_ref = chks(ref_k).dqdv;

    % loss = 10%에 가장 가까운 매칭 체크업 탐색
    sub_j = ICA_TBL(ICA_TBL.cell == c & isfinite(ICA_TBL.ica_ah) & ICA_TBL.loss > 0, :);
    aged_k = NaN;
    aged_loss = NaN;
    if ~isempty(sub_j)
        [~, best_j] = min(abs(sub_j.loss - 10));
        aged_loss = sub_j.loss(best_j);
        target_ah = sub_j.ica_ah(best_j);
        [~, aged_k] = min(abs([chks.ah]' - target_ah));
    end

    ymax = max(y_ref, [], 'omitnan');
    if isfinite(aged_k)
        ymax = max(ymax, max(chks(aged_k).dqdv, [], 'omitnan'));
    end
    ymax = max(ymax * 1.12, 1);

    % 3개 전압 창 음영 표시
    for w = 1:3
        patch(ax, [WINS(w,1) WINS(w,2) WINS(w,2) WINS(w,1)], [0 0 ymax ymax], ...
            win_shades(w, :), 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'HandleVisibility', 'off');
        text(ax, mean(WINS(w, :)), ymax * 0.95, WIN_NAMES(w), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 9, 'Color', [0.35 0.35 0.35], 'FontWeight', 'bold');
    end

    plot(ax, V_GRID, y_ref, '--', 'Color', col, 'LineWidth', 1.6, 'DisplayName', 'Ref');
    if isfinite(aged_k)
        plot(ax, V_GRID, chks(aged_k).dqdv, '-', 'Color', col, 'LineWidth', 2.0, ...
            'DisplayName', sprintf('Aged (%.1f%% loss)', aged_loss));
    end

    xlim(ax, [3.45 4.18]);
    ylim(ax, [0 ymax]);
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'Voltage [V]');
    if p == 1
        ylabel(ax, 'dQ/dV [Ah/V]');
    end
    Tc_val = J.T_C(find(J.cell == c, 1));
    title(ax, sprintf('Cell %d (%s, %d °C) | loss = %.1f%%', c, T_c, round(Tc_val), aged_loss));
    legend(ax, 'Location', 'southwest', 'Box', 'off');
end
exportgraphics(fig1, fullfile(OUT_DIR, 'fig_ica_curves.png'), 'Resolution', 150, 'BackgroundColor', 'white');
close(fig1);
log_msg(log_fid, '저장 -> out\\fig_ica_curves.png\n');

% (2) fig_ica_peaks_vs_loss.png: 3x3 타일 (행=창 W1~W3, 열=height_rel/area_rel/shift_mV)
fig2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 950]);
tl2 = tiledlayout(fig2, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
col_suffix = ["height_rel", "area_rel", "shift_mV"];
col_ylabel = ["Peak height (rel.)", "Window area (rel.)", "Peak shift [mV]"];
ref_lines  = [1, 1, 0];

for w = 1:3
    for mc = 1:3
        ax = nexttile(tl2, (w - 1) * 3 + mc);
        hold(ax, 'on');
        yline(ax, ref_lines(mc), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'HandleVisibility', 'off');

        fname = sprintf('%s_%s', WIN_NAMES(w), col_suffix(mc));
        for ti = 1:numel(T_ORDER)
            T_cur = T_ORDER(ti);
            col = COL_MAP.(T_cur);
            mask = (ICA_TBL.T == T_cur) & isfinite(ICA_TBL.loss) & isfinite(ICA_TBL.(fname));
            scatter(ax, ICA_TBL.loss(mask), ICA_TBL.(fname)(mask), 22, col, 'filled', ...
                'MarkerFaceAlpha', 0.72, 'DisplayName', T_LABEL(ti));
        end

        grid(ax, 'on');
        box(ax, 'on');
        if w == 1
            title(ax, col_ylabel(mc), 'FontWeight', 'bold');
        end
        if mc == 1
            ylabel(ax, sprintf('%s [%.2f–%.2f V]\n%s', WIN_NAMES(w), WINS(w,1), WINS(w,2), col_ylabel(mc)));
        else
            ylabel(ax, col_ylabel(mc));
        end
        if w == 3
            xlabel(ax, 'Capacity loss [%]');
        end
        if w == 1 && mc == 1
            legend(ax, 'Location', 'southwest', 'Box', 'off', 'FontSize', 8.5);
        end
    end
end
exportgraphics(fig2, fullfile(OUT_DIR, 'fig_ica_peaks_vs_loss.png'), 'Resolution', 150, 'BackgroundColor', 'white');
close(fig2);
log_msg(log_fid, '저장 -> out\\fig_ica_peaks_vs_loss.png\n');

% (3) fig_mode_trajectories.png: 1x3 타일 (LLI, LAM_PE, LAM_NE) vs days (선형 x축)
fig3 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1350 430]);
tl3 = tiledlayout(fig3, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for mi = 1:numel(MODES)
    mname = MODES(mi);
    ax = nexttile(tl3, mi);
    hold(ax, 'on');

    for ci = 1:numel(all_cells)
        c = all_cells(ci);
        fld = sprintf('c%d_%s', c, mname);
        if ~isfield(kin_fit_store, fld)
            continue;
        end
        st = kin_fit_store.(fld);
        col = COL_MAP.(st.T);

        % 관측 곡선 (얇은 실선)
        plot(ax, st.days, st.y, '-', 'Color', [col 0.75], 'LineWidth', 1.0, 'HandleVisibility', 'off');

        % 피팅 곡선 (점선): 조밀한 days 격자에서 부드럽게 표시
        [du, ia] = unique(st.days);
        eu = st.efc(ia);
        if numel(du) >= 2
            d_fine = linspace(min(du), max(du), 80)';
            e_fine = interp1(du, eu, d_fine, 'linear');
            y_fit = st.a * sqrt(d_fine) + st.b * e_fine;
            plot(ax, d_fine, y_fit, '--', 'Color', [col 0.85], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        end
    end

    % 범례용 더미 핸들
    if mi == 1
        for ti = 1:numel(T_ORDER)
            plot(ax, NaN, NaN, '-', 'Color', COL_MAP.(T_ORDER(ti)), 'LineWidth', 1.6, 'DisplayName', T_LABEL(ti));
        end
        plot(ax, NaN, NaN, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'DisplayName', 'Data');
        plot(ax, NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'DisplayName', 'Fit (a\surdd + b\cdotEFC)');
        legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8.5);
    end

    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'Days since reference');
    if mi == 1
        ylabel(ax, 'Degradation mode [%]');
    end
    title(ax, strrep(mname, '_', '\_'));
end
exportgraphics(fig3, fullfile(OUT_DIR, 'fig_mode_trajectories.png'), 'Resolution', 150, 'BackgroundColor', 'white');
close(fig3);
log_msg(log_fid, '저장 -> out\\fig_mode_trajectories.png\n');

% (4) fig_link_coeffs.png: 모드별 회귀계수 +- 95% CI (dR0, dRarc 패널, 온도 그룹)
fig4 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1050 450]);
tl4 = tiledlayout(fig4, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
offsets = [-0.22, 0.0, 0.22];

for tgi = 1:numel(TARGETS)
    tname = TARGETS(tgi);
    ax = nexttile(tl4, tgi);
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.1, 'HandleVisibility', 'off');

    for ti = 1:numel(T_ORDER)
        T_cur = T_ORDER(ti);
        col = COL_MAP.(T_cur);
        xs = (1:3)' + offsets(ti);
        ys = nan(3, 1);
        err_lo = nan(3, 1);
        err_hi = nan(3, 1);
        is_clr = false(3, 1);

        for mi = 1:3
            idx_r = find(LINK_TBL.target == tname & LINK_TBL.T == T_cur & LINK_TBL.mode == MODES(mi), 1);
            if ~isempty(idx_r)
                ys(mi)     = LINK_TBL.coef_mOhm_per_pctp(idx_r);
                err_lo(mi) = ys(mi) - LINK_TBL.ci_lo(idx_r);
                err_hi(mi) = LINK_TBL.ci_hi(idx_r) - ys(mi);
                is_clr(mi) = LINK_TBL.clear(idx_r);
            end
        end

        % 에러바 + 전체 표시 (clear=채운 마커, 비clear=빈 마커)
        errorbar(ax, xs, ys, err_lo, err_hi, 'LineStyle', 'none', 'Color', col, ...
            'LineWidth', 1.6, 'CapSize', 6, 'HandleVisibility', 'off');
        if any(is_clr)
            plot(ax, xs(is_clr), ys(is_clr), 'o', 'Color', col, 'MarkerFaceColor', col, ...
                'MarkerSize', 7, 'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
        if any(~is_clr)
            plot(ax, xs(~is_clr), ys(~is_clr), 'o', 'Color', col, 'MarkerFaceColor', 'w', ...
                'MarkerSize', 7, 'LineWidth', 1.4, 'HandleVisibility', 'off');
        end
        % 범례용 더미 핸들
        plot(ax, NaN, NaN, 'o-', 'Color', col, 'MarkerFaceColor', col, ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', T_LABEL(ti));
    end

    if tgi == 1
        plot(ax, NaN, NaN, 'o', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3], ...
            'MarkerSize', 6, 'DisplayName', '95% CI excludes 0 (clear)');
        plot(ax, NaN, NaN, 'o', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', 'w', ...
            'MarkerSize', 6, 'DisplayName', '95% CI includes 0');
        legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 8.5);
    end

    xlim(ax, [0.5 3.5]);
    xticks(ax, 1:3);
    xticklabels(ax, ["LLI", "LAM\_PE", "LAM\_NE"]);
    grid(ax, 'on');
    box(ax, 'on');
    ylabel(ax, 'Coefficient [m\Omega / %p]');
    title(ax, sprintf('%s = c_{LLI}LLI + c_{PE}LAM_{PE} + c_{NE}LAM_{NE}', tname));
end
exportgraphics(fig4, fullfile(OUT_DIR, 'fig_link_coeffs.png'), 'Resolution', 150, 'BackgroundColor', 'white');
close(fig4);
log_msg(log_fid, '저장 -> out\\fig_link_coeffs.png\n');

log_msg(log_fid, '\n=== 전체 분석 완료 (소요 시간: %.1f초) ===\n', toc(t_start));


%% 로컬 함수

function safe_fclose(fid)
    if fid >= 0
        try
            fclose(fid);
        catch
        end
    end
end

function log_msg(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    if fid >= 0
        fprintf(fid, fmt, varargin{:});
    end
end

function v = to_double_vec(col)
    if isnumeric(col)
        v = double(col(:));
    elseif islogical(col)
        v = double(col(:));
    else
        v = str2double(string(col(:)));
    end
end

function dqdv_s = compute_dqdv_smooth(Q_grid, dV, sigma)
    N = numel(Q_grid);
    dqdv_s = nan(N, 1);
    idx_fin = find(isfinite(Q_grid));
    if numel(idx_fin) < 5
        return;
    end
    i1 = idx_fin(1);
    i2 = idx_fin(end);
    q = Q_grid(i1:i2);
    % 내부 NaN이 있을 경우 선형 보간
    if any(~isfinite(q))
        t = (1:numel(q))';
        f = isfinite(q);
        q = interp1(t(f), q(f), t, 'linear');
    end
    m = numel(q);
    dq = zeros(m, 1);
    dq(2:m-1) = (q(3:m) - q(1:m-2)) / (2 * dV);
    dq(1) = (q(2) - q(1)) / dV;
    dq(m) = (q(m) - q(m-1)) / dV;

    % 가우시안 커널 생성 및 끝단 복제 패딩 + conv 스무딩
    hw = ceil(3 * sigma);
    x = (-hw:hw)';
    g = exp(-0.5 * (x / sigma).^2);
    g = g / sum(g);
    dq_pad = [repmat(dq(1), hw, 1); dq; repmat(dq(end), hw, 1)];
    dq_sm = conv(dq_pad, g, 'valid');
    dqdv_s(i1:i2) = dq_sm;
end

function [p_idx, p_prom] = find_local_maxima_prom(y)
    y = y(:);
    N = numel(y);
    peaks = [];
    i = 2;
    while i <= N - 1
        if ~isfinite(y(i)) || ~isfinite(y(i-1)) || ~isfinite(y(i+1))
            i = i + 1;
            continue;
        end
        if y(i) > y(i-1)
            j = i;
            while j < N && isfinite(y(j+1)) && y(j+1) == y(i)
                j = j + 1;
            end
            if j < N && isfinite(y(j+1)) && y(j) > y(j+1)
                peaks(end + 1, 1) = round(0.5 * (i + j)); %#ok<AGROW>
            end
            i = j + 1;
        else
            i = i + 1;
        end
    end

    nP = numel(peaks);
    p_idx = peaks;
    p_prom = zeros(nP, 1);
    if nP == 0
        return;
    end

    fin_all = find(isfinite(y));
    first_f = fin_all(1);
    last_f  = fin_all(end);

    for k = 1:nP
        p = peaks(k);
        yp = y(p);

        % 좌측으로 yp보다 큰 점까지 탐색
        L = first_f;
        for j = p-1:-1:first_f
            if isfinite(y(j)) && y(j) > yp
                L = j;
                break;
            end
        end
        left_min = min(y(L:p), [], 'omitnan');

        % 우측으로 yp보다 큰 점까지 탐색
        R = last_f;
        for j = p+1:last_f
            if isfinite(y(j)) && y(j) > yp
                R = j;
                break;
            end
        end
        right_min = min(y(p:R), [], 'omitnan');

        p_prom(k) = yp - max(left_min, right_min);
    end
end

function r = avg_rank(x)
    x = x(:);
    n = numel(x);
    [xs, ord] = sort(x);
    rs = zeros(n, 1);
    i = 1;
    while i <= n
        j = i;
        while j < n && xs(j + 1) == xs(i)
            j = j + 1;
        end
        rs(i:j) = 0.5 * (i + j);
        i = j + 1;
    end
    r = zeros(n, 1);
    r(ord) = rs;
end

function r = pearson_r(x, y)
    x = x(:);
    y = y(:);
    v = isfinite(x) & isfinite(y);
    x = x(v);
    y = y(v);
    if numel(x) < 2
        r = NaN;
        return;
    end
    dx = x - mean(x);
    dy = y - mean(y);
    denom = sqrt(sum(dx.^2) * sum(dy.^2));
    if denom <= 0
        r = NaN;
    else
        r = sum(dx .* dy) / denom;
    end
end

function r = spearman_r(x, y)
    x = x(:);
    y = y(:);
    v = isfinite(x) & isfinite(y);
    x = x(v);
    y = y(v);
    if numel(x) < 3
        r = NaN;
        return;
    end
    r = pearson_r(avg_rank(x), avg_rank(y));
end

function r = partial_spearman_r(x, y, z)
    x = x(:);
    y = y(:);
    z = z(:);
    v = isfinite(x) & isfinite(y) & isfinite(z);
    x = x(v);
    y = y(v);
    z = z(v);
    n = numel(x);
    if n < 4
        r = NaN;
        return;
    end
    Z = [ones(n, 1), z];
    if rank(Z) < 2
        r = NaN;
        return;
    end
    rx = x - Z * (Z \ x);
    ry = y - Z * (Z \ y);
    r = spearman_r(rx, ry);
end

function val = my_prctile(x, p)
    x = x(isfinite(x));
    n = numel(x);
    if n == 0
        val = NaN;
        return;
    end
    xs = sort(x(:));
    if n == 1
        val = xs(1);
        return;
    end
    q = ((0.5:1:(n - 0.5))' / n) * 100;
    p_c = min(max(p, q(1)), q(end));
    val = interp1(q, xs, p_c, 'linear');
end