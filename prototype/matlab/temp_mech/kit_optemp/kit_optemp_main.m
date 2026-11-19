%% 0 설정
% KIT battery dataset: 분리된 노화 효과(25 degC) vs 순간 온도 효과(Ea) 분석
RAW_DIR      = '..\..\..\..\data\kit_radar\eis_addendum\cell_eisv2_fixed';   % temp_mech로 이동하며 경로 수정
META_CSV     = 'kit_meta.csv';
SOH_CSV      = 'kit_soh.csv';
OUT_DIR      = 'out';
EIS_FEAT_DIR = '..\..';   % temp_mech로 이동하며 경로 수정 (eis_features.m은 matlab 루트)

SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR)
    SCRIPT_DIR = pwd;
end

raw_dir_abs  = resolve_rel_path(SCRIPT_DIR, RAW_DIR);
meta_csv_abs = resolve_rel_path(SCRIPT_DIR, META_CSV);
soh_csv_abs  = resolve_rel_path(SCRIPT_DIR, SOH_CSV);
out_dir_abs  = resolve_rel_path(SCRIPT_DIR, OUT_DIR);
addpath(resolve_rel_path(SCRIPT_DIR, EIS_FEAT_DIR));

if ~exist(out_dir_abs, 'dir')
    mkdir(out_dir_abs);
end

REQ_COLS = {'timestamp_s', 'sd_block_id', 'is_rt', 'soc_nom', 'valid', ...
            't_avg_degC', 'freq_Hz', 'z_amp_mOhm', 'z_ph_deg'};

MIN_FREQ_PTS     = 25;
MAX_PAIR_DT_S    = 2 * 86400;   % 2 days
MIN_PAIR_DTEMP_K = 8.0;         % 8 K
MAX_SOH_EXTRAP_S = 3 * 86400;   % 3 days
R_GAS_KJ         = 8.314e-3;    % kJ/(mol*K)
VALID_T_AGE      = [0, 10, 40];
COMPONENTS       = {'R0', 'R_arc', 'R_lf', 'Re_lowf'};
COMP_LABELS      = {'R_0', 'R_{arc}', 'R_{lf}', 'Re_{lowf}'};
SOH_BIN_NAMES    = {'>=0.98', '0.95-0.98', '0.90-0.95', '<0.90'};

% 메타데이터 및 SOH 로드
meta_tbl = read_csv_clean(meta_csv_abs);
soh_tbl  = read_csv_clean(soh_csv_abs);

meta_tbl.cell  = string(meta_tbl.cell);
meta_tbl.group = string(meta_tbl.group);
soh_tbl.cell   = string(soh_tbl.cell);

meta_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:height(meta_tbl)
    meta_map(char(meta_tbl.cell(i))) = i;
end

%% 1 스펙트럼 추출
file_list = dir(fullfile(raw_dir_abs, 'cell_eisv2_*.csv'));
n_files = numel(file_list);
fprintf('Found %d EIS files in %s\n', n_files, raw_dir_abs);

spec_chunks = cell(n_files, 1);
skip_logs   = {};
total_skip_npts = 0;
total_skip_sign = 0;
n_files_ok      = 0;

for fi = 1:n_files
    fname = file_list(fi).name;
    fpath = fullfile(file_list(fi).folder, fname);
    cid_str = erase(erase(string(fname), 'cell_eisv2_'), '.csv');
    cid_chr = char(cid_str);

    if ~isKey(meta_map, cid_chr)
        skip_logs{end+1} = sprintf('%s\tSKIP_META: cell %s not in kit_meta.csv', fname, cid_chr); %#ok<SAGROW>
        if mod(fi, 20) == 0 || fi == n_files
            fprintf('Progress: %d / %d files processed\n', fi, n_files);
        end
        continue;
    end

    try
        opts = detectImportOptions(fpath, 'FileType', 'text', 'Delimiter', ';');
        missing_cols = setdiff(REQ_COLS, opts.VariableNames);
        if ~isempty(missing_cols)
            error('Missing required columns: %s', strjoin(missing_cols, ', '));
        end
        opts.SelectedVariableNames = REQ_COLS;
        opts = setvartype(opts, REQ_COLS, 'double');
        opts = setvaropts(opts, REQ_COLS, 'TreatAsMissing', {'nan', 'NaN', 'NA', ''});
        Traw = readtable(fpath, opts);

        valid_mask = (Traw.valid == 1) & ...
                     isfinite(Traw.freq_Hz) & ...
                     isfinite(Traw.sd_block_id) & ...
                     isfinite(Traw.z_amp_mOhm) & ...
                     isfinite(Traw.z_ph_deg);

        if ~any(valid_mask)
            skip_logs{end+1} = sprintf('%s\tSKIP_EMPTY: no valid rows', fname); %#ok<SAGROW>
            if mod(fi, 20) == 0 || fi == n_files
                fprintf('Progress: %d / %d files processed\n', fi, n_files);
            end
            continue;
        end

        ts_v   = Traw.timestamp_s(valid_mask);
        blk_v  = Traw.sd_block_id(valid_mask);
        rt_v   = Traw.is_rt(valid_mask);
        soc_v  = Traw.soc_nom(valid_mask);
        tavg_v = Traw.t_avg_degC(valid_mask);
        f_v    = Traw.freq_Hz(valid_mask);
        amp_v  = Traw.z_amp_mOhm(valid_mask);
        ph_v   = Traw.z_ph_deg(valid_mask);

        % 블록별 안정 정렬 (블록 내 원래 순서 유지)
        [blk_s, ord] = sort(blk_v);
        ts_v   = ts_v(ord);
        rt_v   = rt_v(ord);
        soc_v  = soc_v(ord);
        tavg_v = tavg_v(ord);
        f_v    = f_v(ord);
        amp_v  = amp_v(ord);
        ph_v   = ph_v(ord);

        edges = [1; find(diff(blk_s) ~= 0) + 1; numel(blk_s) + 1];
        n_blk = numel(edges) - 1;

        c_block   = zeros(n_blk, 1);
        c_is_rt   = zeros(n_blk, 1);
        c_soc     = zeros(n_blk, 1);
        c_t       = zeros(n_blk, 1);
        c_T_meas  = zeros(n_blk, 1);
        c_R0      = nan(n_blk, 1);
        c_R_arc   = nan(n_blk, 1);
        c_R_lf    = nan(n_blk, 1);
        c_Re_lowf = nan(n_blk, 1);

        ns = 0;
        f_skip_npts = 0;
        f_skip_sign = 0;

        for bi = 1:n_blk
            i1 = edges(bi);
            i2 = edges(bi + 1) - 1;
            npts = i2 - i1 + 1;
            if npts < MIN_FREQ_PTS
                f_skip_npts = f_skip_npts + 1;
                continue;
            end

            fb   = f_v(i1:i2);
            ampb = amp_v(i1:i2);
            phb  = ph_v(i1:i2);

            reb  = ampb .* cos(deg2rad(phb));
            mimb = ampb .* sin(deg2rad(phb));

            % 부호 건전성 검사: 1-100 Hz 구간의 mim 중앙값이 양수여야 함
            mid_mask = (fb >= 1) & (fb <= 100);
            if ~any(mid_mask) || ~(median(mimb(mid_mask), 'omitnan') > 0)
                f_skip_sign = f_skip_sign + 1;
                continue;
            end

            d = eis_features(fb, reb, mimb);
            [~, ilow] = min(fb);
            re_lowf_val = reb(ilow);

            if isfinite(d.R_arc) && isfinite(d.R0)
                r_lf_val = re_lowf_val - (d.R0 + d.R_arc);
            else
                r_lf_val = NaN;
            end

            ns = ns + 1;
            c_block(ns)   = blk_s(i1);
            c_is_rt(ns)   = rt_v(i1);
            c_soc(ns)     = soc_v(i1);
            c_t(ns)       = ts_v(i1);
            c_T_meas(ns)  = mean(tavg_v(i1:i2), 'omitnan');
            c_R0(ns)      = d.R0;
            c_R_arc(ns)   = d.R_arc;
            c_R_lf(ns)    = r_lf_val;
            c_Re_lowf(ns) = re_lowf_val;
        end

        total_skip_npts = total_skip_npts + f_skip_npts;
        total_skip_sign = total_skip_sign + f_skip_sign;

        if ns > 0
            idx = 1:ns;
            spec_chunks{fi} = table( ...
                repmat(cid_str, ns, 1), c_block(idx), c_is_rt(idx), c_soc(idx), ...
                c_t(idx), c_T_meas(idx), c_R0(idx), c_R_arc(idx), c_R_lf(idx), c_Re_lowf(idx), ...
                'VariableNames', {'cell', 'block', 'is_rt', 'soc', 't', 'T_meas', ...
                                  'R0', 'R_arc', 'R_lf', 'Re_lowf'});
            n_files_ok = n_files_ok + 1;
        else
            skip_logs{end+1} = sprintf('%s\tSKIP_QC: 0 valid spectra (npts<25: %d, sign_fail: %d)', ...
                fname, f_skip_npts, f_skip_sign); %#ok<SAGROW>
        end
    catch ME
        skip_logs{end+1} = sprintf('%s\tREAD_ERROR: %s', fname, ME.message); %#ok<SAGROW>
    end

    if mod(fi, 20) == 0 || fi == n_files
        fprintf('Progress: %d / %d files processed\n', fi, n_files);
    end
end

non_empty = ~cellfun(@isempty, spec_chunks);
if any(non_empty)
    spec_tbl = vertcat(spec_chunks{non_empty});
else
    error('No valid spectra extracted from %s.', raw_dir_abs);
end

% skipped.txt 기록
fid_skip = fopen(fullfile(out_dir_abs, 'skipped.txt'), 'w');
if fid_skip ~= -1
    fprintf(fid_skip, '# Skipped files and spectrum QC summary\n');
    fprintf(fid_skip, '# Total files scanned: %d, files with valid spectra: %d, skipped files: %d\n', ...
        n_files, n_files_ok, numel(skip_logs));
    fprintf(fid_skip, '# Spectra skipped for npts < %d: %d\n', MIN_FREQ_PTS, total_skip_npts);
    fprintf(fid_skip, '# Spectra skipped for capacitive sign check (median mim 1-100 Hz <= 0): %d\n', total_skip_sign);
    for si = 1:numel(skip_logs)
        fprintf(fid_skip, '%s\n', skip_logs{si});
    end
    fclose(fid_skip);
end

% T_age 결합 및 T_age별 is_rt==0 스펙트럼 개수 검증 출력
n_spec = height(spec_tbl);
spec_T_age = nan(n_spec, 1);
for i = 1:n_spec
    m_idx = meta_map(char(spec_tbl.cell(i)));
    spec_T_age(i) = meta_tbl.T_age(m_idx);
end
spec_tbl.T_age = spec_T_age;

fprintf('\n=== Spectrum count verification by T_age ===\n');
all_T_ages = sort(unique(meta_tbl.T_age(isfinite(meta_tbl.T_age))));
fprintf('%8s %14s %14s %10s\n', 'T_age', 'n(is_rt==0)', 'n(is_rt==1)', 'n_cells');
for ti = 1:numel(all_T_ages)
    ta = all_T_ages(ti);
    m_ta = (spec_tbl.T_age == ta);
    n_op = sum(m_ta & (spec_tbl.is_rt == 0));
    n_rt = sum(m_ta & (spec_tbl.is_rt == 1));
    nc   = numel(unique(spec_tbl.cell(m_ta)));
    fprintf('%8g %14d %14d %10d\n', ta, n_op, n_rt, nc);
end

%% 2 짝짓기
% SOH 셀별 룩업 맵 구성 (중복 t 평균 및 정렬)
soh_map = build_soh_map(soh_tbl);

% T_age in {0, 10, 40} 셀만 사용
spec_pair_pool = spec_tbl(ismember(spec_tbl.T_age, VALID_T_AGE), :);
pair_cells = unique(spec_pair_pool.cell);
n_pair_cells = numel(pair_cells);
pair_chunks = cell(n_pair_cells, 1);

for ci = 1:n_pair_cells
    cid = pair_cells(ci);
    sub_c = spec_pair_pool(spec_pair_pool.cell == cid, :);
    socs_c = unique(sub_c.soc(isfinite(sub_c.soc)));

    cell_rows = cell(numel(socs_c), 1);
    for si = 1:numel(socs_c)
        sc = socs_c(si);
        op_sub = sub_c(sub_c.soc == sc & sub_c.is_rt == 0, :);
        rt_sub = sub_c(sub_c.soc == sc & sub_c.is_rt == 1, :);
        if isempty(op_sub) || isempty(rt_sub)
            continue;
        end

        nop = height(op_sub);
        p_op_idx = zeros(nop, 1);
        p_rt_idx = zeros(nop, 1);
        np = 0;

        for k = 1:nop
            t_op = op_sub.t(k);
            if ~isfinite(t_op)
                continue;
            end
            [dt_abs, j_best] = min(abs(rt_sub.t - t_op));
            if isempty(dt_abs) || ~(dt_abs < MAX_PAIR_DT_S)
                continue;
            end
            T_op = op_sub.T_meas(k);
            T_rt = rt_sub.T_meas(j_best);
            if ~isfinite(T_op) || ~isfinite(T_rt) || abs(T_op - T_rt) < MIN_PAIR_DTEMP_K
                continue;
            end
            np = np + 1;
            p_op_idx(np) = k;
            p_rt_idx(np) = j_best;
        end

        if np > 0
            io = p_op_idx(1:np);
            ir = p_rt_idx(1:np);
            cell_rows{si} = table( ...
                repmat(cid, np, 1), repmat(sc, np, 1), ...
                op_sub.block(io), rt_sub.block(ir), ...
                op_sub.t(io), rt_sub.t(ir), ...
                (op_sub.t(io) - rt_sub.t(ir)) / 3600, ...
                op_sub.T_meas(io), rt_sub.T_meas(ir), ...
                op_sub.T_meas(io) - rt_sub.T_meas(ir), ...
                op_sub.R0(io), rt_sub.R0(ir), ...
                op_sub.R_arc(io), rt_sub.R_arc(ir), ...
                op_sub.R_lf(io), rt_sub.R_lf(ir), ...
                op_sub.Re_lowf(io), rt_sub.Re_lowf(ir), ...
                'VariableNames', {'cell', 'soc', 'block_op', 'block_rt', ...
                                  't_op', 't_rt', 'dt_hr', ...
                                  'T_meas_op', 'T_meas_rt', 'dT_K', ...
                                  'R0_op', 'R0_rt', 'R_arc_op', 'R_arc_rt', ...
                                  'R_lf_op', 'R_lf_rt', 'Re_lowf_op', 'Re_lowf_rt'});
        end
    end

    ne_s = ~cellfun(@isempty, cell_rows);
    if any(ne_s)
        pair_chunks{ci} = vertcat(cell_rows{ne_s});
    end
end

ne_p = ~cellfun(@isempty, pair_chunks);
if any(ne_p)
    pairs_tbl = vertcat(pair_chunks{ne_p});
else
    error('No valid op-T / RT pairs found.');
end

% RT 스펙트럼 시각에서 SOH, day, efc 선형 보간 (범위 밖 3일 초과 시 NaN)
n_pairs = height(pairs_tbl);
p_soh = nan(n_pairs, 1);
p_day = nan(n_pairs, 1);
p_efc = nan(n_pairs, 1);

u_pcells = unique(pairs_tbl.cell);
for ci = 1:numel(u_pcells)
    cid_chr = char(u_pcells(ci));
    if ~isKey(soh_map, cid_chr)
        continue;
    end
    srec = soh_map(cid_chr);
    m_c = (pairs_tbl.cell == u_pcells(ci));
    tq = pairs_tbl.t_rt(m_c);
    p_soh(m_c) = interp_limited(srec.t, srec.soh, tq, MAX_SOH_EXTRAP_S);
    p_day(m_c) = interp_limited(srec.t, srec.day, tq, MAX_SOH_EXTRAP_S);
    p_efc(m_c) = interp_limited(srec.t, srec.efc, tq, MAX_SOH_EXTRAP_S);
end

pairs_tbl.soh = p_soh;
pairs_tbl.day = p_day;
pairs_tbl.efc = p_efc;

% 메타데이터 열 부착
m_idx_vec = zeros(n_pairs, 1);
for i = 1:n_pairs
    m_idx_vec(i) = meta_map(char(pairs_tbl.cell(i)));
end
pairs_tbl.T_age    = meta_tbl.T_age(m_idx_vec);
pairs_tbl.age_type = meta_tbl.age_type(m_idx_vec);
pairs_tbl.chg_rate = meta_tbl.chg_rate(m_idx_vec);
pairs_tbl.dis_rate = meta_tbl.dis_rate(m_idx_vec);
pairs_tbl.cal_soc  = meta_tbl.cal_soc(m_idx_vec);
pairs_tbl.group    = meta_tbl.group(m_idx_vec);
pairs_tbl.v_min    = meta_tbl.v_min(m_idx_vec);
pairs_tbl.v_max    = meta_tbl.v_max(m_idx_vec);

%% 3 Ea 계산
inv_dT = (1 ./ (pairs_tbl.T_meas_op + 273.15)) - (1 ./ (pairs_tbl.T_meas_rt + 273.15));

for ci = 1:numel(COMPONENTS)
    comp = COMPONENTS{ci};
    x_op = pairs_tbl.([comp '_op']);
    x_rt = pairs_tbl.([comp '_rt']);
    ea   = nan(n_pairs, 1);

    ok = isfinite(x_op) & isfinite(x_rt) & (x_op > 0) & (x_rt > 0) & ...
         isfinite(inv_dT) & (inv_dT ~= 0);
    ea(ok) = R_GAS_KJ .* log(x_op(ok) ./ x_rt(ok)) ./ inv_dT(ok);
    pairs_tbl.(['Ea_' comp]) = ea;
end

writetable(pairs_tbl, fullfile(out_dir_abs, 'optemp_pairs.csv'));

%% 4 요약/추세
% SOH 구간 라벨 부여
soh_bin_col = strings(n_pairs, 1);
soh_bin_col(pairs_tbl.soh >= 0.98) = ">=0.98";
soh_bin_col(pairs_tbl.soh >= 0.95 & pairs_tbl.soh < 0.98) = "0.95-0.98";
soh_bin_col(pairs_tbl.soh >= 0.90 & pairs_tbl.soh < 0.95) = "0.90-0.95";
soh_bin_col(pairs_tbl.soh < 0.90) = "<0.90";

all_socs = sort(unique(pairs_tbl.soc(isfinite(pairs_tbl.soc))));
age_types = [1, 2, 3];

% 4-1. optemp_summary.csv
max_sum_rows = numel(COMPONENTS) * numel(VALID_T_AGE) * numel(age_types) * ...
               numel(SOH_BIN_NAMES) * numel(all_socs);
s_comp     = strings(max_sum_rows, 1);
s_T_age    = nan(max_sum_rows, 1);
s_age_type = nan(max_sum_rows, 1);
s_soh_bin  = strings(max_sum_rows, 1);
s_soc      = nan(max_sum_rows, 1);
s_n_pairs  = zeros(max_sum_rows, 1);
s_n_cells  = zeros(max_sum_rows, 1);
s_med      = nan(max_sum_rows, 1);
s_p25      = nan(max_sum_rows, 1);
s_p75      = nan(max_sum_rows, 1);
sr = 0;

for ci = 1:numel(COMPONENTS)
    comp = COMPONENTS{ci};
    ea_all = pairs_tbl.(['Ea_' comp]);
    for ti = 1:numel(VALID_T_AGE)
        ta = VALID_T_AGE(ti);
        m_ta = (pairs_tbl.T_age == ta);
        for ai = 1:numel(age_types)
            at = age_types(ai);
            m_at = m_ta & (pairs_tbl.age_type == at);
            if ~any(m_at), continue; end
            for bi = 1:numel(SOH_BIN_NAMES)
                bname = string(SOH_BIN_NAMES{bi});
                m_bin = m_at & (soh_bin_col == bname);
                if ~any(m_bin), continue; end
                for si = 1:numel(all_socs)
                    sc = all_socs(si);
                    m_sub = m_bin & (pairs_tbl.soc == sc) & isfinite(ea_all);
                    np_sub = sum(m_sub);
                    if np_sub == 0, continue; end

                    ea_sub = ea_all(m_sub);
                    sr = sr + 1;
                    s_comp(sr)     = string(comp);
                    s_T_age(sr)    = ta;
                    s_age_type(sr) = at;
                    s_soh_bin(sr)  = bname;
                    s_soc(sr)      = sc;
                    s_n_pairs(sr)  = np_sub;
                    s_n_cells(sr)  = numel(unique(pairs_tbl.cell(m_sub)));
                    s_med(sr)      = median(ea_sub);
                    s_p25(sr)      = pctile_base(ea_sub, 25);
                    s_p75(sr)      = pctile_base(ea_sub, 75);
                end
            end
        end
    end
end

idx_s = 1:sr;
summary_tbl = table( ...
    s_comp(idx_s), s_T_age(idx_s), s_age_type(idx_s), s_soh_bin(idx_s), s_soc(idx_s), ...
    s_n_pairs(idx_s), s_n_cells(idx_s), s_med(idx_s), s_p25(idx_s), s_p75(idx_s), ...
    'VariableNames', {'component', 'T_age', 'age_type', 'soh_bin', 'soc', ...
                      'n_pairs', 'n_cells', 'median_Ea', 'p25_Ea', 'p75_Ea'});
writetable(summary_tbl, fullfile(out_dir_abs, 'optemp_summary.csv'));

% 4-2. optemp_trend.csv (셀 단위 집계로 유사반복 방지)
max_tr_rows = numel(COMPONENTS) * numel(VALID_T_AGE) * numel(all_socs);
t_comp      = strings(max_tr_rows, 1);
t_T_age     = nan(max_tr_rows, 1);
t_soc       = nan(max_tr_rows, 1);
t_n_cells   = zeros(max_tr_rows, 1);
t_med_rho   = nan(max_tr_rows, 1);
t_pool_rho  = nan(max_tr_rows, 1);
t_med_slope = nan(max_tr_rows, 1);
t_frac_pos  = nan(max_tr_rows, 1);
t_sign_p    = nan(max_tr_rows, 1);
tr = 0;

for ci = 1:numel(COMPONENTS)
    comp = COMPONENTS{ci};
    ea_all = pairs_tbl.(['Ea_' comp]);
    for ti = 1:numel(VALID_T_AGE)
        ta = VALID_T_AGE(ti);
        for si = 1:numel(all_socs)
            sc = all_socs(si);
            m_sub = (pairs_tbl.T_age == ta) & (pairs_tbl.soc == sc) & ...
                    isfinite(ea_all) & isfinite(pairs_tbl.soh);
            if ~any(m_sub), continue; end

            cells_sub = pairs_tbl.cell(m_sub);
            fade_sub  = 1.0 - pairs_tbl.soh(m_sub);
            ea_sub    = ea_all(m_sub);

            uc = unique(cells_sub);
            cell_slopes = nan(numel(uc), 1);
            cell_rhos   = nan(numel(uc), 1);
            nc_valid = 0;

            for k = 1:numel(uc)
                mc = (cells_sub == uc(k));
                if sum(mc) < 4
                    continue;
                end
                x = fade_sub(mc);
                y = ea_sub(mc);
                xc = x - mean(x);
                denom = sum(xc .^ 2);
                if denom <= 0
                    continue;
                end
                nc_valid = nc_valid + 1;
                cell_slopes(nc_valid) = sum(xc .* (y - mean(y))) / denom;
                cell_rhos(nc_valid)   = calc_spearman(x, y);
            end

            tr = tr + 1;
            t_comp(tr)     = string(comp);
            t_T_age(tr)    = ta;
            t_soc(tr)      = sc;
            t_n_cells(tr)  = nc_valid;
            t_pool_rho(tr) = calc_spearman(fade_sub, ea_sub);

            if nc_valid > 0
                sl = cell_slopes(1:nc_valid);
                rh = cell_rhos(1:nc_valid);
                t_med_rho(tr)   = median(rh, 'omitnan');
                t_med_slope(tr) = median(sl, 'omitnan');
                t_frac_pos(tr)  = mean(sl > 0);
                n_pos = sum(sl > 0);
                n_neg = sum(sl < 0);
                t_sign_p(tr) = binomial_sign_test_p(n_pos, n_neg);
            end
        end
    end
end

idx_t = 1:tr;
trend_tbl = table( ...
    t_comp(idx_t), t_T_age(idx_t), t_soc(idx_t), t_n_cells(idx_t), ...
    t_med_rho(idx_t), t_pool_rho(idx_t), t_med_slope(idx_t), t_frac_pos(idx_t), t_sign_p(idx_t), ...
    'VariableNames', {'component', 'T_age', 'soc', 'n_cells', ...
                      'spearman_rho_cell_median', 'spearman_rho_pooled', ...
                      'median_slope_Ea_vs_fade', 'frac_pos_slope', 'sign_test_p'});
writetable(trend_tbl, fullfile(out_dir_abs, 'optemp_trend.csv'));

% 4-3. 명령 창 요약 출력
fprintf('\n=== Pair counts per T_age ===\n');
fprintf('%8s %10s %10s\n', 'T_age', 'n_pairs', 'n_cells');
for ti = 1:numel(VALID_T_AGE)
    ta = VALID_T_AGE(ti);
    m_ta = (pairs_tbl.T_age == ta);
    fprintf('%8g %10d %10d\n', ta, sum(m_ta), numel(unique(pairs_tbl.cell(m_ta))));
end

fprintf('\n=== SOC 50 Median Ea [kJ/mol]: Fresh (SOH >= 0.98) vs Aged (SOH < 0.95) ===\n');
fprintf('%-10s %6s %16s %10s %16s %10s\n', ...
    'Component', 'T_age', 'Fresh_med_Ea', 'n_fresh', 'Aged_med_Ea', 'n_aged');
m_soc50 = abs(pairs_tbl.soc - 50) < 0.5;
for ci = 1:numel(COMPONENTS)
    comp = COMPONENTS{ci};
    ea_all = pairs_tbl.(['Ea_' comp]);
    for ti = 1:numel(VALID_T_AGE)
        ta = VALID_T_AGE(ti);
        m_base  = m_soc50 & (pairs_tbl.T_age == ta) & isfinite(ea_all) & isfinite(pairs_tbl.soh);
        m_fresh = m_base & (pairs_tbl.soh >= 0.98);
        m_aged  = m_base & (pairs_tbl.soh < 0.95);
        med_f   = median(ea_all(m_fresh), 'omitnan');
        med_a   = median(ea_all(m_aged), 'omitnan');
        fprintf('%-10s %6g %16.2f %10d %16.2f %10d\n', ...
            comp, ta, med_f, sum(m_fresh), med_a, sum(m_aged));
    end
end

%% 5 그림
T_COLORS = [0.00, 0.45, 0.70;   % 0 degC  blue
            0.34, 0.71, 0.91;   % 10 degC light blue
            0.84, 0.37, 0.00];  % 40 degC orange

% 5-1. fig_optemp_Ea_vs_soh.png (2x2 tiles, SOC 50)
fig1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1080, 800]);
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Apparent Activation Energy E_a vs SOH (SOC 50%)', 'FontWeight', 'bold', 'FontSize', 13);

p50 = pairs_tbl(m_soc50 & isfinite(pairs_tbl.soh), :);
soh_edges = [0.00, 0.70:0.02:0.98, 1.05];
leg_handles = gobjects(numel(VALID_T_AGE), 1);
leg_labels  = arrayfun(@(t) sprintf('T_{age} = %d ^\\circC', t), VALID_T_AGE, 'UniformOutput', false);

for ci = 1:numel(COMPONENTS)
    ax = nexttile(tl1);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    comp = COMPONENTS{ci};
    ea_vec = p50.(['Ea_' comp]);
    soh_vec = p50.soh;

    % 산점도
    for ti = 1:numel(VALID_T_AGE)
        ta = VALID_T_AGE(ti);
        col = T_COLORS(ti, :);
        m = (p50.T_age == ta) & isfinite(ea_vec);
        if any(m)
            scatter(ax, soh_vec(m), ea_vec(m), 12, col, 'filled', ...
                'MarkerFaceAlpha', 0.22, 'MarkerEdgeAlpha', 0.22, 'HandleVisibility', 'off');
        end
    end

    % T_age별 구간 중앙값 선
    for ti = 1:numel(VALID_T_AGE)
        ta = VALID_T_AGE(ti);
        col = T_COLORS(ti, :);
        m = (p50.T_age == ta) & isfinite(ea_vec);
        bx = [];
        by = [];
        if any(m)
            sx = soh_vec(m);
            sy = ea_vec(m);
            for bi = 1:(numel(soh_edges) - 1)
                if bi == numel(soh_edges) - 1
                    mb = (sx >= soh_edges(bi)) & (sx <= soh_edges(bi + 1));
                else
                    mb = (sx >= soh_edges(bi)) & (sx < soh_edges(bi + 1));
                end
                if sum(mb) >= 3
                    bx(end+1) = median(sx(mb)); %#ok<SAGROW>
                    by(end+1) = median(sy(mb)); %#ok<SAGROW>
                end
            end
        end
        if ~isempty(bx)
            [bx, ord] = sort(bx, 'descend');
            by = by(ord);
            h = plot(ax, bx, by, '-o', 'Color', col, 'LineWidth', 2.2, ...
                'MarkerSize', 4.5, 'MarkerFaceColor', col, 'DisplayName', leg_labels{ti});
        else
            h = plot(ax, NaN, NaN, '-o', 'Color', col, 'LineWidth', 2.2, ...
                'MarkerSize', 4.5, 'MarkerFaceColor', col, 'DisplayName', leg_labels{ti});
        end
        if ci == 1
            leg_handles(ti) = h;
        end
    end

    set(ax, 'XDir', 'reverse');
    xlabel(ax, 'SOH');
    ylabel(ax, 'E_a [kJ/mol]');
    title(ax, COMP_LABELS{ci}, 'FontWeight', 'bold');

    valid_ea = ea_vec(isfinite(ea_vec));
    if numel(valid_ea) >= 5
        p01 = pctile_base(valid_ea, 1);
        p99 = pctile_base(valid_ea, 99);
        span = p99 - p01;
        if span > 0
            ylim(ax, [min(0, p01 - 0.10 * span), p99 + 0.12 * span]);
        end
    end
    if ~isempty(soh_vec)
        xlim(ax, [max(0.50, min(soh_vec) - 0.01), min(1.03, max(soh_vec) + 0.01)]);
    end
end

lgd = legend(leg_handles, leg_labels, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
exportgraphics(fig1, fullfile(out_dir_abs, 'fig_optemp_Ea_vs_soh.png'), ...
    'Resolution', 150, 'BackgroundColor', 'white');
close(fig1);

% 5-2. fig_optemp_nyquist_example.png (10 degC cycling & 40 degC cycling, SOC 50)
ex_T_list = [10, 40];
fig2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1120, 480]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Example Nyquist Spectra (SOC 50%): Solid = RT, Dashed = op-T', ...
    'FontWeight', 'bold', 'FontSize', 12);

col_fresh = [0.00, 0.45, 0.70];
col_aged  = [0.84, 0.37, 0.00];

for ei = 1:numel(ex_T_list)
    ta = ex_T_list(ei);
    ax = nexttile(tl2);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    m_grp = m_soc50 & (pairs_tbl.T_age == ta) & (pairs_tbl.age_type == 2);
    sub_grp = pairs_tbl(m_grp, :);

    if isempty(sub_grp)
        title(ax, sprintf('No %d degC cycling SOC-50 pairs', ta));
        continue;
    end

    % 해당 그룹에서 SOC-50 쌍 수가 가장 많은 셀 선택
    [g_idx, g_cells] = findgroups(sub_grp.cell);
    counts = accumarray(g_idx, 1);
    [~, best_g] = max(counts);
    sel_cell = g_cells(best_g);

    c_pairs = sub_grp(sub_grp.cell == sel_cell, :);
    c_pairs = sortrows(c_pairs, 't_rt', 'ascend');

    idx_first = 1;
    if height(c_pairs) == 1
        idx_aged = 1;
    else
        soh_rest = c_pairs.soh(2:end);
        if any(isfinite(soh_rest))
            [~, rel_idx] = min(abs(soh_rest - 0.90));
            idx_aged = rel_idx + 1;
        else
            idx_aged = height(c_pairs);
        end
    end

    p_f = c_pairs(idx_first, :);
    p_a = c_pairs(idx_aged, :);

    target_blocks = unique([p_f.block_rt, p_f.block_op, p_a.block_rt, p_a.block_op]);
    cell_csv = fullfile(raw_dir_abs, sprintf('cell_eisv2_%s.csv', char(sel_cell)));
    blk_map = load_blocks_from_csv(cell_csv, target_blocks, REQ_COLS);

    yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');

    % Fresh pair
    if isKey(blk_map, p_f.block_rt)
        sp = blk_map(p_f.block_rt);
        plot(ax, sp.re, sp.mim, '-', 'Color', col_fresh, 'LineWidth', 1.7, ...
            'DisplayName', sprintf('First RT (%.1f °C, SOH=%.3f)', p_f.T_meas_rt, p_f.soh));
    end
    if isKey(blk_map, p_f.block_op)
        sp = blk_map(p_f.block_op);
        plot(ax, sp.re, sp.mim, '--', 'Color', col_fresh, 'LineWidth', 1.7, ...
            'DisplayName', sprintf('First op-T (%.1f °C, SOH=%.3f)', p_f.T_meas_op, p_f.soh));
    end

    % Aged pair
    if idx_aged ~= idx_first
        if isKey(blk_map, p_a.block_rt)
            sp = blk_map(p_a.block_rt);
            plot(ax, sp.re, sp.mim, '-', 'Color', col_aged, 'LineWidth', 1.7, ...
                'DisplayName', sprintf('Aged RT (%.1f °C, SOH=%.3f)', p_a.T_meas_rt, p_a.soh));
        end
        if isKey(blk_map, p_a.block_op)
            sp = blk_map(p_a.block_op);
            plot(ax, sp.re, sp.mim, '--', 'Color', col_aged, 'LineWidth', 1.7, ...
                'DisplayName', sprintf('Aged op-T (%.1f °C, SOH=%.3f)', p_a.T_meas_op, p_a.soh));
        end
    end

    axis(ax, 'equal');
    xlabel(ax, 'Re(Z) [m\Omega]');
    ylabel(ax, '-Im(Z) [m\Omega]');
    title(ax, sprintf('%s | T_age=%d degC (cycling) | SOH %.3f -> %.3f', ...
        char(sel_cell), ta, p_f.soh, p_a.soh), 'Interpreter', 'none', 'FontWeight', 'bold');
    legend(ax, 'Location', 'northwest', 'FontSize', 8);
end

exportgraphics(fig2, fullfile(out_dir_abs, 'fig_optemp_nyquist_example.png'), ...
    'Resolution', 150, 'BackgroundColor', 'white');
close(fig2);

fprintf('\nDone. Outputs written to %s\n', out_dir_abs);

%% Local functions
function p_abs = resolve_rel_path(base_dir, p_in)
if isempty(p_in)
    p_abs = base_dir;
    return;
end
% 절대 경로 판별 (Windows 드라이브 문자/UNC 또는 Unix 루트)
is_abs = startsWith(p_in, filesep) || startsWith(p_in, '/') || ...
         startsWith(p_in, '\\') || (~isempty(regexp(p_in, '^[A-Za-z]:', 'once')));
if is_abs
    p_abs = p_in;
else
    p_norm = strrep(strrep(p_in, '\', filesep), '/', filesep);
    p_abs = fullfile(base_dir, p_norm);
end
end

function T = read_csv_clean(fpath)
opts = detectImportOptions(fpath, 'FileType', 'text', 'Encoding', 'UTF-8');
opts.VariableNamingRule = 'preserve';
T = readtable(fpath, opts);
vnames = T.Properties.VariableNames;
for k = 1:numel(vnames)
    vn = vnames{k};
    vn = erase(vn, char(65279));
    vnames{k} = strtrim(vn);
end
T.Properties.VariableNames = vnames;
end

function soh_map = build_soh_map(soh_tbl)
soh_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
u_cells = unique(soh_tbl.cell);
for i = 1:numel(u_cells)
    cid = u_cells(i);
    sub = soh_tbl(soh_tbl.cell == cid, :);
    valid = isfinite(sub.t);
    if ~any(valid), continue; end
    t_v   = sub.t(valid);
    day_v = sub.day(valid);
    efc_v = sub.efc(valid);
    soh_v = sub.soh(valid);

    [t_u, ~, ic] = unique(t_v, 'sorted');
    if numel(t_u) < numel(t_v)
        day_u = accumarray(ic, day_v, [], @(x) mean(x, 'omitnan'));
        efc_u = accumarray(ic, efc_v, [], @(x) mean(x, 'omitnan'));
        soh_u = accumarray(ic, soh_v, [], @(x) mean(x, 'omitnan'));
    else
        [~, ord] = sort(t_v);
        t_u   = t_v(ord);
        day_u = day_v(ord);
        efc_u = efc_v(ord);
        soh_u = soh_v(ord);
    end
    soh_map(char(cid)) = struct('t', t_u, 'day', day_u, 'efc', efc_u, 'soh', soh_u);
end
end

function yq = interp_limited(t_ref, y_ref, tq, max_extrap_s)
yq = nan(size(tq));
valid = isfinite(t_ref) & isfinite(y_ref);
t_v = t_ref(valid);
y_v = y_ref(valid);
n = numel(t_v);
if n == 0 || isempty(tq)
    return;
end
if n == 1
    in_rng = isfinite(tq) & (abs(tq - t_v(1)) <= max_extrap_s);
    yq(in_rng) = y_v(1);
else
    fin_q = isfinite(tq);
    if any(fin_q)
        vals = interp1(t_v, y_v, tq(fin_q), 'linear', 'extrap');
        oob = (tq(fin_q) < (t_v(1) - max_extrap_s)) | (tq(fin_q) > (t_v(end) + max_extrap_s));
        vals(oob) = NaN;
        yq(fin_q) = vals;
    end
end
end

function p = pctile_base(x, q)
x = x(isfinite(x));
n = numel(x);
if n == 0
    p = nan(size(q));
    return;
end
x = sort(x(:));
if n == 1
    p = repmat(x(1), size(q));
    return;
end
pos = (q / 100) * n + 0.5;
pos = max(1, min(n, pos));
lo = floor(pos);
hi = ceil(pos);
frac = pos - lo;
p = x(lo) + frac .* (x(hi) - x(lo));
end

function r = calc_spearman(x, y)
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if numel(x) < 3 || all(x == x(1)) || all(y == y(1))
    r = NaN;
    return;
end
try
    r = corr(x(:), y(:), 'Type', 'Spearman');
catch
    rx = manual_tiedrank(x(:));
    ry = manual_tiedrank(y(:));
    dx = rx - mean(rx);
    dy = ry - mean(ry);
    denom = sqrt(sum(dx .^ 2) * sum(dy .^ 2));
    if denom > 0
        r = sum(dx .* dy) / denom;
    else
        r = NaN;
    end
end
end

function r = manual_tiedrank(v)
n = numel(v);
[vs, ord] = sort(v(:));
r_sorted = zeros(n, 1);
i = 1;
while i <= n
    j = i;
    while j < n && vs(j + 1) == vs(i)
        j = j + 1;
    end
    r_sorted(i:j) = 0.5 * (i + j);
    i = j + 1;
end
r = zeros(n, 1);
r(ord) = r_sorted;
end

function p = binomial_sign_test_p(n_pos, n_neg)
n_eff = n_pos + n_neg;
if n_eff == 0
    p = NaN;
    return;
end
k_min = min(n_pos, n_neg);
ks = (0:k_min)';
log_pmf = gammaln(n_eff + 1) - gammaln(ks + 1) - gammaln(n_eff - ks + 1) - n_eff * log(2);
p = min(1.0, 2.0 * sum(exp(log_pmf)));
end

function blk_map = load_blocks_from_csv(fpath, target_blocks, req_cols)
blk_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
if ~exist(fpath, 'file')
    return;
end
opts = detectImportOptions(fpath, 'FileType', 'text', 'Delimiter', ';');
opts.SelectedVariableNames = req_cols;
opts = setvartype(opts, req_cols, 'double');
opts = setvaropts(opts, req_cols, 'TreatAsMissing', {'nan', 'NaN', 'NA', ''});
Traw = readtable(fpath, opts);

mask = (Traw.valid == 1) & ...
       ismember(Traw.sd_block_id, target_blocks) & ...
       isfinite(Traw.freq_Hz) & ...
       isfinite(Traw.z_amp_mOhm) & ...
       isfinite(Traw.z_ph_deg);

Tsub = Traw(mask, :);
for k = 1:numel(target_blocks)
    bid = target_blocks(k);
    mb = (Tsub.sd_block_id == bid);
    if ~any(mb), continue; end
    fb = Tsub.freq_Hz(mb);
    ampb = Tsub.z_amp_mOhm(mb);
    phb = Tsub.z_ph_deg(mb);
    [fb, ord] = sort(fb, 'descend');
    ampb = ampb(ord);
    phb  = phb(ord);
    reb  = ampb .* cos(deg2rad(phb));
    mimb = ampb .* sin(deg2rad(phb));
    blk_map(bid) = struct('f', fb, 're', reb, 'mim', mimb);
end
end