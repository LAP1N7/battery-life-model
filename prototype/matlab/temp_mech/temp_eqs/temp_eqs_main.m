%% 0 설정
clear; clc;

out_dir = 'out';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

log_path = fullfile(out_dir, 'log.txt');
log_fid = fopen(log_path, 'w');
if log_fid < 0
    error('로그 파일을 열 수 없습니다: %s', log_path);
end
cleanup_log = onCleanup(@() close_if_open(log_fid));

t_start_all = tic;
log_printf(log_fid, '=== temp_eqs_main 시작 (%s) ===\n', char(datetime('now')));

% 입력 파일 경로
KIT_PATH  = fullfile('..', 'kit_models', 'kit_rt_soc50.csv');
MOH_PATH  = fullfile('..', 'mohtat_modes', 'mohtat_modes_joined.csv');
MREF_PATH = 'mohtat_ref_resistance.csv';

% 물리 상수
R_GAS = 8.314462618;       % J/(mol*K)
T_REF_K = 298.15;          % 25 C in K
X_REF = 1000 / T_REF_K;    % 1000 / 298.15

% 온도별 색상표
TEMP_LIST = [-5, 0, 10, 25, 40, 45];
TEMP_COLORS = [ ...
    0.00, 0.30, 0.60; ...  % -5 C
    0.00, 0.45, 0.70; ...  %  0 C
    0.35, 0.70, 0.90; ...  % 10 C
    0.00, 0.62, 0.45; ...  % 25 C
    0.90, 0.60, 0.00; ...  % 40 C
    0.84, 0.37, 0.00  ...  % 45 C
];

% 데이터 로드
kit_raw = table();
moh_raw = table();
mref_raw = table();

try
    log_printf(log_fid, '[0] 데이터 로딩 중...\n');
    kit_raw  = read_csv_clean(KIT_PATH);
    moh_raw  = read_csv_clean(MOH_PATH);
    mref_raw = read_csv_clean(MREF_PATH);

    kit_raw.cell  = string(kit_raw.cell);
    moh_raw.cell  = string(moh_raw.cell);
    moh_raw.T     = string(moh_raw.T);
    mref_raw.cell = string(mref_raw.cell);
    mref_raw.T    = string(mref_raw.T);

    log_printf(log_fid, '    KIT: %d rows, %d cells\n', height(kit_raw), numel(unique(kit_raw.cell)));
    log_printf(log_fid, '    MOH: %d rows, %d cells\n', height(moh_raw), numel(unique(moh_raw.cell)));
    log_printf(log_fid, '    MREF: %d rows\n', height(mref_raw));
catch ME
    log_error(log_fid, '0 설정/로딩', ME);
    return;
end

% 섹션 간 공유 변수 초기화
K_cells_tbl = table();
K_groups_tbl = table();
K_paired_tbl = table();
S_cells_tbl = table();
S_group_tbl = table();
S_arrh_tbl = table();
S_fit_store = struct();
P_moh_interp_tbl = table();
P_partial_tbl = table();
P_collin_tbl = table();

%% K KIT 동일 용량 손실(equal-loss) 비교
try
    log_printf(log_fid, '\n=== [K] KIT equal-loss 비교 시작 ===\n');

    % KIT 셀별 기준점(최소 day) 계산 및 L, dR0, dRarc 산출
    kit_proc = preprocess_kit(kit_raw);

    % 사이클링(age_type 2,3): L = 5, 10, 20 %
    % 캘린더(age_type 1): L = 5 %
    cyc_levels = [5, 10, 20];
    cal_levels = 5;

    u_cells = unique(kit_proc.cell, 'stable');
    n_u = numel(u_cells);

    rec_cell = strings(0,1);
    rec_age_group = strings(0,1);
    rec_age_type = zeros(0,1);
    rec_T = zeros(0,1);
    rec_proto = strings(0,1);
    rec_level = zeros(0,1);
    rec_R0_ref = zeros(0,1);
    rec_Rarc_ref = zeros(0,1);
    rec_dR0 = zeros(0,1);
    rec_dRarc = zeros(0,1);
    rec_ratio = zeros(0,1);
    rec_rel_dR0 = zeros(0,1);
    rec_rel_dRarc = zeros(0,1);

    for i = 1:n_u
        cid = u_cells(i);
        sub = kit_proc(kit_proc.cell == cid, :);
        sub = sortrows(sub, 'day', 'ascend');
        at = sub.age_type(1);

        if ismember(at, [2, 3])
            ag = "cycling";
            lvls = cyc_levels;
        elseif at == 1
            ag = "calendar";
            lvls = cal_levels;
        else
            continue;
        end

        T_age = sub.T_age(1);
        R0_ref = sub.R0_ref(1);
        Rarc_ref = sub.Rarc_ref(1);
        pkey = make_proto_key(at, sub.chg_rate(1), sub.dis_rate(1), ...
            sub.v_min(1), sub.v_max(1), sub.cal_soc(1));

        Y_interp = interp_cell_cummax(sub.L, [sub.dR0, sub.dRarc], lvls);

        for k = 1:numel(lvls)
            dr0 = Y_interp(k, 1);
            dra = Y_interp(k, 2);
            if ~(isfinite(dr0) && isfinite(dra))
                continue;
            end
            if dr0 > 0
                rat = dra / dr0;
            else
                rat = NaN;
            end
            rec_cell(end+1,1) = cid; %#ok<SAGROW>
            rec_age_group(end+1,1) = ag; %#ok<SAGROW>
            rec_age_type(end+1,1) = at; %#ok<SAGROW>
            rec_T(end+1,1) = T_age; %#ok<SAGROW>
            rec_proto(end+1,1) = pkey; %#ok<SAGROW>
            rec_level(end+1,1) = lvls(k); %#ok<SAGROW>
            rec_R0_ref(end+1,1) = R0_ref; %#ok<SAGROW>
            rec_Rarc_ref(end+1,1) = Rarc_ref; %#ok<SAGROW>
            rec_dR0(end+1,1) = dr0; %#ok<SAGROW>
            rec_dRarc(end+1,1) = dra; %#ok<SAGROW>
            rec_ratio(end+1,1) = rat; %#ok<SAGROW>
            rec_rel_dR0(end+1,1) = dr0 / R0_ref; %#ok<SAGROW>
            rec_rel_dRarc(end+1,1) = dra / Rarc_ref; %#ok<SAGROW>
        end
    end

    K_cells_tbl = table(rec_cell, rec_age_group, rec_age_type, rec_T, rec_proto, ...
        rec_level, rec_R0_ref, rec_Rarc_ref, rec_dR0, rec_dRarc, rec_ratio, ...
        rec_rel_dR0, rec_rel_dRarc, ...
        'VariableNames', {'cell','age_group','age_type','T_age','proto_key', ...
        'level','R0_ref','Rarc_ref','dR0','dRarc','ratio','rel_dR0','rel_dRarc'});

    writetable(K_cells_tbl, fullfile(out_dir, 'K_equal_loss_cells.csv'));

    % 그룹별 요약 (age_group x T_age x level)
    ag_list = ["cycling", "calendar"];
    g_ag = strings(0,1);
    g_T = zeros(0,1);
    g_lvl = zeros(0,1);
    g_n = zeros(0,1);
    g_n_ratio = zeros(0,1);
    g_med_dR0 = zeros(0,1);
    g_med_dRarc = zeros(0,1);
    g_med_ratio = zeros(0,1);
    g_med_rel_dR0 = zeros(0,1);
    g_med_rel_dRarc = zeros(0,1);

    for a = 1:numel(ag_list)
        ag = ag_list(a);
        sub_ag = K_cells_tbl(K_cells_tbl.age_group == ag, :);
        if isempty(sub_ag)
            continue;
        end
        lvls = unique(sub_ag.level);
        Ts = unique(sub_ag.T_age);
        for l_idx = 1:numel(lvls)
            lv = lvls(l_idx);
            for t_idx = 1:numel(Ts)
                Tv = Ts(t_idx);
                sub = sub_ag(sub_ag.level == lv & sub_ag.T_age == Tv, :);
                if isempty(sub)
                    continue;
                end
                g_ag(end+1,1) = ag; %#ok<SAGROW>
                g_T(end+1,1) = Tv; %#ok<SAGROW>
                g_lvl(end+1,1) = lv; %#ok<SAGROW>
                g_n(end+1,1) = height(sub); %#ok<SAGROW>
                g_n_ratio(end+1,1) = sum(isfinite(sub.ratio)); %#ok<SAGROW>
                g_med_dR0(end+1,1) = median(sub.dR0, 'omitnan'); %#ok<SAGROW>
                g_med_dRarc(end+1,1) = median(sub.dRarc, 'omitnan'); %#ok<SAGROW>
                g_med_ratio(end+1,1) = median(sub.ratio, 'omitnan'); %#ok<SAGROW>
                g_med_rel_dR0(end+1,1) = median(sub.rel_dR0, 'omitnan'); %#ok<SAGROW>
                g_med_rel_dRarc(end+1,1) = median(sub.rel_dRarc, 'omitnan'); %#ok<SAGROW>
            end
        end
    end

    K_groups_tbl = table(g_ag, g_T, g_lvl, g_n, g_n_ratio, ...
        g_med_dR0, g_med_dRarc, g_med_ratio, g_med_rel_dR0, g_med_rel_dRarc, ...
        'VariableNames', {'age_group','T_age','level','n','n_ratio', ...
        'med_dR0','med_dRarc','med_ratio','med_rel_dR0','med_rel_dRarc'});

    writetable(K_groups_tbl, fullfile(out_dir, 'K_equal_loss_groups.csv'));

    log_printf(log_fid, '[K] KIT equal-loss 그룹 요약:\n');
    log_printf(log_fid, '%-9s %4s %5s %4s %5s %8s %8s %8s %9s %9s\n', ...
        'age_group','T_C','L(%)','n','n_rat','med_dR0','med_dRa','med_rat','rel_dR0','rel_dRa');
    for i = 1:height(K_groups_tbl)
        log_printf(log_fid, '%-9s %4.0f %5.0f %4d %5d %8.3f %8.3f %8.3f %9.4f %9.4f\n', ...
            K_groups_tbl.age_group(i), K_groups_tbl.T_age(i), K_groups_tbl.level(i), ...
            K_groups_tbl.n(i), K_groups_tbl.n_ratio(i), K_groups_tbl.med_dR0(i), ...
            K_groups_tbl.med_dRarc(i), K_groups_tbl.med_ratio(i), ...
            K_groups_tbl.med_rel_dR0(i), K_groups_tbl.med_rel_dRarc(i));
    end

    % 쌍체(paired) 비교: 사이클링 셀 기준 0 vs 25, 10 vs 25, 40 vs 25
    cyc_tbl = K_cells_tbl(K_cells_tbl.age_group == "cycling", :);
    cmp_temps = [0, 10, 40];
    T_base = 25;

    p_cmp = strings(0,1);
    p_Tcmp = zeros(0,1);
    p_Tref = zeros(0,1);
    p_lvl = zeros(0,1);
    p_n = zeros(0,1);
    p_n_Tcmp_gt_25 = zeros(0,1);
    p_n_cold_gt_warm = zeros(0,1);
    p_med_diff_Tcmp_25 = zeros(0,1);
    p_med_diff_cold_warm = zeros(0,1);
    p_pval = zeros(0,1);

    for c = 1:numel(cmp_temps)
        Tc = cmp_temps(c);
        for l_idx = 1:numel(cyc_levels)
            lv = cyc_levels(l_idx);
            sub_c = cyc_tbl(cyc_tbl.T_age == Tc & cyc_tbl.level == lv & isfinite(cyc_tbl.ratio), :);
            sub_b = cyc_tbl(cyc_tbl.T_age == T_base & cyc_tbl.level == lv & isfinite(cyc_tbl.ratio), :);

            keys_common = intersect(unique(sub_c.proto_key), unique(sub_b.proto_key));
            n_pairs = numel(keys_common);

            if n_pairs > 0
                r_c = zeros(n_pairs, 1);
                r_b = zeros(n_pairs, 1);
                for j = 1:n_pairs
                    kj = keys_common(j);
                    r_c(j) = median(sub_c.ratio(sub_c.proto_key == kj), 'omitnan');
                    r_b(j) = median(sub_b.ratio(sub_b.proto_key == kj), 'omitnan');
                end
                diff_Tcmp_25 = r_c - r_b;
                if Tc < T_base
                    diff_cold_warm = r_c - r_b;
                else
                    diff_cold_warm = r_b - r_c;
                end
                n_Tcmp_gt = sum(diff_Tcmp_25 > 0);
                n_cold_gt = sum(diff_cold_warm > 0);
                n_nz = sum(diff_Tcmp_25 ~= 0);
                pval = my_signtest(sum(diff_Tcmp_25 > 0), n_nz);
                med_d1 = median(diff_Tcmp_25, 'omitnan');
                med_d2 = median(diff_cold_warm, 'omitnan');
            else
                n_Tcmp_gt = 0;
                n_cold_gt = 0;
                med_d1 = NaN;
                med_d2 = NaN;
                pval = NaN;
            end

            p_cmp(end+1,1) = sprintf('%d_vs_%d', Tc, T_base); %#ok<SAGROW>
            p_Tcmp(end+1,1) = Tc; %#ok<SAGROW>
            p_Tref(end+1,1) = T_base; %#ok<SAGROW>
            p_lvl(end+1,1) = lv; %#ok<SAGROW>
            p_n(end+1,1) = n_pairs; %#ok<SAGROW>
            p_n_Tcmp_gt_25(end+1,1) = n_Tcmp_gt; %#ok<SAGROW>
            p_n_cold_gt_warm(end+1,1) = n_cold_gt; %#ok<SAGROW>
            p_med_diff_Tcmp_25(end+1,1) = med_d1; %#ok<SAGROW>
            p_med_diff_cold_warm(end+1,1) = med_d2; %#ok<SAGROW>
            p_pval(end+1,1) = pval; %#ok<SAGROW>
        end
    end

    K_paired_tbl = table(p_cmp, p_Tcmp, p_Tref, p_lvl, p_n, ...
        p_n_Tcmp_gt_25, p_n_cold_gt_warm, p_med_diff_Tcmp_25, ...
        p_med_diff_cold_warm, p_pval, ...
        'VariableNames', {'comparison','T_cmp','T_ref','level','n_pairs', ...
        'n_Tcmp_gt_25','n_colder_gt_warmer','med_diff_Tcmp_minus_25', ...
        'med_diff_colder_minus_warmer','p_sign_two_sided'});

    writetable(K_paired_tbl, fullfile(out_dir, 'K_paired.csv'));

    log_printf(log_fid, '\n[K] Paired sign test (ratio = dRarc/dR0):\n');
    log_printf(log_fid, '%-10s %5s %5s %8s %8s %10s %10s %8s\n', ...
        'pair','L(%)','n','T>25','cold>hot','med(T-25)','med(c-h)','p_sign');
    for i = 1:height(K_paired_tbl)
        log_printf(log_fid, '%-10s %5.0f %5d %8d %8d %10.3f %10.3f %8.4f\n', ...
            K_paired_tbl.comparison(i), K_paired_tbl.level(i), K_paired_tbl.n_pairs(i), ...
            K_paired_tbl.n_Tcmp_gt_25(i), K_paired_tbl.n_colder_gt_warmer(i), ...
            K_paired_tbl.med_diff_Tcmp_minus_25(i), ...
            K_paired_tbl.med_diff_colder_minus_warmer(i), ...
            K_paired_tbl.p_sign_two_sided(i));
    end

catch ME
    log_error(log_fid, 'K equal-loss', ME);
end

%% S 온도의 함수로서의 기울기 s(T)
try
    log_printf(log_fid, '\n=== [S] 온도별 기울기 s(T) 및 Arrhenius 적합 시작 ===\n');

    if ~exist('kit_proc', 'var') || isempty(kit_proc)
        kit_proc = preprocess_kit(kit_raw);
    end

    % 1) 셀별 원점 통과 기울기 s_cell = sum(L*dR)/sum(L^2) (L > 0)
    sc_ds = strings(0,1);
    sc_cell = strings(0,1);
    sc_T = zeros(0,1);
    sc_n_R0 = zeros(0,1);
    sc_n_arc = zeros(0,1);
    sc_maxL = zeros(0,1);
    sc_R0_ref = zeros(0,1);
    sc_Rarc_ref = zeros(0,1);
    sc_s_R0 = zeros(0,1);
    sc_s_arc = zeros(0,1);
    sc_s_R0_norm = zeros(0,1);
    sc_s_arc_norm = zeros(0,1);
    sc_ratio_raw = zeros(0,1);
    sc_ratio_norm = zeros(0,1);

    % KIT: age_type == 2, L > 0 행 >= 3, max L >= 5
    kit_at2 = kit_proc(kit_proc.age_type == 2, :);
    u_kit = unique(kit_at2.cell, 'stable');
    for i = 1:numel(u_kit)
        cid = u_kit(i);
        sub = kit_at2(kit_at2.cell == cid, :);
        ok0 = isfinite(sub.L) & sub.L > 0 & isfinite(sub.dR0);
        oka = isfinite(sub.L) & sub.L > 0 & isfinite(sub.dRarc);
        ok_all = isfinite(sub.L) & sub.L > 0;
        if ~any(ok_all)
            continue;
        end
        maxL = max(sub.L(ok_all));
        if maxL < 5
            continue;
        end
        if sum(ok0) < 3 && sum(oka) < 3
            continue;
        end

        R0_ref = sub.R0_ref(1);
        Rarc_ref = sub.Rarc_ref(1);

        if sum(ok0) >= 3
            s0 = sum(sub.L(ok0) .* sub.dR0(ok0)) / sum(sub.L(ok0).^2);
            s0n = s0 / R0_ref;
        else
            s0 = NaN; s0n = NaN;
        end
        if sum(oka) >= 3
            sa = sum(sub.L(oka) .* sub.dRarc(oka)) / sum(sub.L(oka).^2);
            san = sa / Rarc_ref;
        else
            sa = NaN; san = NaN;
        end

        if isfinite(s0) && s0 > 0 && isfinite(sa) && sa > 0
            rr = sa / s0;
            rn = san / s0n;
        else
            rr = NaN; rn = NaN;
        end

        sc_ds(end+1,1) = "KIT"; %#ok<SAGROW>
        sc_cell(end+1,1) = cid; %#ok<SAGROW>
        sc_T(end+1,1) = sub.T_age(1); %#ok<SAGROW>
        sc_n_R0(end+1,1) = sum(ok0); %#ok<SAGROW>
        sc_n_arc(end+1,1) = sum(oka); %#ok<SAGROW>
        sc_maxL(end+1,1) = maxL; %#ok<SAGROW>
        sc_R0_ref(end+1,1) = R0_ref; %#ok<SAGROW>
        sc_Rarc_ref(end+1,1) = Rarc_ref; %#ok<SAGROW>
        sc_s_R0(end+1,1) = s0; %#ok<SAGROW>
        sc_s_arc(end+1,1) = sa; %#ok<SAGROW>
        sc_s_R0_norm(end+1,1) = s0n; %#ok<SAGROW>
        sc_s_arc_norm(end+1,1) = san; %#ok<SAGROW>
        sc_ratio_raw(end+1,1) = rr; %#ok<SAGROW>
        sc_ratio_norm(end+1,1) = rn; %#ok<SAGROW>
    end

    % Mohtat: 모든 셀, 사용 가능 행(finite dR, loss > 0) >= 3, MREF로 정규화
    u_moh = unique(moh_raw.cell, 'stable');
    for i = 1:numel(u_moh)
        cid = u_moh(i);
        sub = moh_raw(moh_raw.cell == cid, :);
        T_str = sub.T(1);
        T_val = sub.T_C(1);

        m_idx = find(mref_raw.cell == cid & mref_raw.T == T_str, 1);
        if isempty(m_idx)
            m_idx = find(mref_raw.cell == cid, 1);
        end
        if isempty(m_idx)
            continue;
        end
        R0_ref = mref_raw.R0_ref(m_idx);
        Rarc_ref = mref_raw.Rarc_ref(m_idx);

        ok0 = isfinite(sub.loss) & sub.loss > 0 & isfinite(sub.dR0);
        oka = isfinite(sub.loss) & sub.loss > 0 & isfinite(sub.dRarc);
        if sum(ok0) < 3 && sum(oka) < 3
            continue;
        end

        maxL = max(sub.loss(isfinite(sub.loss)));
        if sum(ok0) >= 3
            s0 = sum(sub.loss(ok0) .* sub.dR0(ok0)) / sum(sub.loss(ok0).^2);
            s0n = s0 / R0_ref;
        else
            s0 = NaN; s0n = NaN;
        end
        if sum(oka) >= 3
            sa = sum(sub.loss(oka) .* sub.dRarc(oka)) / sum(sub.loss(oka).^2);
            san = sa / Rarc_ref;
        else
            sa = NaN; san = NaN;
        end

        if isfinite(s0) && s0 > 0 && isfinite(sa) && sa > 0
            rr = sa / s0;
            rn = san / s0n;
        else
            rr = NaN; rn = NaN;
        end

        sc_ds(end+1,1) = "Mohtat"; %#ok<SAGROW>
        sc_cell(end+1,1) = cid; %#ok<SAGROW>
        sc_T(end+1,1) = T_val; %#ok<SAGROW>
        sc_n_R0(end+1,1) = sum(ok0); %#ok<SAGROW>
        sc_n_arc(end+1,1) = sum(oka); %#ok<SAGROW>
        sc_maxL(end+1,1) = maxL; %#ok<SAGROW>
        sc_R0_ref(end+1,1) = R0_ref; %#ok<SAGROW>
        sc_Rarc_ref(end+1,1) = Rarc_ref; %#ok<SAGROW>
        sc_s_R0(end+1,1) = s0; %#ok<SAGROW>
        sc_s_arc(end+1,1) = sa; %#ok<SAGROW>
        sc_s_R0_norm(end+1,1) = s0n; %#ok<SAGROW>
        sc_s_arc_norm(end+1,1) = san; %#ok<SAGROW>
        sc_ratio_raw(end+1,1) = rr; %#ok<SAGROW>
        sc_ratio_norm(end+1,1) = rn; %#ok<SAGROW>
    end

    S_cells_tbl = table(sc_ds, sc_cell, sc_T, sc_n_R0, sc_n_arc, sc_maxL, ...
        sc_R0_ref, sc_Rarc_ref, sc_s_R0, sc_s_arc, sc_s_R0_norm, sc_s_arc_norm, ...
        sc_ratio_raw, sc_ratio_norm, ...
        'VariableNames', {'dataset','cell','T_C','n_rows_R0','n_rows_arc','max_L', ...
        'R0_ref','Rarc_ref','s_R0','s_arc','s_R0_norm','s_arc_norm', ...
        'ratio_raw','ratio_norm'});

    writetable(S_cells_tbl, fullfile(out_dir, 'S_cell_slopes.csv'));

    % 2) 데이터셋 x 온도별 요약 (median, IQR, n)
    ds_list = ["KIT", "Mohtat"];
    sg_ds = strings(0,1);
    sg_T = zeros(0,1);
    sg_n = zeros(0,1);
    sg_med_R0 = zeros(0,1); sg_iqr_R0 = zeros(0,1);
    sg_med_arc = zeros(0,1); sg_iqr_arc = zeros(0,1);
    sg_med_R0n_pct = zeros(0,1); sg_iqr_R0n_pct = zeros(0,1);
    sg_med_arcn_pct = zeros(0,1); sg_iqr_arcn_pct = zeros(0,1);
    sg_med_rat_norm = zeros(0,1); sg_iqr_rat_norm = zeros(0,1);

    for d = 1:numel(ds_list)
        ds = ds_list(d);
        sub_d = S_cells_tbl(S_cells_tbl.dataset == ds, :);
        Ts = sort(unique(sub_d.T_C));
        for t = 1:numel(Ts)
            Tv = Ts(t);
            sub = sub_d(sub_d.T_C == Tv, :);
            sg_ds(end+1,1) = ds; %#ok<SAGROW>
            sg_T(end+1,1) = Tv; %#ok<SAGROW>
            sg_n(end+1,1) = height(sub); %#ok<SAGROW>
            sg_med_R0(end+1,1) = median(sub.s_R0, 'omitnan'); %#ok<SAGROW>
            sg_iqr_R0(end+1,1) = my_iqr(sub.s_R0); %#ok<SAGROW>
            sg_med_arc(end+1,1) = median(sub.s_arc, 'omitnan'); %#ok<SAGROW>
            sg_iqr_arc(end+1,1) = my_iqr(sub.s_arc); %#ok<SAGROW>
            sg_med_R0n_pct(end+1,1) = 100 * median(sub.s_R0_norm, 'omitnan'); %#ok<SAGROW>
            sg_iqr_R0n_pct(end+1,1) = 100 * my_iqr(sub.s_R0_norm); %#ok<SAGROW>
            sg_med_arcn_pct(end+1,1) = 100 * median(sub.s_arc_norm, 'omitnan'); %#ok<SAGROW>
            sg_iqr_arcn_pct(end+1,1) = 100 * my_iqr(sub.s_arc_norm); %#ok<SAGROW>
            sg_med_rat_norm(end+1,1) = median(sub.ratio_norm, 'omitnan'); %#ok<SAGROW>
            sg_iqr_rat_norm(end+1,1) = my_iqr(sub.ratio_norm); %#ok<SAGROW>
        end
    end

    S_group_tbl = table(sg_ds, sg_T, sg_n, sg_med_R0, sg_iqr_R0, ...
        sg_med_arc, sg_iqr_arc, sg_med_R0n_pct, sg_iqr_R0n_pct, ...
        sg_med_arcn_pct, sg_iqr_arcn_pct, sg_med_rat_norm, sg_iqr_rat_norm, ...
        'VariableNames', {'dataset','T_C','n','med_s_R0','iqr_s_R0', ...
        'med_s_arc','iqr_s_arc','med_s_R0_pct_per_pct','iqr_s_R0_pct_per_pct', ...
        'med_s_arc_pct_per_pct','iqr_s_arc_pct_per_pct', ...
        'med_ratio_norm','iqr_ratio_norm'});

    writetable(S_group_tbl, fullfile(out_dir, 'S_group.csv'));

    log_printf(log_fid, '[S] 온도별 정규화 기울기 요약 (%%/% %%, median [IQR]):\n');
    log_printf(log_fid, '%-7s %5s %4s %18s %18s %18s\n', ...
        'dataset','T_C','n','s_R0 (%/%)','s_arc (%/%)','s_arc/s_R0');
    for i = 1:height(S_group_tbl)
        log_printf(log_fid, '%-7s %5.0f %4d %8.3f [%7.3f] %8.3f [%7.3f] %8.3f [%7.3f]\n', ...
            S_group_tbl.dataset(i), S_group_tbl.T_C(i), S_group_tbl.n(i), ...
            S_group_tbl.med_s_R0_pct_per_pct(i), S_group_tbl.iqr_s_R0_pct_per_pct(i), ...
            S_group_tbl.med_s_arc_pct_per_pct(i), S_group_tbl.iqr_s_arc_pct_per_pct(i), ...
            S_group_tbl.med_ratio_norm(i), S_group_tbl.iqr_ratio_norm(i));
    end

    % 3) Arrhenius 적합 (A1, A2) + 셀 부트스트랩 (1000회, rng(1))
    log_printf(log_fid, '\n[S] 참고: Mohtat은 3개 온도(-5/25/45 C)뿐이므로 4-파라미터 A2 모델은 식별 불가 -> Mohtat은 A1만 적합.\n');

    n_boot = 1000;
    targets = {"R0", "arc", "ratio"};
    target_cols = {"s_R0_norm", "s_arc_norm", "ratio_norm"};
    datasets_fit = {"KIT", "Mohtat", "Joint"};

    ar_ds = strings(0,1);
    ar_tgt = strings(0,1);
    ar_mod = strings(0,1);
    ar_n = zeros(0,1);
    ar_k = zeros(0,1);
    ar_sse = zeros(0,1);
    ar_aic = zeros(0,1);
    ar_daic = zeros(0,1);
    ar_sref_K = zeros(0,1);
    ar_sref_M = zeros(0,1);
    ar_a1_K = zeros(0,1);
    ar_a2_K = zeros(0,1);
    ar_a1_M = zeros(0,1);
    ar_a2_M = zeros(0,1);
    ar_Ea = zeros(0,1); ar_Ea_lo = zeros(0,1); ar_Ea_hi = zeros(0,1);
    ar_Eh = zeros(0,1); ar_Eh_lo = zeros(0,1); ar_Eh_hi = zeros(0,1);
    ar_Ec = zeros(0,1); ar_Ec_lo = zeros(0,1); ar_Ec_hi = zeros(0,1);
    ar_eq = strings(0,1);

    for ti = 1:numel(targets)
        tgt = targets{ti};
        col = target_cols{ti};

        for di = 1:numel(datasets_fit)
            ds_name = datasets_fit{di};

            if strcmp(ds_name, 'Joint')
                sub = S_cells_tbl;
            else
                sub = S_cells_tbl(S_cells_tbl.dataset == string(ds_name), :);
            end

            s_val = sub.(col);
            ok = isfinite(s_val) & s_val > 0 & isfinite(sub.T_C);
            sub = sub(ok, :);
            s_val = s_val(ok);

            n_obs = numel(s_val);
            if n_obs < 3
                continue;
            end

            x_val = 1000 ./ (sub.T_C + 273.15);
            y_val = log(s_val);
            ds_id = double(sub.dataset == "Mohtat") + 1; % 1=KIT, 2=Mohtat
            is_joint = strcmp(ds_name, 'Joint');

            % --- A1 적합 ---
            fitA1 = fit_A1(x_val, y_val, ds_id, is_joint, X_REF, R_GAS);

            % 부트스트랩 (dataset x temperature 층화 재표본)
            rng(1);
            strata = build_strata(sub.dataset, sub.T_C);
            boot_Ea = nan(n_boot, 1);
            for b = 1:n_boot
                idx_b = resample_strata(strata);
                fb = fit_A1(x_val(idx_b), y_val(idx_b), ds_id(idx_b), is_joint, X_REF, R_GAS);
                boot_Ea(b) = fb.Ea_kJ;
            end
            ci_Ea = my_quantile(boot_Ea, [0.025, 0.975]);

            [sref_K, sref_M, a1_K, a1_M] = unpack_A1_params(fitA1, ds_name);
            eq_A1 = format_eq_A1(tgt, ds_name, sref_K, sref_M, fitA1.Ea_kJ);

            % 저장
            S_fit_store.(sprintf('%s_%s_A1', ds_name, tgt)) = fitA1;

            ar_ds(end+1,1) = string(ds_name); %#ok<SAGROW>
            ar_tgt(end+1,1) = string(tgt); %#ok<SAGROW>
            ar_mod(end+1,1) = "A1"; %#ok<SAGROW>
            ar_n(end+1,1) = n_obs; %#ok<SAGROW>
            ar_k(end+1,1) = fitA1.k; %#ok<SAGROW>
            ar_sse(end+1,1) = fitA1.SSE; %#ok<SAGROW>
            ar_aic(end+1,1) = fitA1.AIC; %#ok<SAGROW>
            ar_daic(end+1,1) = 0; %#ok<SAGROW>
            ar_sref_K(end+1,1) = sref_K; %#ok<SAGROW>
            ar_sref_M(end+1,1) = sref_M; %#ok<SAGROW>
            ar_a1_K(end+1,1) = a1_K; %#ok<SAGROW>
            ar_a2_K(end+1,1) = NaN; %#ok<SAGROW>
            ar_a1_M(end+1,1) = a1_M; %#ok<SAGROW>
            ar_a2_M(end+1,1) = NaN; %#ok<SAGROW>
            ar_Ea(end+1,1) = fitA1.Ea_kJ; %#ok<SAGROW>
            ar_Ea_lo(end+1,1) = ci_Ea(1); %#ok<SAGROW>
            ar_Ea_hi(end+1,1) = ci_Ea(2); %#ok<SAGROW>
            ar_Eh(end+1,1) = NaN; ar_Eh_lo(end+1,1) = NaN; ar_Eh_hi(end+1,1) = NaN; %#ok<SAGROW>
            ar_Ec(end+1,1) = NaN; ar_Ec_lo(end+1,1) = NaN; ar_Ec_hi(end+1,1) = NaN; %#ok<SAGROW>
            ar_eq(end+1,1) = eq_A1; %#ok<SAGROW>

            % --- A2 적합 (R0, arc 대상; KIT 및 Joint만) ---
            if ~strcmp(tgt, 'ratio') && ~strcmp(ds_name, 'Mohtat')
                fitA2 = fit_A2_fast(x_val, y_val, ds_id, is_joint, X_REF, R_GAS, []);

                % 부트스트랩 (그룹 요약 기반 초고속 동일 목적함수 사용)
                rng(1);
                boot_Eh = nan(n_boot, 1);
                boot_Ec = nan(n_boot, 1);
                for b = 1:n_boot
                    idx_b = resample_strata(strata);
                    fb2 = fit_A2_fast(x_val(idx_b), y_val(idx_b), ds_id(idx_b), ...
                        is_joint, X_REF, R_GAS, fitA2.p_opt);
                    boot_Eh(b) = fb2.Eh_kJ;
                    boot_Ec(b) = fb2.Ec_kJ;
                end
                ci_Eh = my_quantile(boot_Eh, [0.025, 0.975]);
                ci_Ec = my_quantile(boot_Ec, [0.025, 0.975]);

                [sref2_K, sref2_M, a12_K, a22_K, a12_M, a22_M] = unpack_A2_params(fitA2, ds_name);
                eq_A2 = format_eq_A2(tgt, ds_name, a12_K, a22_K, a12_M, a22_M, fitA2.Eh_kJ, fitA2.Ec_kJ);

                S_fit_store.(sprintf('%s_%s_A2', ds_name, tgt)) = fitA2;

                ar_ds(end+1,1) = string(ds_name); %#ok<SAGROW>
                ar_tgt(end+1,1) = string(tgt); %#ok<SAGROW>
                ar_mod(end+1,1) = "A2"; %#ok<SAGROW>
                ar_n(end+1,1) = n_obs; %#ok<SAGROW>
                ar_k(end+1,1) = fitA2.k; %#ok<SAGROW>
                ar_sse(end+1,1) = fitA2.SSE; %#ok<SAGROW>
                ar_aic(end+1,1) = fitA2.AIC; %#ok<SAGROW>
                ar_daic(end+1,1) = fitA2.AIC - fitA1.AIC; %#ok<SAGROW>
                ar_sref_K(end+1,1) = sref2_K; %#ok<SAGROW>
                ar_sref_M(end+1,1) = sref2_M; %#ok<SAGROW>
                ar_a1_K(end+1,1) = a12_K; %#ok<SAGROW>
                ar_a2_K(end+1,1) = a22_K; %#ok<SAGROW>
                ar_a1_M(end+1,1) = a12_M; %#ok<SAGROW>
                ar_a2_M(end+1,1) = a22_M; %#ok<SAGROW>
                ar_Ea(end+1,1) = NaN; ar_Ea_lo(end+1,1) = NaN; ar_Ea_hi(end+1,1) = NaN; %#ok<SAGROW>
                ar_Eh(end+1,1) = fitA2.Eh_kJ; %#ok<SAGROW>
                ar_Eh_lo(end+1,1) = ci_Eh(1); %#ok<SAGROW>
                ar_Eh_hi(end+1,1) = ci_Eh(2); %#ok<SAGROW>
                ar_Ec(end+1,1) = fitA2.Ec_kJ; %#ok<SAGROW>
                ar_Ec_lo(end+1,1) = ci_Ec(1); %#ok<SAGROW>
                ar_Ec_hi(end+1,1) = ci_Ec(2); %#ok<SAGROW>
                ar_eq(end+1,1) = eq_A2; %#ok<SAGROW>
            end
        end
    end

    S_arrh_tbl = table(ar_ds, ar_tgt, ar_mod, ar_n, ar_k, ar_sse, ar_aic, ar_daic, ...
        ar_sref_K, ar_sref_M, ar_a1_K, ar_a2_K, ar_a1_M, ar_a2_M, ...
        ar_Ea, ar_Ea_lo, ar_Ea_hi, ar_Eh, ar_Eh_lo, ar_Eh_hi, ...
        ar_Ec, ar_Ec_lo, ar_Ec_hi, ar_eq, ...
        'VariableNames', {'dataset','target','model','n','k','SSE','AIC','dAIC_vs_A1', ...
        's_ref_KIT','s_ref_MOH','a1_KIT','a2_KIT','a1_MOH','a2_MOH', ...
        'Ea_kJ','Ea_CI_lo','Ea_CI_hi','Eh_kJ','Eh_CI_lo','Eh_CI_hi', ...
        'Ec_kJ','Ec_CI_lo','Ec_CI_hi','equation'});

    writetable(S_arrh_tbl, fullfile(out_dir, 'S_arrhenius.csv'));

    log_printf(log_fid, '\n[S] Arrhenius 적합 결과 요약 (단위: kJ/mol, 양수 Ea = 온도 상승 시 기울기 증가):\n');
    log_printf(log_fid, '%-7s %-6s %-3s %4s %8s %7s | %-28s | %-28s\n', ...
        'dataset','target','mod','n','AIC','dAIC','Ea or Eh [95% CI]','Ec [95% CI]');
    for i = 1:height(S_arrh_tbl)
        if S_arrh_tbl.model(i) == "A1"
            str1 = sprintf('Ea=%6.2f [%6.2f,%6.2f]', ...
                S_arrh_tbl.Ea_kJ(i), S_arrh_tbl.Ea_CI_lo(i), S_arrh_tbl.Ea_CI_hi(i));
            str2 = '-';
        else
            str1 = sprintf('Eh=%6.2f [%6.2f,%6.2f]', ...
                S_arrh_tbl.Eh_kJ(i), S_arrh_tbl.Eh_CI_lo(i), S_arrh_tbl.Eh_CI_hi(i));
            str2 = sprintf('Ec=%6.2f [%6.2f,%6.2f]', ...
                S_arrh_tbl.Ec_kJ(i), S_arrh_tbl.Ec_CI_lo(i), S_arrh_tbl.Ec_CI_hi(i));
        end
        log_printf(log_fid, '%-7s %-6s %-3s %4d %8.2f %7.2f | %-28s | %-28s\n', ...
            S_arrh_tbl.dataset(i), S_arrh_tbl.target(i), S_arrh_tbl.model(i), ...
            S_arrh_tbl.n(i), S_arrh_tbl.AIC(i), S_arrh_tbl.dAIC_vs_A1(i), str1, str2);
    end

    log_printf(log_fid, '\n[S] 최종 방정식 (정규화 기울기 s = s_cell/R_ref, 단위: fraction / %% loss):\n');
    for i = 1:height(S_arrh_tbl)
        log_printf(log_fid, '  [%s | %s | %s] %s\n', ...
            S_arrh_tbl.dataset(i), S_arrh_tbl.target(i), S_arrh_tbl.model(i), ...
            S_arrh_tbl.equation(i));
    end

catch ME
    log_error(log_fid, 'S slope/Arrhenius', ME);
end

%% P Mohtat 동일 온도 내 모드(LAM_PE, LAM_NE) 검정
try
    log_printf(log_fid, '\n=== [P] Mohtat 동일 온도 내 모드 편상관 및 공선성 분석 시작 ===\n');
    log_printf(log_fid, 'n per T is 7 cells at most, treat as exploratory.\n');

    % 1) 셀별 L = 5, 10 % 보간 (cummax loss)
    moh_lvls = [5, 10];
    u_moh = unique(moh_raw.cell, 'stable');

    pm_cell = strings(0,1);
    pm_T = strings(0,1);
    pm_TC = zeros(0,1);
    pm_lvl = zeros(0,1);
    pm_LLI = zeros(0,1);
    pm_LAMPE = zeros(0,1);
    pm_LAMNE = zeros(0,1);
    pm_dR0 = zeros(0,1);
    pm_dRarc = zeros(0,1);
    pm_ratio = zeros(0,1);

    for i = 1:numel(u_moh)
        cid = u_moh(i);
        sub = moh_raw(moh_raw.cell == cid, :);
        sub = sortrows(sub, {'days', 'efc'}, {'ascend', 'ascend'});

        Y_mat = [sub.LLI, sub.LAM_PE, sub.LAM_NE, sub.dR0, sub.dRarc];
        Y_int = interp_cell_cummax(sub.loss, Y_mat, moh_lvls);

        for k = 1:numel(moh_lvls)
            if all(isnan(Y_int(k, :)))
                continue;
            end
            dr0 = Y_int(k, 4);
            dra = Y_int(k, 5);
            if isfinite(dr0) && dr0 > 0 && isfinite(dra)
                rat = dra / dr0;
            else
                rat = NaN;
            end
            pm_cell(end+1,1) = cid; %#ok<SAGROW>
            pm_T(end+1,1) = sub.T(1); %#ok<SAGROW>
            pm_TC(end+1,1) = sub.T_C(1); %#ok<SAGROW>
            pm_lvl(end+1,1) = moh_lvls(k); %#ok<SAGROW>
            pm_LLI(end+1,1) = Y_int(k, 1); %#ok<SAGROW>
            pm_LAMPE(end+1,1) = Y_int(k, 2); %#ok<SAGROW>
            pm_LAMNE(end+1,1) = Y_int(k, 3); %#ok<SAGROW>
            pm_dR0(end+1,1) = dr0; %#ok<SAGROW>
            pm_dRarc(end+1,1) = dra; %#ok<SAGROW>
            pm_ratio(end+1,1) = rat; %#ok<SAGROW>
        end
    end

    P_moh_interp_tbl = table(pm_cell, pm_T, pm_TC, pm_lvl, pm_LLI, pm_LAMPE, pm_LAMNE, ...
        pm_dR0, pm_dRarc, pm_ratio, ...
        'VariableNames', {'cell','T','T_C','level','LLI','LAM_PE','LAM_NE','dR0','dRarc','ratio'});

    % 2) 온도별 및 Pooled(온도 고정효과 = 온도 내 평균 제거) 편스피어만 상관
    T_groups = {"C", "R", "H", "Pooled_FE"};
    T_deg_map = [-5, 25, 45, NaN];
    resp_list = {"dR0", "dRarc"};
    mode_list = {"LAM_PE", "LAM_NE"};

    pp_T = strings(0,1);
    pp_TC = zeros(0,1);
    pp_lvl = zeros(0,1);
    pp_resp = strings(0,1);
    pp_mode = strings(0,1);
    pp_ctrl = strings(0,1);
    pp_rho = zeros(0,1);
    pp_n = zeros(0,1);

    for l_idx = 1:numel(moh_lvls)
        lv = moh_lvls(l_idx);
        sub_l = P_moh_interp_tbl(P_moh_interp_tbl.level == lv, :);

        for g = 1:numel(T_groups)
            tg = T_groups{g};
            is_pooled = strcmp(tg, 'Pooled_FE');
            if is_pooled
                sub_g = sub_l;
            else
                sub_g = sub_l(sub_l.T == string(tg), :);
            end

            for ri = 1:numel(resp_list)
                rname = resp_list{ri};
                for mi = 1:numel(mode_list)
                    mname = mode_list{mi};

                    x = sub_g.(mname);
                    y = sub_g.(rname);
                    z = sub_g.LLI;
                    grp = sub_g.T;

                    ok = isfinite(x) & isfinite(y) & isfinite(z);
                    x = x(ok); y = y(ok); z = z(ok); grp = grp(ok);
                    n_ok = numel(x);

                    if ~is_pooled
                        rho_val = my_partial_spearman(x, y, z);
                    else
                        rho_val = my_partial_spearman_fe(x, y, z, grp);
                    end

                    pp_T(end+1,1) = string(tg); %#ok<SAGROW>
                    pp_TC(end+1,1) = T_deg_map(g); %#ok<SAGROW>
                    pp_lvl(end+1,1) = lv; %#ok<SAGROW>
                    pp_resp(end+1,1) = string(rname); %#ok<SAGROW>
                    pp_mode(end+1,1) = string(mname); %#ok<SAGROW>
                    pp_ctrl(end+1,1) = "LLI"; %#ok<SAGROW>
                    pp_rho(end+1,1) = rho_val; %#ok<SAGROW>
                    pp_n(end+1,1) = n_ok; %#ok<SAGROW>
                end
            end
        end
    end

    P_partial_tbl = table(pp_T, pp_TC, pp_lvl, pp_resp, pp_mode, pp_ctrl, pp_rho, pp_n, ...
        'VariableNames', {'T_group','T_C','level','response','mode','control','rho_partial','n'});

    writetable(P_partial_tbl, fullfile(out_dir, 'P_mode_partial.csv'));

    log_printf(log_fid, '\n[P] 동일 손실(L=5,10%%)에서 LLI 통제 편스피어만 상관 (Exploratory):\n');
    log_printf(log_fid, '%-10s %5s %-6s %-7s %8s %4s\n', ...
        'T_group','L(%)','resp','mode','rho','n');
    for i = 1:height(P_partial_tbl)
        log_printf(log_fid, '%-10s %5.0f %-6s %-7s %8.3f %4d\n', ...
            P_partial_tbl.T_group(i), P_partial_tbl.level(i), ...
            P_partial_tbl.response(i), P_partial_tbl.mode(i), ...
            P_partial_tbl.rho_partial(i), P_partial_tbl.n(i));
    end

    % 3) 전체 check-up(보간 전) 기반 공선성 분석 (온도별)
    col_T = strings(0,1);
    col_TC = zeros(0,1);
    col_n_rows = zeros(0,1);
    col_n_cells = zeros(0,1);
    col_r_LLI_PE = zeros(0,1);
    col_r_LLI_NE = zeros(0,1);
    col_r_PE_NE = zeros(0,1);
    col_rho_LLI_PE = zeros(0,1);
    col_rho_LLI_NE = zeros(0,1);
    col_rho_PE_NE = zeros(0,1);
    col_cond_raw = zeros(0,1);
    col_cond_std = zeros(0,1);

    T_only = {"C", "R", "H"};
    T_only_C = [-5, 25, 45];
    for g = 1:numel(T_only)
        tg = T_only{g};
        sub = moh_raw(moh_raw.T == string(tg), :);
        M = [sub.LLI, sub.LAM_PE, sub.LAM_NE];
        ok = all(isfinite(M), 2);
        M = M(ok, :);
        sub_ok = sub(ok, :);

        if size(M, 1) >= 3
            r12 = my_pearson(M(:,1), M(:,2));
            r13 = my_pearson(M(:,1), M(:,3));
            r23 = my_pearson(M(:,2), M(:,3));
            s12 = my_spearman(M(:,1), M(:,2));
            s13 = my_spearman(M(:,1), M(:,3));
            s23 = my_spearman(M(:,2), M(:,3));
            c_raw = cond(M);
            sd = std(M, 0, 1);
            if all(sd > 0)
                Mz = (M - mean(M, 1)) ./ sd;
                c_std = cond(Mz);
            else
                c_std = NaN;
            end
        else
            r12 = NaN; r13 = NaN; r23 = NaN;
            s12 = NaN; s13 = NaN; s23 = NaN;
            c_raw = NaN; c_std = NaN;
        end

        col_T(end+1,1) = string(tg); %#ok<SAGROW>
        col_TC(end+1,1) = T_only_C(g); %#ok<SAGROW>
        col_n_rows(end+1,1) = size(M, 1); %#ok<SAGROW>
        col_n_cells(end+1,1) = numel(unique(sub_ok.cell)); %#ok<SAGROW>
        col_r_LLI_PE(end+1,1) = r12; %#ok<SAGROW>
        col_r_LLI_NE(end+1,1) = r13; %#ok<SAGROW>
        col_r_PE_NE(end+1,1) = r23; %#ok<SAGROW>
        col_rho_LLI_PE(end+1,1) = s12; %#ok<SAGROW>
        col_rho_LLI_NE(end+1,1) = s13; %#ok<SAGROW>
        col_rho_PE_NE(end+1,1) = s23; %#ok<SAGROW>
        col_cond_raw(end+1,1) = c_raw; %#ok<SAGROW>
        col_cond_std(end+1,1) = c_std; %#ok<SAGROW>
    end

    P_collin_tbl = table(col_T, col_TC, col_n_rows, col_n_cells, ...
        col_r_LLI_PE, col_r_LLI_NE, col_r_PE_NE, ...
        col_rho_LLI_PE, col_rho_LLI_NE, col_rho_PE_NE, ...
        col_cond_raw, col_cond_std, ...
        'VariableNames', {'T','T_C','n_rows','n_cells', ...
        'r_LLI_LAMPE','r_LLI_LAMNE','r_LAMPE_LAMNE', ...
        'rho_LLI_LAMPE','rho_LLI_LAMNE','rho_LAMPE_LAMNE', ...
        'cond_raw','cond_std'});

    writetable(P_collin_tbl, fullfile(out_dir, 'P_collinearity.csv'));

    log_printf(log_fid, '\n[P] 전체 check-up 기준 LLI/LAM_PE/LAM_NE 공선성:\n');
    log_printf(log_fid, '%-4s %4s %5s %5s %8s %8s %8s %9s %9s\n', ...
        'T','T_C','rows','cells','r(L,PE)','r(L,NE)','r(PE,NE)','cond_raw','cond_std');
    for i = 1:height(P_collin_tbl)
        log_printf(log_fid, '%-4s %4.0f %5d %5d %8.3f %8.3f %8.3f %9.2f %9.2f\n', ...
            P_collin_tbl.T(i), P_collin_tbl.T_C(i), P_collin_tbl.n_rows(i), ...
            P_collin_tbl.n_cells(i), P_collin_tbl.r_LLI_LAMPE(i), ...
            P_collin_tbl.r_LLI_LAMNE(i), P_collin_tbl.r_LAMPE_LAMNE(i), ...
            P_collin_tbl.cond_raw(i), P_collin_tbl.cond_std(i));
    end

catch ME
    log_error(log_fid, 'P Mohtat modes', ME);
end

%% 그림 생성 및 저장
try
    log_printf(log_fid, '\n=== [그림] 시각화 생성 시작 ===\n');

    % --- 그림 1: out\fig_K_equal_loss.png ---
    if ~isempty(K_cells_tbl)
        fig1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1250, 420]);
        tl1 = tiledlayout(fig1, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        cyc_cells = K_cells_tbl(K_cells_tbl.age_group == "cycling" & isfinite(K_cells_tbl.ratio), :);
        lvls = [5, 10, 20];
        kit_Ts = [0, 10, 25, 40];
        moh_Ts = [-5, 25, 45];

        rng(42); % 지터 재현성
        for k = 1:3
            lv = lvls(k);
            ax = nexttile(tl1, k);
            hold(ax, 'on');
            box(ax, 'on');
            grid(ax, 'on');
            set(ax, 'Color', 'w', 'GridAlpha', 0.25);

            yline(ax, 1.0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'HandleVisibility', 'off');

            sub_l = cyc_cells(cyc_cells.level == lv, :);
            all_y_panel = sub_l.ratio;

            h_kit_pt = [];
            h_kit_med = [];
            for ti = 1:numel(kit_Ts)
                Tv = kit_Ts(ti);
                col = get_temp_color(Tv, TEMP_LIST, TEMP_COLORS);
                vals = sub_l.ratio(sub_l.T_age == Tv);
                if isempty(vals)
                    continue;
                end
                jit = (rand(size(vals)) - 0.5) * 3.0;
                ht = scatter(ax, Tv + jit, vals, 32, col, 'filled', ...
                    'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', col*0.7, 'LineWidth', 0.5, ...
                    'HandleVisibility', 'off');
                if isempty(h_kit_pt), h_kit_pt = ht; end

                med_v = median(vals, 'omitnan');
                hm = plot(ax, Tv, med_v, 'o', 'MarkerSize', 9, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.3, ...
                    'HandleVisibility', 'off');
                if isempty(h_kit_med), h_kit_med = hm; end
            end

            % Mohtat 5%, 10% 중앙값 다이아몬드 오버레이
            if ismember(lv, [5, 10]) && ~isempty(P_moh_interp_tbl)
                sub_m = P_moh_interp_tbl(P_moh_interp_tbl.level == lv & isfinite(P_moh_interp_tbl.ratio), :);
                for ti = 1:numel(moh_Ts)
                    Tv = moh_Ts(ti);
                    vals_m = sub_m.ratio(sub_m.T_C == Tv);
                    if isempty(vals_m)
                        continue;
                    end
                    med_m = median(vals_m, 'omitnan');
                    all_y_panel = [all_y_panel; med_m]; %#ok<AGROW>
                    col = get_temp_color(Tv, TEMP_LIST, TEMP_COLORS);
                    plot(ax, Tv, med_m, 'd', 'MarkerSize', 10, ...
                        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.4, ...
                        'HandleVisibility', 'off');
                end
            end

            % 범례용 더미 핸들 (첫 번째 패널에만)
            if k == 1
                p1 = scatter(ax, NaN, NaN, 32, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.55);
                p2 = plot(ax, NaN, NaN, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.3 0.3 0.3], ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.3);
                p3 = plot(ax, NaN, NaN, 'd', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.3 0.3], ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.4);
                legend(ax, [p1, p2, p3], {'KIT cells', 'KIT median', 'Mohtat median'}, ...
                    'Location', 'northeast', 'Box', 'on', 'FontSize', 8.5);
            end

            xlim(ax, [-10, 50]);
            xticks(ax, TEMP_LIST);
            xlabel(ax, 'Ageing temperature T_{age} (^{\circ}C)', 'FontSize', 10);
            if k == 1
                ylabel(ax, 'Resistance growth ratio \DeltaR_{arc} / \DeltaR_0', 'FontSize', 10);
            end
            title(ax, sprintf('Capacity loss L = %d%%', lv), 'FontSize', 11, 'FontWeight', 'bold');

            % 극단 이상치에 의한 축 압축 방지 (필요 시 상한 조절 및 주석)
            yf = all_y_panel(isfinite(all_y_panel));
            if ~isempty(yf)
                q95 = my_quantile(yf, 0.95);
                ymax_raw = max(yf);
                ymin_raw = min(0, min(yf));
                if ymax_raw > max(15, 3 * q95) && q95 > 0
                    y_cap = max(5, 2.2 * q95);
                    n_above = sum(yf > y_cap);
                    ylim(ax, [ymin_raw - 0.05*y_cap, y_cap]);
                    text(ax, 0.03, 0.96, sprintf('(%d pt > %.1f clipped)', n_above, y_cap), ...
                        'Units', 'normalized', 'VerticalAlignment', 'top', ...
                        'FontSize', 8, 'Color', [0.4 0.4 0.4]);
                else
                    pad = max(0.2, 0.1 * (ymax_raw - ymin_raw));
                    ylim(ax, [ymin_raw - 0.05*pad, ymax_raw + pad]);
                end
            end
        end

        title(tl1, 'Equal-loss resistance ratio (\DeltaR_{arc}/\DeltaR_0) vs ageing temperature', ...
            'FontSize', 12, 'FontWeight', 'bold');
        fig1_path = fullfile(out_dir, 'fig_K_equal_loss.png');
        exportgraphics(fig1, fig1_path, 'Resolution', 150, 'BackgroundColor', 'white');
        close(fig1);
        log_printf(log_fid, '[그림] 저장 완료: %s\n', fig1_path);
    end

    % --- 그림 2: out\fig_S_arrhenius.png ---
    if ~isempty(S_cells_tbl)
        fig2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1280, 440]);
        tl2 = tiledlayout(fig2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        tgt_names = {"R0", "arc", "ratio"};
        tgt_cols  = {"s_R0_norm", "s_arc_norm", "ratio_norm"};
        tgt_titles = {'R_0 normalised slope', 'R_{arc} normalised slope', 'Slope ratio s_{arc} / s_{R0}'};
        tgt_ylabels = {'ln(s_{R0} / R_{0,ref})', 'ln(s_{arc} / R_{arc,ref})', 'ln(s_{arc,norm} / s_{R0,norm})'};

        x_grid_kit = linspace(1000/(40+273.15), 1000/(0+273.15), 120)';
        x_grid_moh = linspace(1000/(45+273.15), 1000/(-5+273.15), 120)';

        for k = 1:3
            tgt = tgt_names{k};
            col_name = tgt_cols{k};

            ax = nexttile(tl2, k);
            hold(ax, 'on');
            box(ax, 'on');
            grid(ax, 'on');
            set(ax, 'Color', 'w', 'GridAlpha', 0.25);

            % 셀 데이터 플롯 (KIT 원형, Mohtat 삼각형)
            for ti = 1:numel(TEMP_LIST)
                Tv = TEMP_LIST(ti);
                c_rgb = TEMP_COLORS(ti, :);
                xv = 1000 / (Tv + 273.15);

                sub_k = S_cells_tbl(S_cells_tbl.dataset == "KIT" & S_cells_tbl.T_C == Tv, :);
                vk = sub_k.(col_name);
                vk = vk(isfinite(vk) & vk > 0);
                if ~isempty(vk)
                    scatter(ax, repmat(xv, size(vk)), log(vk), 34, c_rgb, 'o', 'filled', ...
                        'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', c_rgb*0.6, 'LineWidth', 0.5, ...
                        'HandleVisibility', 'off');
                end

                sub_m = S_cells_tbl(S_cells_tbl.dataset == "Mohtat" & S_cells_tbl.T_C == Tv, :);
                vm = sub_m.(col_name);
                vm = vm(isfinite(vm) & vm > 0);
                if ~isempty(vm)
                    scatter(ax, repmat(xv, size(vm)), log(vm), 44, c_rgb, '^', 'filled', ...
                        'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'k', 'LineWidth', 0.6, ...
                        'HandleVisibility', 'off');
                end
            end

            % A1 직선 (KIT, Mohtat) 및 A2 곡선 (KIT, 점선)
            fld_k1 = sprintf('KIT_%s_A1', tgt);
            if isfield(S_fit_store, fld_k1)
                fk1 = S_fit_store.(fld_k1);
                yk1 = fk1.a + fk1.b * (x_grid_kit - X_REF);
                plot(ax, x_grid_kit, yk1, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.7, ...
                    'HandleVisibility', 'off');
            end

            fld_m1 = sprintf('Mohtat_%s_A1', tgt);
            if isfield(S_fit_store, fld_m1)
                fm1 = S_fit_store.(fld_m1);
                ym1 = fm1.a + fm1.b * (x_grid_moh - X_REF);
                plot(ax, x_grid_moh, ym1, '-.', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.7, ...
                    'HandleVisibility', 'off');
            end

            fld_k2 = sprintf('KIT_%s_A2', tgt);
            if isfield(S_fit_store, fld_k2)
                fk2 = S_fit_store.(fld_k2);
                yk2 = eval_A2_single(fk2.p_opt, x_grid_kit - X_REF, R_GAS);
                plot(ax, x_grid_kit, yk2, '--', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.8, ...
                    'HandleVisibility', 'off');
            end

            % 온도 가이드 수직선 및 상단 라벨
            xlim(ax, [3.08, 3.80]);
            yl = ylim(ax);
            ylim(ax, [yl(1), yl(2) + 0.14*(yl(2) - yl(1))]);
            yl2 = ylim(ax);
            for ti = 1:numel(TEMP_LIST)
                Tv = TEMP_LIST(ti);
                xv = 1000 / (Tv + 273.15);
                plot(ax, [xv, xv], yl2, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.7, ...
                    'HandleVisibility', 'off');
                text(ax, xv, yl2(2) - 0.02*(yl2(2)-yl2(1)), sprintf('%d°', Tv), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'FontSize', 7.5, 'Color', TEMP_COLORS(ti,:)*0.8);
            end

            if k == 1
                h1 = scatter(ax, NaN, NaN, 34, [0.3 0.3 0.3], 'o', 'filled');
                h2 = scatter(ax, NaN, NaN, 44, [0.3 0.3 0.3], '^', 'filled', 'MarkerEdgeColor', 'k');
                h3 = plot(ax, NaN, NaN, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.7);
                h4 = plot(ax, NaN, NaN, '-.', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.7);
                h5 = plot(ax, NaN, NaN, '--', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.8);
                legend(ax, [h1, h2, h3, h4, h5], ...
                    {'KIT cells', 'Mohtat cells', 'KIT A1', 'Mohtat A1', 'KIT A2'}, ...
                    'Location', 'southwest', 'FontSize', 8, 'Box', 'on');
            end

            xlabel(ax, '1000 / T (K^{-1})', 'FontSize', 10);
            ylabel(ax, tgt_ylabels{k}, 'FontSize', 10);
            title(ax, tgt_titles{k}, 'FontSize', 11, 'FontWeight', 'bold');
        end

        title(tl2, 'Arrhenius analysis of normalised resistance-per-loss slopes', ...
            'FontSize', 12, 'FontWeight', 'bold');
        fig2_path = fullfile(out_dir, 'fig_S_arrhenius.png');
        exportgraphics(fig2, fig2_path, 'Resolution', 150, 'BackgroundColor', 'white');
        close(fig2);
        log_printf(log_fid, '[그림] 저장 완료: %s\n', fig2_path);
    end

catch ME
    log_error(log_fid, '그림 생성', ME);
end

log_printf(log_fid, '\n=== 완료 (총 소요 시간: %.2f s) ===\n', toc(t_start_all));


%% ==================== 로컬 함수 ====================

function close_if_open(fid)
    if fid >= 3
        try
            fclose(fid);
        catch
        end
    end
end

function log_printf(fid, fmt, varargin)
    fprintf(1, fmt, varargin{:});
    if fid >= 3
        fprintf(fid, fmt, varargin{:});
    end
end

function log_error(fid, sec_name, ME)
    log_printf(fid, '[ERROR in %s] %s\n', sec_name, ME.message);
    for k = 1:numel(ME.stack)
        log_printf(fid, '    at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
end

function T = read_csv_clean(fpath)
    opts = detectImportOptions(fpath, 'Delimiter', ',');
    opts.VariableNamingRule = 'preserve';
    T = readtable(fpath, opts);
    % UTF-8 BOM 제거
    vnames = T.Properties.VariableNames;
    for i = 1:numel(vnames)
        vnames{i} = regexprep(vnames{i}, '^\xFEFF', '');
    end
    T.Properties.VariableNames = vnames;
end

function kit_proc = preprocess_kit(kit_raw)
    % 셀별 최소 day 기준점 계산 및 L, dR0, dRarc 산출
    n = height(kit_raw);
    L = nan(n, 1);
    dR0 = nan(n, 1);
    dRarc = nan(n, 1);
    R0_ref_col = nan(n, 1);
    Rarc_ref_col = nan(n, 1);

    u_cells = unique(kit_raw.cell, 'stable');
    for i = 1:numel(u_cells)
        idx = find(kit_raw.cell == u_cells(i));
        days = kit_raw.day(idx);
        [~, ord] = sort(days, 'ascend');
        idx_sorted = idx(ord);

        ref_i = idx_sorted(1);
        soh_ref = kit_raw.soh(ref_i);
        r0_ref = kit_raw.R0(ref_i);
        ra_ref = kit_raw.R_arc(ref_i);

        L(idx) = 100 * (soh_ref - kit_raw.soh(idx)) / soh_ref;
        dR0(idx) = kit_raw.R0(idx) - r0_ref;
        dRarc(idx) = kit_raw.R_arc(idx) - ra_ref;
        R0_ref_col(idx) = r0_ref;
        Rarc_ref_col(idx) = ra_ref;
    end

    kit_proc = kit_raw;
    kit_proc.L = L;
    kit_proc.dR0 = dR0;
    kit_proc.dRarc = dRarc;
    kit_proc.R0_ref = R0_ref_col;
    kit_proc.Rarc_ref = Rarc_ref_col;
end

function key = make_proto_key(age_type, chg_rate, dis_rate, v_min, v_max, cal_soc)
    vals = [double(age_type), double(chg_rate), double(dis_rate), ...
            double(v_min), double(v_max), double(cal_soc)];
    parts = strings(1, numel(vals));
    for i = 1:numel(vals)
        if isnan(vals(i))
            parts(i) = "NaN";
        else
            parts(i) = string(sprintf('%.6g', vals(i)));
        end
    end
    key = strjoin(parts, '|');
end

function Y_out = interp_cell_cummax(L_raw, Y_mat, levels)
    % L을 cummax로 단조 증가화 후 지정 손실 수준에서 선형 보간 (외삽 금지)
    n_lvls = numel(levels);
    q = size(Y_mat, 2);
    Y_out = nan(n_lvls, q);

    okL = isfinite(L_raw);
    if ~any(okL)
        return;
    end
    L_mono = cummax(L_raw(okL));
    Yf = Y_mat(okL, :);

    for j = 1:q
        okY = isfinite(Yf(:, j));
        if ~any(okY)
            continue;
        end
        Lj = L_mono(okY);
        yj = Yf(okY, j);
        [Lu, idx_u] = unique(Lj, 'first');
        yu = yj(idx_u);

        for k = 1:n_lvls
            lv = levels(k);
            if numel(Lu) >= 2 && Lu(1) <= lv && Lu(end) >= lv
                Y_out(k, j) = interp1(Lu, yu, lv, 'linear');
            elseif numel(Lu) == 1 && Lu(1) == lv
                Y_out(k, j) = yu(1);
            end
        end
    end
end

function q = my_quantile(x, p)
    % 베이스 MATLAB 분위수 계산 (선형 보간, NaN 제외)
    xf = sort(x(isfinite(x)));
    if isempty(xf)
        q = nan(size(p));
        return;
    end
    n = numel(xf);
    if n == 1
        q = repmat(xf, size(p));
        return;
    end
    pos = p * n + 0.5;
    pos = min(max(pos, 1), n);
    lo = floor(pos);
    hi = ceil(pos);
    frac = pos - lo;
    xf = xf(:).';   % 행벡터로 맞춰 암시적 확장 방지
    q = xf(lo) + frac .* (xf(hi) - xf(lo));
    q = reshape(q, size(p));
end

function val = my_iqr(x)
    qs = my_quantile(x, [0.25, 0.75]);
    val = qs(2) - qs(1);
end

function pval = my_signtest(k_pos, n_nz)
    % 정확 양측 이항 부호검정 (nchoosek 합 이용)
    if n_nz <= 0
        pval = NaN;
        return;
    end
    k_min = min(k_pos, n_nz - k_pos);
    if 2 * k_min == n_nz
        pval = 1.0;
        return;
    end
    p_tail = 0.0;
    for k = 0:k_min
        if n_nz <= 50
            p_tail = p_tail + nchoosek(n_nz, k) * (0.5 ^ n_nz);
        else
            p_tail = p_tail + exp(gammaln(n_nz + 1) - gammaln(k + 1) - ...
                gammaln(n_nz - k + 1) - n_nz * log(2));
        end
    end
    pval = min(1.0, 2.0 * p_tail);
end

function r = my_tiedrank(x)
    % 평균 동점 순위(average-tie ranks)
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
        avg_rank = 0.5 * (i + j);
        rs(i:j) = avg_rank;
        i = j + 1;
    end
    r = zeros(n, 1);
    r(ord) = rs;
end

function r = my_pearson(x, y)
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 2
        r = NaN;
        return;
    end
    dx = x(:) - mean(x);
    dy = y(:) - mean(y);
    denom = sqrt(sum(dx.^2) * sum(dy.^2));
    if denom == 0
        r = NaN;
    else
        r = sum(dx .* dy) / denom;
    end
end

function rho = my_spearman(x, y)
    ok = isfinite(x) & isfinite(y);
    if sum(ok) < 2
        rho = NaN;
        return;
    end
    rho = my_pearson(my_tiedrank(x(ok)), my_tiedrank(y(ok)));
end

function rho = my_partial_spearman(x, y, z)
    % z 순위로 x, y 순위를 회귀한 잔차 간 Pearson 상관
    n = numel(x);
    if n < 3
        rho = NaN;
        return;
    end
    rx = my_tiedrank(x);
    ry = my_tiedrank(y);
    rz = my_tiedrank(z);
    Zmat = [ones(n, 1), rz];
    ex = rx - Zmat * (Zmat \ rx);
    ey = ry - Zmat * (Zmat \ ry);
    if std(ex) < 1e-12 || std(ey) < 1e-12
        rho = NaN;
    else
        rho = my_pearson(ex, ey);
    end
end

function rho = my_partial_spearman_fe(x, y, z, grp)
    % 온도 고정효과(그룹 내 평균 제거)를 포함한 편스피어만 상관
    n = numel(x);
    if n < 4
        rho = NaN;
        return;
    end
    rx = my_tiedrank(x);
    ry = my_tiedrank(y);
    rz = my_tiedrank(z);

    ug = unique(grp);
    rx_dm = rx; ry_dm = ry; rz_dm = rz;
    for i = 1:numel(ug)
        idx = (grp == ug(i));
        rx_dm(idx) = rx(idx) - mean(rx(idx));
        ry_dm(idx) = ry(idx) - mean(ry(idx));
        rz_dm(idx) = rz(idx) - mean(rz(idx));
    end

    if sum(rz_dm.^2) < 1e-12
        ex = rx_dm;
        ey = ry_dm;
    else
        bx = (rz_dm' * rx_dm) / (rz_dm' * rz_dm);
        by = (rz_dm' * ry_dm) / (rz_dm' * rz_dm);
        ex = rx_dm - rz_dm * bx;
        ey = ry_dm - rz_dm * by;
    end

    if std(ex) < 1e-12 || std(ey) < 1e-12
        rho = NaN;
    else
        rho = my_pearson(ex, ey);
    end
end

function fit = fit_A1(x, y, ds_id, is_joint, x_ref, R_gas)
    % A1: log s = a + b*(x - x_ref), Ea_kJ = -b * R_gas (양수 Ea = 고온에서 기울기 증가)
    dx = x(:) - x_ref;
    y = y(:);
    n = numel(y);

    if ~is_joint
        X = [ones(n, 1), dx];
        beta = X \ y;
        yhat = X * beta;
        fit.a = beta(1);
        fit.a_KIT = beta(1);
        fit.a_MOH = beta(1);
        fit.b = beta(2);
        k = 3; % a, b, sigma^2
    else
        X = [double(ds_id == 1), double(ds_id == 2), dx];
        beta = X \ y;
        yhat = X * beta;
        fit.a = NaN;
        fit.a_KIT = beta(1);
        fit.a_MOH = beta(2);
        fit.b = beta(3);
        k = 4; % a_KIT, a_MOH, b, sigma^2
    end

    res = y - yhat;
    SSE = sum(res.^2);
    fit.Ea_kJ = -fit.b * R_gas;
    fit.SSE = SSE;
    fit.k = k;
    fit.AIC = n * log(max(SSE / n, realmin)) + 2 * k;
end

function log_s = eval_A2_single(p, dx, R_gas)
    % p = [a1, a2, logEh, logEc]
    a1 = p(1);
    a2 = p(2);
    Eh = exp(min(max(p(3), -15), log(500)));
    Ec = exp(min(max(p(4), -15), log(500)));
    u1 = a1 - (Eh / R_gas) * dx;
    u2 = a2 + (Ec / R_gas) * dx;
    m = max(u1, u2);
    log_s = m + log(exp(u1 - m) + exp(u2 - m));
end

function fit = fit_A2_fast(x, y, ds_id, is_joint, x_ref, R_gas, warm_p)
    % A2 두 갈래 Arrhenius 모델 (그룹 축약 항등식을 이용한 고속 정확한 최소제곱)
    dx = x(:) - x_ref;
    y = y(:);
    n = numel(y);

    % 고유 (ds_id, dx) 그룹 축약: sum(y_i - f_g)^2 = SSW + sum n_g*(ybar_g - f_g)^2
    if ~is_joint
        [ux, ~, gidx] = unique(dx);
        ng = numel(ux);
        g_ds = ones(ng, 1);
        g_dx = ux;
    else
        keys = [ds_id(:), dx];
        [uk, ~, gidx] = unique(keys, 'rows');
        ng = size(uk, 1);
        g_ds = uk(:, 1);
        g_dx = uk(:, 2);
    end

    g_n = zeros(ng, 1);
    g_ybar = zeros(ng, 1);
    SSW = 0;
    for g = 1:ng
        yg = y(gidx == g);
        g_n(g) = numel(yg);
        g_ybar(g) = mean(yg);
        SSW = SSW + sum((yg - g_ybar(g)).^2);
    end

    opts = optimset('Display', 'off', 'MaxFunEvals', 1200, 'MaxIter', 800, ...
        'TolX', 1e-5, 'TolFun', 1e-6);

    if ~is_joint
        k = 5; % a1, a2, Eh, Ec, sigma^2
        obj = @(p) obj_A2_single(p, g_dx, g_ybar, g_n, R_gas);
        if isempty(warm_p)
            a0 = mean(y);
            starts = [ ...
                a0 - 0.5, a0 - 0.5, log(20), log(20); ...
                a0 - 0.2, a0 - 1.5, log(40), log(10); ...
                a0 - 1.5, a0 - 0.2, log(10), log(40); ...
                a0 - 0.7, a0 - 0.7, log(60), log(60); ...
                a0 - 0.3, a0 - 2.0, log(15), log(80); ...
                a0 - 2.0, a0 - 0.3, log(80), log(15); ...
                a0 - 0.5, a0 - 0.5, log(5),  log(5);  ...
                a0 - 0.5, a0 - 0.5, log(100), log(25) ...
            ];
        else
            starts = warm_p(:)';
        end
    else
        k = 7; % a1_K, a2_K, a1_M, a2_M, Eh, Ec, sigma^2
        obj = @(p) obj_A2_joint(p, g_ds, g_dx, g_ybar, g_n, R_gas);
        if isempty(warm_p)
            aK = mean(y(ds_id == 1));
            aM = mean(y(ds_id == 2));
            starts = [ ...
                aK - 0.5, aK - 0.5, aM - 0.5, aM - 0.5, log(20), log(20); ...
                aK - 0.2, aK - 1.5, aM - 0.2, aM - 1.5, log(40), log(10); ...
                aK - 1.5, aK - 0.2, aM - 1.5, aM - 0.2, log(10), log(40); ...
                aK - 0.7, aK - 0.7, aM - 0.7, aM - 0.7, log(60), log(60); ...
                aK - 0.3, aK - 2.0, aM - 1.5, aM - 0.3, log(25), log(50); ...
                aK - 2.0, aK - 0.3, aM - 0.3, aM - 2.0, log(50), log(25); ...
                aK - 0.5, aK - 0.5, aM - 0.5, aM - 0.5, log(5),  log(5)   ...
            ];
        else
            starts = warm_p(:)';
        end
    end

    best_f = inf;
    best_p = starts(1, :);
    for s = 1:size(starts, 1)
        [p_try, f_try] = fminsearch(obj, starts(s, :), opts);
        if f_try < best_f
            best_f = f_try;
            best_p = p_try;
        end
    end

    % 최종 SSE는 경계 페널티 없이 정확한 잔차 제곱합으로 계산
    if ~is_joint
        p_cl = best_p;
        p_cl(3) = min(max(p_cl(3), -15), log(500));
        p_cl(4) = min(max(p_cl(4), -15), log(500));
        f_g = eval_A2_single(p_cl, g_dx, R_gas);
        SSE = SSW + sum(g_n .* (g_ybar - f_g).^2);
        fit.p_opt = p_cl;
        fit.a1_KIT = p_cl(1);
        fit.a2_KIT = p_cl(2);
        fit.a1_MOH = NaN;
        fit.a2_MOH = NaN;
        fit.Eh_kJ = exp(p_cl(3));
        fit.Ec_kJ = exp(p_cl(4));
    else
        p_cl = best_p;
        p_cl(5) = min(max(p_cl(5), -15), log(500));
        p_cl(6) = min(max(p_cl(6), -15), log(500));
        f_g = zeros(ng, 1);
        idxK = (g_ds == 1);
        idxM = (g_ds == 2);
        f_g(idxK) = eval_A2_single([p_cl(1), p_cl(2), p_cl(5), p_cl(6)], g_dx(idxK), R_gas);
        f_g(idxM) = eval_A2_single([p_cl(3), p_cl(4), p_cl(5), p_cl(6)], g_dx(idxM), R_gas);
        SSE = SSW + sum(g_n .* (g_ybar - f_g).^2);
        fit.p_opt = p_cl;
        fit.a1_KIT = p_cl(1);
        fit.a2_KIT = p_cl(2);
        fit.a1_MOH = p_cl(3);
        fit.a2_MOH = p_cl(4);
        fit.Eh_kJ = exp(p_cl(5));
        fit.Ec_kJ = exp(p_cl(6));
    end

    fit.SSE = SSE;
    fit.k = k;
    fit.AIC = n * log(max(SSE / n, realmin)) + 2 * k;
end

function val = obj_A2_single(p, g_dx, g_ybar, g_n, R_gas)
    pen = 0;
    if p(3) < -15, pen = pen + 1e3 * (p(3) + 15)^2; end
    if p(3) > log(500), pen = pen + 1e3 * (p(3) - log(500))^2; end
    if p(4) < -15, pen = pen + 1e3 * (p(4) + 15)^2; end
    if p(4) > log(500), pen = pen + 1e3 * (p(4) - log(500))^2; end
    if abs(p(1)) > 50, pen = pen + 1e2 * (abs(p(1)) - 50)^2; end
    if abs(p(2)) > 50, pen = pen + 1e2 * (abs(p(2)) - 50)^2; end

    f_g = eval_A2_single(p, g_dx, R_gas);
    val = sum(g_n .* (g_ybar - f_g).^2) + pen;
end

function val = obj_A2_joint(p, g_ds, g_dx, g_ybar, g_n, R_gas)
    pen = 0;
    if p(5) < -15, pen = pen + 1e3 * (p(5) + 15)^2; end
    if p(5) > log(500), pen = pen + 1e3 * (p(5) - log(500))^2; end
    if p(6) < -15, pen = pen + 1e3 * (p(6) + 15)^2; end
    if p(6) > log(500), pen = pen + 1e3 * (p(6) - log(500))^2; end
    for j = 1:4
        if abs(p(j)) > 50, pen = pen + 1e2 * (abs(p(j)) - 50)^2; end
    end

    f_g = zeros(numel(g_dx), 1);
    idxK = (g_ds == 1);
    idxM = (g_ds == 2);
    f_g(idxK) = eval_A2_single([p(1), p(2), p(5), p(6)], g_dx(idxK), R_gas);
    f_g(idxM) = eval_A2_single([p(3), p(4), p(5), p(6)], g_dx(idxM), R_gas);
    val = sum(g_n .* (g_ybar - f_g).^2) + pen;
end

function strata = build_strata(ds_col, T_col)
    keys = ds_col + "|" + string(T_col);
    uk = unique(keys, 'stable');
    strata = cell(numel(uk), 1);
    for i = 1:numel(uk)
        strata{i} = find(keys == uk(i));
    end
end

function idx = resample_strata(strata)
    n_s = numel(strata);
    idx = [];
    for i = 1:n_s
        s_idx = strata{i};
        m = numel(s_idx);
        pick = s_idx(randi(m, m, 1));
        idx = [idx; pick]; %#ok<AGROW>
    end
end

function [sref_K, sref_M, a1_K, a1_M] = unpack_A1_params(fitA1, ds_name)
    if strcmp(ds_name, 'KIT')
        a1_K = fitA1.a_KIT; a1_M = NaN;
        sref_K = exp(a1_K); sref_M = NaN;
    elseif strcmp(ds_name, 'Mohtat')
        a1_K = NaN; a1_M = fitA1.a_MOH;
        sref_K = NaN; sref_M = exp(a1_M);
    else
        a1_K = fitA1.a_KIT; a1_M = fitA1.a_MOH;
        sref_K = exp(a1_K); sref_M = exp(a1_M);
    end
end

function [sref_K, sref_M, a1_K, a2_K, a1_M, a2_M] = unpack_A2_params(fitA2, ds_name)
    if strcmp(ds_name, 'KIT')
        a1_K = fitA2.a1_KIT; a2_K = fitA2.a2_KIT;
        a1_M = NaN; a2_M = NaN;
        sref_K = exp(a1_K) + exp(a2_K);
        sref_M = NaN;
    else
        a1_K = fitA2.a1_KIT; a2_K = fitA2.a2_KIT;
        a1_M = fitA2.a1_MOH; a2_M = fitA2.a2_MOH;
        sref_K = exp(a1_K) + exp(a2_K);
        sref_M = exp(a1_M) + exp(a2_M);
    end
end

function eq_str = format_eq_A1(tgt, ds_name, sref_K, sref_M, Ea_kJ)
    if strcmp(ds_name, 'KIT')
        sref_txt = sprintf('s_ref = %.4g', sref_K);
    elseif strcmp(ds_name, 'Mohtat')
        sref_txt = sprintf('s_ref = %.4g', sref_M);
    else
        sref_txt = sprintf('s_ref(KIT) = %.4g, s_ref(MOH) = %.4g', sref_K, sref_M);
    end
    eq_str = string(sprintf('s_%s(T) = s_ref*exp[-Ea/R*(1/T - 1/298.15)], %s, Ea = %.2f kJ/mol', ...
        tgt, sref_txt, Ea_kJ));
end

function eq_str = format_eq_A2(tgt, ds_name, a1_K, a2_K, a1_M, a2_M, Eh_kJ, Ec_kJ)
    if strcmp(ds_name, 'KIT')
        amp_txt = sprintf('s_h = %.4g, s_c = %.4g', exp(a1_K), exp(a2_K));
    else
        amp_txt = sprintf('KIT(s_h=%.4g,s_c=%.4g), MOH(s_h=%.4g,s_c=%.4g)', ...
            exp(a1_K), exp(a2_K), exp(a1_M), exp(a2_M));
    end
    eq_str = string(sprintf('s_%s(T) = s_h*exp[-Eh/R*(1/T - 1/298.15)] + s_c*exp[+Ec/R*(1/T - 1/298.15)], %s, Eh = %.2f kJ/mol, Ec = %.2f kJ/mol', ...
        tgt, amp_txt, Eh_kJ, Ec_kJ));
end

function col = get_temp_color(Tv, temp_list, temp_colors)
    idx = find(temp_list == Tv, 1);
    if isempty(idx)
        col = [0.4, 0.4, 0.4];
    else
        col = temp_colors(idx, :);
    end
end