%% 0 설정
% KIT 배터리 저항 열화 분석 (R0, R_arc): (A) 아크 분해, (B) LOTO 모델 비교, (C) WLTP 검증
clearvars;
close all;
clc;
t_start_all = tic;

% 입력 경로 및 출력 폴더
path_rt    = 'kit_rt_soc50.csv';
path_pairs = fullfile('..', 'kit_optemp', 'out', 'optemp_pairs.csv');
out_dir    = 'out';
if ~isfolder(out_dir)
    mkdir(out_dir);
end
log_file = fullfile(out_dir, 'fit_log.txt');
fid_log = fopen(log_file, 'w');
if fid_log > 0
    fprintf(fid_log, '=== KIT Models Fit Log (%s) ===\n', char(datetime('now')));
    fclose(fid_log);
end

% 물리 상수 및 가정 파라미터
Rg      = 8.314e-3;   % kJ / (mol * K)
Tref    = 298.15;     % K (25 degC)
I_WLTP  = 1.0;        % WLTP 유효 충전 C-rate 가정값
SOC_CYC = 50.0;       % 사이클/WLTP 유효 SOC (%)

% 난수 시드 및 최적화 옵션
rng(1, 'twister');
n_boot       = 200;   % (A) 셀 단위 부트스트랩 반복 횟수
n_multistart = 8;     % (B) 멀티스타트 횟수 (>= 6)

% lsqnonlin 사용 가능 여부 확인 (Optimization Toolbox 없으면 fminsearch로 대체)
use_lsq = false;
if exist('lsqnonlin', 'file') > 0
    try
        opt_test = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 5);
        lsqnonlin(@(p) p, 0.1, -1, 1, opt_test);
        use_lsq = true;
    catch ME
        log_msg(log_file, sprintf('lsqnonlin check failed (%s); falling back to fminsearch.', ME.message));
        use_lsq = false;
    end
else
    log_msg(log_file, 'lsqnonlin not found; using fminsearch fallback with bounded transforms.');
end

if ~use_lsq
    n_multistart = 6; % fminsearch 폴백 시 속도 확보 (>= 6 유지)
end

fprintf('[0 설정] 시작 | Solver: %s | Multistarts: %d | Bootstrap: %d\n', ...
    ternary(use_lsq, 'lsqnonlin (analytical Jacobian)', 'fminsearch (bounded transform)'), ...
    n_multistart, n_boot);

%% 1 전처리
fprintf('[1 전처리] 데이터 로딩 및 전처리 중...\n');

if ~isfile(path_rt)
    error('입력 파일을 찾을 수 없습니다: %s', path_rt);
end
if ~isfile(path_pairs)
    error('입력 파일을 찾을 수 없습니다: %s', path_pairs);
end

raw_rt = readtable(path_rt, 'TextType', 'string');
raw_rt = sanitize_colnames(raw_rt);
raw_rt = ensure_numeric_cols(raw_rt, {'day','efc','soh','R0','R_arc','T_age','age_type','chg_rate','dis_rate','cal_soc'});

% RT 데이터 셀별 전처리
T_rt = preprocess_rt_data(raw_rt, Rg, Tref, I_WLTP, SOC_CYC);
fprintf('  RT 데이터: 총 %d행, %d셀 (cal=%d, cyc=%d, WLTP=%d)\n', ...
    height(T_rt), numel(unique(T_rt.cell)), ...
    numel(unique(T_rt.cell(T_rt.age_type == 1))), ...
    numel(unique(T_rt.cell(T_rt.age_type == 2))), ...
    numel(unique(T_rt.cell(T_rt.age_type == 3))));

% 쌍(paired) 데이터 로딩 및 SOC==50 필터링
raw_pairs = readtable(path_pairs, 'TextType', 'string');
raw_pairs = sanitize_colnames(raw_pairs);
raw_pairs = ensure_numeric_cols(raw_pairs, {'soc','T_meas_op','T_meas_rt','R0_op','R0_rt','R_arc_op','R_arc_rt','soh','day','efc','T_age','age_type'});

T_pairs_arc = preprocess_pairs_data(raw_pairs, 'R_arc', Rg);
T_pairs_r0  = preprocess_pairs_data(raw_pairs, 'R0', Rg);
fprintf('  Paired 데이터 (SOC=50, T_age in {0,10,40}): R_arc=%d행 (%d셀), R0=%d행 (%d셀)\n', ...
    height(T_pairs_arc), numel(unique(T_pairs_arc.cell)), ...
    height(T_pairs_r0),  numel(unique(T_pairs_r0.cell)));

%% 2 (A) 아크 분해
fprintf('\n[2 (A) 아크 분해] 모델 적합 및 부트스트랩 실행 중...\n');

A_rows = {};
A_fits = struct();
A_comp_names = {'R_arc', 'R0'};
A_tables     = {T_pairs_arc, T_pairs_r0};
T_age_A      = [0, 10, 40];
lb_E         = -50;
ub_E         = 150;

for ic = 1:2
    cname = A_comp_names{ic};
    Tp    = A_tables{ic};
    fprintf('  - Component: %s (n=%d, cells=%d)\n', cname, height(Tp), numel(unique(Tp.cell)));

    % A0, A1, A2 적합
    variants = {'A0', 'A1', 'A2'};
    for iv = 1:numel(variants)
        vname = variants{iv};
        switch vname
            case 'A0'
                pnames = {'E_ct'};
                lb = lb_E; ub = ub_E;
                p0_list = [40; 10; 70];
                do_boot = ~use_lsq;
            case 'A1'
                pnames = {'E_ct', 'E_new'};
                lb = [lb_E, lb_E]; ub = [ub_E, ub_E];
                p0_list = [50, 15; 20, 5; 70, 30; 30, -5];
                do_boot = true; % A1은 항상 셀 단위 부트스트랩 수행
            case 'A2'
                pnames = {'E_ct', 'E_new_T0', 'E_new_T10', 'E_new_T40'};
                lb = repmat(lb_E, 1, 4); ub = repmat(ub_E, 1, 4);
                p0_list = [50, 10, 20, 30; 20, 5, 10, 15; 70, 20, 30, 40; 35, -5, 15, 25];
                do_boot = true; % A2는 항상 셀 단위 부트스트랩 수행
        end

        k = numel(pnames);
        n = height(Tp);

        % 본 적합
        best_p   = nan(1, k);
        best_rss = inf;
        best_J   = [];
        fit_ok   = false;

        for ist = 1:size(p0_list, 1)
            try
                obj_fun = @(p) res_A_model(p, Tp, vname);
                [p_hat, rss_hat, J_hat] = fit_nonlinear_lsq(obj_fun, p0_list(ist, :), lb, ub, use_lsq, true);
                if isfinite(rss_hat) && rss_hat < best_rss
                    best_p   = p_hat;
                    best_rss = rss_hat;
                    best_J   = J_hat;
                    fit_ok   = true;
                end
            catch ME
                log_msg(log_file, sprintf('[A %s %s start %d] %s', cname, vname, ist, ME.message));
            end
        end

        if ~fit_ok
            log_msg(log_file, sprintf('[A %s %s] All starts failed.', cname, vname));
            ci_lo = nan(1, k); ci_hi = nan(1, k);
            aic = NaN; bic = NaN;
        else
            aic = n * log(max(best_rss, realmin) / n) + 2 * k;
            bic = n * log(max(best_rss, realmin) / n) + k * log(n);

            % 95% 신뢰구간 계산 (Jacobian 또는 셀 단위 부트스트랩)
            if do_boot
                rng(100 + ic * 10 + iv, 'twister');
                boot_P = run_A_bootstrap(Tp, vname, best_p, lb, ub, n_boot, use_lsq, log_file);
                ci_lo = nan(1, k); ci_hi = nan(1, k);
                for ip = 1:k
                    ci_lo(ip) = my_prctile(boot_P(:, ip), 2.5);
                    ci_hi(ip) = my_prctile(boot_P(:, ip), 97.5);
                end
            else
                [ci_lo, ci_hi] = jacobian_ci95(best_p, best_rss, best_J, n, k);
            end
        end

        % 적합 결과 저장 (그림용)
        A_fits.(cname).(vname).p     = best_p;
        A_fits.(cname).(vname).RSS   = best_rss;
        A_fits.(cname).(vname).n     = n;
        A_fits.(cname).(vname).k     = k;
        A_fits.(cname).(vname).AIC   = aic;
        A_fits.(cname).(vname).BIC   = bic;
        A_fits.(cname).(vname).ci_lo = ci_lo;
        A_fits.(cname).(vname).ci_hi = ci_hi;

        for ip = 1:k
            A_rows(end+1, :) = {vname, cname, pnames{ip}, best_p(ip), ci_lo(ip), ci_hi(ip), best_rss, n, k, aic, bic}; %#ok<SAGROW>
        end
    end

    % Fresh-only 체크 (SOH >= 0.98, T_age별 단일 Ea)
    for it = 1:numel(T_age_A)
        tval = T_age_A(it);
        sub_f = Tp(Tp.T_age == tval & Tp.soh >= 0.98, :);
        vname = sprintf('fresh_T%d', tval);
        nf = height(sub_f);
        if nf < 2
            log_msg(log_file, sprintf('[A %s %s] Too few fresh pairs (n=%d).', cname, vname, nf));
            A_rows(end+1, :) = {vname, cname, 'E_fresh', NaN, NaN, NaN, NaN, nf, 1, NaN, NaN}; %#ok<SAGROW>
            continue;
        end
        try
            obj_f = @(p) res_A_fresh(p, sub_f);
            [p_f, rss_f, J_f] = fit_nonlinear_lsq(obj_f, 40, lb_E, ub_E, use_lsq, true);
            aic_f = nf * log(max(rss_f, realmin) / nf) + 2;
            bic_f = nf * log(max(rss_f, realmin) / nf) + log(nf);
            if use_lsq && ~isempty(J_f)
                [lo_f, hi_f] = jacobian_ci95(p_f, rss_f, J_f, nf, 1);
            else
                rng(300 + ic * 10 + it, 'twister');
                boot_f = run_A_bootstrap(sub_f, 'fresh', p_f, lb_E, ub_E, n_boot, use_lsq, log_file);
                lo_f = my_prctile(boot_f(:, 1), 2.5);
                hi_f = my_prctile(boot_f(:, 1), 97.5);
            end
            A_rows(end+1, :) = {vname, cname, 'E_fresh', p_f(1), lo_f(1), hi_f(1), rss_f, nf, 1, aic_f, bic_f}; %#ok<SAGROW>
        catch ME
            log_msg(log_file, sprintf('[A %s %s] Fit failed: %s', cname, vname, ME.message));
            A_rows(end+1, :) = {vname, cname, 'E_fresh', NaN, NaN, NaN, NaN, nf, 1, NaN, NaN}; %#ok<SAGROW>
        end
    end
end

T_A_out = cell2table(A_rows, 'VariableNames', ...
    {'variant','component','param','value','ci_lo','ci_hi','RSS','n','k','AIC','BIC'});
writetable(T_A_out, fullfile(out_dir, 'A_arc_decomposition.csv'));

fprintf('\n--- [표 A] 아크 분해 결과 (out/A_arc_decomposition.csv) ---\n');
fprintf('%-10s %-6s %-10s %8s [%8s, %8s] %9s %6s %9s %9s\n', ...
    'Variant', 'Comp', 'Param', 'Value', 'CI_lo', 'CI_hi', 'RSS', 'n', 'AIC', 'BIC');
for r = 1:height(T_A_out)
    fprintf('%-10s %-6s %-10s %8.2f [%8.2f, %8.2f] %9.3f %6d %9.1f %9.1f\n', ...
        T_A_out.variant{r}, T_A_out.component{r}, T_A_out.param{r}, ...
        T_A_out.value(r), T_A_out.ci_lo(r), T_A_out.ci_hi(r), ...
        T_A_out.RSS(r), T_A_out.n(r), T_A_out.AIC(r), T_A_out.BIC(r));
end

%% 3 (B) LOTO
fprintf('\n[3 (B) LOTO] 온도 제외 교차검증 및 전체 적합 실행 중...\n');

T_rt_12 = T_rt(T_rt.age_type == 1 | T_rt.age_type == 2, :);
B_comp_names = {'dR0', 'dRarc'};
B_models     = {'M1', 'M2', 'M3', 'M4'};
T_ho_list    = [0, 10, 25, 40];

B_loto_rows   = {};
B_param_rows  = {};
B_insample    = struct();

for ic = 1:numel(B_comp_names)
    cname = B_comp_names{ic};
    T_sub = T_rt_12(isfinite(T_rt_12.(cname)), :);
    fprintf('  - Component: %s (age_type 1+2: %d행, %d셀)\n', ...
        cname, height(T_sub), numel(unique(T_sub.cell)));

    % 3.1 전체 온도 인샘플 적합 (순서 M1 -> M3 -> M2 -> M4로 하여 M3 해를 M2/M4 초기값에 활용)
    fit_order = {'M1', 'M3', 'M2', 'M4'};
    insample_p = struct();
    for im = 1:numel(fit_order)
        mname = fit_order{im};
        rng(1, 'twister');
        seed_p = [];
        if strcmp(mname, 'M2') && isfield(insample_p, 'M3') && all(isfinite(insample_p.M3))
            p3 = insample_p.M3;
            seed_p = [p3(1), p3(2), p3(3), 0.5, p3(4), p3(5), p3(6), 1.0];
        elseif strcmp(mname, 'M4') && isfield(insample_p, 'M3') && all(isfinite(insample_p.M3))
            p3 = insample_p.M3;
            seed_p = [p3(1), p3(2), p3(3), p3(4), p3(5), p3(6), -3.0, 45.0];
        end

        [p_best, ~, ok] = fit_B_model_multistart(T_sub, cname, mname, n_multistart, use_lsq, seed_p, log_file, sprintf('Insample_%s_%s', cname, mname));
        insample_p.(mname) = p_best;
        spec = get_B_model_spec(mname);
        k = numel(spec.names);
        n_in = height(T_sub);

        if ok
            pred_in = eval_B_model(p_best, T_sub, mname);
            err_in  = pred_in - T_sub.(cname);
            rss_in  = sum(err_in .^ 2);
            rmse_in = sqrt(rss_in / n_in);
            aic_in  = n_in * log(max(rss_in, realmin) / n_in) + 2 * k;
            bic_in  = n_in * log(max(rss_in, realmin) / n_in) + k * log(n_in);
        else
            rss_in = NaN; rmse_in = NaN; aic_in = NaN; bic_in = NaN;
        end

        B_insample.(cname).(mname).p    = p_best;
        B_insample.(cname).(mname).RMSE = rmse_in;
        B_insample.(cname).(mname).AIC  = aic_in;
        B_insample.(cname).(mname).BIC  = bic_in;
        B_insample.(cname).(mname).ok   = ok;
    end

    % M1..M4 순서로 인샘플 파라미터 저장
    for im = 1:numel(B_models)
        mname = B_models{im};
        spec  = get_B_model_spec(mname);
        info  = B_insample.(cname).(mname);
        n_in  = height(T_sub);
        k     = numel(spec.names);
        rss_in = info.RMSE^2 * n_in;
        for ip = 1:k
            B_param_rows(end+1, :) = {cname, mname, spec.names{ip}, info.p(ip), info.RMSE, rss_in, n_in, k, info.AIC, info.BIC}; %#ok<SAGROW>
        end
    end

    % 3.2 LOTO 교차검증 (각 held-out 온도별)
    for it = 1:numel(T_ho_list)
        t_ho = T_ho_list(it);
        T_tr = T_sub(T_sub.T_age ~= t_ho, :);
        T_te = T_sub(T_sub.T_age == t_ho, :);

        loto_p = struct();
        for im = 1:numel(fit_order)
            mname = fit_order{im};
            rng(1, 'twister');
            seed_p = [];
            if strcmp(mname, 'M2') && isfield(loto_p, 'M3') && all(isfinite(loto_p.M3))
                p3 = loto_p.M3;
                seed_p = [p3(1), p3(2), p3(3), 0.5, p3(4), p3(5), p3(6), 1.0];
            elseif strcmp(mname, 'M4') && isfield(loto_p, 'M3') && all(isfinite(loto_p.M3))
                p3 = loto_p.M3;
                seed_p = [p3(1), p3(2), p3(3), p3(4), p3(5), p3(6), -3.0, 45.0];
            end
            [p_fold, ~, ok_fold] = fit_B_model_multistart(T_tr, cname, mname, n_multistart, use_lsq, seed_p, log_file, sprintf('LOTO_%s_%s_T%d', cname, mname, t_ho));
            loto_p.(mname) = p_fold;
            loto_p.([mname '_ok']) = ok_fold;
        end

        for im = 1:numel(B_models)
            mname = B_models{im};
            if loto_p.([mname '_ok'])
                pred_te = eval_B_model(loto_p.(mname), T_te, mname);
                metrics = eval_pred_metrics(T_te, cname, pred_te);
            else
                metrics = struct('RMSE', NaN, 'median_ratio', NaN, 'frac_over', NaN, 'n_rows', height(T_te), 'n_cells_last', 0);
            end
            B_loto_rows(end+1, :) = {cname, mname, t_ho, metrics.RMSE, metrics.median_ratio, metrics.frac_over, metrics.n_rows, metrics.n_cells_last}; %#ok<SAGROW>
        end
        fprintf('    LOTO held-out T=%2d degC 완료\n', t_ho);
    end
end

T_B_loto = cell2table(B_loto_rows, 'VariableNames', ...
    {'component','model','T_heldout','RMSE','median_ratio','frac_over','n_rows','n_cells_last'});
T_B_params = cell2table(B_param_rows, 'VariableNames', ...
    {'component','model','param','value','RMSE','RSS','n','k','AIC','BIC'});

writetable(T_B_loto,   fullfile(out_dir, 'B_loto.csv'));
writetable(T_B_params, fullfile(out_dir, 'B_insample_params.csv'));

% 요약표 출력: 행=모델, 열=held-out T, 셀="RMSE | median ratio"
fprintf('\n--- [표 B] LOTO 결과: RMSE (mOhm) | median(pred/meas) ---\n');
for ic = 1:numel(B_comp_names)
    cname = B_comp_names{ic};
    fprintf('\n[Component: %s]\n', cname);
    fprintf('%-6s | %-18s | %-18s | %-18s | %-18s | %-16s\n', ...
        'Model', 'T=0 degC', 'T=10 degC', 'T=25 degC', 'T=40 degC', 'In-sample RMSE/AIC');
    fprintf('%s\n', repmat('-', 1, 108));
    for im = 1:numel(B_models)
        mname = B_models{im};
        cells_str = cell(1, numel(T_ho_list));
        for it = 1:numel(T_ho_list)
            row = T_B_loto(strcmp(T_B_loto.component, cname) & strcmp(T_B_loto.model, mname) & T_B_loto.T_heldout == T_ho_list(it), :);
            cells_str{it} = sprintf('%6.3f | %6.3f', row.RMSE(1), row.median_ratio(1));
        end
        info = B_insample.(cname).(mname);
        fprintf('%-6s | %-18s | %-18s | %-18s | %-18s | %6.3f / %8.1f\n', ...
            mname, cells_str{1}, cells_str{2}, cells_str{3}, cells_str{4}, info.RMSE, info.AIC);
    end
end

%% 4 (C) WLTP
fprintf('\n[4 (C) WLTP] 정조건(age_type 1+2) 모델로 WLTP(age_type 3) 예측 중...\n');

T_wltp_all = T_rt(T_rt.age_type == 3, :);
T_wltp_temps = sort(unique(T_wltp_all.T_age));
C_rows = {};

for ic = 1:numel(B_comp_names)
    cname = B_comp_names{ic};
    T_w = T_wltp_all(isfinite(T_wltp_all.(cname)), :);

    for im = 1:numel(B_models)
        mname = B_models{im};
        info  = B_insample.(cname).(mname);

        % 기본 I_WLTP = 1.0 예측
        if info.ok && height(T_w) > 0
            pred_w = eval_B_model(info.p, T_w, mname);
            m_all  = eval_pred_metrics(T_w, cname, pred_w);
        else
            pred_w = nan(height(T_w), 1);
            m_all  = struct('RMSE', NaN, 'median_ratio', NaN, 'frac_over', NaN, 'n_rows', height(T_w), 'n_cells_last', 0);
        end
        C_rows(end+1, :) = {cname, mname, I_WLTP, 'ALL', m_all.RMSE, m_all.median_ratio, m_all.frac_over, m_all.n_rows, m_all.n_cells_last}; %#ok<SAGROW>

        % T_age별 성능
        for it = 1:numel(T_wltp_temps)
            tval = T_wltp_temps(it);
            mask_t = (T_w.T_age == tval);
            if any(mask_t) && info.ok
                m_t = eval_pred_metrics(T_w(mask_t, :), cname, pred_w(mask_t));
            else
                m_t = struct('RMSE', NaN, 'median_ratio', NaN, 'frac_over', NaN, 'n_rows', sum(mask_t), 'n_cells_last', 0);
            end
            C_rows(end+1, :) = {cname, mname, I_WLTP, sprintf('%g', tval), m_t.RMSE, m_t.median_ratio, m_t.frac_over, m_t.n_rows, m_t.n_cells_last}; %#ok<SAGROW>
        end
    end
end

% 민감도 분석: I_WLTP = 0.5, 2.0 (M2~M4)
I_sens_vals = [0.5, 2.0];
sens_summary = {};
for ic = 1:numel(B_comp_names)
    cname = B_comp_names{ic};
    T_w = T_wltp_all(isfinite(T_wltp_all.(cname)), :);
    for im = 2:4
        mname = B_models{im};
        info  = B_insample.(cname).(mname);
        med_vals = nan(1, 3);
        % nominal (1.0)
        row_nom = T_w;
        if info.ok && height(row_nom) > 0
            m_nom = eval_pred_metrics(row_nom, cname, eval_B_model(info.p, row_nom, mname));
            med_vals(2) = m_nom.median_ratio;
        end
        for isv = 1:numel(I_sens_vals)
            ival = I_sens_vals(isv);
            T_ws = T_w;
            T_ws.I(:) = ival;
            if info.ok && height(T_ws) > 0
                pred_s = eval_B_model(info.p, T_ws, mname);
                m_s    = eval_pred_metrics(T_ws, cname, pred_s);
            else
                pred_s = nan(height(T_ws), 1);
                m_s    = struct('RMSE', NaN, 'median_ratio', NaN, 'frac_over', NaN, 'n_rows', height(T_ws), 'n_cells_last', 0);
            end
            idx_store = ternary(isv == 1, 1, 3);
            med_vals(idx_store) = m_s.median_ratio;

            C_rows(end+1, :) = {cname, mname, ival, 'ALL', m_s.RMSE, m_s.median_ratio, m_s.frac_over, m_s.n_rows, m_s.n_cells_last}; %#ok<SAGROW>
            for it = 1:numel(T_wltp_temps)
                tval = T_wltp_temps(it);
                mask_t = (T_ws.T_age == tval);
                if any(mask_t) && info.ok
                    m_st = eval_pred_metrics(T_ws(mask_t, :), cname, pred_s(mask_t));
                else
                    m_st = struct('RMSE', NaN, 'median_ratio', NaN, 'frac_over', NaN, 'n_rows', sum(mask_t), 'n_cells_last', 0);
                end
                C_rows(end+1, :) = {cname, mname, ival, sprintf('%g', tval), m_st.RMSE, m_st.median_ratio, m_st.frac_over, m_st.n_rows, m_st.n_cells_last}; %#ok<SAGROW>
            end
        end
        sens_summary(end+1, :) = {cname, mname, min(med_vals, [], 'omitnan'), max(med_vals, [], 'omitnan'), med_vals(1), med_vals(2), med_vals(3)}; %#ok<SAGROW>
    end
end

T_C_out = cell2table(C_rows, 'VariableNames', ...
    {'component','model','I_WLTP','T_age','RMSE','median_ratio','frac_over','n_rows','n_cells_last'});
writetable(T_C_out, fullfile(out_dir, 'C_wltp.csv'));

fprintf('\n--- [표 C] WLTP 검증 결과 (I_WLTP = %.1f, out/C_wltp.csv) ---\n', I_WLTP);
T_C_nom = T_C_out(abs(T_C_out.I_WLTP - I_WLTP) < 1e-9, :);
fprintf('%-6s %-5s %-6s %9s %11s %10s %7s %7s\n', ...
    'Comp', 'Model', 'T_age', 'RMSE', 'MedRatio', 'FracOver', 'n_rows', 'n_cell');
for r = 1:height(T_C_nom)
    fprintf('%-6s %-5s %-6s %9.3f %11.3f %10.3f %7d %7d\n', ...
        T_C_nom.component{r}, T_C_nom.model{r}, T_C_nom.T_age{r}, ...
        T_C_nom.RMSE(r), T_C_nom.median_ratio(r), T_C_nom.frac_over(r), ...
        T_C_nom.n_rows(r), T_C_nom.n_cells_last(r));
end

fprintf('\n[참고] I_WLTP = %.1f C는 메타데이터 부재로 인한 유효 충전 C-rate 가정값입니다.\n', I_WLTP);
fprintf('[민감도] I_WLTP = 0.5 ~ 2.0 C 변화 시 전체 WLTP median(pred/meas) 범위:\n');
for r = 1:size(sens_summary, 1)
    fprintf('  %s %s: [%.3f, %.3f]  (I=0.5: %.3f | I=1.0: %.3f | I=2.0: %.3f)\n', ...
        sens_summary{r,1}, sens_summary{r,2}, sens_summary{r,3}, sens_summary{r,4}, ...
        sens_summary{r,5}, sens_summary{r,6}, sens_summary{r,7});
end

%% 5 그림
fprintf('\n[5 그림] 결과 그림 생성 및 저장 중...\n');

% 5.1 Fig A: 겉보기 Ea vs SOH (실측 binned median, A1/A2 모델 예측 곡선)
try
    figA = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 11, 4.8], 'Color', 'w');
    tlA  = tiledlayout(figA, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    colors_T = [0.00, 0.45, 0.74;   % 0 degC (파랑)
                0.20, 0.63, 0.17;   % 10 degC (초록)
                0.85, 0.33, 0.10];  % 40 degC (주황/빨강)
    soh_edges = 0.80:0.02:1.005;

    for ic = 1:2
        cname = A_comp_names{ic};
        Tp    = A_tables{ic};
        ax    = nexttile(tlA);
        hold(ax, 'on');
        grid(ax, 'on');

        % 데이터 및 A1, A2 모델 기반 겉보기 Ea 계산
        invdT_K = (1 ./ Tp.T_op) - (1 ./ Tp.T_rt);
        Ea_obs  = Rg .* log(Tp.R_op ./ Tp.R_rt) ./ invdT_K;

        pA1 = A_fits.(cname).A1.p;
        pA2 = A_fits.(cname).A2.p;
        hasA1 = all(isfinite(pA1));
        hasA2 = all(isfinite(pA2));

        if hasA1
            log_Rpred_A1 = eval_A_logRpred(pA1, Tp, 'A1');
            Ea_A1 = Rg .* (log_Rpred_A1 - log(Tp.R_rt)) ./ invdT_K;
        else
            Ea_A1 = nan(height(Tp), 1);
        end
        if hasA2
            log_Rpred_A2 = eval_A_logRpred(pA2, Tp, 'A2');
            Ea_A2 = Rg .* (log_Rpred_A2 - log(Tp.R_rt)) ./ invdT_K;
        else
            Ea_A2 = nan(height(Tp), 1);
        end

        h_leg = gobjects(0);
        l_leg = {};

        for it = 1:numel(T_age_A)
            tval = T_age_A(it);
            col  = colors_T(it, :);
            mT   = (Tp.T_age == tval) & isfinite(Ea_obs);

            [bx, by_obs] = binned_median_xy(Tp.soh(mT), Ea_obs(mT), soh_edges, 2);
            [~,  by_a1]  = binned_median_xy(Tp.soh(mT), Ea_A1(mT),  soh_edges, 2);
            [~,  by_a2]  = binned_median_xy(Tp.soh(mT), Ea_A2(mT),  soh_edges, 2);

            h1 = plot(ax, bx, by_obs, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 5.5);
            h_leg(end+1) = h1; %#ok<SAGROW>
            l_leg{end+1} = sprintf('Data T=%d^\\circC', tval); %#ok<SAGROW>

            if hasA1
                h2 = plot(ax, bx, by_a1, '--', 'Color', col, 'LineWidth', 1.4);
                if it == 1
                    h_leg(end+1) = h2; %#ok<SAGROW>
                    l_leg{end+1} = 'A1 (shared E_{new})'; %#ok<SAGROW>
                end
            end
            if hasA2
                h3 = plot(ax, bx, by_a2, '-', 'Color', col, 'LineWidth', 2.0);
                if it == 1
                    h_leg(end+1) = h3; %#ok<SAGROW>
                    l_leg{end+1} = 'A2 (T-specific E_{new})'; %#ok<SAGROW>
                end
            end
        end

        set(ax, 'XDir', 'reverse');
        xlim(ax, [0.795, 1.005]);
        xlabel(ax, 'SOH (-)');
        ylabel(ax, 'Apparent E_a (kJ/mol)');
        if strcmp(cname, 'R_arc')
            title(ax, 'R_{arc} (mid-frequency arc)');
        else
            title(ax, 'R_0 (high-frequency intercept)');
        end
        if ic == 1 && ~isempty(h_leg)
            legend(ax, h_leg, l_leg, 'Location', 'best', 'FontSize', 8);
        end
    end

    title(tlA, '(A) Apparent Activation Energy vs. SOH: Data vs. Arc Decomposition Models');
    exportgraphics(figA, fullfile(out_dir, 'fig_A_decomp.png'), 'Resolution', 150, 'BackgroundColor', 'white');
    close(figA);
    fprintf('  저장 완료: %s\n', fullfile(out_dir, 'fig_A_decomp.png'));
catch ME
    log_msg(log_file, sprintf('[Fig A] Failed: %s', ME.message));
    fprintf('  [경고] Fig A 생성 실패: %s\n', ME.message);
end

% 5.2 Fig B: LOTO held-out T별 median(pred/meas) 그룹 막대 그래프 (log y축, y=1 기준선)
try
    figB = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 10.5, 4.6], 'Color', 'w');
    tlB  = tiledlayout(figB, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    bar_colors = [0.45, 0.45, 0.45;   % M1
                  0.00, 0.45, 0.74;   % M2
                  0.93, 0.69, 0.13;   % M3
                  0.85, 0.33, 0.10];  % M4

    for ic = 1:numel(B_comp_names)
        cname = B_comp_names{ic};
        ax    = nexttile(tlB);

        M_ratio = nan(numel(T_ho_list), numel(B_models));
        for it = 1:numel(T_ho_list)
            for im = 1:numel(B_models)
                row = T_B_loto(strcmp(T_B_loto.component, cname) & ...
                               strcmp(T_B_loto.model, B_models{im}) & ...
                               T_B_loto.T_heldout == T_ho_list(it), :);
                if ~isempty(row) && isfinite(row.median_ratio(1))
                    M_ratio(it, im) = max(row.median_ratio(1), 1e-3);
                end
            end
        end

        bh = bar(ax, 1:numel(T_ho_list), M_ratio, 'grouped', 'BaseValue', 1.0);
        for im = 1:numel(B_models)
            bh(im).FaceColor = bar_colors(im, :);
        end
        set(ax, 'YScale', 'log');
        hold(ax, 'on');
        grid(ax, 'on');
        yline(ax, 1.0, 'k--', 'LineWidth', 1.2);

        valid_v = M_ratio(isfinite(M_ratio) & M_ratio > 0);
        if ~isempty(valid_v)
            ymin = min(0.2, min(valid_v) * 0.7);
            ymax = max(5.0, max(valid_v) * 1.4);
            ylim(ax, [ymin, ymax]);
        end

        xticks(ax, 1:numel(T_ho_list));
        xticklabels(ax, {'0^\circC', '10^\circC', '25^\circC', '40^\circC'});
        xlabel(ax, 'Held-out Ageing Temperature');
        ylabel(ax, 'Median (pred / meas) at last check-up');
        if strcmp(cname, 'dR0')
            title(ax, '\DeltaR_0');
        else
            title(ax, '\DeltaR_{arc}');
        end
        if ic == 1
            legend(ax, bh, B_models, 'Location', 'best', 'FontSize', 9);
        end
    end

    title(tlB, '(B) Leave-One-Temperature-Out (LOTO): Median Prediction Ratio');
    exportgraphics(figB, fullfile(out_dir, 'fig_B_loto.png'), 'Resolution', 150, 'BackgroundColor', 'white');
    close(figB);
    fprintf('  저장 완료: %s\n', fullfile(out_dir, 'fig_B_loto.png'));
catch ME
    log_msg(log_file, sprintf('[Fig B] Failed: %s', ME.message));
    fprintf('  [경고] Fig B 생성 실패: %s\n', ME.message);
end

fprintf('\n전체 분석 완료 (총 소요시간: %.1f초).\n', toc(t_start_all));

%% ========================================================================
%  Local Functions
% =========================================================================

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function log_msg(log_file, msg)
fid = fopen(log_file, 'a');
if fid > 0
    fprintf(fid, '[%s] %s\n', char(datetime('now')), msg);
    fclose(fid);
end
end

function T = sanitize_colnames(T)
vnames = T.Properties.VariableNames;
for i = 1:numel(vnames)
    vn = vnames{i};
    if endsWith(vn, 'cell') && ~strcmp(vn, 'cell')
        vnames{i} = 'cell';
    end
end
T.Properties.VariableNames = vnames;
if ismember('cell', T.Properties.VariableNames)
    T.cell = string(T.cell);
end
end

function T = ensure_numeric_cols(T, cols)
for i = 1:numel(cols)
    c = cols{i};
    if ismember(c, T.Properties.VariableNames)
        col = T.(c);
        if ~isnumeric(col)
            T.(c) = str2double(string(col));
        else
            T.(c) = double(col);
        end
    end
end
end

function T_out = preprocess_rt_data(raw, Rg, Tref, I_WLTP, SOC_CYC)
cells = unique(raw.cell, 'stable');
keep_blocks = cell(numel(cells), 1);
n_keep = 0;

for i = 1:numel(cells)
    sub = raw(raw.cell == cells(i), :);
    sub = sub(isfinite(sub.day) & isfinite(sub.soh) & isfinite(sub.T_age), :);
    if isempty(sub)
        continue;
    end
    sub = sortrows(sub, 'day');

    idx_ge99 = find(sub.soh >= 0.99);
    if isempty(idx_ge99)
        continue;
    end
    idx_bol = idx_ge99(1:min(2, numel(idx_ge99)));

    R0_bol   = mean(sub.R0(idx_bol), 'omitnan');
    Rarc_bol = mean(sub.R_arc(idx_bol), 'omitnan');
    day_bol  = sub.day(idx_bol(1));
    efc_bol  = sub.efc(idx_bol(1));

    % 첫 soh < 0.80 행까지(포함) 유지
    idx_lt80 = find(sub.soh < 0.80, 1, 'first');
    if ~isempty(idx_lt80)
        sub = sub(1:idx_lt80, :);
    end

    % 기준점 이동 및 증분 계산
    sub.dR0   = sub.R0    - R0_bol;
    sub.dRarc = sub.R_arc - Rarc_bol;
    sub.day   = max(sub.day - day_bol, 0);
    sub.efc   = max(sub.efc - efc_bol, 0);
    sub.efc(sub.age_type == 1) = 0;

    % 스트레스 변수
    TK = sub.T_age + 273.15;
    sub.TK  = TK;
    sub.x   = (1 ./ TK - 1 / Tref) / Rg;

    SOC = repmat(SOC_CYC, height(sub), 1);
    is_cal = (sub.age_type == 1);
    SOC(is_cal) = sub.cal_soc(is_cal);
    sub.SOC = SOC;

    I_eff = zeros(height(sub), 1);
    is_cyc  = (sub.age_type == 2);
    is_wltp = (sub.age_type == 3);
    I_eff(is_cyc)  = sub.chg_rate(is_cyc);
    I_eff(is_wltp) = I_WLTP;
    sub.I = I_eff;

    n_keep = n_keep + 1;
    keep_blocks{n_keep} = sub;
end

if n_keep == 0
    error('전처리 후 유효한 RT 셀 데이터가 없습니다.');
end
T_out = vertcat(keep_blocks{1:n_keep});
end

function T_out = preprocess_pairs_data(raw, comp, Rg)
if strcmp(comp, 'R_arc')
    col_op = 'R_arc_op';
    col_rt = 'R_arc_rt';
else
    col_op = 'R0_op';
    col_rt = 'R0_rt';
end

mask = (raw.soc == 50) & ismember(raw.T_age, [0, 10, 40]) & ...
       isfinite(raw.(col_op)) & isfinite(raw.(col_rt)) & ...
       (raw.(col_op) > 0) & (raw.(col_rt) > 0) & ...
       isfinite(raw.T_meas_op) & isfinite(raw.T_meas_rt) & ...
       isfinite(raw.soh) & isfinite(raw.day);
raw = raw(mask, :);

cells = unique(raw.cell, 'stable');
keep_blocks = cell(numel(cells), 1);
n_keep = 0;

for i = 1:numel(cells)
    sub = raw(raw.cell == cells(i), :);
    sub = sortrows(sub, 'day');

    idx_ge99 = find(sub.soh >= 0.99, 1, 'first');
    if isempty(idx_ge99)
        continue;
    end
    R_b = sub.(col_rt)(idx_ge99);
    if ~isfinite(R_b) || R_b <= 0
        continue;
    end

    % 첫 soh < 0.80 행까지(포함) 유지
    idx_lt80 = find(sub.soh < 0.80, 1, 'first');
    if ~isempty(idx_lt80)
        sub = sub(1:idx_lt80, :);
    end

    T_op  = sub.T_meas_op + 273.15;
    T_rt  = sub.T_meas_rt + 273.15;
    invdT = (1 ./ T_op - 1 ./ T_rt) / Rg;

    valid_T = isfinite(invdT) & (abs(invdT) > 1e-9);
    sub   = sub(valid_T, :);
    T_op  = T_op(valid_T);
    T_rt  = T_rt(valid_T);
    invdT = invdT(valid_T);
    if isempty(sub)
        continue;
    end

    R_op = sub.(col_op);
    R_rt = sub.(col_rt);
    dR   = R_rt - R_b;

    blk = table(sub.cell, sub.day, sub.soh, sub.T_age, sub.age_type, ...
                T_op, T_rt, invdT, R_op, R_rt, repmat(R_b, height(sub), 1), dR, ...
                'VariableNames', {'cell','day','soh','T_age','age_type','T_op','T_rt','invdT','R_op','R_rt','R_b','dR'});
    n_keep = n_keep + 1;
    keep_blocks{n_keep} = blk;
end

if n_keep == 0
    error('전처리 후 유효한 Paired(%s) 데이터가 없습니다.', comp);
end
T_out = vertcat(keep_blocks{1:n_keep});
end

function [log_Rpred, w_ct, w_new] = eval_A_logRpred(p, Tp, variant)
invdT = Tp.invdT;
R_b   = Tp.R_b;
dR    = Tp.dR;
n     = height(Tp);

E_ct = p(1);
switch variant
    case 'A0'
        E_new_vec = repmat(E_ct, n, 1);
    case 'A1'
        E_new_vec = repmat(p(2), n, 1);
    case 'A2'
        E_new_vec = zeros(n, 1);
        E_new_vec(Tp.T_age == 0)  = p(2);
        E_new_vec(Tp.T_age == 10) = p(3);
        E_new_vec(Tp.T_age == 40) = p(4);
end

a = log(R_b) + E_ct .* invdT;
log_Rpred = a;
w_new = zeros(n, 1);

pos = (dR > 0);
if any(pos)
    b_pos = log(dR(pos)) + E_new_vec(pos) .* invdT(pos);
    a_pos = a(pos);
    m_pos = max(a_pos, b_pos);
    log_Rpred(pos) = m_pos + log(exp(a_pos - m_pos) + exp(b_pos - m_pos));
    diff_ab = min(max(a_pos - b_pos, -50), 50);
    w_new(pos) = 1 ./ (1 + exp(diff_ab));
end
w_ct = 1 - w_new;
end

function [r, J] = res_A_model(p, Tp, variant)
[log_Rpred, w_ct, w_new] = eval_A_logRpred(p, Tp, variant);
r = log(Tp.R_op) - log_Rpred;

if nargout > 1
    invdT = Tp.invdT;
    switch variant
        case 'A0'
            J = -invdT;
        case 'A1'
            J = [-(w_ct .* invdT), -(w_new .* invdT)];
        case 'A2'
            m0  = double(Tp.T_age == 0);
            m10 = double(Tp.T_age == 10);
            m40 = double(Tp.T_age == 40);
            J = [-(w_ct .* invdT), -(m0 .* w_new .* invdT), -(m10 .* w_new .* invdT), -(m40 .* w_new .* invdT)];
    end
end
end

function [r, J] = res_A_fresh(p, Tp)
r = log(Tp.R_op ./ Tp.R_rt) - p(1) .* Tp.invdT;
if nargout > 1
    J = -Tp.invdT;
end
end

function boot_P = run_A_bootstrap(Tp, variant, p_init, lb, ub, n_boot, use_lsq, log_file)
k = numel(p_init);
boot_P = nan(n_boot, k);

% 온도군별 셀 리스트 사전 인덱싱 (층화 셀 리샘플링)
t_groups = unique(Tp.T_age);
grp_cells = cell(numel(t_groups), 1);
grp_rows  = cell(numel(t_groups), 1);
for ig = 1:numel(t_groups)
    m_g = (Tp.T_age == t_groups(ig));
    c_g = unique(Tp.cell(m_g), 'stable');
    grp_cells{ig} = c_g;
    r_g = cell(numel(c_g), 1);
    for ic = 1:numel(c_g)
        r_g{ic} = find(m_g & (Tp.cell == c_g(ic)));
    end
    grp_rows{ig} = r_g;
end

for b = 1:n_boot
    idx_list = {};
    for ig = 1:numel(t_groups)
        nc = numel(grp_cells{ig});
        draw = randi(nc, nc, 1);
        for d = 1:nc
            idx_list{end+1} = grp_rows{ig}{draw(d)}; %#ok<AGROW>
        end
    end
    idx_all = vertcat(idx_list{:});
    Tb = Tp(idx_all, :);
    try
        if strcmp(variant, 'fresh')
            obj_b = @(p) res_A_fresh(p, Tb);
        else
            obj_b = @(p) res_A_model(p, Tb, variant);
        end
        [p_b, ~, ~] = fit_nonlinear_lsq(obj_b, p_init, lb, ub, use_lsq, false);
        boot_P(b, :) = p_b;
    catch ME
        if b <= 3
            log_msg(log_file, sprintf('[Bootstrap %s b=%d] %s', variant, b, ME.message));
        end
    end
end
end

function [p_hat, rss, J_out] = fit_nonlinear_lsq(obj_fun, p0, lb, ub, use_lsq, need_J)
p0 = min(max(p0(:)', lb), ub);
J_out = [];

if use_lsq
    opts = optimoptions('lsqnonlin', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', true, ...
        'MaxIterations', 300, ...
        'MaxFunctionEvaluations', 1500, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);
    [p_hat, rss, ~, ~, ~, ~, J_sp] = lsqnonlin(obj_fun, p0, lb, ub, opts);
    p_hat = p_hat(:)';
    if need_J
        J_out = full(J_sp);
    end
else
    % fminsearch 폴백: 로짓 변환으로 경계 제약 처리
    u0 = box_to_unconstrained(p0, lb, ub);
    sse_fun = @(u) eval_sse_bounded(u, obj_fun, lb, ub);
    opts = optimset('Display', 'off', 'MaxFunEvals', 4000, 'MaxIter', 2000, 'TolX', 1e-6, 'TolFun', 1e-6);
    [u_hat, rss] = fminsearch(sse_fun, u0, opts);
    p_hat = unconstrained_to_box(u_hat, lb, ub);
end
end

function sse = eval_sse_bounded(u, obj_fun, lb, ub)
p = unconstrained_to_box(u, lb, ub);
r = obj_fun(p);
if any(~isfinite(r))
    sse = 1e30;
else
    sse = sum(r .^ 2);
end
end

function u = box_to_unconstrained(p, lb, ub)
frac = (p - lb) ./ (ub - lb);
frac = min(max(frac, 1e-6), 1 - 1e-6);
u = log(frac ./ (1 - frac));
end

function p = unconstrained_to_box(u, lb, ub)
uc = min(max(u, -35), 35);
sig = 1 ./ (1 + exp(-uc));
p = lb + (ub - lb) .* sig;
end

function [ci_lo, ci_hi] = jacobian_ci95(p, rss, J, n, k)
if isempty(J) || any(~isfinite(J(:))) || n <= k
    ci_lo = nan(1, k);
    ci_hi = nan(1, k);
    return;
end
df = max(n - k, 1);
s2 = rss / df;
covp = s2 * pinv(J' * J);
se = sqrt(max(diag(covp)', 0));
% 95% 양측 t-분포 임계값 (Hill/Cornish-Fisher 근사, base MATLAB 전용)
z = 1.959963984540054;
tval = z + (z^3 + z)/(4*df) + (5*z^5 + 16*z^3 + 3*z)/(96*df^2) + (3*z^7 + 19*z^5 + 17*z^3 - 15*z)/(384*df^3);
ci_lo = p - tval .* se;
ci_hi = p + tval .* se;
end

function q = my_prctile(x, pct)
x = x(isfinite(x));
if isempty(x)
    q = NaN;
    return;
end
x = sort(x(:));
n = numel(x);
if n == 1
    q = x(1);
    return;
end
pos = 1 + (n - 1) * (pct / 100);
lo = floor(pos);
hi = ceil(pos);
q = x(lo) + (pos - lo) * (x(hi) - x(lo));
end

function spec = get_B_model_spec(mname)
switch mname
    case 'M1'
        spec.names = {'a', 'Ea', 'z'};
        spec.lb    = [-30, -150, 0.1];
        spec.ub    = [ 20,  150, 2.0];
        spec.p_nom = [-4, 40, 0.5];
        spec.s_lo  = [-10, -20, 0.25];
        spec.s_hi  = [  0,  90, 1.30];
    case 'M2'
        spec.names = {'a1', 'E1', 's', 'z', 'a2', 'E2', 'm', 'w'};
        spec.lb    = [-30, -150, -6, 0.1, -30, -150, -1, 0.1];
        spec.ub    = [ 20,  150,  6, 2.0,  20,  150,  3, 2.0];
        spec.p_nom = [-4, 40, 0.5, 0.5, -6, 20, 0.5, 1.0];
        spec.s_lo  = [-10, -20, -2, 0.25, -12, -40, -0.5, 0.3];
        spec.s_hi  = [  0,  90,  2, 1.30,  -1,  80,  2.0, 1.6];
    case 'M3'
        spec.names = {'a1', 'E1', 's', 'a2', 'E2', 'm'};
        spec.lb    = [-30, -150, -6, -30, -150, -1];
        spec.ub    = [ 20,  150,  6,  20,  150,  3];
        spec.p_nom = [-4, 40, 0.5, -6, 20, 0.5];
        spec.s_lo  = [-10, -20, -2, -12, -40, -0.5];
        spec.s_hi  = [  0,  90,  2,  -1,  80,  2.0];
    case 'M4'
        spec.names = {'a1', 'E1', 's', 'a2', 'E2', 'm', 'c', 'Ec'};
        spec.lb    = [-30, -150, -6, -30, -150, -1, -10,   0];
        spec.ub    = [ 20,  150,  6,  20,  150,  3,  10, 150];
        spec.p_nom = [-4, 40, 0.5, -6, 30, 0.5, -3, 50];
        spec.s_lo  = [-10, -20, -2, -12, -20, -0.5, -6,  10];
        spec.s_hi  = [  0,  90,  2,  -1,  90,  2.0,  3, 110];
end
end

function [pred, J] = eval_B_model(p, T, mname)
x     = T.x;
day   = max(T.day, 0);
efc   = max(T.efc, 0);
soc_n = (T.SOC - 50) / 50;
I_eff = max(T.I, 1e-3);
n     = height(T);

switch mname
    case 'M1'
        a = p(1); Ea = p(2); z = p(3);
        u1_raw = a - Ea .* x;
        u1 = min(max(u1_raw, -50), 50);
        day_z = day .^ z;
        f1 = exp(u1) .* day_z;
        pred = f1;
        if nargout > 1
            du1 = double(u1_raw > -50 & u1_raw < 50);
            log_day = zeros(n, 1);
            md = (day > 0);
            log_day(md) = log(day(md));
            J = [f1 .* du1, -x .* f1 .* du1, f1 .* log_day];
        end

    case 'M2'
        a1 = p(1); E1 = p(2); s = p(3); z = p(4);
        a2 = p(5); E2 = p(6); m = p(7); w = p(8);

        u1_raw = a1 - E1 .* x + s .* soc_n;
        u1 = min(max(u1_raw, -50), 50);
        f1 = exp(u1) .* (day .^ z);

        u2_raw = a2 - E2 .* x;
        u2 = min(max(u2_raw, -50), 50);
        f2 = exp(u2) .* (I_eff .^ m) .* (efc .^ w);
        pred = f1 + f2;

        if nargout > 1
            du1 = double(u1_raw > -50 & u1_raw < 50);
            du2 = double(u2_raw > -50 & u2_raw < 50);
            log_day = zeros(n, 1); md = (day > 0); log_day(md) = log(day(md));
            log_efc = zeros(n, 1); me = (efc > 0); log_efc(me) = log(efc(me));
            log_I   = log(I_eff);
            J = [f1 .* du1, -x .* f1 .* du1, soc_n .* f1 .* du1, f1 .* log_day, ...
                 f2 .* du2, -x .* f2 .* du2, f2 .* log_I, f2 .* log_efc];
        end

    case 'M3'
        a1 = p(1); E1 = p(2); s = p(3);
        a2 = p(4); E2 = p(5); m = p(6);

        u1_raw = a1 - E1 .* x + s .* soc_n;
        u1 = min(max(u1_raw, -50), 50);
        f1 = exp(u1) .* sqrt(day);

        u2_raw = a2 - E2 .* x;
        u2 = min(max(u2_raw, -50), 50);
        f2 = exp(u2) .* (I_eff .^ m) .* efc;
        pred = f1 + f2;

        if nargout > 1
            du1 = double(u1_raw > -50 & u1_raw < 50);
            du2 = double(u2_raw > -50 & u2_raw < 50);
            log_I = log(I_eff);
            J = [f1 .* du1, -x .* f1 .* du1, soc_n .* f1 .* du1, ...
                 f2 .* du2, -x .* f2 .* du2, f2 .* log_I];
        end

    case 'M4'
        a1 = p(1); E1 = p(2); s = p(3);
        a2 = p(4); E2 = p(5); m = p(6); c = p(7); Ec = p(8);

        u1_raw = a1 - E1 .* x + s .* soc_n;
        u1 = min(max(u1_raw, -50), 50);
        f1 = exp(u1) .* sqrt(day);

        u2a_raw = a2 - E2 .* x;
        u2a = min(max(u2a_raw, -50), 50);
        u2b_raw = a2 + c + Ec .* x;
        u2b = min(max(u2b_raw, -50), 50);

        Im_efc = (I_eff .^ m) .* efc;
        f2a = exp(u2a) .* Im_efc;
        f2b = exp(u2b) .* Im_efc;
        f2  = f2a + f2b;
        pred = f1 + f2;

        if nargout > 1
            du1  = double(u1_raw  > -50 & u1_raw  < 50);
            du2a = double(u2a_raw > -50 & u2a_raw < 50);
            du2b = double(u2b_raw > -50 & u2b_raw < 50);
            log_I = log(I_eff);
            J = [f1 .* du1, -x .* f1 .* du1, soc_n .* f1 .* du1, ...
                 f2a .* du2a + f2b .* du2b, -x .* f2a .* du2a, f2 .* log_I, ...
                 f2b .* du2b, x .* f2b .* du2b];
        end
end
end

function [r, J] = res_B_weighted(p, T, y, w_cell, mname)
if nargout > 1
    [pred, Jp] = eval_B_model(p, T, mname);
    r = w_cell .* (pred - y);
    J = w_cell .* Jp;
else
    pred = eval_B_model(p, T, mname);
    r = w_cell .* (pred - y);
end
end

function [p_best, rss_best, ok] = fit_B_model_multistart(T_tr, cname, mname, n_starts, use_lsq, seed_p, log_file, tag)
spec = get_B_model_spec(mname);
k = numel(spec.names);
y = T_tr.(cname);

% 셀별 가중치 1/sqrt(n_rows of that cell)
[~, ~, g] = unique(T_tr.cell);
counts = accumarray(g, 1);
w_cell = 1 ./ sqrt(counts(g));

% 초기값 생성 (기준값 + 시드값 + 범위 내 랜덤 멀티스타트)
starts = zeros(n_starts, k);
starts(1, :) = spec.p_nom;
idx_next = 2;
if ~isempty(seed_p) && numel(seed_p) == k && all(isfinite(seed_p))
    starts(idx_next, :) = min(max(seed_p, spec.lb), spec.ub);
    idx_next = idx_next + 1;
end
for s = idx_next:n_starts
    if s == n_starts
        % 마지막 스타트는 전체 범위에서 넓게 샘플링
        starts(s, :) = spec.lb + rand(1, k) .* (spec.ub - spec.lb);
    else
        starts(s, :) = spec.s_lo + rand(1, k) .* (spec.s_hi - spec.s_lo);
    end
end

p_best   = nan(1, k);
rss_best = inf;
ok       = false;

obj_fun = @(p) res_B_weighted(p, T_tr, y, w_cell, mname);

for s = 1:n_starts
    try
        [p_hat, rss_hat, ~] = fit_nonlinear_lsq(obj_fun, starts(s, :), spec.lb, spec.ub, use_lsq, false);
        if isfinite(rss_hat) && rss_hat < rss_best
            p_best   = p_hat;
            rss_best = rss_hat;
            ok       = true;
        end
    catch ME
        log_msg(log_file, sprintf('[%s start %d] %s', tag, s, ME.message));
    end
end

if ~ok
    log_msg(log_file, sprintf('[%s] All %d starts failed.', tag, n_starts));
end
end

function metrics = eval_pred_metrics(T_te, cname, pred)
y = T_te.(cname);
valid = isfinite(y) & isfinite(pred);
y_v = y(valid);
p_v = pred(valid);

if isempty(y_v)
    metrics = struct('RMSE', NaN, 'median_ratio', NaN, 'frac_over', NaN, 'n_rows', 0, 'n_cells_last', 0);
    return;
end

rmse = sqrt(mean((p_v - y_v) .^ 2));

% 셀별 마지막 행 추출 (day 최대)
T_v = T_te(valid, :);
cells = unique(T_v.cell, 'stable');
meas_last = nan(numel(cells), 1);
pred_last = nan(numel(cells), 1);
for i = 1:numel(cells)
    idx = find(T_v.cell == cells(i));
    [~, imax] = max(T_v.day(idx));
    ilast = idx(imax);
    meas_last(i) = y_v(ilast);
    pred_last(i) = p_v(ilast);
end

mask_pos = (meas_last > 0.05) & isfinite(meas_last) & isfinite(pred_last);
if any(mask_pos)
    ratios = pred_last(mask_pos) ./ meas_last(mask_pos);
    med_ratio = median(ratios, 'omitnan');
    frac_over = mean(pred_last(mask_pos) > meas_last(mask_pos));
    n_cells_eval = sum(mask_pos);
else
    med_ratio = NaN;
    frac_over = NaN;
    n_cells_eval = 0;
end

metrics = struct('RMSE', rmse, 'median_ratio', med_ratio, 'frac_over', frac_over, ...
                 'n_rows', numel(y_v), 'n_cells_last', n_cells_eval);
end

function [bx, by] = binned_median_xy(x, y, edges, min_pts)
nb = numel(edges) - 1;
bx = nan(nb, 1);
by = nan(nb, 1);
for i = 1:nb
    if i < nb
        m = (x >= edges(i)) & (x < edges(i+1)) & isfinite(x) & isfinite(y);
    else
        m = (x >= edges(i)) & (x <= edges(i+1)) & isfinite(x) & isfinite(y);
    end
    if sum(m) >= min_pts
        bx(i) = median(x(m), 'omitnan');
        by(i) = median(y(m), 'omitnan');
    end
end
valid = isfinite(bx) & isfinite(by);
bx = bx(valid);
by = by(valid);
[bx, ord] = sort(bx, 'descend');
by = by(ord);
end