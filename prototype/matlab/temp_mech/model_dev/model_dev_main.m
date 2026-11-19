%% 0 설정
clearvars;
close all;

t_total = tic;
if ~exist('out', 'dir')
    mkdir('out');
end
logFile = fullfile('out', 'log.txt');
fid = fopen(logFile, 'w');
if fid ~= -1
    fclose(fid);
end

logmsg(logFile, sprintf('=== model_dev_main 시작: %s ===', char(datetime('now'))));

try
    t_sec = tic;
    JOINED = fullfile('..', 'mohtat_modes', 'mohtat_modes_joined.csv');
    TEFF = 'cell_temp.csv';

    % 색상 정의
    COL = struct();
    COL.C = [0, 0.45, 0.70];
    COL.R = [0, 0.62, 0.45];
    COL.H = [0.84, 0.37, 0];
    T_LIST = ["C", "R", "H"];

    % 조인 데이터 로드
    optsJ = detectImportOptions(JOINED);
    optsJ = setvartype(optsJ, {'T', 'dis', 'ref_not_bol'}, 'string');
    tbl_orig = readtable(JOINED, optsJ);
    tbl_orig.T = strtrim(string(tbl_orig.T));
    tbl_orig.dis = strtrim(string(tbl_orig.dis));

    % 유효 온도 데이터 로드 및 병합
    tbl_teff = readtable(TEFF);
    if ~ismember('T_eff_mean', tbl_teff.Properties.VariableNames)
        error('cell_temp.csv에 T_eff_mean 열이 없습니다.');
    end
    [tf_map, loc_map] = ismember(tbl_orig.cell, tbl_teff.cell);
    if ~all(tf_map)
        error('일부 셀의 T_eff_mean 값을 cell_temp.csv에서 찾을 수 없습니다.');
    end
    tbl_orig.T_eff = double(tbl_teff.T_eff_mean(loc_map));

    % 프로토콜 키 생성 및 정렬 (cell, ah 순)
    tbl_orig.proto = compose('%.6g|%s|%.6g', double(tbl_orig.chg), tbl_orig.dis, double(tbl_orig.dod));
    tbl_orig = sortrows(tbl_orig, {'cell', 'ah'});

    n_all_cells = numel(unique(tbl_orig.cell));
    logmsg(logFile, sprintf('[0 설정] 데이터 로드 완료: %d행, %d셀 (%.2f s)', ...
        height(tbl_orig), n_all_cells, toc(t_sec)));
catch ME
    logmsg(logFile, sprintf('[0 설정 오류] %s', getReport(ME, 'extended', 'hyperlinks', 'off')));
    return;
end

%% D1 기준점 민감도 분석 (Reference-point sensitivity)
try
    t_sec = tic;
    logmsg(logFile, '--- [D1] 기준점 민감도 분석 시작 ---');

    % 1) shifted 데이터셋 생성
    h_rows = tbl_orig(tbl_orig.T == "H", :);
    h_cells = unique(h_rows.cell, 'stable');
    h_proto = strings(numel(h_cells), 1);
    h_Aref = nan(numel(h_cells), 1);
    for i = 1:numel(h_cells)
        subH = h_rows(h_rows.cell == h_cells(i), :);
        h_proto(i) = subH.proto(1);
        idx0 = find(abs(subH.loss) < 1e-12, 1, 'first');
        if isempty(idx0)
            idx0 = 1;
        end
        h_Aref(i) = subH.ah(idx0);
    end

    all_cells = unique(tbl_orig.cell, 'stable');
    shifted_parts = cell(numel(all_cells), 1);
    keep_cnt = 0;
    dropped_cells = [];

    for i = 1:numel(all_cells)
        cid = all_cells(i);
        sub = tbl_orig(tbl_orig.cell == cid, :);
        sub = sortrows(sub, 'ah');
        t_val = sub.T(1);

        if t_val == "H"
            keep_cnt = keep_cnt + 1;
            shifted_parts{keep_cnt} = sub;
        else
            p_key = sub.proto(1);
            idx_h = find(h_proto == p_key, 1, 'first');
            if isempty(idx_h)
                logmsg(logFile, sprintf('  [D1 경고] 셀 %d (%s, %s): 동일 프로토콜 H 셀 없음 -> 제외', cid, t_val, p_key));
                dropped_cells(end+1) = cid; %#ok<SAGROW>
                continue;
            end
            A_ref = h_Aref(idx_h);
            idx_ref = find(sub.ah >= A_ref, 1, 'first');
            if isempty(idx_ref)
                logmsg(logFile, sprintf('  [D1 제외] 셀 %d (%s, %s): A_ref=%.2f Ah 이상 체크업 없음 (max ah=%.2f) -> 제외', ...
                    cid, t_val, p_key, A_ref, max(sub.ah)));
                dropped_cells(end+1) = cid; %#ok<SAGROW>
                continue;
            end

            loss_ref = sub.loss(idx_ref);
            LLI_ref = sub.LLI(idx_ref);
            LAM_PE_ref = sub.LAM_PE(idx_ref);
            LAM_NE_ref = sub.LAM_NE(idx_ref);
            dR0_ref = sub.dR0(idx_ref);
            dRarc_ref = sub.dRarc(idx_ref);
            days_ref = sub.days(idx_ref);

            sub2 = sub(idx_ref:end, :);
            sub2.loss = 100 * (1 - (1 - sub2.loss / 100) ./ (1 - loss_ref / 100));
            sub2.LLI = sub2.LLI - LLI_ref;
            sub2.LAM_PE = sub2.LAM_PE - LAM_PE_ref;
            sub2.LAM_NE = sub2.LAM_NE - LAM_NE_ref;
            sub2.dR0 = sub2.dR0 - dR0_ref;
            sub2.dRarc = sub2.dRarc - dRarc_ref;
            sub2.days = sub2.days - days_ref;

            % 새 기준 행의 수치 오차 제거
            sub2.loss(1) = 0;
            sub2.LLI(1) = 0;
            sub2.LAM_PE(1) = 0;
            sub2.LAM_NE(1) = 0;
            sub2.dR0(1) = 0;
            sub2.dRarc(1) = 0;
            sub2.days(1) = 0;

            keep_cnt = keep_cnt + 1;
            shifted_parts{keep_cnt} = sub2;
            logmsg(logFile, sprintf('  [D1 이동] 셀 %d (%s): A_ref=%.2f -> 새 기준 ah=%.2f (loss_ref=%.2f%%), 잔여행=%d', ...
                cid, t_val, A_ref, sub.ah(idx_ref), loss_ref, height(sub2)));
        end
    end

    tbl_shifted = vertcat(shifted_parts{1:keep_cnt});
    logmsg(logFile, sprintf('  [D1] shifted 데이터셋 완료: %d행, %d셀 (제외 셀 수: %d)', ...
        height(tbl_shifted), numel(unique(tbl_shifted.cell)), numel(dropped_cells)));

    % 2) M4 및 E2 적합 + 셀 부트스트랩 (B=1000, rng(1))
    B_D1 = 1000;
    rng(1);

    versions = ["orig", "shifted"];
    datasets = {tbl_orig, tbl_shifted};
    targets_M4 = ["dR0", "dRarc"];
    modes_E2 = ["LLI", "LAM_PE", "LAM_NE"];
    all_qty_D1 = [targets_M4, modes_E2];

    res_ver = strings(0,1);
    res_qty = strings(0,1);
    res_T = strings(0,1);
    res_est = zeros(0,1);
    res_lo = zeros(0,1);
    res_hi = zeros(0,1);
    res_nc = zeros(0,1);
    res_nr = zeros(0,1);

    dif_ver = strings(0,1);
    dif_qty = strings(0,1);
    dif_pair = strings(0,1);
    dif_est = zeros(0,1);
    dif_lo = zeros(0,1);
    dif_hi = zeros(0,1);
    dif_ex0 = false(0,1);

    for vIdx = 1:2
        vName = versions(vIdx);
        ds = datasets{vIdx};

        for qIdx = 1:numel(all_qty_D1)
            qName = all_qty_D1(qIdx);
            isM4 = any(qName == targets_M4);

            est_T = nan(1, 3);
            boot_T = nan(B_D1, 3);

            for tIdx = 1:3
                tName = T_LIST(tIdx);
                if isM4
                    ycol = double(ds.(char(qName)));
                    mask = (ds.T == tName) & (ds.loss > 0) & isfinite(ds.loss) & isfinite(ycol);
                    subT = ds(mask, :);
                    yT = ycol(mask);
                    [est_val, boot_vals, n_c, n_r] = fit_M4_and_boot(subT.cell, subT.loss, yT, B_D1);
                else
                    mcol = double(ds.(char(qName)));
                    mask = (ds.T == tName) & (mcol > 0.05) & isfinite(mcol) & (ds.days > 0) & isfinite(ds.days);
                    subT = ds(mask, :);
                    mT = mcol(mask);
                    [est_val, boot_vals, n_c, n_r] = fit_E2_and_boot(subT.cell, subT.days, mT, B_D1);
                end

                est_T(tIdx) = est_val;
                boot_T(:, tIdx) = boot_vals;

                ci_lo = my_prctile(boot_vals, 2.5);
                ci_hi = my_prctile(boot_vals, 97.5);

                res_ver(end+1,1) = vName; %#ok<SAGROW>
                res_qty(end+1,1) = qName; %#ok<SAGROW>
                res_T(end+1,1) = tName; %#ok<SAGROW>
                res_est(end+1,1) = est_val; %#ok<SAGROW>
                res_lo(end+1,1) = ci_lo; %#ok<SAGROW>
                res_hi(end+1,1) = ci_hi; %#ok<SAGROW>
                res_nc(end+1,1) = n_c; %#ok<SAGROW>
                res_nr(end+1,1) = n_r; %#ok<SAGROW>
            end

            % 온도간 차이 분포 (C-H, R-H)
            pairs = {"C-H", 1, 3; "R-H", 2, 3};
            for pIdx = 1:2
                pLabel = string(pairs{pIdx, 1});
                i1 = pairs{pIdx, 2};
                i2 = pairs{pIdx, 3};
                d_est = est_T(i1) - est_T(i2);
                d_boot = boot_T(:, i1) - boot_T(:, i2);
                d_lo = my_prctile(d_boot, 2.5);
                d_hi = my_prctile(d_boot, 97.5);
                ex0 = (d_lo > 0) || (d_hi < 0);

                dif_ver(end+1,1) = vName; %#ok<SAGROW>
                dif_qty(end+1,1) = qName; %#ok<SAGROW>
                dif_pair(end+1,1) = pLabel; %#ok<SAGROW>
                dif_est(end+1,1) = d_est; %#ok<SAGROW>
                dif_lo(end+1,1) = d_lo; %#ok<SAGROW>
                dif_hi(end+1,1) = d_hi; %#ok<SAGROW>
                dif_ex0(end+1,1) = ex0; %#ok<SAGROW>
            end
        end
    end

    tbl_D1_res = table(res_ver, res_qty, res_T, res_est, res_lo, res_hi, res_nc, res_nr, ...
        'VariableNames', {'version', 'target_mode', 'T', 'estimate', 'ci_lo', 'ci_hi', 'n_cells', 'n_rows'});
    tbl_D1_dif = table(dif_ver, dif_qty, dif_pair, dif_est, dif_lo, dif_hi, dif_ex0, ...
        'VariableNames', {'version', 'quantity', 'pair', 'diff', 'ci_lo', 'ci_hi', 'excludes0'});

    writetable(tbl_D1_res, fullfile('out', 'D1_ref_sensitivity.csv'));
    writetable(tbl_D1_dif, fullfile('out', 'D1_differences.csv'));

    logmsg(logFile, sprintf('[D1] 완료 및 저장 (%.2f s)', toc(t_sec)));
catch ME
    logmsg(logFile, sprintf('[D1 오류] %s', getReport(ME, 'extended', 'hyperlinks', 'off')));
end

%% D2 s에 대한 구조화 모델 (Structured models for s)
try
    t_sec = tic;
    logmsg(logFile, '--- [D2] 구조화 모델 적합 시작 ---');

    targets_D2 = ["dR0", "dRarc"];
    model_names = ["S0", "S1", "S2", "S3", "S4"];
    B_D2 = 500;
    rng(2);

    d2_target = strings(0,1);
    d2_model = strings(0,1);
    d2_params = strings(0,1);
    d2_n = zeros(0,1);
    d2_k = zeros(0,1);
    d2_rss = zeros(0,1);
    d2_aic = zeros(0,1);
    d2_r2 = zeros(0,1);
    d2_loco_rmse = zeros(0,1);
    d2_ratio_med = zeros(0,1);

    ci_target = strings(0,1);
    ci_model = strings(0,1);
    ci_param = strings(0,1);
    ci_est = zeros(0,1);
    ci_lo = zeros(0,1);
    ci_hi = zeros(0,1);
    ci_ex0 = false(0,1);

    % 그림용 저장 구조체
    D2_plot = struct();

    for tIdx = 1:numel(targets_D2)
        tName = targets_D2(tIdx);
        ycol = double(tbl_orig.(char(tName)));
        mask = (tbl_orig.loss > 0) & isfinite(tbl_orig.loss) & isfinite(ycol);
        sub = tbl_orig(mask, :);
        y = ycol(mask);
        L = double(sub.loss);
        T_idx = zeros(height(sub), 1);
        T_idx(sub.T == "C") = 1;
        T_idx(sub.T == "R") = 2;
        T_idx(sub.T == "H") = 3;
        T_amb = double(sub.T_C);
        T_eff = double(sub.T_eff);
        log_chg = log(double(sub.chg));
        cells_vec = sub.cell;
        ah_vec = double(sub.ah);

        u_cells = unique(cells_vec, 'stable');
        n_cells = numel(u_cells);
        n_obs = numel(y);
        SST = sum((y - mean(y)).^2);

        % 셀 인덱스 매핑 (1..n_cells)
        [~, cell_id] = ismember(cells_vec, u_cells);
        cell_T = zeros(n_cells, 1);
        for c = 1:n_cells
            cell_T(c) = T_idx(find(cell_id == c, 1, 'first'));
        end

        % --- 전체 데이터 적합 ---
        % S0 (= M4)
        [s0_pT, s0_sCell, s0_rss] = fit_D2_S0(cell_id, T_idx, L, y, n_cells);
        k_S0 = n_cells + 3;
        aic_S0 = n_obs * log(s0_rss / n_obs) + 2 * k_S0;
        r2_S0 = 1 - s0_rss / SST;
        pstr_S0 = sprintf('p_C=%.4f; p_R=%.4f; p_H=%.4f; s_cell(%d)', s0_pT(1), s0_pT(2), s0_pT(3), n_cells);

        % S1
        [th_S1, rss_S1] = fit_D2_S1(T_idx, L, y, s0_pT);
        k_S1 = 6;
        aic_S1 = n_obs * log(rss_S1 / n_obs) + 2 * k_S1;
        r2_S1 = 1 - rss_S1 / SST;
        pstr_S1 = sprintf('a_C=%.4f; a_R=%.4f; a_H=%.4f; p_C=%.4f; p_R=%.4f; p_H=%.4f', ...
            th_S1(1), th_S1(2), th_S1(3), th_S1(4), th_S1(5), th_S1(6));

        % S2
        [th_S2, rss_S2] = fit_D2_S2(T_idx, log_chg, L, y, th_S1, s0_pT, s0_sCell, cell_id, n_cells);
        k_S2 = 7;
        aic_S2 = n_obs * log(rss_S2 / n_obs) + 2 * k_S2;
        r2_S2 = 1 - rss_S2 / SST;
        pstr_S2 = sprintf('a_C=%.4f; a_R=%.4f; a_H=%.4f; alpha=%.4f; p_C=%.4f; p_R=%.4f; p_H=%.4f', ...
            th_S2(1), th_S2(2), th_S2(3), th_S2(4), th_S2(5), th_S2(6), th_S2(7));

        % S3 (T_amb)
        [th_S3, rss_S3] = fit_D2_S34(T_amb, log_chg, L, y, th_S2, th_S1, s0_pT, s0_sCell, cell_id, n_cells, []);
        k_S3 = 5;
        aic_S3 = n_obs * log(rss_S3 / n_obs) + 2 * k_S3;
        r2_S3 = 1 - rss_S3 / SST;
        pstr_S3 = sprintf('a0=%.4f; beta=%.5f; alpha=%.4f; b0=%.4f; b1=%.5f', ...
            th_S3(1), th_S3(2), th_S3(3), th_S3(4), th_S3(5));

        % S4 (T_eff)
        [th_S4, rss_S4] = fit_D2_S34(T_eff, log_chg, L, y, th_S2, th_S1, s0_pT, s0_sCell, cell_id, n_cells, th_S3);
        k_S4 = 5;
        aic_S4 = n_obs * log(rss_S4 / n_obs) + 2 * k_S4;
        r2_S4 = 1 - rss_S4 / SST;
        pstr_S4 = sprintf('a0=%.4f; beta=%.5f; alpha=%.4f; b0=%.4f; b1=%.5f', ...
            th_S4(1), th_S4(2), th_S4(3), th_S4(4), th_S4(5));

        % --- LOCO (S1-S4) ---
        loco_rmse = nan(1, 5);
        ratio_med = nan(1, 5);
        loco_preds_all = nan(n_obs, 4); % S1..S4

        for mIdx = 1:4
            mName = model_names(mIdx + 1);
            res_all = nan(n_obs, 1);
            last_ratios = nan(n_cells, 1);

            for c = 1:n_cells
                tr = (cell_id ~= c);
                te = (cell_id == c);

                switch mName
                    case "S1"
                        th_c = fit_D2_S1_fast(T_idx(tr), L(tr), y(tr), th_S1(4:6));
                        yhat_c = predict_D2("S1", th_c, T_idx(te), T_amb(te), T_eff(te), log_chg(te), L(te));
                    case "S2"
                        th_c = fit_D2_S2_fast(T_idx(tr), log_chg(tr), L(tr), y(tr), ones(sum(tr),1), th_S2);
                        yhat_c = predict_D2("S2", th_c, T_idx(te), T_amb(te), T_eff(te), log_chg(te), L(te));
                    case "S3"
                        th_c = fit_D2_S34_fast(T_amb(tr), log_chg(tr), L(tr), y(tr), ones(sum(tr),1), th_S3);
                        yhat_c = predict_D2("S3", th_c, T_idx(te), T_amb(te), T_eff(te), log_chg(te), L(te));
                    case "S4"
                        th_c = fit_D2_S34_fast(T_eff(tr), log_chg(tr), L(tr), y(tr), ones(sum(tr),1), th_S4);
                        yhat_c = predict_D2("S4", th_c, T_idx(te), T_amb(te), T_eff(te), log_chg(te), L(te));
                end

                loco_preds_all(te, mIdx) = yhat_c;
                res_all(te) = y(te) - yhat_c;

                % 마지막 체크업에서의 pred/meas
                te_ah = ah_vec(te);
                te_y = y(te);
                [~, idx_last] = max(te_ah);
                if te_y(idx_last) > 0 && isfinite(yhat_c(idx_last))
                    last_ratios(c) = yhat_c(idx_last) / te_y(idx_last);
                end
            end

            loco_rmse(mIdx + 1) = sqrt(mean(res_all.^2));
            ratio_med(mIdx + 1) = median(last_ratios, 'omitnan');
        end

        % 결과 테이블 누적
        pstrs = [string(pstr_S0); string(pstr_S1); string(pstr_S2); string(pstr_S3); string(pstr_S4)];
        ks = [k_S0; k_S1; k_S2; k_S3; k_S4];
        rsss = [s0_rss; rss_S1; rss_S2; rss_S3; rss_S4];
        aics = [aic_S0; aic_S1; aic_S2; aic_S3; aic_S4];
        r2s = [r2_S0; r2_S1; r2_S2; r2_S3; r2_S4];

        for mIdx = 1:5
            d2_target(end+1,1) = tName; %#ok<SAGROW>
            d2_model(end+1,1) = model_names(mIdx); %#ok<SAGROW>
            d2_params(end+1,1) = pstrs(mIdx); %#ok<SAGROW>
            d2_n(end+1,1) = n_obs; %#ok<SAGROW>
            d2_k(end+1,1) = ks(mIdx); %#ok<SAGROW>
            d2_rss(end+1,1) = rsss(mIdx); %#ok<SAGROW>
            d2_aic(end+1,1) = aics(mIdx); %#ok<SAGROW>
            d2_r2(end+1,1) = r2s(mIdx); %#ok<SAGROW>
            d2_loco_rmse(end+1,1) = loco_rmse(mIdx); %#ok<SAGROW>
            d2_ratio_med(end+1,1) = ratio_med(mIdx); %#ok<SAGROW>
        end

        % 사전 등록된 점검 출력 (S2, S3, S4 vs S1 <= 0.90, S4 vs S3)
        r_S2_S1 = loco_rmse(3) / loco_rmse(2);
        r_S3_S1 = loco_rmse(4) / loco_rmse(2);
        r_S4_S1 = loco_rmse(5) / loco_rmse(2);
        r_S4_S3 = loco_rmse(5) / loco_rmse(4);
        logmsg(logFile, sprintf('  [D2 사전등록 점검 - %s]', tName));
        logmsg(logFile, sprintf('    LOCO RMSE: S1=%.4f, S2=%.4f, S3=%.4f, S4=%.4f', ...
            loco_rmse(2), loco_rmse(3), loco_rmse(4), loco_rmse(5)));
        logmsg(logFile, sprintf('    S2/S1 = %.4f (<=0.90: %d) | S3/S1 = %.4f (<=0.90: %d) | S4/S1 = %.4f (<=0.90: %d) | S4/S3 = %.4f (S4<S3: %d)', ...
            r_S2_S1, r_S2_S1 <= 0.90, r_S3_S1, r_S3_S1 <= 0.90, r_S4_S1, r_S4_S1 <= 0.90, r_S4_S3, r_S4_S3 < 1.0));

        % 그림용 데이터 저장
        [~, best_idx] = min(loco_rmse(2:5));
        D2_plot.(char(tName)).loco_rmse = loco_rmse(2:5);
        D2_plot.(char(tName)).aic = aics(2:5)';
        D2_plot.(char(tName)).best_model = model_names(best_idx + 1);
        D2_plot.(char(tName)).y_meas = y;
        D2_plot.(char(tName)).y_pred = loco_preds_all(:, best_idx);
        D2_plot.(char(tName)).T = sub.T;

        % --- 셀 부트스트랩 (500, rng(2), 온도별 층화 리샘플링) ---
        boot_S2_alpha = nan(B_D2, 1);
        boot_S3_par = nan(B_D2, 3); % alpha, beta, b1
        boot_S4_par = nan(B_D2, 3); % alpha, beta, b1

        cells_by_T = {find(cell_T == 1), find(cell_T == 2), find(cell_T == 3)};

        for b = 1:B_D2
            w_cell = zeros(n_cells, 1);
            for ti = 1:3
                pool = cells_by_T{ti};
                draws = pool(randi(numel(pool), numel(pool), 1));
                for d = 1:numel(draws)
                    w_cell(draws(d)) = w_cell(draws(d)) + 1;
                end
            end
            w_row = w_cell(cell_id);

            th2_b = fit_D2_S2_fast(T_idx, log_chg, L, y, w_row, th_S2);
            boot_S2_alpha(b) = th2_b(4);

            th3_b = fit_D2_S34_fast(T_amb, log_chg, L, y, w_row, th_S3);
            boot_S3_par(b, :) = [th3_b(3), th3_b(2), th3_b(5)];

            th4_b = fit_D2_S34_fast(T_eff, log_chg, L, y, w_row, th_S4);
            boot_S4_par(b, :) = [th4_b(3), th4_b(2), th4_b(5)];
        end

        % CI 기록
        ci_entries = {
            "S2", "alpha", th_S2(4), boot_S2_alpha;
            "S3", "alpha", th_S3(3), boot_S3_par(:,1);
            "S3", "beta",  th_S3(2), boot_S3_par(:,2);
            "S3", "b1",    th_S3(5), boot_S3_par(:,3);
            "S4", "alpha", th_S4(3), boot_S4_par(:,1);
            "S4", "beta",  th_S4(2), boot_S4_par(:,2);
            "S4", "b1",    th_S4(5), boot_S4_par(:,3)
        };
        for e = 1:size(ci_entries, 1)
            lo_e = my_prctile(ci_entries{e,4}, 2.5);
            hi_e = my_prctile(ci_entries{e,4}, 97.5);
            ci_target(end+1,1) = tName; %#ok<SAGROW>
            ci_model(end+1,1) = ci_entries{e,1}; %#ok<SAGROW>
            ci_param(end+1,1) = ci_entries{e,2}; %#ok<SAGROW>
            ci_est(end+1,1) = ci_entries{e,3}; %#ok<SAGROW>
            ci_lo(end+1,1) = lo_e; %#ok<SAGROW>
            ci_hi(end+1,1) = hi_e; %#ok<SAGROW>
            ci_ex0(end+1,1) = (lo_e > 0) || (hi_e < 0); %#ok<SAGROW>
        end
    end

    tbl_D2_models = table(d2_target, d2_model, d2_params, d2_k, d2_n, d2_rss, d2_aic, d2_r2, d2_loco_rmse, d2_ratio_med, ...
        'VariableNames', {'target', 'model', 'params', 'k', 'n', 'RSS', 'AIC', 'R2', 'LOCO_RMSE', 'ratio_median'});
    tbl_D2_ci = table(ci_target, ci_model, ci_param, ci_est, ci_lo, ci_hi, ci_ex0, ...
        'VariableNames', {'target', 'model', 'param', 'estimate', 'ci_lo', 'ci_hi', 'excludes0'});

    writetable(tbl_D2_models, fullfile('out', 'D2_models.csv'));
    writetable(tbl_D2_ci, fullfile('out', 'D2_params_ci.csv'));

    logmsg(logFile, sprintf('[D2] 완료 및 저장 (%.2f s)', toc(t_sec)));
catch ME
    logmsg(logFile, sprintf('[D2 오류] %s', getReport(ME, 'extended', 'hyperlinks', 'off')));
end

%% D3 가속시험 전이 (hot -> room/cold)
try
    t_sec = tic;
    logmsg(logFile, '--- [D3] 가속시험 전이 분석 시작 ---');

    versions = ["orig", "shifted"];
    datasets = {tbl_orig, tbl_shifted};
    modes_D3 = ["LLI", "LAM_PE", "LAM_NE"];
    res_D3 = ["dR0", "dRarc"];
    qtys_D3 = [modes_D3, res_D3];
    levels_D3 = [5, 10];
    T_eval = ["C", "R"];

    % 셀별 결과
    c_ver = strings(0,1);
    c_T = strings(0,1);
    c_cell = zeros(0,1);
    c_qty = strings(0,1);
    c_lvl = zeros(0,1);
    c_meas = zeros(0,1);
    c_pred = zeros(0,1);
    c_rat = zeros(0,1);

    % 요약 결과
    s_ver = strings(0,1);
    s_T = strings(0,1);
    s_qty = strings(0,1);
    s_lvl = zeros(0,1);
    s_n = zeros(0,1);
    s_med = zeros(0,1);
    s_iqr = zeros(0,1);
    s_q25 = zeros(0,1);
    s_q75 = zeros(0,1);
    s_rep = false(0,1);

    for vIdx = 1:2
        vName = versions(vIdx);
        ds = datasets{vIdx};
        subH = ds(ds.T == "H", :);

        % H 셀 적합
        H_fits = struct();
        for mIdx = 1:numel(modes_D3)
            q = modes_D3(mIdx);
            mv = double(subH.(char(q)));
            maskH = (subH.loss > 0) & isfinite(subH.loss) & isfinite(mv);
            LH = double(subH.loss(maskH));
            mH = mv(maskH);
            kH = sum(LH .* mH) / sum(LH .^ 2);
            H_fits.(char(q)) = struct('type', 'linear', 'k', kH);
            logmsg(logFile, sprintf('  [D3 H적합 %s] %s: k=%.4f', vName, q, kH));
        end
        for rIdx = 1:numel(res_D3)
            q = res_D3(rIdx);
            rv = double(subH.(char(q)));
            maskH = (subH.loss > 0) & isfinite(subH.loss) & isfinite(rv);
            LH = double(subH.loss(maskH));
            rH = rv(maskH);
            [sH, pH] = fit_pooled_power(LH, rH);
            H_fits.(char(q)) = struct('type', 'power', 's', sH, 'p', pH);
            logmsg(logFile, sprintf('  [D3 H적합 %s] %s: s=%.5f, p=%.4f', vName, q, sH, pH));
        end

        % C, R 셀 평가
        for tIdx = 1:numel(T_eval)
            tName = T_eval(tIdx);
            subT = ds(ds.T == tName, :);
            u_cells = unique(subT.cell, 'stable');

            for qIdx = 1:numel(qtys_D3)
                q = qtys_D3(qIdx);
                fitObj = H_fits.(char(q));

                for lIdx = 1:numel(levels_D3)
                    lvl = levels_D3(lIdx);
                    if strcmp(fitObj.type, 'linear')
                        pred_val = fitObj.k * lvl;
                    else
                        pred_val = fitObj.s * (lvl ^ fitObj.p);
                    end

                    ratios_grp = zeros(0,1);
                    for ci = 1:numel(u_cells)
                        cid = u_cells(ci);
                        sc = subT(subT.cell == cid, :);
                        sc = sortrows(sc, 'ah');
                        Lc = double(sc.loss);
                        yc = double(sc.(char(q)));
                        valid = isfinite(Lc) & isfinite(yc);
                        Lc = Lc(valid);
                        yc = yc(valid);
                        if numel(Lc) < 2
                            continue;
                        end
                        L_mon = cummax(Lc);
                        if max(L_mon) < lvl || min(L_mon) > lvl
                            continue;
                        end
                        [L_u, ia] = unique(L_mon, 'stable');
                        y_u = yc(ia);
                        if numel(L_u) < 2
                            continue;
                        end
                        meas_val = interp1(L_u, y_u, lvl, 'linear');
                        if ~isfinite(meas_val) || meas_val <= 0
                            continue;
                        end
                        rat = pred_val / meas_val;
                        ratios_grp(end+1,1) = rat; %#ok<SAGROW>

                        c_ver(end+1,1) = vName; %#ok<SAGROW>
                        c_T(end+1,1) = tName; %#ok<SAGROW>
                        c_cell(end+1,1) = cid; %#ok<SAGROW>
                        c_qty(end+1,1) = q; %#ok<SAGROW>
                        c_lvl(end+1,1) = lvl; %#ok<SAGROW>
                        c_meas(end+1,1) = meas_val; %#ok<SAGROW>
                        c_pred(end+1,1) = pred_val; %#ok<SAGROW>
                        c_rat(end+1,1) = rat; %#ok<SAGROW>
                    end

                    n_g = numel(ratios_grp);
                    if n_g > 0
                        med_g = median(ratios_grp);
                        q25_g = my_prctile(ratios_grp, 25);
                        q75_g = my_prctile(ratios_grp, 75);
                        iqr_g = q75_g - q25_g;
                        rep_g = (med_g >= 0.8) && (med_g <= 1.25);
                    else
                        med_g = NaN;
                        q25_g = NaN;
                        q75_g = NaN;
                        iqr_g = NaN;
                        rep_g = false;
                    end

                    s_ver(end+1,1) = vName; %#ok<SAGROW>
                    s_T(end+1,1) = tName; %#ok<SAGROW>
                    s_qty(end+1,1) = q; %#ok<SAGROW>
                    s_lvl(end+1,1) = lvl; %#ok<SAGROW>
                    s_n(end+1,1) = n_g; %#ok<SAGROW>
                    s_med(end+1,1) = med_g; %#ok<SAGROW>
                    s_iqr(end+1,1) = iqr_g; %#ok<SAGROW>
                    s_q25(end+1,1) = q25_g; %#ok<SAGROW>
                    s_q75(end+1,1) = q75_g; %#ok<SAGROW>
                    s_rep(end+1,1) = rep_g; %#ok<SAGROW>
                end
            end
        end
    end

    tbl_D3_cells = table(c_ver, c_T, c_cell, c_qty, c_lvl, c_meas, c_pred, c_rat, ...
        'VariableNames', {'version', 'T', 'cell', 'quantity', 'level', 'meas', 'pred', 'ratio'});
    tbl_D3_sum = table(s_ver, s_T, s_qty, s_lvl, s_n, s_med, s_iqr, s_q25, s_q75, s_rep, ...
        'VariableNames', {'version', 'T', 'quantity', 'level', 'n', 'median_ratio', 'IQR', 'q25', 'q75', 'reproduced'});

    writetable(tbl_D3_sum, fullfile('out', 'D3_transfer.csv'));
    writetable(tbl_D3_cells, fullfile('out', 'D3_transfer_cells.csv'));

    logmsg(logFile, sprintf('[D3] 완료 및 저장 (%.2f s)', toc(t_sec)));
catch ME
    logmsg(logFile, sprintf('[D3 오류] %s', getReport(ME, 'extended', 'hyperlinks', 'off')));
end

%% 그림 생성 (Figures)
try
    t_sec = tic;
    logmsg(logFile, '--- [그림] 생성 시작 ---');

    % --- Fig D1: 2x1 p_T (R0, arc) & z_T (3 modes) ---
    fig1 = figure('Color', 'w', 'Position', [100, 100, 960, 720], 'Visible', 'off');
    tl1 = tiledlayout(fig1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 상단: p_T
    ax1 = nexttile(tl1);
    hold(ax1, 'on');
    grid(ax1, 'on');
    p_qtys = ["dR0", "dRarc"];
    p_labels = ["\DeltaR_0", "\DeltaR_{arc}"];
    base_x = [1, 2, 3, 4.6, 5.6, 6.6];
    x_ticks_1 = base_x;
    x_lbls_1 = ["C", "R", "H", "C", "R", "H"];

    for qIdx = 1:2
        for tIdx = 1:3
            x0 = base_x((qIdx - 1)*3 + tIdx);
            tName = T_LIST(tIdx);
            c_rgb = COL.(char(tName));
            for vIdx = 1:2
                vName = versions(vIdx);
                row = tbl_D1_res(tbl_D1_res.version == vName & tbl_D1_res.target_mode == p_qtys(qIdx) & tbl_D1_res.T == tName, :);
                if isempty(row), continue; end
                dx = -0.14 + 0.28 * (vIdx - 1);
                est = row.estimate(1);
                neg = est - row.ci_lo(1);
                pos = row.ci_hi(1) - est;
                if vIdx == 1
                    errorbar(ax1, x0 + dx, est, neg, pos, 'o', 'Color', c_rgb, ...
                        'MarkerFaceColor', c_rgb, 'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 5);
                else
                    errorbar(ax1, x0 + dx, est, neg, pos, 's', 'Color', c_rgb, ...
                        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 5);
                end
            end
        end
    end
    xline(ax1, 3.8, '--', 'Color', [0.5 0.5 0.5]);
    xticks(ax1, x_ticks_1);
    xticklabels(ax1, x_lbls_1);
    xlim(ax1, [0.3, 7.3]);
    yl1 = ylim(ax1);
    text(ax1, 2.0, yl1(2) - 0.06*(yl1(2)-yl1(1)), p_labels(1), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    text(ax1, 5.6, yl1(2) - 0.06*(yl1(2)-yl1(1)), p_labels(2), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    ylabel(ax1, 'p_T (M4 exponent)');
    title(ax1, 'D1 Reference-point sensitivity: Resistance exponent p_T (orig vs shifted)');

    % 범례용 더미 핸들
    h_leg = gobjects(5,1);
    h_leg(1) = plot(ax1, nan, nan, 'o', 'Color', COL.C, 'MarkerFaceColor', COL.C, 'LineWidth', 1.4);
    h_leg(2) = plot(ax1, nan, nan, 'o', 'Color', COL.R, 'MarkerFaceColor', COL.R, 'LineWidth', 1.4);
    h_leg(3) = plot(ax1, nan, nan, 'o', 'Color', COL.H, 'MarkerFaceColor', COL.H, 'LineWidth', 1.4);
    h_leg(4) = plot(ax1, nan, nan, 'o', 'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2], 'LineWidth', 1.4);
    h_leg(5) = plot(ax1, nan, nan, 's', 'Color', [0.2 0.2 0.2], 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
    legend(ax1, h_leg, {'C (-5 °C)', 'R (25 °C)', 'H (45 °C)', 'orig (filled)', 'shifted (open)'}, ...
        'Location', 'northeast', 'Orientation', 'horizontal', 'NumColumns', 5);

    % 하단: z_T
    ax2 = nexttile(tl1);
    hold(ax2, 'on');
    grid(ax2, 'on');
    z_qtys = ["LLI", "LAM_PE", "LAM_NE"];
    z_labels = ["LLI", "LAM_{PE}", "LAM_{NE}"];
    base_x2 = [1, 2, 3, 4.6, 5.6, 6.6, 8.2, 9.2, 10.2];
    x_lbls_2 = repmat(["C", "R", "H"], 1, 3);

    for qIdx = 1:3
        for tIdx = 1:3
            x0 = base_x2((qIdx - 1)*3 + tIdx);
            tName = T_LIST(tIdx);
            c_rgb = COL.(char(tName));
            for vIdx = 1:2
                vName = versions(vIdx);
                row = tbl_D1_res(tbl_D1_res.version == vName & tbl_D1_res.target_mode == z_qtys(qIdx) & tbl_D1_res.T == tName, :);
                if isempty(row), continue; end
                dx = -0.14 + 0.28 * (vIdx - 1);
                est = row.estimate(1);
                neg = est - row.ci_lo(1);
                pos = row.ci_hi(1) - est;
                if vIdx == 1
                    errorbar(ax2, x0 + dx, est, neg, pos, 'o', 'Color', c_rgb, ...
                        'MarkerFaceColor', c_rgb, 'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 5);
                else
                    errorbar(ax2, x0 + dx, est, neg, pos, 's', 'Color', c_rgb, ...
                        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 5);
                end
            end
        end
    end
    xline(ax2, 3.8, '--', 'Color', [0.5 0.5 0.5]);
    xline(ax2, 7.4, '--', 'Color', [0.5 0.5 0.5]);
    xticks(ax2, base_x2);
    xticklabels(ax2, x_lbls_2);
    xlim(ax2, [0.3, 10.9]);
    yl2 = ylim(ax2);
    text(ax2, 2.0, yl2(2) - 0.06*(yl2(2)-yl2(1)), z_labels(1), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    text(ax2, 5.6, yl2(2) - 0.06*(yl2(2)-yl2(1)), z_labels(2), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    text(ax2, 9.2, yl2(2) - 0.06*(yl2(2)-yl2(1)), z_labels(3), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    ylabel(ax2, 'z_T (E2 time exponent)');
    title(ax2, 'D1 Reference-point sensitivity: Mode time exponent z_T (orig vs shifted)');

    exportgraphics(fig1, fullfile('out', 'fig_D1_ref.png'), 'Resolution', 150, 'BackgroundColor', 'white');
    close(fig1);

    % --- Fig D2: 2x2 (열: R0, arc / 행: LOCO RMSE & AIC 막대, 최적 LOCO 예측 vs 실측 산점도) ---
    fig2 = figure('Color', 'w', 'Position', [100, 100, 1100, 760], 'Visible', 'off');
    tl2 = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    t_names_D2 = ["dR0", "dRarc"];
    t_titles_D2 = ["\DeltaR_0", "\DeltaR_{arc}"];

    for colIdx = 1:2
        tName = t_names_D2(colIdx);
        pd = D2_plot.(char(tName));

        % 상단: LOCO RMSE & AIC 막대 그래프
        ax_top = nexttile(tl2, colIdx);
        hold(ax_top, 'on');
        grid(ax_top, 'on');
        yyaxis(ax_top, 'left');
        b1 = bar(ax_top, (1:4) - 0.18, pd.loco_rmse, 0.34, 'FaceColor', [0.20 0.45 0.72]);
        ylabel(ax_top, 'LOCO RMSE (m\Omega)');
        ax_top.YColor = [0.20 0.45 0.72];

        yyaxis(ax_top, 'right');
        b2 = bar(ax_top, (1:4) + 0.18, pd.aic, 0.34, 'FaceColor', [0.85 0.52 0.18]);
        ylabel(ax_top, 'AIC');
        ax_top.YColor = [0.75 0.40 0.10];

        xticks(ax_top, 1:4);
        xticklabels(ax_top, {'S1', 'S2', 'S3', 'S4'});
        xlim(ax_top, [0.4, 4.6]);
        title(ax_top, sprintf('%s: LOCO RMSE & AIC (S1–S4)', t_titles_D2(colIdx)));
        legend(ax_top, [b1, b2], {'LOCO RMSE', 'AIC'}, 'Location', 'best');

        % 하단: 최적 LOCO 모델의 예측 vs 실측 산점도
        ax_bot = nexttile(tl2, colIdx + 2);
        hold(ax_bot, 'on');
        grid(ax_bot, 'on');
        for ti = 1:3
            tn = T_LIST(ti);
            mk = (pd.T == tn);
            scatter(ax_bot, pd.y_meas(mk), pd.y_pred(mk), 28, COL.(char(tn)), 'filled', ...
                'MarkerFaceAlpha', 0.75, 'DisplayName', char(tn));
        end
        all_v = [pd.y_meas(:); pd.y_pred(:)];
        all_v = all_v(isfinite(all_v));
        lo_v = min(0, min(all_v));
        hi_v = max(all_v) * 1.05;
        plot(ax_bot, [lo_v, hi_v], [lo_v, hi_v], '--k', 'LineWidth', 1.1, 'DisplayName', '1:1');
        xlim(ax_bot, [lo_v, hi_v]);
        ylim(ax_bot, [lo_v, hi_v]);
        xlabel(ax_bot, sprintf('Measured %s (m\\Omega)', t_titles_D2(colIdx)));
        ylabel(ax_bot, sprintf('LOCO Predicted %s (m\\Omega)', t_titles_D2(colIdx)));
        title(ax_bot, sprintf('%s: Best LOCO (%s) Pred vs Meas', t_titles_D2(colIdx), pd.best_model));
        legend(ax_bot, 'Location', 'northwest');
    end

    exportgraphics(fig2, fullfile('out', 'fig_D2_models.png'), 'Resolution', 150, 'BackgroundColor', 'white');
    close(fig2);

    % --- Fig D3: 1x2 (L=5, L=10) median pred/meas (log scale) with IQR bars & 0.8-1.25 band ---
    fig3 = figure('Color', 'w', 'Position', [100, 100, 1120, 480], 'Visible', 'off');
    tl3 = tiledlayout(fig3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    q_labels_D3 = {'LLI', 'LAM_{PE}', 'LAM_{NE}', '\DeltaR_0', '\DeltaR_{arc}'};
    offsets_D3 = [-0.21, -0.07, 0.07, 0.21]; % C-orig, C-shifted, R-orig, R-shifted

    for lIdx = 1:2
        lvl = levels_D3(lIdx);
        ax = nexttile(tl3);
        hold(ax, 'on');
        grid(ax, 'on');

        % 0.8 - 1.25 재현 대역 음영
        fill(ax, [0.4, 5.6, 5.6, 0.4], [0.8, 0.8, 1.25, 1.25], [0.85, 0.94, 0.85], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.65);
        yline(ax, 1.0, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);

        y_vals_collect = [0.8; 1.25];
        for qIdx = 1:numel(qtys_D3)
            q = qtys_D3(qIdx);
            combo = 0;
            for ti = 1:2
                tName = T_eval(ti);
                c_rgb = COL.(char(tName));
                for vIdx = 1:2
                    vName = versions(vIdx);
                    combo = combo + 1;
                    row = tbl_D3_sum(tbl_D3_sum.version == vName & tbl_D3_sum.T == tName & ...
                        tbl_D3_sum.quantity == q & tbl_D3_sum.level == lvl, :);
                    if isempty(row) || row.n(1) == 0 || ~isfinite(row.median_ratio(1))
                        continue;
                    end
                    med_v = row.median_ratio(1);
                    neg_v = max(0, med_v - row.q25(1));
                    pos_v = max(0, row.q75(1) - med_v);
                    xpos = qIdx + offsets_D3(combo);

                    if vIdx == 1
                        errorbar(ax, xpos, med_v, neg_v, pos_v, 'o', 'Color', c_rgb, ...
                            'MarkerFaceColor', c_rgb, 'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 4);
                    else
                        errorbar(ax, xpos, med_v, neg_v, pos_v, 's', 'Color', c_rgb, ...
                            'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.4, 'CapSize', 4);
                    end
                    y_vals_collect = [y_vals_collect; row.q25(1); row.q75(1); med_v]; %#ok<AGROW>
                end
            end
        end

        set(ax, 'YScale', 'log');
        xticks(ax, 1:5);
        xticklabels(ax, q_labels_D3);
        xlim(ax, [0.4, 5.6]);
        y_pos = y_vals_collect(isfinite(y_vals_collect) & y_vals_collect > 0);
        ylim(ax, [min(y_pos)*0.75, max(y_pos)*1.35]);
        ylabel(ax, 'Predicted / Measured (log scale)');
        title(ax, sprintf('D3 Transfer (H \\rightarrow C, R) at L = %d%%', lvl));

        if lIdx == 1
            hl = gobjects(5,1);
            hl(1) = plot(ax, nan, nan, 'o', 'Color', COL.C, 'MarkerFaceColor', COL.C, 'LineWidth', 1.4);
            hl(2) = plot(ax, nan, nan, 's', 'Color', COL.C, 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
            hl(3) = plot(ax, nan, nan, 'o', 'Color', COL.R, 'MarkerFaceColor', COL.R, 'LineWidth', 1.4);
            hl(4) = plot(ax, nan, nan, 's', 'Color', COL.R, 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
            hl(5) = fill(ax, nan, nan, [0.85, 0.94, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.65);
            legend(ax, hl, {'C orig', 'C shifted', 'R orig', 'R shifted', '0.8–1.25 band'}, 'Location', 'best');
        end
    end

    exportgraphics(fig3, fullfile('out', 'fig_D3_transfer.png'), 'Resolution', 150, 'BackgroundColor', 'white');
    close(fig3);

    logmsg(logFile, sprintf('[그림] 3개 그림 저장 완료 (%.2f s)', toc(t_sec)));
catch ME
    logmsg(logFile, sprintf('[그림 오류] %s', getReport(ME, 'extended', 'hyperlinks', 'off')));
end

logmsg(logFile, sprintf('=== 전체 완료: 총 소요시간 %.2f s ===', toc(t_total)));

%% ==================== 로컬 함수들 ====================

function logmsg(logFile, msg)
    fprintf('%s\n', msg);
    fid = fopen(logFile, 'a');
    if fid ~= -1
        fprintf(fid, '%s\n', msg);
        fclose(fid);
    end
end

function val = my_prctile(x, p)
    x = x(isfinite(x));
    n = numel(x);
    if n == 0
        val = NaN;
        return;
    end
    if n == 1
        val = x(1);
        return;
    end
    xs = sort(x(:));
    q = 100 * ((1:n)' - 0.5) / n;
    pc = min(max(p, q(1)), q(end));
    val = interp1(q, xs, pc, 'linear');
end

function [p_hat, boot_p, n_cells, n_rows] = fit_M4_and_boot(cell_vec, L_vec, y_vec, B)
    n_rows = numel(y_vec);
    u_cells = unique(cell_vec, 'stable');
    n_cells = numel(u_cells);
    if n_cells == 0 || n_rows == 0
        p_hat = NaN;
        boot_p = nan(B, 1);
        return;
    end

    L_c = cell(n_cells, 1);
    y_c = cell(n_cells, 1);
    for i = 1:n_cells
        mk = (cell_vec == u_cells(i));
        L_c{i} = double(L_vec(mk));
        y_c{i} = double(y_vec(mk));
    end

    % 초기값 탐색용 그리드 사전 계산 (n_cells x n_grid)
    p_grid = linspace(0.2, 4.0, 20);
    rss_grid = zeros(n_cells, numel(p_grid));
    for g = 1:numel(p_grid)
        pg = p_grid(g);
        for i = 1:n_cells
            xi = L_c{i} .^ pg;
            si = sum(xi .* y_c{i}) / sum(xi .^ 2);
            rss_grid(i, g) = sum((y_c{i} - si * xi) .^ 2);
        end
    end

    opts = optimset('Display', 'off', 'TolX', 1e-5, 'TolFun', 1e-6, 'MaxFunEvals', 800, 'MaxIter', 400);

    % 점추정
    w1 = ones(n_cells, 1);
    [~, idx0] = min(w1' * rss_grid);
    p0 = p_grid(idx0);
    p_hat = fminsearch(@(p) m4_weighted_rss(p, L_c, y_c, w1), p0, opts);
    p_hat = min(max(p_hat, 0.2), 4.0);

    % 셀 부트스트랩
    boot_p = nan(B, 1);
    for b = 1:B
        draws = randi(n_cells, n_cells, 1);
        wb = zeros(n_cells, 1);
        for d = 1:n_cells
            wb(draws(d)) = wb(draws(d)) + 1;
        end
        [~, idxb] = min(wb' * rss_grid);
        pb = fminsearch(@(p) m4_weighted_rss(p, L_c, y_c, wb), p_grid(idxb), opts);
        boot_p(b) = min(max(pb, 0.2), 4.0);
    end
end

function rss = m4_weighted_rss(p, L_c, y_c, w)
    if p < 0.2 || p > 4.0 || ~isfinite(p)
        dist = max(0.2 - p, 0) + max(p - 4.0, 0);
        rss = 1e30 * (1 + dist^2);
        return;
    end
    rss = 0;
    for i = 1:numel(w)
        wi = w(i);
        if wi == 0, continue; end
        xi = L_c{i} .^ p;
        si = sum(xi .* y_c{i}) / sum(xi .^ 2);
        r = y_c{i} - si * xi;
        rss = rss + wi * sum(r .^ 2);
    end
end

function [z_hat, boot_z, n_cells, n_rows] = fit_E2_and_boot(cell_vec, days_vec, m_vec, B)
    n_rows = numel(m_vec);
    u_cells = unique(cell_vec, 'stable');
    n_cells = numel(u_cells);
    if n_cells == 0 || n_rows == 0
        z_hat = NaN;
        boot_z = nan(B, 1);
        return;
    end

    Sxx = zeros(n_cells, 1);
    Sxy = zeros(n_cells, 1);
    for i = 1:n_cells
        mk = (cell_vec == u_cells(i));
        xi = log(double(days_vec(mk)));
        yi = log(double(m_vec(mk)));
        dxi = xi - mean(xi);
        dyi = yi - mean(yi);
        Sxx(i) = sum(dxi .^ 2);
        Sxy(i) = sum(dxi .* dyi);
    end

    denom = sum(Sxx);
    if denom > 0
        z_hat = sum(Sxy) / denom;
    else
        z_hat = NaN;
    end

    boot_z = nan(B, 1);
    for b = 1:B
        draws = randi(n_cells, n_cells, 1);
        wb = zeros(n_cells, 1);
        for d = 1:n_cells
            wb(draws(d)) = wb(draws(d)) + 1;
        end
        db = sum(wb .* Sxx);
        if db > 0
            boot_z(b) = sum(wb .* Sxy) / db;
        end
    end
end

function [p_T, s_cell, total_rss] = fit_D2_S0(cell_id, T_idx, L, y, n_cells)
    p_T = zeros(1, 3);
    s_cell = zeros(n_cells, 1);
    total_rss = 0;
    opts = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-7, 'MaxFunEvals', 2000, 'MaxIter', 1000);
    p_grid = linspace(0.1, 4.0, 25);

    for ti = 1:3
        mkT = (T_idx == ti);
        c_list = unique(cell_id(mkT), 'stable');
        Lc = cell(numel(c_list), 1);
        yc = cell(numel(c_list), 1);
        for i = 1:numel(c_list)
            mkc = (cell_id == c_list(i));
            Lc{i} = L(mkc);
            yc{i} = y(mkc);
        end

        best_r = inf;
        p0 = 1.0;
        for g = 1:numel(p_grid)
            rg = s0_sub_rss(p_grid(g), Lc, yc);
            if rg < best_r
                best_r = rg;
                p0 = p_grid(g);
            end
        end
        [p_opt, rss_opt] = fminsearch(@(p) s0_sub_rss(p, Lc, yc), p0, opts);
        p_T(ti) = p_opt;
        total_rss = total_rss + rss_opt;
        for i = 1:numel(c_list)
            xi = Lc{i} .^ p_opt;
            s_cell(c_list(i)) = sum(xi .* yc{i}) / sum(xi .^ 2);
        end
    end
end

function rss = s0_sub_rss(p, Lc, yc)
    if p < 0.05 || p > 6.0 || ~isfinite(p)
        rss = 1e30 * (1 + max(0.05 - p, 0)^2 + max(p - 6.0, 0)^2);
        return;
    end
    rss = 0;
    for i = 1:numel(Lc)
        xi = Lc{i} .^ p;
        si = sum(xi .* yc{i}) / sum(xi .^ 2);
        r = yc{i} - si * xi;
        rss = rss + sum(r .^ 2);
    end
end

function [th_opt, rss_opt] = fit_D2_S1(T_idx, L, y, p_init)
    th0 = fit_D2_S1_fast(T_idx, L, y, p_init);
    opts = optimset('Display', 'off', 'TolX', 1e-7, 'TolFun', 1e-7, 'MaxFunEvals', 20000, 'MaxIter', 10000);
    [th_opt, rss_opt] = fminsearch(@(th) rss_D2_S1_full(th, T_idx, L, y), th0, opts);
end

function th = fit_D2_S1_fast(T_idx, L, y, p_init)
    a_T = zeros(1, 3);
    p_T = zeros(1, 3);
    opts = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-7, 'MaxFunEvals', 1000, 'MaxIter', 500);
    p_grid = linspace(0.1, 4.0, 20);
    for ti = 1:3
        mk = (T_idx == ti);
        Lt = L(mk);
        yt = y(mk);
        cands = unique([p_init(ti), p_grid]);
        best_r = inf;
        p0 = p_init(ti);
        for g = 1:numel(cands)
            rg = pooled_p_rss(cands(g), Lt, yt);
            if rg < best_r
                best_r = rg;
                p0 = cands(g);
            end
        end
        p_opt = fminsearch(@(p) pooled_p_rss(p, Lt, yt), p0, opts);
        x = Lt .^ p_opt;
        s_opt = max(1e-12, sum(x .* yt) / sum(x .^ 2));
        a_T(ti) = log(s_opt);
        p_T(ti) = p_opt;
    end
    th = [a_T, p_T];
end

function rss = pooled_p_rss(p, L, y)
    if p < 0.05 || p > 6.0 || ~isfinite(p)
        rss = 1e30 * (1 + max(0.05 - p, 0)^2 + max(p - 6.0, 0)^2);
        return;
    end
    x = L .^ p;
    s = max(1e-12, sum(x .* y) / sum(x .^ 2));
    r = y - s * x;
    rss = sum(r .^ 2);
end

function rss = rss_D2_S1_full(th, T_idx, L, y)
    a_T = th(1:3);
    p_T = th(4:6);
    if any(abs(a_T) > 30) || any(p_T < 0.05) || any(p_T > 6.0) || any(~isfinite(th))
        rss = 1e30;
        return;
    end
    yhat = exp(a_T(T_idx)') .* (L .^ (p_T(T_idx)'));  % 수정: 전치 연산자 우선순위로 행벡터가 되어 NxN 브로드캐스트 발생, 괄호 추가
    r = y - yhat;
    rss = sum(r .^ 2);
end

function [th_opt, rss_opt] = fit_D2_S2(T_idx, log_chg, L, y, th_S1, s0_pT, s0_sCell, cell_id, n_cells)
    % 다중 시작점 생성 (S1, S0/M4 회귀, alpha 섭동)
    starts4 = zeros(4, 4); % [alpha, p_C, p_R, p_H]
    starts4(1, :) = [0.0, th_S1(4:6)];

    % M4 s_cell 로그 회귀로 초기값 추정
    c_T = zeros(n_cells, 1);
    c_lchg = zeros(n_cells, 1);
    for c = 1:n_cells
        idx = find(cell_id == c, 1, 'first');
        c_T(c) = T_idx(idx);
        c_lchg(c) = log_chg(idx);
    end
    pos_s = s0_sCell > 0;
    if sum(pos_s) >= 5
        Xreg = [double(c_T(pos_s)==1), double(c_T(pos_s)==2), double(c_T(pos_s)==3), c_lchg(pos_s)];
        b_reg = Xreg \ log(s0_sCell(pos_s));
        starts4(2, :) = [b_reg(4), s0_pT];
    else
        starts4(2, :) = [0.0, s0_pT];
    end
    starts4(3, :) = [0.4, th_S1(4:6)];
    starts4(4, :) = [-0.4, th_S1(4:6)];

    w1 = ones(numel(y), 1);
    opts4 = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxFunEvals', 8000, 'MaxIter', 4000);
    opts7 = optimset('Display', 'off', 'TolX', 1e-7, 'TolFun', 1e-7, 'MaxFunEvals', 25000, 'MaxIter', 12000);

    rss_opt = inf;
    th_opt = [th_S1(1:3), 0, th_S1(4:6)];

    for s = 1:size(starts4, 1)
        u_opt = fminsearch(@(u) rss_D2_S2_prof(u, T_idx, log_chg, L, y, w1), starts4(s, :), opts4);
        [~, a_T_0] = rss_D2_S2_prof(u_opt, T_idx, log_chg, L, y, w1);
        th_init = [a_T_0, u_opt(1), u_opt(2:4)];
        [th_cand, rss_cand] = fminsearch(@(th) rss_D2_S2_full(th, T_idx, log_chg, L, y, w1), th_init, opts7);
        if rss_cand < rss_opt
            rss_opt = rss_cand;
            th_opt = th_cand;
        end
    end
end

function th = fit_D2_S2_fast(T_idx, log_chg, L, y, w, th_warm)
    opts = optimset('Display', 'off', 'TolX', 1e-5, 'TolFun', 1e-5, 'MaxFunEvals', 2500, 'MaxIter', 1200);
    u0 = [th_warm(4), th_warm(5:7)];
    u_opt = fminsearch(@(u) rss_D2_S2_prof(u, T_idx, log_chg, L, y, w), u0, opts);
    [~, a_T] = rss_D2_S2_prof(u_opt, T_idx, log_chg, L, y, w);
    th = [a_T, u_opt(1), u_opt(2:4)];
end

function [rss, a_T] = rss_D2_S2_prof(u, T_idx, log_chg, L, y, w)
    alpha = u(1);
    p_T = u(2:4);
    a_T = zeros(1, 3);
    if abs(alpha) > 10 || any(p_T < 0.05) || any(p_T > 6.0) || any(~isfinite(u))
        rss = 1e30;
        return;
    end
    rss = 0;
    for ti = 1:3
        mk = (T_idx == ti) & (w > 0);
        if ~any(mk)
            a_T(ti) = 0;
            continue;
        end
        z = exp(alpha * log_chg(mk)) .* (L(mk) .^ p_T(ti));
        wt = w(mk);
        yt = y(mk);
        s_val = max(1e-12, sum(wt .* z .* yt) / sum(wt .* (z .^ 2)));
        a_T(ti) = log(s_val);
        r = yt - s_val * z;
        rss = rss + sum(wt .* (r .^ 2));
    end
end

function rss = rss_D2_S2_full(th, T_idx, log_chg, L, y, w)
    a_T = th(1:3);
    alpha = th(4);
    p_T = th(5:7);
    if any(abs(a_T) > 30) || abs(alpha) > 10 || any(p_T < 0.05) || any(p_T > 6.0) || any(~isfinite(th))
        rss = 1e30;
        return;
    end
    ls = a_T(T_idx)' + alpha * log_chg;
    if any(abs(ls) > 30)
        rss = 1e30;
        return;
    end
    yhat = exp(ls) .* (L .^ (p_T(T_idx)'));  % 수정: 전치 연산자 우선순위로 행벡터가 되어 NxN 브로드캐스트 발생, 괄호 추가
    r = y - yhat;
    rss = sum(w .* (r .^ 2));
end

function [th_opt, rss_opt] = fit_D2_S34(T_var, log_chg, L, y, th_S2, th_S1, s0_pT, s0_sCell, cell_id, n_cells, th_ref)
    % 스케일 변환 Ts = (T - 20) / 25 로 수치 안정성 확보
    T_mu = 20.0;
    T_sc = 25.0;
    Ts = (T_var - T_mu) / T_sc;

    T_nom = [-5; 25; 45];
    Ts_nom = (T_nom - T_mu) / T_sc;

    % 시작점 1: S2 파라미터 선형 근사
    ba2 = [ones(3,1), Ts_nom] \ th_S2(1:3)';
    bp2 = [ones(3,1), Ts_nom] \ th_S2(5:7)';
    st1 = [ba2(2), th_S2(4), bp2(1), bp2(2)];

    % 시작점 2: S1 선형 근사
    ba1 = [ones(3,1), Ts_nom] \ th_S1(1:3)';
    bp1 = [ones(3,1), Ts_nom] \ th_S1(4:6)';
    st2 = [ba1(2), 0.0, bp1(1), bp1(2)];

    % 시작점 3: M4 s_cell 회귀
    c_Ts = zeros(n_cells, 1);
    c_lchg = zeros(n_cells, 1);
    for c = 1:n_cells
        idx = find(cell_id == c, 1, 'first');
        c_Ts(c) = Ts(idx);
        c_lchg(c) = log_chg(idx);
    end
    pos_s = s0_sCell > 0;
    if sum(pos_s) >= 4
        breg = [ones(sum(pos_s),1), c_Ts(pos_s), c_lchg(pos_s)] \ log(s0_sCell(pos_s));
        bp0 = [ones(3,1), Ts_nom] \ s0_pT';
        st3 = [breg(2), breg(3), bp0(1), bp0(2)];
    else
        st3 = st1;
    end

    starts4 = [st1; st2; st3; st1 + [0.3, 0.2, 0.0, 0.2]; st1 - [0.3, 0.2, 0.0, 0.2]];
    if ~isempty(th_ref)
        % 이전 모델(S3) 해를 스케일 좌표로 변환하여 추가
        beta_s = th_ref(2) * T_sc;
        b1_s = th_ref(5) * T_sc;
        b0_s = th_ref(4) + th_ref(5) * T_mu;
        starts4 = [starts4; beta_s, th_ref(3), b0_s, b1_s];
    end

    w1 = ones(numel(y), 1);
    opts4 = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxFunEvals', 8000, 'MaxIter', 4000);
    opts5 = optimset('Display', 'off', 'TolX', 1e-7, 'TolFun', 1e-7, 'MaxFunEvals', 25000, 'MaxIter', 12000);

    rss_opt = inf;
    th_s_opt = [0, st1];

    for s = 1:size(starts4, 1)
        u_opt = fminsearch(@(u) rss_D2_S34_prof(u, Ts, log_chg, L, y, w1), starts4(s, :), opts4);
        [~, a0_s] = rss_D2_S34_prof(u_opt, Ts, log_chg, L, y, w1);
        th_s_init = [a0_s, u_opt(1), u_opt(2), u_opt(3), u_opt(4)];
        [th_s_cand, rss_cand] = fminsearch(@(th) rss_D2_S34_full(th, Ts, log_chg, L, y, w1), th_s_init, opts5);
        if rss_cand < rss_opt
            rss_opt = rss_cand;
            th_s_opt = th_s_cand;
        end
    end

    % 원래 온도 스케일로 역변환: [a0, beta, alpha, b0, b1]
    beta = th_s_opt(2) / T_sc;
    a0 = th_s_opt(1) - T_mu * beta;
    alpha = th_s_opt(3);
    b1 = th_s_opt(5) / T_sc;
    b0 = th_s_opt(4) - T_mu * b1;
    th_opt = [a0, beta, alpha, b0, b1];
end

function th = fit_D2_S34_fast(T_var, log_chg, L, y, w, th_warm)
    T_mu = 20.0;
    T_sc = 25.0;
    Ts = (T_var - T_mu) / T_sc;

    beta_s0 = th_warm(2) * T_sc;
    b1_s0 = th_warm(5) * T_sc;
    b0_s0 = th_warm(4) + th_warm(5) * T_mu;
    u0 = [beta_s0, th_warm(3), b0_s0, b1_s0];

    opts = optimset('Display', 'off', 'TolX', 1e-5, 'TolFun', 1e-5, 'MaxFunEvals', 2500, 'MaxIter', 1200);
    u_opt = fminsearch(@(u) rss_D2_S34_prof(u, Ts, log_chg, L, y, w), u0, opts);
    [~, a0_s] = rss_D2_S34_prof(u_opt, Ts, log_chg, L, y, w);

    beta = u_opt(1) / T_sc;
    a0 = a0_s - T_mu * beta;
    alpha = u_opt(2);
    b1 = u_opt(4) / T_sc;
    b0 = u_opt(3) - T_mu * b1;
    th = [a0, beta, alpha, b0, b1];
end

function [rss, a0_s] = rss_D2_S34_prof(u, Ts, log_chg, L, y, w)
    beta_s = u(1);
    alpha = u(2);
    b0_s = u(3);
    b1_s = u(4);
    a0_s = 0;

    if abs(beta_s) > 20 || abs(alpha) > 10 || abs(b0_s) > 10 || abs(b1_s) > 10 || any(~isfinite(u))
        rss = 1e30;
        return;
    end
    p = b0_s + b1_s * Ts;
    if any(p < 0.05) || any(p > 6.0)
        rss = 1e30;
        return;
    end
    ls_rel = beta_s * Ts + alpha * log_chg;
    if any(abs(ls_rel) > 25)
        rss = 1e30;
        return;
    end
    z = exp(ls_rel) .* (L .^ p);
    s0 = max(1e-12, sum(w .* z .* y) / sum(w .* (z .^ 2)));
    a0_s = log(s0);
    r = y - s0 * z;
    rss = sum(w .* (r .^ 2));
end

function rss = rss_D2_S34_full(th_s, Ts, log_chg, L, y, w)
    a0_s = th_s(1);
    beta_s = th_s(2);
    alpha = th_s(3);
    b0_s = th_s(4);
    b1_s = th_s(5);

    if abs(a0_s) > 30 || abs(beta_s) > 20 || abs(alpha) > 10 || abs(b0_s) > 10 || abs(b1_s) > 10 || any(~isfinite(th_s))
        rss = 1e30;
        return;
    end
    p = b0_s + b1_s * Ts;
    if any(p < 0.05) || any(p > 6.0)
        rss = 1e30;
        return;
    end
    ls = a0_s + beta_s * Ts + alpha * log_chg;
    if any(abs(ls) > 30)
        rss = 1e30;
        return;
    end
    yhat = exp(ls) .* (L .^ p);
    r = y - yhat;
    rss = sum(w .* (r .^ 2));
end

function yhat = predict_D2(mName, th, T_idx, T_amb, T_eff, log_chg, L)
    switch mName
        case "S1"
            a_T = th(1:3);
            p_T = th(4:6);
            yhat = exp(a_T(T_idx)') .* (L .^ (p_T(T_idx)'));  % 수정: 전치 연산자 우선순위로 행벡터가 되어 NxN 브로드캐스트 발생, 괄호 추가
        case "S2"
            a_T = th(1:3);
            alpha = th(4);
            p_T = th(5:7);
            yhat = exp(a_T(T_idx)' + alpha * log_chg) .* (L .^ (p_T(T_idx)'));  % 수정: 전치 연산자 우선순위로 행벡터가 되어 NxN 브로드캐스트 발생, 괄호 추가
        case "S3"
            ls = th(1) + th(2) * T_amb + th(3) * log_chg;
            p = th(4) + th(5) * T_amb;
            yhat = exp(ls) .* (L .^ p);
        case "S4"
            ls = th(1) + th(2) * T_eff + th(3) * log_chg;
            p = th(4) + th(5) * T_eff;
            yhat = exp(ls) .* (L .^ p);
    end
end

function [s_opt, p_opt] = fit_pooled_power(L, y)
    p_grid = linspace(0.2, 4.0, 25);
    best_r = inf;
    p0 = 1.0;
    for g = 1:numel(p_grid)
        rg = pooled_p_rss_bounded(p_grid(g), L, y);
        if rg < best_r
            best_r = rg;
            p0 = p_grid(g);
        end
    end
    opts = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-7, 'MaxFunEvals', 1000, 'MaxIter', 500);
    p_opt = fminsearch(@(p) pooled_p_rss_bounded(p, L, y), p0, opts);
    p_opt = min(max(p_opt, 0.2), 4.0);
    x = L .^ p_opt;
    s_opt = sum(x .* y) / sum(x .^ 2);
end

function rss = pooled_p_rss_bounded(p, L, y)
    if p < 0.2 || p > 4.0 || ~isfinite(p)
        dist = max(0.2 - p, 0) + max(p - 4.0, 0);
        rss = 1e30 * (1 + dist^2);
        return;
    end
    x = L .^ p;
    s = sum(x .* y) / sum(x .^ 2);
    r = y - s * x;
    rss = sum(r .^ 2);
end