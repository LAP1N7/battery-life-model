%% 0 설정
clearvars;
close all;

tStart = tic;
outDir = 'out';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

logFile = fullfile(outDir, 'log.txt');
if exist(logFile, 'file')
    delete(logFile);
end
diary(logFile);

try
    fprintf('=== mohtat_eqs_main 시작: %s ===\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

    JOINED = fullfile('..', 'mohtat_modes', 'mohtat_modes_joined.csv');
    if ~exist(JOINED, 'file')
        error('입력 파일을 찾을 수 없습니다: %s', JOINED);
    end

    % 온도 순서 및 색상 정의
    T_list = ["C", "R", "H"];
    T_deg  = [-5, 25, 45];
    T_col  = [0.00 0.45 0.70; ...
              0.00 0.62 0.45; ...
              0.84 0.37 0.00];

    modes   = ["LLI", "LAM_PE", "LAM_NE"];
    targets = ["dR0", "dRarc"];

    % 데이터 불러오기 (UTF-8 및 BOM 방어)
    opts = detectImportOptions(JOINED, 'Encoding', 'UTF-8', 'VariableNamingRule', 'preserve');
    vnames0 = opts.VariableNames;
    vnamesClean = regexprep(vnames0, '^\xFEFF', '');
    if ~ismember('cell', vnamesClean) && endsWith(vnamesClean{1}, 'cell')
        vnamesClean{1} = 'cell';
    end
    opts.VariableNames = vnamesClean;

    strCols = {'T', 'dis'};
    for k = 1:numel(strCols)
        if ismember(strCols{k}, opts.VariableNames)
            opts = setvartype(opts, strCols{k}, 'string');
        end
    end
    numCols = {'T_C', 'chg', 'dod', 'ah', 'days', 'efc', 'loss', 'LLI', 'LAM_PE', 'LAM_NE', 'dR0', 'dRarc'};
    for k = 1:numel(numCols)
        if ismember(numCols{k}, opts.VariableNames)
            opts = setvartype(opts, numCols{k}, 'double');
        end
    end

    raw = readtable(JOINED, opts);
    raw.T   = strtrim(string(raw.T));
    raw.dis = strtrim(string(raw.dis));

    % 기본 필터: loss > 0, days > 0
    maskBase = isfinite(raw.loss) & raw.loss > 0 & isfinite(raw.days) & raw.days > 0;
    dfBase = raw(maskBase, :);

    % 저항 방정식용 필터: 기본 필터 + dR0, dRarc 유한값
    maskRes = maskBase & isfinite(raw.dR0) & isfinite(raw.dRarc);
    dfRes = raw(maskRes, :);

    fprintf('입력 파일: %s\n', JOINED);
    fprintf('전체 행: %d | 기본 유효 행(loss>0, days>0): %d | 저항 유효 행: %d\n', ...
        height(raw), height(dfBase), height(dfRes));

    for iT = 1:numel(T_list)
        tKey = T_list(iT);
        subB = dfBase(dfBase.T == tKey, :);
        subR = dfRes(dfRes.T == tKey, :);
        uB = unique(subB.cell);
        uR = unique(subR.cell);
        nGe3 = 0;
        for ic = 1:numel(uR)
            if sum(subR.cell == uR(ic)) >= 3
                nGe3 = nGe3 + 1;
            end
        end
        fprintf('  T=%s (%+3d C): 기본 %3d행(%d셀), 저항 %3d행(%d셀, >=3행 %d셀)\n', ...
            tKey, T_deg(iT), height(subB), numel(uB), height(subR), numel(uR), nGe3);
    end

    %% E1 열화모드 분배 vs 용량손실 (온도별)
    fprintf('\n==================== E1: 모드 분배 vs 용량손실 ====================\n');
    nE1 = numel(T_list) * numel(modes);
    e1_T       = strings(nE1, 1);
    e1_mode    = strings(nE1, 1);
    e1_nrows   = zeros(nE1, 1);
    e1_ncells  = zeros(nE1, 1);
    e1_k_lin   = nan(nE1, 1);
    e1_R2_lin  = nan(nE1, 1);
    e1_RMSE_lin= nan(nE1, 1);
    e1_AIC_lin = nan(nE1, 1);
    e1_k_pow   = nan(nE1, 1);
    e1_q_pow   = nan(nE1, 1);
    e1_R2_pow  = nan(nE1, 1);
    e1_RMSE_pow= nan(nE1, 1);
    e1_AIC_pow = nan(nE1, 1);
    e1_dAIC    = nan(nE1, 1);

    idx = 0;
    for iT = 1:numel(T_list)
        tKey = T_list(iT);
        subT = dfBase(dfBase.T == tKey, :);
        for iM = 1:numel(modes)
            idx = idx + 1;
            mName = modes(iM);
            mVal = subT.(mName);
            valid = isfinite(mVal) & isfinite(subT.loss) & subT.loss > 0;
            L = subT.loss(valid);
            y = mVal(valid);
            cIds = subT.cell(valid);
            n = numel(y);

            e1_T(idx)      = tKey;
            e1_mode(idx)   = mName;
            e1_nrows(idx)  = n;
            e1_ncells(idx) = numel(unique(cIds));

            if n < 2
                continue;
            end

            % (a) 원점 통과 선형: m = k * L
            k_a = (L' * y) / (L' * L);
            yhat_a = k_a * L;
            rss_a = sum((y - yhat_a).^2);
            e1_k_lin(idx)    = k_a;
            e1_R2_lin(idx)   = calc_r2(y, yhat_a);
            e1_RMSE_lin(idx) = sqrt(rss_a / n);
            e1_AIC_lin(idx)  = calc_aic(rss_a, n, 1);

            % (b) 거듭제곱: m = k * L^q (로그 파라미터 fminsearch)
            [k_b, q_b, rss_b] = fit_e1_power(L, y, k_a);
            yhat_b = k_b * (L .^ q_b);
            e1_k_pow(idx)    = k_b;
            e1_q_pow(idx)    = q_b;
            e1_R2_pow(idx)   = calc_r2(y, yhat_b);
            e1_RMSE_pow(idx) = sqrt(rss_b / n);
            e1_AIC_pow(idx)  = calc_aic(rss_b, n, 2);
            e1_dAIC(idx)     = e1_AIC_pow(idx) - e1_AIC_lin(idx);
        end
    end

    T_E1 = table(e1_T, e1_mode, e1_nrows, e1_ncells, ...
        e1_k_lin, e1_R2_lin, e1_RMSE_lin, e1_AIC_lin, ...
        e1_k_pow, e1_q_pow, e1_R2_pow, e1_RMSE_pow, e1_AIC_pow, e1_dAIC, ...
        'VariableNames', {'T','mode','n_rows','n_cells', ...
        'k_lin','R2_lin','RMSE_lin','AIC_lin', ...
        'k_pow','q_pow','R2_pow','RMSE_pow','AIC_pow','dAIC_pow_minus_lin'});

    writetable(T_E1, fullfile(outDir, 'E1_modes_vs_loss.csv'));

    fprintf('%-2s %-7s %4s %3s | %7s %6s %6s %7s | %7s %6s %6s %6s %7s %7s\n', ...
        'T', 'Mode', 'n', 'nc', 'k_lin', 'R2_l', 'RMSE_l', 'AIC_l', ...
        'k_pow', 'q_pow', 'R2_p', 'RMSE_p', 'AIC_p', 'dAIC');
    for i = 1:height(T_E1)
        fprintf('%-2s %-7s %4d %3d | %7.4f %6.3f %6.3f %7.1f | %7.4f %6.3f %6.3f %6.3f %7.1f %+7.1f\n', ...
            T_E1.T(i), T_E1.mode(i), T_E1.n_rows(i), T_E1.n_cells(i), ...
            T_E1.k_lin(i), T_E1.R2_lin(i), T_E1.RMSE_lin(i), T_E1.AIC_lin(i), ...
            T_E1.k_pow(i), T_E1.q_pow(i), T_E1.R2_pow(i), T_E1.RMSE_pow(i), T_E1.AIC_pow(i), T_E1.dAIC_pow_minus_lin(i));
    end

    %% E2 열화모드 시간 형상 (온도별 고정효과)
    fprintf('\n==================== E2: 모드 시간 형상 (days & efc) ====================\n');
    nE2 = numel(T_list) * numel(modes);
    e2_T        = strings(nE2, 1);
    e2_mode     = strings(nE2, 1);
    e2_nrows_d  = zeros(nE2, 1);
    e2_ncells_d = zeros(nE2, 1);
    e2_z_d      = nan(nE2, 1);
    e2_Amed_d   = nan(nE2, 1);
    e2_Amin_d   = nan(nE2, 1);
    e2_Amax_d   = nan(nE2, 1);
    e2_R2_d     = nan(nE2, 1);

    e2_nrows_e  = zeros(nE2, 1);
    e2_ncells_e = zeros(nE2, 1);
    e2_z_e      = nan(nE2, 1);
    e2_Amed_e   = nan(nE2, 1);
    e2_Amin_e   = nan(nE2, 1);
    e2_Amax_e   = nan(nE2, 1);
    e2_R2_e     = nan(nE2, 1);

    idx = 0;
    for iT = 1:numel(T_list)
        tKey = T_list(iT);
        subT = dfBase(dfBase.T == tKey, :);
        for iM = 1:numel(modes)
            idx = idx + 1;
            mName = modes(iM);
            mVal = subT.(mName);
            e2_T(idx)    = tKey;
            e2_mode(idx) = mName;

            % days 기준 (m > 0.05, 셀당 >= 3행)
            [resD] = fit_e2_fe(subT.cell, subT.days, mVal, 0.05, 3);
            e2_nrows_d(idx)  = resD.n_rows;
            e2_ncells_d(idx) = resD.n_cells;
            e2_z_d(idx)      = resD.z;
            e2_Amed_d(idx)   = resD.A_med;
            e2_Amin_d(idx)   = resD.A_min;
            e2_Amax_d(idx)   = resD.A_max;
            e2_R2_d(idx)     = resD.R2_orig;

            % efc 기준 (m > 0.05, efc > 0, 셀당 >= 3행)
            [resE] = fit_e2_fe(subT.cell, subT.efc, mVal, 0.05, 3);
            e2_nrows_e(idx)  = resE.n_rows;
            e2_ncells_e(idx) = resE.n_cells;
            e2_z_e(idx)      = resE.z;
            e2_Amed_e(idx)   = resE.A_med;
            e2_Amin_e(idx)   = resE.A_min;
            e2_Amax_e(idx)   = resE.A_max;
            e2_R2_e(idx)     = resE.R2_orig;
        end
    end

    T_E2 = table(e2_T, e2_mode, ...
        e2_nrows_d, e2_ncells_d, e2_z_d, e2_Amed_d, e2_Amin_d, e2_Amax_d, e2_R2_d, ...
        e2_nrows_e, e2_ncells_e, e2_z_e, e2_Amed_e, e2_Amin_e, e2_Amax_e, e2_R2_e, ...
        'VariableNames', {'T','mode', ...
        'n_rows_days','n_cells_days','z_T','A_med','A_min','A_max','R2_orig_days', ...
        'n_rows_efc','n_cells_efc','z_T_efc','A_efc_med','A_efc_min','A_efc_max','R2_orig_efc'});

    writetable(T_E2, fullfile(outDir, 'E2_mode_time.csv'));

    fprintf('%-2s %-7s | %4s %2s %6s %8s [%7s,%7s] %6s | %6s %8s [%7s,%7s] %6s\n', ...
        'T', 'Mode', 'n', 'nc', 'z_day', 'A_med', 'A_min', 'A_max', 'R2_d', ...
        'z_efc', 'Ae_med', 'Ae_min', 'Ae_max', 'R2_e');
    for i = 1:height(T_E2)
        fprintf('%-2s %-7s | %4d %2d %6.3f %8.4f [%7.4f,%7.4f] %6.3f | %6.3f %8.4f [%7.4f,%7.4f] %6.3f\n', ...
            T_E2.T(i), T_E2.mode(i), T_E2.n_rows_days(i), T_E2.n_cells_days(i), ...
            T_E2.z_T(i), T_E2.A_med(i), T_E2.A_min(i), T_E2.A_max(i), T_E2.R2_orig_days(i), ...
            T_E2.z_T_efc(i), T_E2.A_efc_med(i), T_E2.A_efc_min(i), T_E2.A_efc_max(i), T_E2.R2_orig_efc(i));
    end

    %% E3 저항 증가 vs 용량손실 (온도별, 주 결과)
    fprintf('\n==================== E3: 저항 증가 vs 용량손실 ====================\n');
    nBoot = 1000;
    rng(1);

    nE3 = numel(T_list) * numel(targets) * 4;
    e3_T       = strings(nE3, 1);
    e3_target  = strings(nE3, 1);
    e3_model   = strings(nE3, 1);
    e3_nrows   = zeros(nE3, 1);
    e3_ncells  = zeros(nE3, 1);
    e3_s       = nan(nE3, 1);
    e3_s_lo    = nan(nE3, 1);
    e3_s_hi    = nan(nE3, 1);
    e3_s_med   = nan(nE3, 1);
    e3_s_min   = nan(nE3, 1);
    e3_s_max   = nan(nE3, 1);
    e3_p       = nan(nE3, 1);
    e3_p_lo    = nan(nE3, 1);
    e3_p_hi    = nan(nE3, 1);
    e3_R2      = nan(nE3, 1);
    e3_RMSE    = nan(nE3, 1);
    e3_AIC     = nan(nE3, 1);

    % 셀별 기울기 저장용
    cs_cell = [];
    cs_T    = strings(0, 1);
    cs_tgt  = strings(0, 1);
    cs_sM3  = [];
    cs_sM4  = [];

    % 분산 분해 저장용
    nVD = numel(T_list) * numel(targets);
    vd_T     = strings(nVD, 1);
    vd_tgt   = strings(nVD, 1);
    vd_nP    = zeros(nVD, 1);
    vd_nC    = zeros(nVD, 1);
    vd_R2_M1 = nan(nVD, 1);
    vd_R2_M1c= nan(nVD, 1);
    vd_dM3M1 = nan(nVD, 1);
    vd_dM4M3 = nan(nVD, 1);
    vd_R2_M4 = nan(nVD, 1);

    rowE3 = 0;
    rowVD = 0;

    for iT = 1:numel(T_list)
        tKey = T_list(iT);
        subPool = dfRes(dfRes.T == tKey, :);

        % 셀별 적합용 (>= 3행 셀만 사용)
        [uCellsP, ~, cidP] = unique(subPool.cell);
        cntP = accumarray(cidP, 1);
        keepCells = uCellsP(cntP >= 3);
        subCell = subPool(ismember(subPool.cell, keepCells), :);
        [uCellsC, ~, cidC] = unique(subCell.cell);

        nP = height(subPool);
        nCellsP = numel(uCellsP);
        nC = height(subCell);
        nCellsC = numel(uCellsC);

        L_pool = subPool.loss;
        L_cell = subCell.loss;
        G_cell = make_group_matrix(cidC, nCellsC);
        sameSet = (nCellsP == nCellsC) && isequal(uCellsP, uCellsC);

        for iY = 1:numel(targets)
            yName = targets(iY);
            y_pool = subPool.(yName);
            y_cell = subCell.(yName);

            % --- M1: 풀링 선형 dR = s * L ---
            s_M1 = (L_pool' * y_pool) / (L_pool' * L_pool);
            yhat_M1 = s_M1 * L_pool;
            rss_M1 = sum((y_pool - yhat_M1).^2);
            R2_M1 = calc_r2(y_pool, yhat_M1);
            RMSE_M1 = sqrt(rss_M1 / nP);
            AIC_M1 = calc_aic(rss_M1, nP, 1);

            % --- M2: 풀링 거듭제곱 dR = s * L^p ---
            wUnitP = ones(nCellsP, 1);
            [s_M2, p_M2, rss_M2] = fit_e3_m2(L_pool, y_pool, cidP, wUnitP, [], true);
            yhat_M2 = s_M2 * (L_pool .^ p_M2);
            R2_M2 = calc_r2(y_pool, yhat_M2);
            RMSE_M2 = sqrt(rss_M2 / nP);
            AIC_M2 = calc_aic(rss_M2, nP, 2);

            % --- M3: 셀별 선형 dR = s_cell * L ---
            numM3 = G_cell * (L_cell .* y_cell);
            denM3 = G_cell * (L_cell .* L_cell);
            sVec_M3 = numM3 ./ denM3;
            yhat_M3 = sVec_M3(cidC) .* L_cell;
            rss_M3 = sum((y_cell - yhat_M3).^2);
            R2_M3 = calc_r2(y_cell, yhat_M3);
            RMSE_M3 = sqrt(rss_M3 / nC);
            AIC_M3 = calc_aic(rss_M3, nC, nCellsC);

            % --- M4: 셀별 기울기 + 공통 지수 dR = s_cell * L^p_T ---
            wUnitC = ones(nCellsC, 1);
            [sVec_M4, p_M4, rss_M4] = fit_e3_m4(L_cell, y_cell, cidC, G_cell, wUnitC, [], true);
            yhat_M4 = sVec_M4(cidC) .* (L_cell .^ p_M4);
            R2_M4 = calc_r2(y_cell, yhat_M4);
            RMSE_M4 = sqrt(rss_M4 / nC);
            AIC_M4 = calc_aic(rss_M4, nC, nCellsC + 1);

            % --- 셀 클러스터 부트스트랩 (1000회) ---
            boot_sM1 = nan(nBoot, 1);
            boot_sM2 = nan(nBoot, 1);
            boot_pM2 = nan(nBoot, 1);
            boot_pM4 = nan(nBoot, 1);

            for b = 1:nBoot
                drawP = randi(nCellsP, nCellsP, 1);
                wP = accumarray(drawP, 1, [nCellsP, 1]);
                rwP = wP(cidP);

                % M1 bootstrap
                wL = rwP .* L_pool;
                boot_sM1(b) = (wL' * y_pool) / (wL' * L_pool);

                % M2 bootstrap
                [sb2, pb2] = fit_e3_m2(L_pool, y_pool, cidP, wP, p_M2, false);
                boot_sM2(b) = sb2;
                boot_pM2(b) = pb2;

                % M4 bootstrap
                if sameSet
                    wC = wP;
                else
                    drawC = randi(nCellsC, nCellsC, 1);
                    wC = accumarray(drawC, 1, [nCellsC, 1]);
                end
                [~, pb4] = fit_e3_m4(L_cell, y_cell, cidC, G_cell, wC, p_M4, false);
                boot_pM4(b) = pb4;
            end

            ci_sM1 = my_percentile(boot_sM1, [2.5, 97.5]);
            ci_sM2 = my_percentile(boot_sM2, [2.5, 97.5]);
            ci_pM2 = my_percentile(boot_pM2, [2.5, 97.5]);
            ci_pM4 = my_percentile(boot_pM4, [2.5, 97.5]);

            % 결과 테이블 기록 (M1..M4)
            % M1
            rowE3 = rowE3 + 1;
            e3_T(rowE3) = tKey; e3_target(rowE3) = yName; e3_model(rowE3) = "M1";
            e3_nrows(rowE3) = nP; e3_ncells(rowE3) = nCellsP;
            e3_s(rowE3) = s_M1; e3_s_lo(rowE3) = ci_sM1(1); e3_s_hi(rowE3) = ci_sM1(2);
            e3_p(rowE3) = 1.0;
            e3_R2(rowE3) = R2_M1; e3_RMSE(rowE3) = RMSE_M1; e3_AIC(rowE3) = AIC_M1;

            % M2
            rowE3 = rowE3 + 1;
            e3_T(rowE3) = tKey; e3_target(rowE3) = yName; e3_model(rowE3) = "M2";
            e3_nrows(rowE3) = nP; e3_ncells(rowE3) = nCellsP;
            e3_s(rowE3) = s_M2; e3_s_lo(rowE3) = ci_sM2(1); e3_s_hi(rowE3) = ci_sM2(2);
            e3_p(rowE3) = p_M2; e3_p_lo(rowE3) = ci_pM2(1); e3_p_hi(rowE3) = ci_pM2(2);
            e3_R2(rowE3) = R2_M2; e3_RMSE(rowE3) = RMSE_M2; e3_AIC(rowE3) = AIC_M2;

            % M3
            rowE3 = rowE3 + 1;
            e3_T(rowE3) = tKey; e3_target(rowE3) = yName; e3_model(rowE3) = "M3";
            e3_nrows(rowE3) = nC; e3_ncells(rowE3) = nCellsC;
            e3_s_med(rowE3) = median(sVec_M3); e3_s_min(rowE3) = min(sVec_M3); e3_s_max(rowE3) = max(sVec_M3);
            e3_p(rowE3) = 1.0;
            e3_R2(rowE3) = R2_M3; e3_RMSE(rowE3) = RMSE_M3; e3_AIC(rowE3) = AIC_M3;

            % M4
            rowE3 = rowE3 + 1;
            e3_T(rowE3) = tKey; e3_target(rowE3) = yName; e3_model(rowE3) = "M4";
            e3_nrows(rowE3) = nC; e3_ncells(rowE3) = nCellsC;
            e3_s_med(rowE3) = median(sVec_M4); e3_s_min(rowE3) = min(sVec_M4); e3_s_max(rowE3) = max(sVec_M4);
            e3_p(rowE3) = p_M4; e3_p_lo(rowE3) = ci_pM4(1); e3_p_hi(rowE3) = ci_pM4(2);
            e3_R2(rowE3) = R2_M4; e3_RMSE(rowE3) = RMSE_M4; e3_AIC(rowE3) = AIC_M4;

            % 셀별 기울기 누적
            cs_cell = [cs_cell; uCellsC(:)]; %#ok<AGROW>
            cs_T    = [cs_T; repmat(tKey, nCellsC, 1)]; %#ok<AGROW>
            cs_tgt  = [cs_tgt; repmat(yName, nCellsC, 1)]; %#ok<AGROW>
            cs_sM3  = [cs_sM3; sVec_M3(:)]; %#ok<AGROW>
            cs_sM4  = [cs_sM4; sVec_M4(:)]; %#ok<AGROW>

            % 분산 분해 기록
            rowVD = rowVD + 1;
            s_M1c = (L_cell' * y_cell) / (L_cell' * L_cell);
            R2_M1c = calc_r2(y_cell, s_M1c * L_cell);
            vd_T(rowVD)     = tKey;
            vd_tgt(rowVD)   = yName;
            vd_nP(rowVD)    = nP;
            vd_nC(rowVD)    = nC;
            vd_R2_M1(rowVD) = R2_M1;
            vd_R2_M1c(rowVD)= R2_M1c;
            vd_dM3M1(rowVD) = R2_M3 - R2_M1;
            vd_dM4M3(rowVD) = R2_M4 - R2_M3;
            vd_R2_M4(rowVD) = R2_M4;
        end
    end

    T_E3 = table(e3_T, e3_target, e3_model, e3_nrows, e3_ncells, ...
        e3_s, e3_s_lo, e3_s_hi, e3_s_med, e3_s_min, e3_s_max, ...
        e3_p, e3_p_lo, e3_p_hi, e3_R2, e3_RMSE, e3_AIC, ...
        'VariableNames', {'T','target','model','n_rows','n_cells', ...
        's','s_ci_lo','s_ci_hi','s_cell_med','s_cell_min','s_cell_max', ...
        'p','p_ci_lo','p_ci_hi','R2','RMSE','AIC'});

    T_E3_slopes = table(cs_cell, cs_T, cs_tgt, cs_sM3, cs_sM4, ...
        'VariableNames', {'cell','T','target','s_cell_M3','s_cell_M4'});

    writetable(T_E3, fullfile(outDir, 'E3_resistance_models.csv'));
    writetable(T_E3_slopes, fullfile(outDir, 'E3_cell_slopes.csv'));

    fprintf('%-2s %-5s %-2s %3s %2s | %22s | %21s | %6s %6s %7s\n', ...
        'T', 'Tgt', 'Md', 'n', 'nc', 's (95%CI or med[rng])', 'p (95%CI)', 'R2', 'RMSE', 'AIC');
    for i = 1:height(T_E3)
        if T_E3.model(i) == "M1" || T_E3.model(i) == "M2"
            sStr = sprintf('%.4f [%.4f,%.4f]', T_E3.s(i), T_E3.s_ci_lo(i), T_E3.s_ci_hi(i));
        else
            sStr = sprintf('%.4f [%.4f,%.4f]', T_E3.s_cell_med(i), T_E3.s_cell_min(i), T_E3.s_cell_max(i));
        end
        if isfinite(T_E3.p_ci_lo(i))
            pStr = sprintf('%.3f [%.3f,%.3f]', T_E3.p(i), T_E3.p_ci_lo(i), T_E3.p_ci_hi(i));
        else
            pStr = sprintf('%.3f (fixed)', T_E3.p(i));
        end
        fprintf('%-2s %-5s %-2s %3d %2d | %-22s | %-21s | %6.3f %6.3f %7.1f\n', ...
            T_E3.T(i), T_E3.target(i), T_E3.model(i), T_E3.n_rows(i), T_E3.n_cells(i), ...
            sStr, pStr, T_E3.R2(i), T_E3.RMSE(i), T_E3.AIC(i));
    end

    fprintf('\n--- E3 분산 분해 (Variance Decomposition) ---\n');
    fprintf('%-2s %-5s %4s %4s | %7s %11s %11s %7s\n', ...
        'T', 'Tgt', 'nP', 'nC', 'R2_M1', 'Gain(M3-M1)', 'Gain(M4-M3)', 'R2_M4');
    for i = 1:nVD
        fprintf('%-2s %-5s %4d %4d | %7.4f %+11.4f %+11.4f %7.4f', ...
            vd_T(i), vd_tgt(i), vd_nP(i), vd_nC(i), vd_R2_M1(i), vd_dM3M1(i), vd_dM4M3(i), vd_R2_M4(i));
        if vd_nP(i) ~= vd_nC(i)
            fprintf('  (M1 on >=3-row cells: %.4f)', vd_R2_M1c(i));
        end
        fprintf('\n');
    end

    %% E4 셀별 기울기의 프로토콜 의존성
    fprintf('\n==================== E4: 프로토콜 의존성 (M3 s_cell) ====================\n');
    subR0  = T_E3_slopes(T_E3_slopes.target == "dR0", :);
    subArc = T_E3_slopes(T_E3_slopes.target == "dRarc", :);

    nE4 = height(subR0);
    e4_cell = subR0.cell;
    e4_T    = subR0.T;
    e4_chg  = nan(nE4, 1);
    e4_dis  = strings(nE4, 1);
    e4_dod  = nan(nE4, 1);
    e4_sR0  = subR0.s_cell_M3;
    e4_sArc = nan(nE4, 1);

    for i = 1:nE4
        cId = e4_cell(i);
        tKey = e4_T(i);
        mRow = find(dfRes.cell == cId & dfRes.T == tKey, 1, 'first');
        if ~isempty(mRow)
            e4_chg(i) = dfRes.chg(mRow);
            e4_dis(i) = dfRes.dis(mRow);
            e4_dod(i) = dfRes.dod(mRow);
        end
        aRow = find(subArc.cell == cId & subArc.T == tKey, 1, 'first');
        if ~isempty(aRow)
            e4_sArc(i) = subArc.s_cell_M3(aRow);
        end
    end
    e4_ratio = e4_sArc ./ e4_sR0;

    T_E4 = table(e4_cell, e4_T, e4_chg, e4_dis, e4_dod, e4_sR0, e4_sArc, e4_ratio, ...
        'VariableNames', {'cell','T','chg','dis','dod','s_cell_R0','s_cell_arc','ratio_arc_R0'});
    writetable(T_E4, fullfile(outDir, 'E4_protocol.csv'));

    fprintf('%-5s %-2s %5s %-6s %5s | %9s %10s %10s\n', ...
        'cell', 'T', 'chg', 'dis', 'dod', 's_cell_R0', 's_cell_arc', 's_arc/s_R0');
    for i = 1:height(T_E4)
        fprintf('%-5s %-2s %5.2f %-6s %5.1f | %9.4f %10.4f %10.3f\n', ...
            format_cell_id(T_E4.cell(i)), T_E4.T(i), T_E4.chg(i), T_E4.dis(i), T_E4.dod(i), ...
            T_E4.s_cell_R0(i), T_E4.s_cell_arc(i), T_E4.ratio_arc_R0(i));
    end

    fprintf('\n--- E4 Spearman 상관 (s_cell vs chg) ---\n');
    fprintf('%-4s | %3s %10s | %3s %10s\n', 'Grp', 'n', 'rho(dR0)', 'n', 'rho(dRarc)');
    grpList = [T_list, "ALL"];
    for ig = 1:numel(grpList)
        g = grpList(ig);
        if g == "ALL"
            mask = true(height(T_E4), 1);
        else
            mask = (T_E4.T == g);
        end
        [rhoR0, nR0]   = my_spearman(T_E4.s_cell_R0(mask), T_E4.chg(mask));
        [rhoArc, nArc] = my_spearman(T_E4.s_cell_arc(mask), T_E4.chg(mask));
        fprintf('%-4s | %3d %+10.4f | %3d %+10.4f\n', g, nR0, rhoR0, nArc, rhoArc);
    end

    %% E5 합성 방정식 검증
    fprintf('\n==================== E5: 합성 방정식 dR(t) 검증 ====================\n');
    e5_scope = strings(0, 1);
    e5_T     = strings(0, 1);
    e5_tgt   = strings(0, 1);
    e5_cell  = strings(0, 1);
    e5_nrows = [];
    e5_A     = [];
    e5_z     = [];
    e5_kLLI  = [];
    e5_sM2   = [];
    e5_pM2   = [];
    e5_Ccomp = [];
    e5_exp   = [];
    e5_medR  = [];
    e5_minR  = [];
    e5_maxR  = [];
    e5_R2    = [];
    e5_RMSE  = [];

    fprintf('%-2s %-5s | %8s %6s %6s %7s %6s | dR(t) = %9s * t^%-6s | %4s %2s | %8s [%6s,%6s] %7s\n', ...
        'T', 'Tgt', 'A_LLI', 'z_LLI', 'k_LLI', 's_M2', 'p_M2', 'C_comp', 'q_comp', ...
        'n', 'nc', 'med(p/m)', 'cMin', 'cMax', 'R2_T');

    for iT = 1:numel(T_list)
        tKey = T_list(iT);

        % E2 LLI 시간법칙 파라미터
        rE2 = find(T_E2.T == tKey & T_E2.mode == "LLI", 1, 'first');
        A_lli = T_E2.A_med(rE2);
        z_lli = T_E2.z_T(rE2);

        % E1(a) LLI 분배 기울기
        rE1 = find(T_E1.T == tKey & T_E1.mode == "LLI", 1, 'first');
        k_lli = T_E1.k_lin(rE1);

        % >= 3행 저항 셀 선택
        subT = dfRes(dfRes.T == tKey, :);
        [uC, ~, cid] = unique(subT.cell);
        cnt = accumarray(cid, 1);
        keepC = uC(cnt >= 3);
        subUse = subT(ismember(subT.cell, keepC), :);
        uUse = unique(subUse.cell);

        for iY = 1:numel(targets)
            yName = targets(iY);
            rE3 = find(T_E3.T == tKey & T_E3.target == yName & T_E3.model == "M2", 1, 'first');
            s_m2 = T_E3.s(rE3);
            p_m2 = T_E3.p(rE3);

            C_comp = s_m2 * ((A_lli / k_lli) ^ p_m2);
            q_comp = z_lli * p_m2;

            predAll = [];
            measAll = [];
            cellMedRatios = nan(numel(uUse), 1);

            for ic = 1:numel(uUse)
                cId = uUse(ic);
                sc = subUse(subUse.cell == cId, :);
                tDays = sc.days;
                meas = sc.(yName);
                pred = s_m2 * (((A_lli * (tDays .^ z_lli)) / k_lli) .^ p_m2);

                pos = isfinite(meas) & meas > 0 & isfinite(pred);
                if any(pos)
                    rat = pred(pos) ./ meas(pos);
                    cMedR = median(rat);
                    cMinR = min(rat);
                    cMaxR = max(rat);
                else
                    cMedR = NaN; cMinR = NaN; cMaxR = NaN;
                end
                cR2 = calc_r2(meas, pred);
                cRMSE = sqrt(mean((meas - pred).^2));
                cellMedRatios(ic) = cMedR;

                predAll = [predAll; pred]; %#ok<AGROW>
                measAll = [measAll; meas]; %#ok<AGROW>

                e5_scope(end+1,1) = "per_cell"; %#ok<AGROW>
                e5_T(end+1,1)     = tKey; %#ok<AGROW>
                e5_tgt(end+1,1)   = yName; %#ok<AGROW>
                e5_cell(end+1,1)  = string(format_cell_id(cId)); %#ok<AGROW>
                e5_nrows(end+1,1) = height(sc); %#ok<AGROW>
                e5_A(end+1,1)     = A_lli; %#ok<AGROW>
                e5_z(end+1,1)     = z_lli; %#ok<AGROW>
                e5_kLLI(end+1,1)  = k_lli; %#ok<AGROW>
                e5_sM2(end+1,1)   = s_m2; %#ok<AGROW>
                e5_pM2(end+1,1)   = p_m2; %#ok<AGROW>
                e5_Ccomp(end+1,1) = C_comp; %#ok<AGROW>
                e5_exp(end+1,1)   = q_comp; %#ok<AGROW>
                e5_medR(end+1,1)  = cMedR; %#ok<AGROW>
                e5_minR(end+1,1)  = cMinR; %#ok<AGROW>
                e5_maxR(end+1,1)  = cMaxR; %#ok<AGROW>
                e5_R2(end+1,1)    = cR2; %#ok<AGROW>
                e5_RMSE(end+1,1)  = cRMSE; %#ok<AGROW>
            end

            posAll = isfinite(measAll) & measAll > 0 & isfinite(predAll);
            if any(posAll)
                tMedR = median(predAll(posAll) ./ measAll(posAll));
            else
                tMedR = NaN;
            end
            validCR = cellMedRatios(isfinite(cellMedRatios));
            if ~isempty(validCR)
                cRmin = min(validCR);
                cRmax = max(validCR);
            else
                cRmin = NaN; cRmax = NaN;
            end
            tR2   = calc_r2(measAll, predAll);
            tRMSE = sqrt(mean((measAll - predAll).^2));

            e5_scope(end+1,1) = "per_T"; %#ok<AGROW>
            e5_T(end+1,1)     = tKey; %#ok<AGROW>
            e5_tgt(end+1,1)   = yName; %#ok<AGROW>
            e5_cell(end+1,1)  = "ALL"; %#ok<AGROW>
            e5_nrows(end+1,1) = numel(measAll); %#ok<AGROW>
            e5_A(end+1,1)     = A_lli; %#ok<AGROW>
            e5_z(end+1,1)     = z_lli; %#ok<AGROW>
            e5_kLLI(end+1,1)  = k_lli; %#ok<AGROW>
            e5_sM2(end+1,1)   = s_m2; %#ok<AGROW>
            e5_pM2(end+1,1)   = p_m2; %#ok<AGROW>
            e5_Ccomp(end+1,1) = C_comp; %#ok<AGROW>
            e5_exp(end+1,1)   = q_comp; %#ok<AGROW>
            e5_medR(end+1,1)  = tMedR; %#ok<AGROW>
            e5_minR(end+1,1)  = cRmin; %#ok<AGROW>
            e5_maxR(end+1,1)  = cRmax; %#ok<AGROW>
            e5_R2(end+1,1)    = tR2; %#ok<AGROW>
            e5_RMSE(end+1,1)  = tRMSE; %#ok<AGROW>

            fprintf('%-2s %-5s | %8.4f %6.3f %6.3f %7.4f %6.3f | dR(t) = %9.4g * t^%-6.3f | %4d %2d | %8.3f [%6.2f,%6.2f] %7.3f\n', ...
                tKey, yName, A_lli, z_lli, k_lli, s_m2, p_m2, C_comp, q_comp, ...
                numel(measAll), numel(uUse), tMedR, cRmin, cRmax, tR2);
        end
    end

    T_E5 = table(e5_scope, e5_T, e5_tgt, e5_cell, e5_nrows, ...
        e5_A, e5_z, e5_kLLI, e5_sM2, e5_pM2, e5_Ccomp, e5_exp, ...
        e5_medR, e5_minR, e5_maxR, e5_R2, e5_RMSE, ...
        'VariableNames', {'scope','T','target','cell','n_rows', ...
        'A_LLI_med','z_LLI','k_LLI','s_M2','p_M2','C_comp','q_comp', ...
        'med_ratio_pred_meas','ratio_min','ratio_max','R2','RMSE'});
    writetable(T_E5, fullfile(outDir, 'E5_composite.csv'));

    %% 그림
    fprintf('\n==================== 그림 생성 ====================\n');

    % --- 그림 1: E3 저항 모델 적합 (2x3) ---
    fig1 = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1260 720]);
    tl1 = tiledlayout(fig1, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    yLabelsE3 = ["\DeltaR_0 (m\Omega)", "\DeltaR_{arc} (m\Omega)"];

    for iY = 1:numel(targets)
        yName = targets(iY);
        for iT = 1:numel(T_list)
            tKey = T_list(iT);
            ax = nexttile(tl1);
            hold(ax, 'on');
            grid(ax, 'on');
            box(ax, 'on');

            subP = dfRes(dfRes.T == tKey, :);
            uC = unique(subP.cell);
            cMap = lines(max(numel(uC), 7));

            for ic = 1:numel(uC)
                sc = subP(subP.cell == uC(ic), :);
                [Ls, ord] = sort(sc.loss);
                ys = sc.(yName)(ord);
                plot(ax, Ls, ys, '-o', 'Color', cMap(ic, :), ...
                    'LineWidth', 0.85, 'MarkerSize', 4, ...
                    'MarkerFaceColor', cMap(ic, :), 'HandleVisibility', 'off');
            end

            rM1 = find(T_E3.T == tKey & T_E3.target == yName & T_E3.model == "M1", 1, 'first');
            rM2 = find(T_E3.T == tKey & T_E3.target == yName & T_E3.model == "M2", 1, 'first');
            rM4 = find(T_E3.T == tKey & T_E3.target == yName & T_E3.model == "M4", 1, 'first');

            Lmax = max(subP.loss);
            Lgrid = linspace(0, Lmax, 200)';
            hM1 = plot(ax, Lgrid, T_E3.s(rM1) * Lgrid, '--k', 'LineWidth', 1.4, 'DisplayName', 'M1 pooled linear');
            hM2 = plot(ax, Lgrid, T_E3.s(rM2) * (Lgrid .^ T_E3.p(rM2)), '-k', 'LineWidth', 1.9, 'DisplayName', 'M2 pooled power');

            title(ax, {sprintf('%s (%d^{\\circ}C) | s=%.3g, p=%.2f', tKey, T_deg(iT), T_E3.s(rM2), T_E3.p(rM2)), ...
                       sprintf('R^2 M1/M2/M4 = %.2f / %.2f / %.2f', T_E3.R2(rM1), T_E3.R2(rM2), T_E3.R2(rM4))}, ...
                  'FontSize', 10);

            if iY == 2
                xlabel(ax, 'Capacity loss L (%)', 'FontSize', 10);
            end
            if iT == 1
                ylabel(ax, yLabelsE3(iY), 'FontSize', 10);
            end
            if iY == 1 && iT == 1
                legend(ax, [hM1, hM2], 'Location', 'northwest', 'FontSize', 8.5);
            end

            xlim(ax, [0, Lmax * 1.05]);
            ymin = min(0, min(subP.(yName)) * 1.05);
            ymax = max(subP.(yName)) * 1.08;
            if isfinite(ymin) && isfinite(ymax) && ymax > ymin
                ylim(ax, [ymin, ymax]);
            end
        end
    end

    figE3Path = fullfile(outDir, 'fig_E3_fits.png');
    exportgraphics(fig1, figE3Path, 'Resolution', 150, 'BackgroundColor', 'white');
    close(fig1);
    fprintf('저장 완료: %s\n', figE3Path);

    % --- 그림 2: E1 열화모드 vs 용량손실 (1x3) ---
    fig2 = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1260 420]);
    tl2 = tiledlayout(fig2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    modeLabels = ["LLI (%)", "LAM_{PE} (%)", "LAM_{NE} (%)"];

    for iM = 1:numel(modes)
        mName = modes(iM);
        ax = nexttile(tl2);
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax, 'on');

        hLines = gobjects(numel(T_list), 1);
        legTxt = strings(numel(T_list), 1);
        LmaxAll = 0;

        for iT = 1:numel(T_list)
            tKey = T_list(iT);
            subT = dfBase(dfBase.T == tKey, :);
            mVal = subT.(mName);
            valid = isfinite(mVal) & isfinite(subT.loss) & subT.loss > 0;
            L = subT.loss(valid);
            y = mVal(valid);
            if isempty(L)
                continue;
            end
            LmaxT = max(L);
            LmaxAll = max(LmaxAll, LmaxT);

            scatter(ax, L, y, 22, T_col(iT, :), 'filled', ...
                'MarkerFaceAlpha', 0.65, 'HandleVisibility', 'off');

            rE1 = find(T_E1.T == tKey & T_E1.mode == mName, 1, 'first');
            k_lin = T_E1.k_lin(rE1);
            R2_lin = T_E1.R2_lin(rE1);

            Lgrid = linspace(0, LmaxT, 100)';
            hLines(iT) = plot(ax, Lgrid, k_lin * Lgrid, '-', 'Color', T_col(iT, :), 'LineWidth', 2.0);
            legTxt(iT) = sprintf('%s (%d^{\\circ}C): k=%.2f, R^2=%.2f', tKey, T_deg(iT), k_lin, R2_lin);
        end

        xlabel(ax, 'Capacity loss L (%)', 'FontSize', 10);
        ylabel(ax, modeLabels(iM), 'FontSize', 10);
        title(ax, strrep(mName, '_', '\_'), 'FontSize', 11);
        if LmaxAll > 0
            xlim(ax, [0, LmaxAll * 1.05]);
        end
        validH = isgraphics(hLines);
        if any(validH)
            legend(ax, hLines(validH), legTxt(validH), 'Location', 'northwest', 'FontSize', 8.5);
        end
    end

    figE1Path = fullfile(outDir, 'fig_E1_modes.png');
    exportgraphics(fig2, figE1Path, 'Resolution', 150, 'BackgroundColor', 'white');
    close(fig2);
    fprintf('저장 완료: %s\n', figE1Path);

    fprintf('\n=== 전체 완료 (소요 시간: %.2f초) ===\n', toc(tStart));
    diary off;

catch ME
    fprintf(2, '\n[오류 발생] %s\n', ME.message);
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    diary off;
    rethrow(ME);
end

%% 로컬 함수

function r2 = calc_r2(y, yhat)
    y = y(:);
    yhat = yhat(:);
    ssTot = sum((y - mean(y)).^2);
    if ssTot <= 0 || ~isfinite(ssTot)
        r2 = NaN;
        return;
    end
    ssRes = sum((y - yhat).^2);
    r2 = 1 - ssRes / ssTot;
end

function aic = calc_aic(rss, n, npar)
    if n <= 0 || ~isfinite(rss)
        aic = NaN;
        return;
    end
    aic = n * log(max(rss / n, realmin)) + 2 * npar;
end

function [k_best, q_best, rss_best] = fit_e1_power(L, y, k_lin)
    opt = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10, ...
                   'MaxIter', 3000, 'MaxFunEvals', 5000);
    k0 = max(k_lin, 1e-6);
    qStarts = [1.0, 0.7, 1.3];
    rss_best = inf;
    k_best = k0;
    q_best = 1.0;

    for iq = 1:numel(qStarts)
        th0 = [log(k0), log(qStarts(iq))];
        [thHat, fval] = fminsearch(@(th) e1_pow_rss(th, L, y), th0, opt);
        if fval < rss_best
            rss_best = fval;
            k_best = exp(thHat(1));
            q_best = exp(thHat(2));
        end
    end
end

function rss = e1_pow_rss(th, L, y)
    if any(~isfinite(th)) || th(1) < -30 || th(1) > 30 || th(2) < -5 || th(2) > 3
        rss = 1e30;
        return;
    end
    k = exp(th(1));
    q = exp(th(2));
    yhat = k * (L .^ q);
    if any(~isfinite(yhat))
        rss = 1e30;
        return;
    end
    rss = sum((y - yhat).^2);
end

function res = fit_e2_fe(cellCol, xCol, mCol, mThresh, minRowsPerCell)
    res = struct('n_rows', 0, 'n_cells', 0, 'z', NaN, ...
                 'A_med', NaN, 'A_min', NaN, 'A_max', NaN, 'R2_orig', NaN);

    valid = isfinite(mCol) & mCol > mThresh & isfinite(xCol) & xCol > 0;
    cSub = cellCol(valid);
    xSub = xCol(valid);
    mSub = mCol(valid);

    if isempty(mSub)
        return;
    end

    [uC, ~, cid0] = unique(cSub);
    cnt = accumarray(cid0, 1);
    keepC = uC(cnt >= minRowsPerCell);
    if isempty(keepC)
        return;
    end

    maskKeep = ismember(cSub, keepC);
    cUse = cSub(maskKeep);
    xUse = xSub(maskKeep);
    mUse = mSub(maskKeep);
    [~, ~, cid] = unique(cUse);
    nCells = numel(keepC);
    nRows = numel(mUse);

    lx = log(xUse);
    ly = log(mUse);

    % 고정효과(셀별 절편) OLS: 셀 내 평균 제거
    G = make_group_matrix(cid, nCells);
    nPerCell = sum(G, 2);
    meanLx = (G * lx) ./ nPerCell;
    meanLy = (G * ly) ./ nPerCell;

    lxDm = lx - meanLx(cid);
    lyDm = ly - meanLy(cid);
    den = lxDm' * lxDm;
    if den <= 0
        return;
    end

    z = (lxDm' * lyDm) / den;
    logA = meanLy - z * meanLx;
    A_cell = exp(logA);

    mHat = A_cell(cid) .* (xUse .^ z);

    res.n_rows  = nRows;
    res.n_cells = nCells;
    res.z       = z;
    res.A_med   = median(A_cell);
    res.A_min   = min(A_cell);
    res.A_max   = max(A_cell);
    res.R2_orig = calc_r2(mUse, mHat);
end

function G = make_group_matrix(cid, nGroups)
    n = numel(cid);
    G = zeros(nGroups, n);
    linIdx = cid(:) + (0:n-1)' * nGroups;
    G(linIdx) = 1;
end

function [s_hat, p_hat, rss_hat] = fit_e3_m2(L, y, cid, wCell, pInit, tightTol)
    rw = wCell(cid);
    if tightTol
        opt = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-12, ...
                       'MaxIter', 2000, 'MaxFunEvals', 3000);
        pGrid = [0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0];
    else
        opt = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-8, ...
                       'MaxIter', 500, 'MaxFunEvals', 800);
        pGrid = [pInit, 0.5, 1.0, 1.5, 2.0, 3.0];
    end

    bestP0 = pGrid(1);
    bestR0 = inf;
    for k = 1:numel(pGrid)
        rTmp = e3_m2_rss(pGrid(k), L, y, rw);
        if rTmp < bestR0
            bestR0 = rTmp;
            bestP0 = pGrid(k);
        end
    end

    [p_hat, rss_hat] = fminsearch(@(p) e3_m2_rss(p, L, y, rw), bestP0, opt);
    u = L .^ p_hat;
    wu = rw .* u;
    s_hat = (wu' * y) / (wu' * u);
end

function rss = e3_m2_rss(p, L, y, rw)
    if ~isfinite(p) || p <= 0.01 || p >= 15
        rss = 1e30;
        return;
    end
    u = L .^ p;
    wu = rw .* u;
    den = wu' * u;
    if ~isfinite(den) || den <= 0
        rss = 1e30;
        return;
    end
    s = (wu' * y) / den;
    res = y - s * u;
    rss = (rw .* res)' * res;
    if ~isfinite(rss)
        rss = 1e30;
    end
end

function [sVec, p_hat, rss_hat] = fit_e3_m4(L, y, cid, G, wCell, pInit, tightTol)
    if tightTol
        opt = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-12, ...
                       'MaxIter', 2000, 'MaxFunEvals', 3000);
        pGrid = [0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0];
    else
        opt = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-8, ...
                       'MaxIter', 500, 'MaxFunEvals', 800);
        pGrid = [pInit, 0.5, 1.0, 1.5, 2.0, 3.0];
    end

    bestP0 = pGrid(1);
    bestR0 = inf;
    for k = 1:numel(pGrid)
        rTmp = e3_m4_rss(pGrid(k), L, y, cid, G, wCell);
        if rTmp < bestR0
            bestR0 = rTmp;
            bestP0 = pGrid(k);
        end
    end

    [p_hat, rss_hat] = fminsearch(@(p) e3_m4_rss(p, L, y, cid, G, wCell), bestP0, opt);
    u = L .^ p_hat;
    sVec = (G * (u .* y)) ./ (G * (u .* u));
end

function rss = e3_m4_rss(p, L, y, cid, G, wCell)
    if ~isfinite(p) || p <= 0.01 || p >= 15
        rss = 1e30;
        return;
    end
    u = L .^ p;
    den = G * (u .* u);
    if any(~isfinite(den)) || any(den <= 0)
        rss = 1e30;
        return;
    end
    sVec = (G * (u .* y)) ./ den;
    res = y - sVec(cid) .* u;
    r2Cell = G * (res .^ 2);
    rss = wCell' * r2Cell;
    if ~isfinite(rss)
        rss = 1e30;
    end
end

function q = my_percentile(x, pct)
    x = x(isfinite(x));
    n = numel(x);
    q = nan(size(pct));
    if n == 0
        return;
    end
    xs = sort(x(:));
    if n == 1
        q(:) = xs(1);
        return;
    end
    for k = 1:numel(pct)
        pos = (pct(k) / 100) * n + 0.5;
        if pos <= 1
            q(k) = xs(1);
        elseif pos >= n
            q(k) = xs(n);
        else
            lo = floor(pos);
            hi = ceil(pos);
            frac = pos - lo;
            q(k) = xs(lo) + frac * (xs(hi) - xs(lo));
        end
    end
end

function r = my_tiedrank(v)
    v = v(:);
    n = numel(v);
    [vs, ord] = sort(v);
    rSorted = zeros(n, 1);
    i = 1;
    while i <= n
        j = i;
        while j < n && vs(j + 1) == vs(i)
            j = j + 1;
        end
        rSorted(i:j) = 0.5 * (i + j);
        i = j + 1;
    end
    r = zeros(n, 1);
    r(ord) = rSorted;
end

function r = my_pearson(x, y)
    x = x(:);
    y = y(:);
    dx = x - mean(x);
    dy = y - mean(y);
    den = sqrt((dx' * dx) * (dy' * dy));
    if den <= 0 || ~isfinite(den)
        r = NaN;
    else
        r = (dx' * dy) / den;
    end
end

function [rho, n] = my_spearman(x, y)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    n = numel(x);
    if n < 2
        rho = NaN;
        return;
    end
    rx = my_tiedrank(x);
    ry = my_tiedrank(y);
    rho = my_pearson(rx, ry);
end

function s = format_cell_id(c)
    if isnumeric(c)
        s = sprintf('%g', c);
    else
        s = char(string(c));
    end
end