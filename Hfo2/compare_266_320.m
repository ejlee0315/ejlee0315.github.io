function cmb = compare_266_320(res266, res320, cfg)
% compare_266_320  같은 (h, P) 그리드에서 두 파장 동시 성능 분석
%   ★ 기본 mode = 'joint': 두 파장이 같은 물리적 구조를 공유하는 가정으로
%     각 (h, P) 에서 phase bin 마다 (T266·T320) 최대 구조 1개를 골라
%     joint coverage / joint mean T 를 계산.
%
%   >> compare_266_320                       % 가장 간단 (auto-load + joint)
%   >> cmb = compare_266_320(res266, res320);
%   >> cmb = compare_266_320('search_results_266_320.mat');
%   >> cmb = compare_266_320(res266, res320, struct('mode','independent'))  % 옛날 방식
%
%   cfg 옵션:
%     .mode       : 'joint'(기본) | 'independent'
%                     joint       = 같은 구조를 양 파장 동시에 사용
%                     independent = 각 파장이 자기 best 구조를 따로 선택
%     .primary    : 'lam266'(기본) | 'lam320'  — phase bin 기준 파장 (joint only)
%     .T_min      : 양 파장 모두 이 이상 (기본 0.5, joint only)
%     .N_bin      : phase bin 수 (기본 64, joint only)
%     .joint_obj  : 'T_prod'(기본) | 'T_min' | 'T_sum'  (bin 내 선택 기준)
%     .save_joint_txt : true(기본) — best (h,P) 의 joint 라이브러리 txt 저장
%
%     .score_mode : 'worst'(기본) | 'product' | 'mean' | 'weighted'
%     .w_cov      : coverage 가중치 (기본 2)
%     .w_T        : mean_T 가중치 (기본 1)
%     .plot       : 요약 플롯 (기본 true)
%     .topN       : 상위 N개 출력 (기본 10)

    %% 입력 파싱 -----------------------------------------------------------
    if nargin == 0
        default_mat = 'search_results_266_320.mat';
        if ~exist(default_mat, 'file')
            error('%s 가 현재 폴더에 없습니다. 먼저 run_optimize_266_320 실행.', default_mat);
        end
        fprintf('[compare_266_320] auto-loading %s ...\n', default_mat);
        S = load(default_mat); res266 = S.res266; res320 = S.res320;
        cfg = struct();
    elseif nargin == 1 && (ischar(res266) || isstring(res266))
        S = load(res266); res266 = S.res266; res320 = S.res320;
        cfg = struct();
    elseif nargin < 3
        cfg = struct();
    end
    cfg = set_defaults_cmp(cfg);

    %% (h,P) 그리드 정렬 --------------------------------------------------
    t266 = res266.table; t320 = res320.table;
    if height(t266) ~= height(t320)
        error('266/320 (h,P) 그리드 크기 다름 (%d vs %d)', height(t266), height(t320));
    end
    key266 = t266.h_m*1e12 + t266.P_m*1e6;
    key320 = t320.h_m*1e12 + t320.P_m*1e6;
    [~, i266] = sort(key266); [~, i320] = sort(key320);
    t266 = t266(i266,:); t320 = t320(i320,:);
    if ~all(abs(t266.h_m - t320.h_m) < 1e-15) || ...
       ~all(abs(t266.P_m - t320.P_m) < 1e-15)
        error('(h,P) 그리드 불일치. run_optimize_266_320 다시 실행.');
    end
    files266 = res266.all(i266, :);
    files320 = res320.all(i320, :);

    %% mode 분기 ----------------------------------------------------------
    h_m = t266.h_m;  P_m = t266.P_m;
    if strcmpi(cfg.mode, 'joint')
        fprintf('[compare_266_320] mode = JOINT (양 파장 같은 구조 공유)\n');
        [c266, T266vec, c320, T320vec] = joint_analysis( ...
            files266(:,1), files320(:,1), cfg);
    else
        fprintf('[compare_266_320] mode = INDEPENDENT (파장별 따로 최적)\n');
        c266 = t266.coverage; T266vec = t266.meanT;
        c320 = t320.coverage; T320vec = t320.meanT;
    end

    %% 결합 score ---------------------------------------------------------
    switch cfg.score_mode
        case 'worst',    score = min(c266,c320).^cfg.w_cov .* min(T266vec,T320vec).^cfg.w_T;
        case 'product',  score = (c266.*c320).^(cfg.w_cov/2) .* (T266vec.*T320vec).^(cfg.w_T/2);
        case 'mean',     score = ((c266+c320)/2).^cfg.w_cov .* ((T266vec+T320vec)/2).^cfg.w_T;
        case 'weighted', score = 0.5*(c266.^cfg.w_cov.*T266vec.^cfg.w_T) + ...
                                 0.5*(c320.^cfg.w_cov.*T320vec.^cfg.w_T);
        otherwise, error('Unknown score_mode: %s', cfg.score_mode);
    end

    tbl = table(h_m, P_m, c266, T266vec, c320, T320vec, score, ...
        'VariableNames', {'h_m','P_m','cov266','T266','cov320','T320','score'});

    %% 정렬 + 출력 -------------------------------------------------------
    [~, order] = sort(score, 'descend');
    tbl = tbl(order,:);
    files266 = files266(order,:); files320 = files320(order,:);

    topN = min(cfg.topN, height(tbl));
    fprintf('\n========== %s combined best (score=%s) ==========\n', ...
        upper(cfg.mode), cfg.score_mode);
    fprintf(' rank |  h   |  P   | cov266  T266  | cov320  T320  | score\n');
    fprintf(' -----+------+------+---------------+---------------+-------\n');
    for k = 1:topN
        fprintf(' %3d  | %4d | %4d | %5.1f%% %5.3f | %5.1f%% %5.3f | %.4f\n', ...
            k, round(tbl.h_m(k)*1e9), round(tbl.P_m(k)*1e9), ...
            tbl.cov266(k)*100, tbl.T266(k), tbl.cov320(k)*100, tbl.T320(k), tbl.score(k));
    end
    fprintf('=====================================================================\n');

    best = struct( ...
        'h',tbl.h_m(1), 'P',tbl.P_m(1), ...
        'cov266',tbl.cov266(1), 'T266',tbl.T266(1), ...
        'cov320',tbl.cov320(1), 'T320',tbl.T320(1), ...
        'score',tbl.score(1), ...
        'mat_266',files266{1,1}, 'txt_266',files266{1,2}, ...
        'mat_320',files320{1,1}, 'txt_320',files320{1,2}, ...
        'joint_txt','');

    fprintf('\n>> BEST (%s) : h=%d nm  P=%d nm\n', ...
        cfg.mode, round(best.h*1e9), round(best.P*1e9));
    fprintf('   266nm : cov=%.1f%%  meanT=%.3f\n', best.cov266*100, best.T266);
    fprintf('   320nm : cov=%.1f%%  meanT=%.3f\n', best.cov320*100, best.T320);

    %% joint best 의 라이브러리 txt 저장 ------------------------------------
    if strcmpi(cfg.mode,'joint') && cfg.save_joint_txt
        best.joint_txt = save_joint_txt(best.mat_266, best.mat_320, cfg);
        fprintf('   joint lib : %s\n', best.joint_txt);
    end

    %% 단독 best 참고 ----------------------------------------------------
    fprintf('\n(참고) 단독 best (independent)\n');
    fprintf('   266 only : h=%d P=%d  cov=%.1f%% T=%.3f\n', ...
        round(res266.best.h*1e9), round(res266.best.P*1e9), ...
        res266.best.coverage*100, res266.best.mean_T);
    fprintf('   320 only : h=%d P=%d  cov=%.1f%% T=%.3f\n\n', ...
        round(res320.best.h*1e9), round(res320.best.P*1e9), ...
        res320.best.coverage*100, res320.best.mean_T);

    cmb = struct();
    cmb.table=tbl; cmb.best=best;
    cmb.best_266=res266.best; cmb.best_320=res320.best;
    cmb.cfg=cfg;

    if cfg.plot, plot_compare(tbl, best, cfg); end
end

% ====================================================================
function [cov266, T266vec, cov320, T320vec] = joint_analysis(mat_paths_266, mat_paths_320, cfg)
% 각 (h, P) 에서 두 파장 라이브러리 합쳐 joint bin 선택 후 (cov, mean T) 반환
    N = numel(mat_paths_266);
    cov266 = zeros(N,1); T266vec = zeros(N,1);
    cov320 = zeros(N,1); T320vec = zeros(N,1);

    for k = 1:N
        S1 = load(mat_paths_266{k}); L266 = S1.lib;
        S2 = load(mat_paths_320{k}); L320 = S2.lib;
        if numel(L266.t) ~= numel(L320.t)
            warning('[%d] 구조 개수 불일치 (266=%d, 320=%d)', k, numel(L266.t), numel(L320.t));
            continue;
        end

        T1 = L266.T(:); T2 = L320.T(:);
        phi1 = L266.phi(:); phi2 = L320.phi(:);

        valid = (T1 >= cfg.T_min) & (T2 >= cfg.T_min);

        switch cfg.joint_obj
            case 'T_prod', score = T1 .* T2;
            case 'T_min',  score = min(T1, T2);
            case 'T_sum',  score = (T1 + T2) / 2;
            otherwise,     score = T1 .* T2;
        end

        if strcmpi(cfg.primary, 'lam266')
            phi_prim = phi1;
        else
            phi_prim = phi2;
        end

        edges = linspace(-pi, pi, cfg.N_bin+1);
        sel = nan(cfg.N_bin, 1);
        for b = 1:cfg.N_bin
            m = valid & (phi_prim >= edges(b)) & (phi_prim < edges(b+1));
            if any(m)
                idxs = find(m);
                [~, kk] = max(score(idxs));
                sel(b) = idxs(kk);
            end
        end

        picked = ~isnan(sel);
        idx = sel(picked);
        cov266(k) = mean(picked);   % primary cov (266 if cfg.primary='lam266')
        cov320(k) = cov266(k);      % joint 라 같은 bin 그룹
        if any(picked)
            T266vec(k) = mean(T1(idx));
            T320vec(k) = mean(T2(idx));
        end
    end
end

% ====================================================================
function joint_txt_path = save_joint_txt(mat266, mat320, cfg)
% best (h,P) 의 joint library 를 8-col txt 로 저장
    S1 = load(mat266); L266 = S1.lib;
    S2 = load(mat320); L320 = S2.lib;
    T1 = L266.T(:); T2 = L320.T(:);
    phi1 = L266.phi(:); phi2 = L320.phi(:);
    valid = (T1 >= cfg.T_min) & (T2 >= cfg.T_min);

    switch cfg.joint_obj
        case 'T_prod', score = T1 .* T2;
        case 'T_min',  score = min(T1, T2);
        case 'T_sum',  score = (T1 + T2) / 2;
        otherwise,     score = T1 .* T2;
    end

    if strcmpi(cfg.primary, 'lam266')
        phi_prim = phi1; tag = '266';
    else
        phi_prim = phi2; tag = '320';
    end

    edges = linspace(-pi, pi, cfg.N_bin+1);
    M = [];
    for b = 1:cfg.N_bin
        m = valid & (phi_prim >= edges(b)) & (phi_prim < edges(b+1));
        if ~any(m), continue; end
        idxs = find(m);
        [~, kk] = max(score(idxs));
        i = idxs(kk);
        switch L266.shape{i}
            case 'cyl',    sid = 1; r_eff = L266.p1(i);
            case 'sq',     sid = 2; r_eff = sqrt(L266.p1(i)^2/pi);
            case 'cross',  sid = 3; r_eff = sqrt((2*L266.p1(i)*L266.p2(i)-L266.p2(i)^2)/pi);
            case 'xcross', sid = 4; r_eff = sqrt((2*L266.p1(i)*L266.p2(i)-L266.p2(i)^2)/pi);
        end
        p2v = L266.p2(i); if isnan(p2v), p2v = 0; end
        M(end+1,:) = [r_eff, phi1(i), sqrt(T1(i)), phi2(i), sqrt(T2(i)), ...
                      sid, L266.p1(i), p2v]; %#ok<AGROW>
    end

    [~,fname,~] = fileparts(mat266);
    tok = regexp(fname,'_h(\d+)nm_P(\d+)nm','tokens','once');
    h_nm = tok{1}; P_nm = tok{2};
    joint_txt_path = sprintf('joint_h%snm_P%snm_primary%s.txt', h_nm, P_nm, tag);
    fid = fopen(joint_txt_path,'w');
    fprintf(fid,'%% joint library (same structure used at 266 & 320 nm)\n');
    fprintf(fid,'%% h=%s nm  P=%s nm  primary=%s  obj=%s  T_min=%.2f\n', ...
        h_nm, P_nm, tag, cfg.joint_obj, cfg.T_min);
    fprintf(fid,'%% col1=r_eff[m] col2=phi266 col3=|t266| col4=phi320 col5=|t320| col6=shape(1cyl/2sq/3cross/4xcross) col7=p1[m] col8=p2[m]\n');
    fprintf(fid,'%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%d\t%.12e\t%.12e\n', M.');
    fclose(fid);
end

% ====================================================================
function plot_compare(tbl, best, cfg)
    figure('Position', [80 80 1200 850], 'Color', 'w');

    h_vec = unique(tbl.h_m*1e9);
    P_vec = unique(tbl.P_m*1e9);
    [HH, ~] = meshgrid(h_vec, P_vec);

    C266 = nan(size(HH)); T266m = nan(size(HH));
    C320 = nan(size(HH)); T320m = nan(size(HH));
    SCO  = nan(size(HH));
    for k = 1:height(tbl)
        ih = find(abs(h_vec - tbl.h_m(k)*1e9) < 1e-6);
        ip = find(abs(P_vec - tbl.P_m(k)*1e9) < 1e-6);
        C266(ip,ih)=tbl.cov266(k); T266m(ip,ih)=tbl.T266(k);
        C320(ip,ih)=tbl.cov320(k); T320m(ip,ih)=tbl.T320(k);
        SCO(ip,ih) =tbl.score(k);
    end

    hx = best.h*1e9; py = best.P*1e9;

    subplot(2,3,1); imagesc(h_vec,P_vec,C266); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('266 cov'); hold on;
    plot(hx,py,'kp','MarkerSize',18,'MarkerFaceColor','y','LineWidth',1.3);

    subplot(2,3,2); imagesc(h_vec,P_vec,T266m); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('266 mean T'); hold on;
    plot(hx,py,'kp','MarkerSize',18,'MarkerFaceColor','y','LineWidth',1.3);

    subplot(2,3,3); imagesc(h_vec,P_vec,SCO); axis xy; colorbar;
    xlabel('h (nm)'); ylabel('P (nm)');
    title(sprintf('combined score (%s)',cfg.score_mode)); hold on;
    plot(hx,py,'kp','MarkerSize',18,'MarkerFaceColor','y','LineWidth',1.3);

    subplot(2,3,4); imagesc(h_vec,P_vec,C320); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('320 cov'); hold on;
    plot(hx,py,'kp','MarkerSize',18,'MarkerFaceColor','y','LineWidth',1.3);

    subplot(2,3,5); imagesc(h_vec,P_vec,T320m); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('320 mean T'); hold on;
    plot(hx,py,'kp','MarkerSize',18,'MarkerFaceColor','y','LineWidth',1.3);

    subplot(2,3,6); hold on;
    scatter(tbl.T266,tbl.T320,60,tbl.score,'filled');
    plot(best.T266,best.T320,'kp','MarkerSize',18,'MarkerFaceColor','y','LineWidth',1.3);
    plot([0 1],[0 1],'k:'); xlabel('T 266'); ylabel('T 320');
    xlim([0 1]); ylim([0 1]); colorbar; grid on; box on;
    title('T vs T (star=BEST)');
    sgtitle(sprintf('%s mode | BEST h=%d P=%d', upper(cfg.mode), ...
        round(best.h*1e9), round(best.P*1e9)),'FontWeight','bold');
end

% ====================================================================
function cfg = set_defaults_cmp(cfg)
    if ~isfield(cfg,'mode')         || isempty(cfg.mode),         cfg.mode = 'joint'; end
    if ~isfield(cfg,'primary')      || isempty(cfg.primary),      cfg.primary = 'lam266'; end
    if ~isfield(cfg,'T_min')        || isempty(cfg.T_min),        cfg.T_min = 0.5; end
    if ~isfield(cfg,'N_bin')        || isempty(cfg.N_bin),        cfg.N_bin = 64; end
    if ~isfield(cfg,'joint_obj')    || isempty(cfg.joint_obj),    cfg.joint_obj = 'T_prod'; end
    if ~isfield(cfg,'save_joint_txt')|| isempty(cfg.save_joint_txt), cfg.save_joint_txt = true; end
    if ~isfield(cfg,'score_mode')   || isempty(cfg.score_mode),   cfg.score_mode = 'worst'; end
    if ~isfield(cfg,'w_cov')        || isempty(cfg.w_cov),        cfg.w_cov = 2; end
    if ~isfield(cfg,'w_T')          || isempty(cfg.w_T),          cfg.w_T = 1; end
    if ~isfield(cfg,'plot')         || isempty(cfg.plot),         cfg.plot = true; end
    if ~isfield(cfg,'topN')         || isempty(cfg.topN),         cfg.topN = 10; end
end
