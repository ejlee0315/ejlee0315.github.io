function cmb = compare_266_320(res266, res320, cfg)
% compare_266_320  같은 (h, P) 그리드에서 두 파장 동시 성능 분석
%
%   >> compare_266_320                                    % 가장 간단!
%      -> 현재 폴더의 search_results_266_320.mat 자동 로드
%
%   >> cmb = compare_266_320(res266, res320);             % 워크스페이스 변수 직접 전달
%   >> cmb = compare_266_320('path/to/results.mat');      % 파일 경로 지정
%   >> cmb = compare_266_320(res266, res320, cfg);        % 옵션 지정
%
%   cfg 옵션:
%     .score_mode : 'worst'(기본) | 'product' | 'mean' | 'weighted'
%     .w_cov      : coverage 가중치 (기본 2)
%     .w_T        : mean_T 가중치 (기본 1)
%     .plot       : 요약 플롯 (기본 true)
%     .topN       : 상위 N개 출력 (기본 10)
%
%   score_mode 정의 (coverage, T 모두 [0,1] 범위):
%     'worst'    : min(cov266, cov320)^w_cov * min(T266, T320)^w_T   ← 보수적
%     'product'  : (cov266*cov320)^(w_cov/2) * (T266*T320)^(w_T/2)
%     'mean'     : ((cov266+cov320)/2)^w_cov * ((T266+T320)/2)^w_T
%     'weighted' : 0.5*(cov266^w_cov*T266^w_T) + 0.5*(cov320^w_cov*T320^w_T)
%
%   OUTPUT cmb:
%     .table     : (h, P, cov266, T266, cov320, T320, score) 전체 행
%     .best      : 종합 최고 (h, P, txt_path_266, txt_path_320, ...)
%     .best_266  : 266 단독 기준 best  (참고용)
%     .best_320  : 320 단독 기준 best  (참고용)

    %% 입력 파싱 -----------------------------------------------------------
    %  인자 없이 호출 → 현재 폴더에서 자동 로드
    if nargin == 0
        default_mat = 'search_results_266_320.mat';
        if ~exist(default_mat, 'file')
            error(['%s 가 현재 폴더에 없습니다. 먼저 run_optimize_266_320 을 ' ...
                   '실행하거나, cmb=compare_266_320(<.mat 경로>) 로 호출하세요.'], default_mat);
        end
        fprintf('[compare_266_320] auto-loading %s ...\n', default_mat);
        S = load(default_mat);
        res266 = S.res266; res320 = S.res320;
        cfg = struct();
    elseif nargin == 1 && (ischar(res266) || isstring(res266))
        S = load(res266); res266 = S.res266; res320 = S.res320;
        cfg = struct();
    elseif nargin < 3
        cfg = struct();
    end
    cfg = set_defaults_cmp(cfg);

    t266 = res266.table;   % columns: h_m, P_m, coverage, meanT, score
    t320 = res320.table;

    %% 같은 (h, P) 결합 --------------------------------------------------
    if height(t266) ~= height(t320)
        error('266/320 의 (h,P) 그리드 크기가 다름 (%d vs %d)', ...
            height(t266), height(t320));
    end

    % 정렬 (h, P) 기준 — 같은 순서라고 가정하지만 혹시 모르니 key 맞춤
    key266 = t266.h_m*1e12 + t266.P_m*1e6;   % 고유 key
    key320 = t320.h_m*1e12 + t320.P_m*1e6;
    [~, i266] = sort(key266);
    [~, i320] = sort(key320);
    t266 = t266(i266,:); t320 = t320(i320,:);

    if ~all(abs(t266.h_m - t320.h_m) < 1e-15) || ...
       ~all(abs(t266.P_m - t320.P_m) < 1e-15)
        error('(h,P) 그리드가 일치하지 않음. run_optimize_266_320 을 다시 실행하세요.');
    end

    % all_files 재정렬해서 txt_path 매칭
    files266 = res266.all(i266, :);
    files320 = res320.all(i320, :);

    %% 결합 테이블 --------------------------------------------------------
    N = height(t266);
    h_m  = t266.h_m;   P_m = t266.P_m;
    c266 = t266.coverage;  T266 = t266.meanT;
    c320 = t320.coverage;  T320 = t320.meanT;

    switch cfg.score_mode
        case 'worst'
            score = min(c266,c320).^cfg.w_cov .* min(T266,T320).^cfg.w_T;
        case 'product'
            score = (c266.*c320).^(cfg.w_cov/2) .* (T266.*T320).^(cfg.w_T/2);
        case 'mean'
            score = ((c266+c320)/2).^cfg.w_cov .* ((T266+T320)/2).^cfg.w_T;
        case 'weighted'
            score = 0.5*(c266.^cfg.w_cov .* T266.^cfg.w_T) + ...
                    0.5*(c320.^cfg.w_cov .* T320.^cfg.w_T);
        otherwise
            error('Unknown score_mode: %s', cfg.score_mode);
    end

    tbl = table(h_m, P_m, c266, T266, c320, T320, score, ...
        'VariableNames', {'h_m','P_m','cov266','T266','cov320','T320','score'});

    %% 상위 N 정렬 -------------------------------------------------------
    [~, order] = sort(score, 'descend');
    tbl = tbl(order, :);
    files266 = files266(order, :);
    files320 = files320(order, :);

    topN = min(cfg.topN, N);
    fprintf('\n========== combined best (score_mode=%s, w_cov=%.1f, w_T=%.1f) ==========\n', ...
        cfg.score_mode, cfg.w_cov, cfg.w_T);
    fprintf(' rank |  h   |  P   | cov266  T266  | cov320  T320  | score\n');
    fprintf(' -----+------+------+---------------+---------------+-------\n');
    for k = 1:topN
        fprintf(' %3d  | %4d | %4d | %5.1f%% %5.3f | %5.1f%% %5.3f | %.4f\n', ...
            k, round(tbl.h_m(k)*1e9), round(tbl.P_m(k)*1e9), ...
            tbl.cov266(k)*100, tbl.T266(k), ...
            tbl.cov320(k)*100, tbl.T320(k), ...
            tbl.score(k));
    end
    fprintf('=====================================================================\n\n');

    %% 최종 best 요약 ----------------------------------------------------
    best = struct();
    best.h          = tbl.h_m(1);
    best.P          = tbl.P_m(1);
    best.cov266     = tbl.cov266(1);
    best.T266       = tbl.T266(1);
    best.cov320     = tbl.cov320(1);
    best.T320       = tbl.T320(1);
    best.score      = tbl.score(1);
    best.mat_266    = files266{1, 1};
    best.txt_266    = files266{1, 2};
    best.mat_320    = files320{1, 1};
    best.txt_320    = files320{1, 2};

    fprintf('>> BEST (combined) : h=%d nm, P=%d nm\n', ...
        round(best.h*1e9), round(best.P*1e9));
    fprintf('   266nm : coverage=%.1f%%  meanT=%.3f  lib=%s\n', ...
        best.cov266*100, best.T266, best.txt_266);
    fprintf('   320nm : coverage=%.1f%%  meanT=%.3f  lib=%s\n', ...
        best.cov320*100, best.T320, best.txt_320);

    % 개별 파장 best (참고)
    best_266 = res266.best;
    best_320 = res320.best;
    fprintf('\n(참고) 단독 best\n');
    fprintf('   266 only : h=%d P=%d  cov=%.1f%% T=%.3f\n', ...
        round(best_266.h*1e9), round(best_266.P*1e9), ...
        best_266.coverage*100, best_266.mean_T);
    fprintf('   320 only : h=%d P=%d  cov=%.1f%% T=%.3f\n\n', ...
        round(best_320.h*1e9), round(best_320.P*1e9), ...
        best_320.coverage*100, best_320.mean_T);

    %% 결과 구조체 -------------------------------------------------------
    cmb = struct();
    cmb.table    = tbl;
    cmb.best     = best;
    cmb.best_266 = best_266;
    cmb.best_320 = best_320;
    cmb.cfg      = cfg;

    %% 플롯 --------------------------------------------------------------
    if cfg.plot
        plot_compare(tbl, best, cfg);
    end
end

% ========================================================================
function plot_compare(tbl, best, cfg)
    figure('Position', [80 80 1200 850], 'Color', 'w');

    h_vec = unique(tbl.h_m*1e9);
    P_vec = unique(tbl.P_m*1e9);
    [HH, PP] = meshgrid(h_vec, P_vec);

    % 2D 재배열
    C266 = nan(size(HH)); T266m = nan(size(HH));
    C320 = nan(size(HH)); T320m = nan(size(HH));
    SCO  = nan(size(HH));
    for k = 1:height(tbl)
        ih = find(abs(h_vec - tbl.h_m(k)*1e9) < 1e-6);
        ip = find(abs(P_vec - tbl.P_m(k)*1e9) < 1e-6);
        C266(ip, ih) = tbl.cov266(k);
        T266m(ip, ih) = tbl.T266(k);
        C320(ip, ih) = tbl.cov320(k);
        T320m(ip, ih) = tbl.T320(k);
        SCO(ip, ih)  = tbl.score(k);
    end

    hp_mark = @() plot(best.h*1e9, best.P*1e9, 'kp', ...
        'MarkerSize', 18, 'MarkerFaceColor', 'y', 'LineWidth', 1.3);

    subplot(2,3,1); imagesc(h_vec, P_vec, C266); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('266 nm coverage'); hold on; hp_mark();

    subplot(2,3,2); imagesc(h_vec, P_vec, T266m); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('266 nm mean T'); hold on; hp_mark();

    subplot(2,3,3); imagesc(h_vec, P_vec, SCO); axis xy; colorbar;
    xlabel('h (nm)'); ylabel('P (nm)');
    title(sprintf('combined score (%s)', cfg.score_mode)); hold on; hp_mark();

    subplot(2,3,4); imagesc(h_vec, P_vec, C320); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('320 nm coverage'); hold on; hp_mark();

    subplot(2,3,5); imagesc(h_vec, P_vec, T320m); axis xy; colorbar; caxis([0 1]);
    xlabel('h (nm)'); ylabel('P (nm)'); title('320 nm mean T'); hold on; hp_mark();

    subplot(2,3,6); hold on;
    scatter(tbl.T266, tbl.T320, 60, tbl.score, 'filled');
    plot(best.T266, best.T320, 'kp', 'MarkerSize', 18, ...
         'MarkerFaceColor', 'y', 'LineWidth', 1.3);
    plot([0 1], [0 1], 'k:', 'LineWidth', 0.8);
    xlabel('T 266 nm'); ylabel('T 320 nm');
    xlim([0 1]); ylim([0 1]); colorbar; grid on; box on;
    title('T vs T (stars = BEST)');
    sgtitle(sprintf('266 vs 320 nm comparison  |  BEST h=%d P=%d', ...
        round(best.h*1e9), round(best.P*1e9)), 'FontWeight', 'bold');
end

% ========================================================================
function cfg = set_defaults_cmp(cfg)
    if ~isfield(cfg,'score_mode') || isempty(cfg.score_mode), cfg.score_mode = 'worst'; end
    if ~isfield(cfg,'w_cov')      || isempty(cfg.w_cov),      cfg.w_cov = 2; end
    if ~isfield(cfg,'w_T')        || isempty(cfg.w_T),        cfg.w_T   = 1; end
    if ~isfield(cfg,'plot')       || isempty(cfg.plot),       cfg.plot  = true; end
    if ~isfield(cfg,'topN')       || isempty(cfg.topN),       cfg.topN  = 10; end
end
