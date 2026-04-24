function res = A_Rect_search(cfg)
% A_Rect_search  (h, P) grid search 외부 루프 + 형상 sweep + bin 최적화
%
%   res = A_Rect_search(cfg)
%
%   cfg 필드 (모두 기본값 있음):
%     .lam0          [m]    설계 파장 (예: 266e-9)
%     .h_list        [m]    높이 후보 (예: [400 500 600 700]*nm)
%     .P_list        [m]    주기 후보 (예: [220 240 260 280]*nm). 없으면
%                            자동: P_list = lam0/(2*NA_max)*ratio_list
%     .NA_max               (P_list 자동 생성 시 사용, 기본 0.6)
%     .ratio_list           Nyquist 대비 ratio (기본 [0.6 0.7 0.8 0.9 1.0])
%     .nSiO2, .n_air        (기본 1.46, 1.0)
%     .nx, .ny, .eta        fmm harmonics, eta (기본 6, 6, 45)
%     .addpath_fmm          fmm 경로
%     .xlsx_path            HfO2 엑셀 경로
%     .base_layer_thickness 얇은 바닥 레진층 (기본 100*nm; 없애려면 0)
%
%     .opt           A_Rect_opt 의 cfg (예: .N_bin, .T_min, .prefer)
%
%     .score_fun     @(opt) 점수함수 (default: opt.coverage*opt.mean_T)
%     .save_dir      결과 저장 디렉토리 (기본 'opt_<lam>nm')
%     .verbose       (true)
%
%   OUTPUT res:
%     .table   : (h, P, coverage, mean_T, score)
%     .best    : 최고 score의 (lib, opt, h, P, txt_path, mat_path)
%     .all     : 모든 (h,P) 조합의 lib, opt 저장 파일명 리스트

    if nargin < 1, cfg = struct(); end
    cfg = set_defaults_search(cfg);

    if ~exist(cfg.save_dir, 'dir'); mkdir(cfg.save_dir); end

    H = cfg.h_list(:); Pset = cfg.P_list(:);
    nH = numel(H); nP = numel(Pset);
    nTot = nH * nP;

    if cfg.verbose
        fprintf('\n========================================================\n');
        fprintf(' A_Rect_search  lam0=%.1f nm  (h x P grid: %d x %d = %d)\n', ...
            cfg.lam0*1e9, nH, nP, nTot);
        fprintf('   h_list = %s nm\n', mat2str(round(H*1e9)));
        fprintf('   P_list = %s nm\n', mat2str(round(Pset*1e9)));
        fprintf('========================================================\n');
    end

    table_rows = zeros(nTot, 5);
    all_files = cell(nTot, 2);
    all_libs  = cell(nTot, 1);
    all_opts  = cell(nTot, 1);
    k = 0; t_all = tic;

    for ih = 1:nH
        for ip = 1:nP
            k = k + 1;
            h = H(ih); P = Pset(ip);
            tag = sprintf('lam%dnm_h%dnm_P%dnm', round(cfg.lam0*1e9), ...
                          round(h*1e9), round(P*1e9));
            mat_path = fullfile(cfg.save_dir, ['lib_' tag '.mat']);
            txt_path = fullfile(cfg.save_dir, ['lib_' tag '_opt.txt']);

            if cfg.verbose
                fprintf('\n[%d/%d] h=%dnm  P=%dnm  ->  %s\n', k, nTot, ...
                    round(h*1e9), round(P*1e9), tag);
            end

            % --- 라이브러리 sweep ---
            lib_cfg = build_lib_cfg(cfg, h, P, mat_path);
            lib = A_Rect_lib(lib_cfg);

            % --- 위상 bin 최적화 ---
            opt_cfg = cfg.opt;
            opt_cfg.txt_path = txt_path;
            opt_cfg.plot     = false;        % 그리드 도는 동안엔 plot 끔
            opt = A_Rect_opt(lib, opt_cfg);

            score = cfg.score_fun(opt);
            table_rows(k, :) = [h, P, opt.coverage, opt.mean_T, score];
            all_files(k, :) = {mat_path, txt_path};
            all_libs{k}  = lib;
            all_opts{k}  = opt;

            if cfg.verbose
                fprintf('  -> coverage=%.1f%%  mean_T=%.3f  score=%.4f\n', ...
                    opt.coverage*100, opt.mean_T, score);
            end
        end
    end

    % ==========================================================
    %  BEST 선정 규칙 (lexicographic):
    %    1순위: coverage 최대 (full 2pi cover 를 가장 잘 한 것)
    %    2순위: mean_T 최대 (그 중에서 투과율 가장 높은 것)
    % ==========================================================
    covs = table_rows(:,3);  Ts = table_rows(:,4);
    max_cov = max(covs);
    tie_idx = find(covs >= max_cov - 1e-9);
    [~, jj] = max(Ts(tie_idx));
    ib = tie_idx(jj);
    best = struct( ...
        'score',    table_rows(ib, 5), ...
        'h',        table_rows(ib, 1), ...
        'P',        table_rows(ib, 2), ...
        'coverage', covs(ib), ...
        'mean_T',   Ts(ib), ...
        'lib',      all_libs{ib}, ...
        'opt',      all_opts{ib}, ...
        'mat_path', all_files{ib, 1}, ...
        'txt_path', all_files{ib, 2});
    best.tag = sprintf('lam%dnm_h%dnm_P%dnm', round(cfg.lam0*1e9), ...
                       round(best.h*1e9), round(best.P*1e9));

    res = struct();
    res.table = array2table(table_rows, ...
        'VariableNames', {'h_m','P_m','coverage','meanT','score'});
    res.best = best;
    res.all  = all_files;
    res.cfg  = cfg;

    save(fullfile(cfg.save_dir, sprintf('search_summary_lam%dnm.mat', ...
        round(cfg.lam0*1e9))), 'res');

    if cfg.verbose
        fprintf('\n========================================================\n');
        fprintf(' DONE  (%.1f s)\n', toc(t_all));
        fprintf(' 선정 규칙: coverage 최대 -> 동률시 mean_T 최대\n');
        fprintf('\n  h[nm] | P[nm] | cov(%%) | meanT\n');
        fprintf('  -------+-------+--------+------\n');
        for k2 = 1:size(table_rows,1)
            marker = ''; if k2 == ib, marker = '  <-- BEST'; end
            fprintf('  %5d | %5d | %6.1f | %.3f%s\n', ...
                round(table_rows(k2,1)*1e9), round(table_rows(k2,2)*1e9), ...
                table_rows(k2,3)*100, table_rows(k2,4), marker);
        end
        fprintf('\n BEST  h=%d nm  P=%d nm  coverage=%.1f%%  meanT=%.3f\n', ...
            round(best.h*1e9), round(best.P*1e9), ...
            best.opt.coverage*100, best.opt.mean_T);
        fprintf(' best lib mat : %s\n', best.mat_path);
        fprintf(' best lib txt : %s\n', best.txt_path);
        fprintf('========================================================\n');

        % 베스트 결과는 plot 한번 그려줌
        opt_cfg = cfg.opt; opt_cfg.plot = true;
        opt_cfg.txt_path = best.txt_path;
        A_Rect_opt(best.lib, opt_cfg);
    end
end

% =====================================================================
function lib_cfg = build_lib_cfg(cfg, h, P, mat_path)
    nm = 1e-9;
    lib_cfg = struct();
    lib_cfg.lam0 = cfg.lam0;
    lib_cfg.P    = P;
    lib_cfg.h    = h;
    lib_cfg.nSiO2 = cfg.nSiO2;
    lib_cfg.n_air = cfg.n_air;
    lib_cfg.addpath_fmm = cfg.addpath_fmm;
    lib_cfg.xlsx_path   = cfg.xlsx_path;
    lib_cfg.nx = cfg.nx; lib_cfg.ny = cfg.ny; lib_cfg.eta = cfg.eta;
    lib_cfg.base_layer_thickness = cfg.base_layer_thickness;
    lib_cfg.save_mat = mat_path;

    % ===== 형상 sweep 범위 (새 규칙) =====
    %   원기둥 r : 40 ~ (P/2 - 40) nm
    %   사각기둥 s: 80 ~ (P  -  80) nm
    %   십자가 L : 80 ~ (P  -  80) nm
    %   십자가 W : 80 ~ (P  -  80) nm   (W < L 만 사용; W==L 은 사각과 동일)
    r_min   = 40*nm;   r_max = P/2 - 40*nm;
    sq_min  = 80*nm;   sq_max = P - 80*nm;
    crL_min = 80*nm;   crL_max = P - 80*nm;
    crW_min = 80*nm;   crW_max = P - 80*nm;

    if r_max   < r_min,   warning('[P=%dnm] cyl range 무효 (r_max<r_min)',   round(P*1e9)); end
    if sq_max  < sq_min,  warning('[P=%dnm] sq range 무효',  round(P*1e9)); end
    if crL_max < crL_min, warning('[P=%dnm] cross L range 무효', round(P*1e9)); end

    lib_cfg.cyl.r   = linspace(r_min,   max(r_min, r_max),   cfg.nR);
    lib_cfg.sq.s    = linspace(sq_min,  max(sq_min, sq_max), cfg.nS);
    lib_cfg.cross.L = linspace(crL_min, max(crL_min, crL_max), cfg.nL);
    lib_cfg.cross.W = linspace(crW_min, max(crW_min, crW_max), cfg.nW);
end

% =====================================================================
function cfg = set_defaults_search(cfg)
    nm = 1e-9;

    cfg = sd(cfg, 'lam0', 266*nm);
    cfg = sd(cfg, 'h_list', [400 500 600 700 800]*nm);
    cfg = sd(cfg, 'P_list', [200 210 220 230 240 250]*nm);

    cfg = sd(cfg, 'nSiO2', 1.46);
    cfg = sd(cfg, 'n_air', 1.0);
    cfg = sd(cfg, 'addpath_fmm', 'G:\내 드라이브\Reserach_Source\A. Simulation\fmm');
    cfg = sd(cfg, 'xlsx_path', '');
    cfg = sd(cfg, 'nx', 6); cfg = sd(cfg, 'ny', 6); cfg = sd(cfg, 'eta', 45);
    cfg = sd(cfg, 'base_layer_thickness', 100*nm);
    cfg = sd(cfg, 'verbose', true);
    cfg = sd(cfg, 'save_dir', sprintf('opt_%dnm', round(cfg.lam0*1e9)));

    cfg = sd(cfg, 'nR', 30);   % cylinder radius 샘플 수
    cfg = sd(cfg, 'nS', 30);   % square side 샘플 수
    cfg = sd(cfg, 'nL', 10);   % cross length 샘플 수
    cfg = sd(cfg, 'nW', 6);    % cross width 샘플 수

    if ~isfield(cfg, 'opt') || isempty(cfg.opt)
        cfg.opt = struct('N_bin', 64, 'T_min', 0.5, ...
                         'prefer', {{'cyl','sq','cross'}}, 'allow_wrap', true);
    end

    if ~isfield(cfg, 'score_fun') || isempty(cfg.score_fun)
        % 기본: coverage * mean_T  (둘 다 [0,1] 범위)
        cfg.score_fun = @(o) o.coverage * o.mean_T;
    end
end

function cfg = sd(cfg, f, v)
    if ~isfield(cfg, f) || isempty(cfg.(f)); cfg.(f) = v; end
end
