function opt = A_Rect_opt(lib, opt_cfg)
% A_Rect_opt  full 2pi 위상 커버리지 + 최대 투과율이 되도록 라이브러리 최적화
%
%   opt = A_Rect_opt(lib)                % lib 구조체 직접 입력
%   opt = A_Rect_opt('lib_266nm.mat')    % .mat 파일 경로 입력
%   opt = A_Rect_opt(lib, opt_cfg)
%
%   opt_cfg 필드:
%     .N_bin      : 위상 bin 수 (기본 64)
%     .T_min      : 허용 최소 투과율 (기본 0.5)
%     .prefer     : 동점일 때 선호 구조 order (기본 {'cyl','sq','cross'})
%     .txt_path   : 출력 .txt 경로 (기본: lib_<lam>nm_opt.txt)
%     .plot       : 플롯 여부 (기본 true)
%     .allow_wrap : phase bin이 비었을 때 가까운 bin 차용 허용 (기본 true)
%
%   OUTPUT 구조체 opt:
%     .phase_target     : bin 중심 [-pi, pi]
%     .phase_atom       : 실제 선택된 구조의 위상
%     .T_atom           : 선택된 구조의 투과율
%     .shape_id         : 1=cyl, 2=sq, 3=cross
%     .p1, .p2          : 파라미터 [m]
%     .complex_t        : 복소 투과계수
%     .phase_error      : (target - atom) wrapped
%     .coverage         : 2pi 커버 여부
%     .mean_T           : 선택된 구조들의 평균 T
%
%   또한 TXT 라이브러리 저장. TXT 컬럼:
%     col1 = effective_radius [m]   (렌즈 생성 코드용: sqrt(area/pi))
%     col2 = phase  [rad]
%     col3 = |t|    (amplitude)
%     col4 = shape_id  (1=cyl / 2=sq / 3=cross / 4=xcross)
%     col5 = p1 [m]
%     col6 = p2 [m]  (cross, xcross 만 의미있음)

    if ischar(lib) || isstring(lib)
        S = load(lib); lib = S.lib;
    end

    if nargin < 2, opt_cfg = struct(); end
    opt_cfg = set_defaults_opt(opt_cfg, lib);

    N = numel(lib.t);
    shape_id_all = zeros(N, 1);
    for i = 1:N
        switch lib.shape{i}
            case 'cyl',    shape_id_all(i) = 1;
            case 'sq',     shape_id_all(i) = 2;
            case 'cross',  shape_id_all(i) = 3;
            case 'xcross', shape_id_all(i) = 4;
        end
    end

    % bin 경계
    N_bin = opt_cfg.N_bin;
    edges = linspace(-pi, pi, N_bin+1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    phi = lib.phi(:);
    T   = lib.T(:);
    t   = lib.t(:);
    p1  = lib.p1(:);
    p2  = lib.p2(:);

    % 최소 T 필터
    valid = T >= opt_cfg.T_min;
    fprintf('[A_Rect_opt] 총 %d 구조 중 T>=%.2f 만족: %d개\n', ...
        N, opt_cfg.T_min, sum(valid));

    % 각 bin마다 최대 T 구조 선택
    sel = nan(N_bin, 1);
    for b = 1:N_bin
        in_bin = (phi >= edges(b)) & (phi < edges(b+1)) & valid;
        if ~any(in_bin)
            if opt_cfg.allow_wrap
                % 모든 valid 구조 중 bin 중심과 위상차가 가장 가까운 것 중 최대 T
                d = wrap_diff(phi, centers(b));
                cand = valid & (abs(d) < pi/N_bin * opt_cfg.wrap_range);
                if any(cand)
                    idxs = find(cand);
                    [~, kk] = max(T(idxs));
                    sel(b) = idxs(kk);
                end
            end
            continue;
        end
        idxs = find(in_bin);
        % 우선 최대 T, 동점이면 prefer 순서대로
        Tb = T(idxs);
        [Tmax, ~] = max(Tb);
        tie = idxs(Tb > Tmax - 1e-6);
        if numel(tie) == 1
            sel(b) = tie;
        else
            prefer_ids = opt_cfg.prefer_ids;
            for sid = prefer_ids
                found = tie(shape_id_all(tie) == sid);
                if ~isempty(found)
                    sel(b) = found(1); break;
                end
            end
            if isnan(sel(b)), sel(b) = tie(1); end
        end
    end

    picked = ~isnan(sel);
    coverage = mean(picked);
    mean_T = mean(T(sel(picked)));

    fprintf('[A_Rect_opt] 2pi 커버리지: %.1f%% (%d/%d bins), 평균 T = %.3f\n', ...
        coverage*100, sum(picked), N_bin, mean_T);

    % 결과 패키지
    opt = struct();
    opt.phase_target = centers(:);
    opt.phase_atom   = nan(N_bin,1);
    opt.T_atom       = nan(N_bin,1);
    opt.shape_id     = nan(N_bin,1);
    opt.p1           = nan(N_bin,1);
    opt.p2           = nan(N_bin,1);
    opt.complex_t    = nan(N_bin,1) + 1i*nan(N_bin,1);
    for b = 1:N_bin
        if ~picked(b), continue; end
        k = sel(b);
        opt.phase_atom(b) = phi(k);
        opt.T_atom(b)     = T(k);
        opt.shape_id(b)   = shape_id_all(k);
        opt.p1(b)         = p1(k);
        opt.p2(b)         = p2(k);
        opt.complex_t(b)  = t(k);
    end
    opt.phase_error = wrap_diff(opt.phase_atom, opt.phase_target);
    opt.coverage = coverage;
    opt.mean_T = mean_T;
    opt.cfg_lib = lib.cfg;
    opt.cfg_opt = opt_cfg;

    % TXT 저장: effective radius (렌즈 GDS 생성용 circ 가정 호환)
    r_eff = eff_radius(opt.shape_id, opt.p1, opt.p2);
    M = [r_eff, opt.phase_atom, abs(opt.complex_t), opt.shape_id, opt.p1, opt.p2];
    M = M(all(~isnan(M),2), :);
    fid = fopen(opt_cfg.txt_path, 'w');
    fprintf(fid, '%% col1=r_eff[m] col2=phase[rad] col3=|t| col4=shape(1=cyl,2=sq,3=cross) col5=p1[m] col6=p2[m]\n');
    fprintf(fid, '%.12e\t%.12e\t%.12e\t%d\t%.12e\t%.12e\n', M.');
    fclose(fid);
    fprintf('[A_Rect_opt] TXT 저장: %s  (%d 엔트리)\n', opt_cfg.txt_path, size(M,1));

    if opt_cfg.plot
        plot_result(lib, opt);
    end
end

% =====================================================================
function d = wrap_diff(a, b)
% (a - b) 를 [-pi, pi] 로 래핑
    d = mod(a - b + pi, 2*pi) - pi;
end

function r = eff_radius(sid, p1, p2)
% shape별 단면적에서 환산한 유효 반경
    r = nan(size(sid));
    for i = 1:numel(sid)
        switch sid(i)
            case 1  % cyl
                r(i) = p1(i);
            case 2  % square side = p1
                A = p1(i)^2;
                r(i) = sqrt(A/pi);
            case 3  % cross (+) : L=p1, W=p2, area = 2*L*W - W^2
                A = 2*p1(i)*p2(i) - p2(i)^2;
                r(i) = sqrt(A/pi);
            case 4  % xcross (×, 45도 회전): 면적은 + 와 동일 (회전만)
                A = 2*p1(i)*p2(i) - p2(i)^2;
                r(i) = sqrt(A/pi);
            otherwise
                r(i) = NaN;
        end
    end
end

function plot_result(lib, opt)
    figure('Position',[100 100 1100 800], 'Color','w');

    % (1) 산점도: 전체 라이브러리 (phase, T), 모양별 컬러
    subplot(2,2,1); hold on;
    cmap = [0.1 0.4 0.8; 0.85 0.3 0.1; 0.2 0.6 0.2; 0.6 0.3 0.8];
    names = {'cyl','sq','cross','xcross'};
    for sid = 1:4
        idx = strcmp(lib.shape, names{sid});
        scatter(lib.phi(idx), lib.T(idx), 18, cmap(sid,:), 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', names{sid});
    end
    % 선택된 구조 강조
    ok = ~isnan(opt.shape_id);
    scatter(opt.phase_atom(ok), opt.T_atom(ok), 70, 'k', 'o', 'LineWidth', 1.2, ...
        'DisplayName','selected');
    xlabel('phase (rad)'); ylabel('|t|^2'); xlim([-pi pi]); ylim([0 1]);
    legend('Location','best'); grid on; box on;
    title(sprintf('lam0=%.0f nm  lib(%d)',lib.cfg.lam0*1e9, numel(lib.t)));

    % (2) target phase vs. selected atom phase
    subplot(2,2,2); hold on;
    plot([-pi pi], [-pi pi], 'k--', 'DisplayName','ideal');
    scatter(opt.phase_target, opt.phase_atom, 40, opt.T_atom, 'filled', ...
        'DisplayName','selected');
    colormap(gca, 'parula'); caxis([0 1]); cb=colorbar; cb.Label.String='|t|^2';
    xlabel('target phase'); ylabel('atom phase');
    xlim([-pi pi]); ylim([-pi pi]); grid on; box on;
    title(sprintf('coverage=%.1f%%, mean T=%.3f', opt.coverage*100, opt.mean_T));

    % (3) phase error histogram
    subplot(2,2,3);
    histogram(opt.phase_error(~isnan(opt.phase_error)), 30);
    xlabel('phase error (rad)'); ylabel('#bins');
    grid on; box on; title('phase error');

    % (4) shape 별 분포 (선택된 것만)
    subplot(2,2,4);
    shape_counts = histcounts(opt.shape_id(~isnan(opt.shape_id)), 0.5:1:4.5);
    bar(shape_counts); set(gca,'XTickLabel',{'cyl','sq','cross','xcross'});
    ylabel('# selected atoms'); grid on; box on;
    title('shape usage');
end

% =====================================================================
function cfg = set_defaults_opt(cfg, lib)
    if ~isfield(cfg,'N_bin')    || isempty(cfg.N_bin),    cfg.N_bin = 64; end
    if ~isfield(cfg,'T_min')    || isempty(cfg.T_min),    cfg.T_min = 0.5; end
    if ~isfield(cfg,'prefer')   || isempty(cfg.prefer),   cfg.prefer = {'cyl','sq','cross','xcross'}; end
    if ~isfield(cfg,'plot')     || isempty(cfg.plot),     cfg.plot = true; end
    if ~isfield(cfg,'allow_wrap')|| isempty(cfg.allow_wrap), cfg.allow_wrap = true; end
    if ~isfield(cfg,'wrap_range')|| isempty(cfg.wrap_range), cfg.wrap_range = 2; end
    if ~isfield(cfg,'txt_path') || isempty(cfg.txt_path)
        cfg.txt_path = sprintf('lib_%dnm_opt.txt', round(lib.cfg.lam0*1e9));
    end
    ids = zeros(1, numel(cfg.prefer));
    for i = 1:numel(cfg.prefer)
        switch cfg.prefer{i}
            case 'cyl',    ids(i) = 1;
            case 'sq',     ids(i) = 2;
            case 'cross',  ids(i) = 3;
            case 'xcross', ids(i) = 4;
        end
    end
    cfg.prefer_ids = ids;
end
