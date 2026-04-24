function result = A_RECT_dualband(cfg)
% A_RECT_dualband  joint 라이브러리로 한 장의 물리 렌즈 만들고 266/320 동시 시뮬
%
%   한 번의 구조 선택(primary 파장 phase 타깃) -> 같은 구조 맵으로 양쪽 파장 시뮬.
%   출력: 한 개의 GDS + (FWHM, 효율, 위상맵, 초점 이미지) × 2파장
%
%   cfg:
%     .lib_txt       : 8-col joint 라이브러리 경로 (필수)
%     .D             : 렌즈 지름 [m]
%     .P             : 주기 [m]
%     .ni            : 입사 매질 굴절률 (기본 1)
%     .primary_wl    : 'lam266'(기본) | 'lam320'  (구조 선택 기준 파장)
%     .NA_primary    : primary 파장에서의 NA (기본 0.05)
%     .NA_secondary  : secondary 파장에서의 NA (기본 0.2)  -> f_secondary 결정
%                       없으면 f_primary 와 동일한 거리로 propagate
%     .gds_name      : GDS 파일 basename
%     .nc_poly       : 원기둥 폴리곤 꼭지점 수 (기본 24)
%     .layers        : [cyl, sq, cross, xcross] 레이어 번호 (기본 [10 20 30 40])
%     .inc_angle     : 입사각 [deg] (기본 0)
%     .view_um       : 초점 이미지 view 범위 (기본 60)
%     .plot          : 기본 true
%
%   OUTPUT result:
%     .primary : {lam_nm, NA, f, FWHM, efficiency, I, ...}
%     .secondary: {lam_nm, NA, f, FWHM, efficiency, I, ...}
%     .shape_id_map, .p1_map, .p2_map      (단일 구조 맵)
%     .gds_path                             (한 개만)

    nm = 1e-9; um = 1e-6;
    if nargin < 1, cfg = struct(); end
    cfg = dualband_defaults(cfg);

    %% --------------------- 파장 / 초점거리 결정 ---------------------
    lam266 = 266*nm; lam320 = 320*nm;
    if strcmpi(cfg.primary_wl,'lam266')
        lam_p = lam266;  lam_s = lam320;  tag_p = '266'; tag_s = '320';
    else
        lam_p = lam320;  lam_s = lam266;  tag_p = '320'; tag_s = '266';
    end
    kp = 2*pi/lam_p;  ks = 2*pi/lam_s;

    theta_p = asin(cfg.NA_primary/cfg.ni);
    f_p = cfg.D/2/tan(theta_p);

    if isempty(cfg.NA_secondary)
        f_s = f_p;
    else
        theta_s = asin(cfg.NA_secondary/cfg.ni);
        f_s = cfg.D/2/tan(theta_s);
    end

    fprintf('\n===== A_RECT_dualband =====\n');
    fprintf('  lib       : %s\n', cfg.lib_txt);
    fprintf('  D / P     : %.1f um / %.1f nm\n', cfg.D/um, cfg.P/nm);
    fprintf('  primary   : lam = %.0f nm, NA = %.3f, f = %.2f um\n', lam_p/nm, cfg.NA_primary, f_p/um);
    fprintf('  secondary : lam = %.0f nm, NA = %.3f, f = %.2f um\n', lam_s/nm, cfg.NA_secondary, f_s/um);

    %% --------------------- grid + primary phase map ---------------------
    x_range = -cfg.D/2 + cfg.P/2 : cfg.P : cfg.D/2 - cfg.P/2;
    y_range = x_range;
    [XX, YY] = meshgrid(x_range, y_range);
    R2 = XX.^2 + YY.^2;
    mask_aperture = sqrt(R2) <= cfg.D/2;
    phase_map_p = mod(-cfg.ni*kp*(R2/(2*f_p)), 2*pi) - pi;
    phase_map_p(~mask_aperture) = 0;

    %% --------------------- joint library 로드 (8-col) ---------------------
    M = readmatrix(cfg.lib_txt,'CommentStyle','%');
    M = M(all(~isnan(M),2), :);
    if size(M,2) ~= 8
        error('joint library 는 8-col 이어야 합니다 (현재 %d col).', size(M,2));
    end
    r_lib   = M(:,1);
    phi266L = M(:,2); a266L = M(:,3);
    phi320L = M(:,4); a320L = M(:,5);
    sidL    = M(:,6);
    p1L     = M(:,7);  p2L = M(:,8);

    if strcmpi(cfg.primary_wl,'lam266')
        phi_p_lib = phi266L;  a_p_lib = a266L;
        phi_s_lib = phi320L;  a_s_lib = a320L;
    else
        phi_p_lib = phi320L;  a_p_lib = a320L;
        phi_s_lib = phi266L;  a_s_lib = a266L;
    end
    phi_p_lib_wr = mod(phi_p_lib + pi, 2*pi) - pi;

    fprintf('  library   : %d entries (cyl=%d, sq=%d, cross=%d, xcross=%d)\n', ...
        numel(sidL), sum(sidL==1), sum(sidL==2), sum(sidL==3), sum(sidL==4));

    %% --------------------- 각 픽셀에 구조 1개 선택 (primary phase 기준) ---------------------
    [Ny, Nx] = size(phase_map_p);
    shape_id_map = zeros(Ny,Nx);
    p1_map       = zeros(Ny,Nx);
    p2_map       = zeros(Ny,Nx);
    phi_err_p    = zeros(Ny,Nx);
    T_p_map      = zeros(Ny,Nx);
    T_s_map      = zeros(Ny,Nx);
    complex_p    = zeros(Ny,Nx);
    complex_s    = zeros(Ny,Nx);

    for idx = 1:(Ny*Nx)
        if ~mask_aperture(idx), continue; end
        target = phase_map_p(idx);
        d = mod(target - phi_p_lib_wr + pi, 2*pi) - pi;
        [~, k] = min(abs(d));
        shape_id_map(idx) = sidL(k);
        p1_map(idx)       = p1L(k);
        p2_map(idx)       = p2L(k);
        phi_err_p(idx)    = d(k);
        T_p_map(idx)      = a_p_lib(k)^2;
        T_s_map(idx)      = a_s_lib(k)^2;
        complex_p(idx)    = a_p_lib(k) * exp(1i*phi_p_lib(k));
        complex_s(idx)    = a_s_lib(k) * exp(1i*phi_s_lib(k));
    end

    %% --------------------- 회절 시뮬 (양 파장) ---------------------
    res_p = run_focus_sim(complex_p, x_range, y_range, cfg.D, cfg.P, lam_p, f_p, ...
                          mask_aperture, cfg, tag_p);
    res_s = run_focus_sim(complex_s, x_range, y_range, cfg.D, cfg.P, lam_s, f_s, ...
                          mask_aperture, cfg, tag_s);

    %% --------------------- 플롯 ---------------------
    if cfg.plot
        plot_dualband(x_range, y_range, phase_map_p, phi_err_p, T_p_map, T_s_map, ...
                      shape_id_map, res_p, res_s, cfg, lam_p, lam_s, tag_p, tag_s);
    end

    %% --------------------- GDS (한 장만) ---------------------
    gds_path = write_gds_one(XX, YY, shape_id_map, p1_map, p2_map, mask_aperture, cfg);
    fprintf('  GDS saved : %s\n', gds_path);

    %% --------------------- 리턴 ---------------------
    result = struct();
    result.primary   = struct('lam_nm',lam_p/nm,'NA',cfg.NA_primary,  'f',f_p, ...
                              'FWHM_um',res_p.FWHM,'E_focus_pct',res_p.E_focus_pct);
    result.secondary = struct('lam_nm',lam_s/nm,'NA',cfg.NA_secondary,'f',f_s, ...
                              'FWHM_um',res_s.FWHM,'E_focus_pct',res_s.E_focus_pct);
    result.shape_id_map = shape_id_map;
    result.p1_map = p1_map;  result.p2_map = p2_map;
    result.gds_path = gds_path;

    fprintf('\n>> SUMMARY (단일 물리 렌즈, 양 파장 시뮬)\n');
    fprintf('   %s nm : FWHM=%.3f um,  focus eff(3xFWHM)=%.2f%%\n', ...
        tag_p, res_p.FWHM, res_p.E_focus_pct);
    fprintf('   %s nm : FWHM=%.3f um,  focus eff(3xFWHM)=%.2f%%\n', ...
        tag_s, res_s.FWHM, res_s.E_focus_pct);
end

% ====================================================================
function res = run_focus_sim(u0, x_range, y_range, D, P, lam, f, mask_aperture, cfg, tag)
    um = 1e-6;
    u2 = RS_2D_Diffraction_FW(D, P, lam, f, u0);
    cen = round(length(u2)/2);
    half = round(cfg.view_um/2);
    si = max(1, cen-half);  ei = min(length(u2), cen+half);
    I = abs(u2(si:ei, si:ei)).^2;
    x_v = x_range(si:ei);  y_v = y_range(si:ei);

    [~, mi] = max(I(:));
    [rc, cc] = ind2sub(size(I), mi);
    prof = I(rc,:); prof = prof / max(prof);
    x_um = x_v/um;
    ab = prof >= 0.5; ed = diff(ab);
    ri = find(ed==1, 1, 'first'); fi = find(ed==-1, 1, 'last');
    if ~isempty(ri) && ~isempty(fi)
        x1 = interp1(prof(ri:ri+1), x_um(ri:ri+1), 0.5);
        x2 = interp1(prof(fi:fi+1), x_um(fi:fi+1), 0.5);
        FWHM = abs(x2 - x1);
    else
        FWHM = NaN;
    end

    [XU, YU] = meshgrid(x_v/um, y_v/um);
    x0 = XU(rc,cc); y0 = YU(rc,cc);
    R_um = sqrt((XU-x0).^2 + (YU-y0).^2);
    dA = (P/um)^2;
    if ~isnan(FWHM)
        mask3 = R_um <= 3*FWHM;
        E_focus = sum(I(mask3), 'all') * dA;
    else
        E_focus = NaN;
    end
    E_inc = sum(mask_aperture(:)) * (P/um)^2;
    E_pct = E_focus / E_inc * 100;

    res = struct('I',I, 'x_v',x_v, 'y_v',y_v, 'FWHM',FWHM, ...
                 'E_focus_pct',E_pct, 'prof',prof, 'x_um',x_um, 'tag',tag);
end

% ====================================================================
function plot_dualband(x_range, y_range, phase_map, phase_err, Tp, Ts, sid, ...
                        res_p, res_s, cfg, lam_p, lam_s, tag_p, tag_s)
    um = 1e-6;
    figure('Position',[60 60 1400 900],'Color','w');

    subplot(3,3,1); imagesc(x_range/um, y_range/um, phase_map); axis image xy;
    title(sprintf('target phase (%s nm)', tag_p)); colorbar; caxis([-pi pi]);

    subplot(3,3,2); imagesc(x_range/um, y_range/um, phase_err); axis image xy;
    title('phase error (primary)'); colorbar; caxis([-pi/4 pi/4]);

    subplot(3,3,3); imagesc(x_range/um, y_range/um, sid); axis image xy;
    title('shape map  (1cyl 2sq 3cross 4xcross)');
    colormap(gca,[1 1 1; 0.1 0.4 0.9; 0.9 0.4 0.1; 0.2 0.7 0.2; 0.6 0.3 0.8]);
    caxis([0 4]); colorbar;

    subplot(3,3,4); imagesc(x_range/um, y_range/um, Tp); axis image xy;
    title(sprintf('T map @ %s nm', tag_p)); colorbar; caxis([0 1]);

    subplot(3,3,5); imagesc(x_range/um, y_range/um, Ts); axis image xy;
    title(sprintf('T map @ %s nm (same structures!)', tag_s)); colorbar; caxis([0 1]);

    subplot(3,3,6); hold on;
    plot(res_p.x_um, res_p.prof, 'r-', 'LineWidth', 1.8, 'DisplayName', [tag_p ' nm']);
    plot(res_s.x_um, res_s.prof, 'b--','LineWidth', 1.8, 'DisplayName', [tag_s ' nm']);
    yline(0.5, 'k:'); grid on; box on;
    xlabel('x (\mum)'); ylabel('normalized I');
    legend('Location','best'); title('line profiles (center)');

    subplot(3,3,7); imagesc(res_p.x_v/um, res_p.y_v/um, res_p.I); axis image xy;
    colormap(gca,'hot'); colorbar;
    title(sprintf('|E|^2 @ %s nm  FWHM=%.2fum  \\eta=%.1f%%', ...
        tag_p, res_p.FWHM, res_p.E_focus_pct));

    subplot(3,3,8); imagesc(res_s.x_v/um, res_s.y_v/um, res_s.I); axis image xy;
    colormap(gca,'hot'); colorbar;
    title(sprintf('|E|^2 @ %s nm  FWHM=%.2fum  \\eta=%.1f%%', ...
        tag_s, res_s.FWHM, res_s.E_focus_pct));

    subplot(3,3,9);
    cnts = [sum(sid(:)==1) sum(sid(:)==2) sum(sid(:)==3) sum(sid(:)==4)];
    b = bar(cnts); b.FaceColor='flat';
    b.CData = [0.1 0.4 0.9; 0.9 0.4 0.1; 0.2 0.7 0.2; 0.6 0.3 0.8];
    set(gca,'XTickLabel',{'cyl','sq','cross','xcross'});
    grid on; box on; title('shape usage');

    sgtitle(sprintf('DUAL-BAND  primary=%s  P=%d h=(lib embeds)  |  %s→FWHM %.2fum, %s→FWHM %.2fum', ...
        tag_p, round(cfg.P*1e9), tag_p, res_p.FWHM, tag_s, res_s.FWHM), 'FontWeight','bold');
end

% ====================================================================
function gds_path = write_gds_one(XX, YY, sid, p1, p2, mask, cfg)
    gds_path = [cfg.gds_name '.txt'];
    fid = fopen(gds_path, 'w');
    fprintf(fid,'HEADER 3;\r\n');
    fprintf(fid,'BGNLIB;\r\n');
    fprintf(fid,['LIBNAME ', cfg.gds_name, ';\r\n']);
    fprintf(fid,'UNITS 1.000000e+000 1.000000e-009;\r\n');
    fprintf(fid,'BGNSTR;\r\n');
    fprintf(fid,['STRNAME ', cfg.gds_name, ';\r\n']);
    idxs = find(mask(:));
    for k = 1:numel(idxs)
        ii = idxs(k);
        s = sid(ii);
        if s == 0, continue; end
        cx = XX(ii); cy = YY(ii);
        switch s
            case 1
                r = p1(ii); if r<=0, continue; end
                th = linspace(0, 2*pi, cfg.nc_poly+1); th(end) = [];
                xy = [cx + r*cos(th); cy + r*sin(th)];
                write_poly(fid, xy, cfg.layers(1));
            case 2
                sz = p1(ii); if sz<=0, continue; end
                hx = sz/2;
                xy = [cx-hx cx+hx cx+hx cx-hx; cy-hx cy-hx cy+hx cy+hx];
                write_poly(fid, xy, cfg.layers(2));
            case 3
                L = p1(ii); W = p2(ii);
                if L<=0 || W<=0, continue; end
                hL = L/2; hW = W/2;
                xy = [hW hW hL hL hW hW -hW -hW -hL -hL -hW -hW;
                      hL hW hW -hW -hW -hL -hL -hW -hW hW hW hL];
                xy = xy + [cx; cy];
                write_poly(fid, xy, cfg.layers(3));
            case 4
                L = p1(ii); W = p2(ii);
                if L<=0 || W<=0, continue; end
                hL = L/2; hW = W/2;
                xy0 = [hW hW hL hL hW hW -hW -hW -hL -hL -hW -hW;
                       hL hW hW -hW -hW -hL -hL -hW -hW hW hW hL];
                R45 = [cosd(45) -sind(45); sind(45) cosd(45)];
                xy = R45*xy0 + [cx; cy];
                write_poly(fid, xy, cfg.layers(4));
        end
    end
    fprintf(fid,'ENDSTR\r\n'); fprintf(fid,'ENDLIB\r\n'); fclose(fid);
end

function write_poly(fid, xy_m, layer)
    xy_nm = round(xy_m*1e9);
    fprintf(fid,'BOUNDARY\r\n');
    fprintf(fid,'LAYER %d;\r\n', layer);
    fprintf(fid,'DATATYPE %d;\r\n', layer);
    fprintf(fid,'XY\r\n');
    fprintf(fid,'%d\t:\t%d\r\n', xy_nm);
    fprintf(fid,'ENDEL\r\n');
end

% ====================================================================
function cfg = dualband_defaults(cfg)
    nm = 1e-9; um = 1e-6;
    if ~isfield(cfg,'lib_txt') || isempty(cfg.lib_txt)
        error('cfg.lib_txt (8-col joint library) 필수');
    end
    if ~isfield(cfg,'D') || isempty(cfg.D),  cfg.D = 114*um; end
    if ~isfield(cfg,'P') || isempty(cfg.P),  cfg.P = 180*nm; end
    if ~isfield(cfg,'ni'),          cfg.ni = 1; end
    if ~isfield(cfg,'primary_wl'),  cfg.primary_wl = 'lam266'; end
    if ~isfield(cfg,'NA_primary'),  cfg.NA_primary = 0.05; end
    if ~isfield(cfg,'NA_secondary'),cfg.NA_secondary = 0.2; end
    if ~isfield(cfg,'gds_name') || isempty(cfg.gds_name)
        cfg.gds_name = sprintf('metalens_dualband_primary%s', cfg.primary_wl);
    end
    if ~isfield(cfg,'nc_poly'),     cfg.nc_poly = 24; end
    if ~isfield(cfg,'layers'),      cfg.layers = [10 20 30 40]; end
    if ~isfield(cfg,'inc_angle'),   cfg.inc_angle = 0; end
    if ~isfield(cfg,'view_um'),     cfg.view_um = 60; end
    if ~isfield(cfg,'plot'),        cfg.plot = true; end
end
