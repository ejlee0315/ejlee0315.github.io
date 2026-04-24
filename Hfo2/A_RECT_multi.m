function result = A_RECT_multi(cfg)
% A_RECT_multi  원기둥 / 사각기둥 / 십자가 혼용 메타렌즈 + GDS 생성
%
%   result = A_RECT_multi(cfg)
%
%   cfg 필드 (기본값은 내부 set_defaults 에서):
%     .lam0         설계 파장 [m]              (예: 266e-9 / 320e-9)
%     .D            렌즈 지름 [m]              (예: 114e-6)
%     .P            단위 셀 주기 [m]
%     .NA           목표 NA                   (예: 0.6)
%     .ni           입사 매질 굴절률           (기본 1)
%     .lib_txt      A_Rect_opt.m 가 만든 라이브러리 .txt 경로
%                   (col1=r_eff[m] col2=phase[rad] col3=|t|
%                    col4=shape(1=cyl,2=sq,3=cross,4=xcross) col5=p1[m] col6=p2[m])
%     .gds_name     GDS 파일 basename (확장자 제외)
%     .inc_angle    입사 각도 [deg]            (기본 0)
%     .nc_poly      원기둥 폴리곤 꼭짓점 수    (기본 24)
%     .layers       [cyl_layer, sq_layer, cross_layer, xcross_layer]  (기본 [10, 20, 30, 40])
%     .plot         플롯 여부                 (기본 true)
%     .view_um      diffraction view window (µm, 대략) (기본 45)
%     .factor       구조 축소량 [m]           (기본 0, A_RECT.m 원본과 동일)
%
%   OUTPUT 구조체 result:
%     .phase_map, .phase_err_map, .T_map, .shape_id_map
%     .f            초점거리 [m]
%     .E_focus_pct  3*FWHM 원 내부 에너지 / 입사 에너지 [%]
%     .FWHM_um
%     .gds_path
%
%   EXAMPLE
%     cfg.lam0 = 266e-9; cfg.D = 114e-6; cfg.P = 260e-9;
%     cfg.NA = 0.6; cfg.lib_txt = 'lib_266nm_opt.txt';
%     cfg.gds_name = 'metalens_266nm';
%     A_RECT_multi(cfg);

    if nargin < 1, cfg = struct(); end
    cfg = set_defaults(cfg);

    nm = 1e-9; um = 1e-6;

    %% --------------------- 기본 계산 ------------------------------------
    k0 = 2*pi/cfg.lam0;
    theta = asin(cfg.NA/cfg.ni);
    f = cfg.D/2 / tan(theta);
    DOF = cfg.lam0 / (1 - cos(asin(cfg.NA))) / um;

    if cfg.P > cfg.lam0/2/cfg.NA
        warning('P=%.1f nm > Nyquist(%.1f nm). Aliasing 가능.', ...
            cfg.P/nm, cfg.lam0/2/cfg.NA/nm);
    end

    fprintf('\n===== A_RECT_multi =====\n');
    fprintf('  lam0   = %.1f nm\n', cfg.lam0/nm);
    fprintf('  D      = %.2f um,  P = %.1f nm\n', cfg.D/um, cfg.P/nm);
    fprintf('  NA     = %.2f  ->  f = %.2f um,  DOF = %.2f um\n', cfg.NA, f/um, DOF);
    fprintf('  lib    = %s\n', cfg.lib_txt);

    %% --------------------- 격자 + 위상맵 ---------------------------------
    x_range = -cfg.D/2 + cfg.P/2 : cfg.P : cfg.D/2 - cfg.P/2;
    y_range = x_range;
    [XX, YY] = meshgrid(x_range, y_range);
    R2       = XX.^2 + YY.^2;

    % 얇은 렌즈 (paraxial) 위상
    phase_map = mod(-cfg.ni*k0.*(R2/(2*f)), 2*pi) - pi;
    mask_aperture = (sqrt(R2) <= cfg.D/2);
    phase_map(~mask_aperture) = 0;

    %% --------------------- 라이브러리 로드 --------------------------------
    lib = load_lib_txt(cfg.lib_txt);
    % lib.phase: col2, lib.T: col3^2 (|t|^2),  lib.shape_id, lib.p1, lib.p2
    fprintf('  library entries = %d  (cyl=%d, sq=%d, cross=%d)\n', ...
        numel(lib.phase), sum(lib.shape_id==1), ...
        sum(lib.shape_id==2), sum(lib.shape_id==3));

    %% --------------------- 각 픽셀 룩업 -----------------------------------
    [Ny, Nx] = size(phase_map);
    N_total = Ny * Nx;
    shape_id_map = zeros(Ny, Nx);
    p1_map       = zeros(Ny, Nx);
    p2_map       = zeros(Ny, Nx);
    phase_err_map= zeros(Ny, Nx);
    T_map        = zeros(Ny, Nx);
    complex_map  = zeros(Ny, Nx);

    lib_phase_wrapped = mod(lib.phase + pi, 2*pi) - pi;
    for idx = 1:N_total
        if ~mask_aperture(idx), continue; end
        phi = phase_map(idx);
        d = mod(phi - lib_phase_wrapped + pi, 2*pi) - pi;
        [~, kmin] = min(abs(d));
        shape_id_map(idx) = lib.shape_id(kmin);
        p1_map(idx)       = max(lib.p1(kmin) - cfg.factor/2, 0);
        p2_map(idx)       = max(lib.p2(kmin) - cfg.factor/2, 0);
        phase_err_map(idx)= d(kmin);
        T_map(idx)        = lib.amp(kmin)^2;
        complex_map(idx)  = lib.amp(kmin) * exp(1i*lib.phase(kmin));
    end

    %% --------------------- 회절 시뮬레이션 --------------------------------
    u0 = complex_map .* exp(1i*k0*sind(cfg.inc_angle));
    u2 = RS_2D_Diffraction_FW(cfg.D, cfg.P, cfg.lam0, f, u0);

    center_idx = round(length(u2)/2);
    half_view  = round(cfg.view_um/2);
    s_i = max(1, center_idx - half_view);
    e_i = min(length(u2), center_idx + half_view);
    I_o = abs(u2(s_i:e_i, s_i:e_i)).^2;
    x_view = x_range(s_i:e_i);  y_view = y_range(s_i:e_i);

    % FWHM & 3×FWHM circle energy
    [~, maxIdx] = max(I_o(:));
    [row_c, col_c] = ind2sub(size(I_o), maxIdx);
    profile = I_o(row_c, :); profile = profile / max(profile);
    x_um = x_view/um;
    ab = profile >= 0.5; ed = diff(ab);
    ri = find(ed==1,1,'first'); fi = find(ed==-1,1,'last');
    if ~isempty(ri) && ~isempty(fi)
        x1 = interp1(profile(ri:ri+1), x_um(ri:ri+1), 0.5);
        x2 = interp1(profile(fi:fi+1), x_um(fi:fi+1), 0.5);
        FWHM = abs(x2 - x1);
    else
        FWHM = NaN;
    end

    [XX_um, YY_um] = meshgrid(x_view/um, y_view/um);
    x0 = XX_um(row_c,col_c); y0 = YY_um(row_c,col_c);
    R_um = sqrt((XX_um-x0).^2 + (YY_um-y0).^2);
    dA = (cfg.P/um)^2;
    mask3f = R_um <= 3*FWHM;
    E_focus = sum(I_o(mask3f),'all') * dA;
    E_inc   = sum(mask_aperture(:)) * (cfg.P/um)^2;
    E_focus_pct = E_focus / E_inc * 100;

    fprintf('  FWHM         = %.3f um\n', FWHM);
    fprintf('  focus energy = %.2f %% (3xFWHM / incident)\n', E_focus_pct);

    %% --------------------- 플롯 -----------------------------------------
    if cfg.plot
        figure('Position',[80 80 1300 800],'Color','w');
        subplot(2,3,1);
        imagesc(x_range/um, y_range/um, phase_map); axis image xy;
        title('target phase'); colorbar; caxis([-pi pi]);

        subplot(2,3,2);
        imagesc(x_range/um, y_range/um, phase_err_map); axis image xy;
        title('phase error'); colorbar; caxis([-pi/8 pi/8]);

        subplot(2,3,3);
        imagesc(x_range/um, y_range/um, T_map); axis image xy;
        title('transmission |t|^2'); colorbar; caxis([0 1]);

        subplot(2,3,4);
        imagesc(x_range/um, y_range/um, shape_id_map); axis image xy;
        title('shape (1=cyl,2=sq,3=cross,4=xcross)');
        colormap(gca,[1 1 1; 0.1 0.4 0.9; 0.9 0.4 0.1; 0.2 0.7 0.2; 0.6 0.3 0.8]);
        caxis([0 4]); colorbar;

        subplot(2,3,5);
        imagesc(x_view/um, y_view/um, I_o); axis image xy;
        colormap(gca,'hot'); colorbar;
        title(sprintf('|E|^2 @ focus, FWHM=%.2fum, \\eta_{3F}=%.1f%%', ...
            FWHM, E_focus_pct));

        subplot(2,3,6);
        plot(x_um, profile, 'k-', 'LineWidth', 1.8); grid on; box on;
        xlabel('x (\mum)'); ylabel('norm. intensity');
        title('line profile');
    end

    %% --------------------- GDS 출력 --------------------------------------
    gds_path = write_gds(XX, YY, shape_id_map, p1_map, p2_map, ...
        mask_aperture, cfg);
    fprintf('  GDS saved: %s\n', gds_path);

    %% --------------------- 리턴 ------------------------------------------
    result = struct();
    result.phase_map     = phase_map;
    result.phase_err_map = phase_err_map;
    result.T_map         = T_map;
    result.shape_id_map  = shape_id_map;
    result.p1_map        = p1_map;
    result.p2_map        = p2_map;
    result.f             = f;
    result.FWHM_um       = FWHM;
    result.E_focus_pct   = E_focus_pct;
    result.gds_path      = gds_path;
    result.cfg           = cfg;
end

% =====================================================================
function lib = load_lib_txt(path)
    M = readmatrix(path, 'CommentStyle','%');
    M = M(all(~isnan(M),2), :);
    lib = struct();
    lib.r_eff    = M(:,1);
    lib.phase    = M(:,2);
    lib.amp      = M(:,3);
    lib.shape_id = M(:,4);
    lib.p1       = M(:,5);
    lib.p2       = M(:,6);
end

% =====================================================================
function gds_path = write_gds(XX, YY, sid, p1, p2, mask, cfg)
    gds_path = [cfg.gds_name '.txt'];   % GDSII 는 바이너리라 여기서는 텍스트 포맷 사용
    fid = fopen(gds_path, 'w');

    fprintf(fid, 'HEADER 3;\r\n');
    fprintf(fid, 'BGNLIB;\r\n');
    fprintf(fid, ['LIBNAME ', cfg.gds_name, ';\r\n']);
    fprintf(fid, 'UNITS 1.000000e+000 1.000000e-009;\r\n');
    fprintf(fid, 'BGNSTR;\r\n');
    fprintf(fid, ['STRNAME ', cfg.gds_name, ';\r\n']);

    idxs = find(mask(:));
    for k = 1:numel(idxs)
        ii = idxs(k);
        s = sid(ii);
        if s == 0, continue; end
        cx = XX(ii); cy = YY(ii);
        switch s
            case 1  % cylinder
                r = p1(ii); if r<=0, continue; end
                th = linspace(0, 2*pi, cfg.nc_poly+1); th(end) = [];
                xy = [cx + r*cos(th); cy + r*sin(th)];
                write_poly(fid, xy, cfg.layers(1));
            case 2  % square
                sz = p1(ii); if sz<=0, continue; end
                hx = sz/2;
                xy = [cx-hx cx+hx cx+hx cx-hx;
                      cy-hx cy-hx cy+hx cy+hx];
                write_poly(fid, xy, cfg.layers(2));
            case 3  % cross (plus-shape, 12 corners)
                L = p1(ii); W = p2(ii);
                if L<=0 || W<=0, continue; end
                hL = L/2; hW = W/2;
                % 시계방향 12점 (윗팔 우측상단에서 시작)
                xy = [hW hW hL hL hW hW -hW -hW -hL -hL -hW -hW;
                      hL hW hW -hW -hW -hL -hL -hW -hW hW hW hL];
                xy = xy + [cx; cy];
                write_poly(fid, xy, cfg.layers(3));
            case 4  % xcross (×, 45도 회전된 12점 폴리곤)
                L = p1(ii); W = p2(ii);
                if L<=0 || W<=0, continue; end
                hL = L/2; hW = W/2;
                xy0 = [hW hW hL hL hW hW -hW -hW -hL -hL -hW -hW;
                       hL hW hW -hW -hW -hL -hL -hW -hW hW hW hL];
                R45 = [cosd(45) -sind(45); sind(45) cosd(45)];
                xy = R45 * xy0 + [cx; cy];
                write_poly(fid, xy, cfg.layers(4));
        end
    end

    fprintf(fid, 'ENDSTR\r\n');
    fprintf(fid, 'ENDLIB\r\n');
    fclose(fid);
end

function write_poly(fid, xy_m, layer)
    xy_nm = round(xy_m * 1e9);
    fprintf(fid, 'BOUNDARY\r\n');
    fprintf(fid, 'LAYER %d;\r\n', layer);
    fprintf(fid, 'DATATYPE %d;\r\n', layer);
    fprintf(fid, 'XY\r\n');
    fprintf(fid, '%d\t:\t%d\r\n', xy_nm);
    fprintf(fid, 'ENDEL\r\n');
end

% =====================================================================
function cfg = set_defaults(cfg)
    nm = 1e-9; um = 1e-6;
    cfg = sfe(cfg, 'lam0', 266*nm);
    cfg = sfe(cfg, 'D',    114*um);
    cfg = sfe(cfg, 'P',    260*nm);
    cfg = sfe(cfg, 'NA',   0.6);
    cfg = sfe(cfg, 'ni',   1);
    cfg = sfe(cfg, 'lib_txt',   sprintf('lib_%dnm_opt.txt', round(cfg.lam0*1e9)));
    cfg = sfe(cfg, 'gds_name',  sprintf('metalens_%dnm', round(cfg.lam0*1e9)));
    cfg = sfe(cfg, 'inc_angle', 0);
    cfg = sfe(cfg, 'nc_poly',   24);
    cfg = sfe(cfg, 'layers',    [10, 20, 30, 40]);
    cfg = sfe(cfg, 'plot',      true);
    cfg = sfe(cfg, 'view_um',   60);
    cfg = sfe(cfg, 'factor',    0);
end

function cfg = sfe(cfg, f, v)
    if ~isfield(cfg,f) || isempty(cfg.(f)), cfg.(f) = v; end
end
