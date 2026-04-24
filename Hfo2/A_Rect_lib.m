function lib = A_Rect_lib(cfg)
% A_Rect_lib  원기둥/사각기둥/십자가/45°십자가 4종 구조의 라이브러리 스윕
%
%   lib = A_Rect_lib(cfg)
%
%   cfg 필드 (없으면 기본값 사용):
%     .lam0        [m]     설계 파장 (예: 266e-9, 320e-9)
%     .P           [m]     주기
%     .h           [m]     필러 높이
%     .nSiO2               기판 굴절률
%     .n_air               상부 매질
%     .addpath_fmm         fmm 라이브러리 경로 (string)
%     .xlsx_path           HfO2 엑셀 경로 (기본값 로컬 파일)
%     .nx, .ny             fmm harmonics 수 (기본 6)
%     .eta                 fmm eta (기본 45)
%     .save_mat            저장할 .mat 파일 경로 (기본: lib_<lam>nm.mat)
%
%     .cyl.r        [m]    원기둥 반경 sweep (벡터)
%     .sq.s         [m]    사각기둥 한 변 sweep
%     .cross.L      [m]    십자가 외곽 길이 sweep
%     .cross.W      [m]    십자가 팔 너비 sweep
%     .xcross.L     [m]    45°회전 십자가 외곽 길이 sweep
%     .xcross.W     [m]    45°회전 십자가 팔 너비 sweep
%
%   OUTPUT (lib 구조체)
%     .shape : 'cyl' | 'sq' | 'cross' | 'xcross' (각 엔트리)
%     .p1, .p2 : 파라미터 (cyl: r, p2=NaN / sq: s, p2=NaN / cross, xcross: L, W)
%     .t  : 복소 투과계수
%     .T  : 전력 투과율
%     .phi: 위상 [-pi, pi]
%     .cfg: 사용된 cfg 복사본
%
%   또한 .mat 파일에 lib 저장.
%
%   NOTE: 이 함수는 fmm (Fourier Modal Method) 객체가 PATH에 있어야 작동합니다.
%         실제 시뮬레이션은 fmm 프로젝트(G:\내 드라이브\Reserach_Source\A. Simulation\fmm)가 필요.

    if nargin < 1, cfg = struct(); end
    cfg = set_defaults(cfg);

    if ~isempty(cfg.addpath_fmm) && exist(cfg.addpath_fmm, 'dir')
        addpath(cfg.addpath_fmm);
    end

    [n_mat, k_mat, nc] = HfO2_load_nk(cfg.lam0, cfg.xlsx_path);
    fprintf('[A_Rect_lib] lam0=%.1f nm  ->  n=%.4f, k=%.4f  (nc=%.4f%+.4fi)\n', ...
        cfg.lam0*1e9, n_mat, k_mat, real(nc), imag(nc));
    cfg.n_mat = n_mat; cfg.k_mat = k_mat; cfg.nc = nc;

    % 준비: 네 타입의 파라미터 그리드
    jobs = {};
    for i = 1:numel(cfg.cyl.r)
        jobs{end+1} = struct('shape','cyl','p1',cfg.cyl.r(i),'p2',NaN); %#ok<AGROW>
    end
    for i = 1:numel(cfg.sq.s)
        jobs{end+1} = struct('shape','sq','p1',cfg.sq.s(i),'p2',NaN); %#ok<AGROW>
    end
    n_cross = 0;
    for i = 1:numel(cfg.cross.L)
        for j = 1:numel(cfg.cross.W)
            if cfg.cross.W(j) >= cfg.cross.L(i), continue; end  % 팔너비 < 전체길이
            jobs{end+1} = struct('shape','cross','p1',cfg.cross.L(i),'p2',cfg.cross.W(j)); %#ok<AGROW>
            n_cross = n_cross + 1;
        end
    end
    % 45° 회전 십자가 (x-cross): 회전 후 bbox = (L+W)/sqrt(2) <= (P - gap_min)
    %                             → (L+W) <= (P - gap_min)*sqrt(2)
    %   gap_min = 60 nm (최소 인접 pillar 간격)
    P_diag_limit = (cfg.P - 60e-9) * sqrt(2);
    n_xcross = 0;
    for i = 1:numel(cfg.xcross.L)
        for j = 1:numel(cfg.xcross.W)
            L = cfg.xcross.L(i); W = cfg.xcross.W(j);
            if W >= L, continue; end
            if (L + W) > P_diag_limit, continue; end  % 회전 후 cell 밖으로 삐져나감 방지
            jobs{end+1} = struct('shape','xcross','p1',L,'p2',W); %#ok<AGROW>
            n_xcross = n_xcross + 1;
        end
    end
    N = numel(jobs);
    fprintf('[A_Rect_lib] 총 %d 개 구조 스윕 (cyl=%d, sq=%d, cross=%d, xcross=%d)\n', ...
        N, numel(cfg.cyl.r), numel(cfg.sq.s), n_cross, n_xcross);

    shape = cell(N,1); p1 = zeros(N,1); p2 = zeros(N,1);
    t = zeros(N,1); T = zeros(N,1); phi = zeros(N,1);

    t0_start = tic;
    for idx = 1:N
        job = jobs{idx};
        [t(idx), T(idx)] = run_fmm(cfg, job);
        phi(idx) = angle(t(idx));
        shape{idx} = job.shape;
        p1(idx) = job.p1;
        p2(idx) = job.p2;
        if mod(idx, max(1,floor(N/20))) == 0
            fprintf('  [%3d/%3d] %s p1=%.1fnm p2=%.1fnm  T=%.3f phi=%+.2f rad\n', ...
                idx, N, job.shape, job.p1*1e9, job.p2*1e9, T(idx), phi(idx));
        end
    end
    fprintf('[A_Rect_lib] 완료. 경과 %.1f s\n', toc(t0_start));

    lib = struct('shape',{shape},'p1',p1,'p2',p2,'t',t,'T',T,'phi',phi,'cfg',cfg);

    if ~isempty(cfg.save_mat)
        save(cfg.save_mat, 'lib');
        fprintf('[A_Rect_lib] 저장: %s\n', cfg.save_mat);
    end
end

% =====================================================================
function [t, T] = run_fmm(cfg, job)
% 단일 구조에 대해 fmm 실행하여 (t, T) 반환
    c = fmm;
    c.setopt('verbose', false, 'nvm', true, 'basis', 'lr');
    c.set('lam0', cfg.lam0, 'ax', cfg.P, 'ay', cfg.P, ...
          'nx', cfg.nx, 'ny', cfg.ny, 'eta', cfg.eta, ...
          'n1', cfg.nSiO2, 'n2', cfg.n_air);

    % 얇은 레진 바닥층 (원본 코드와 동일한 관행). 없애고 싶으면 제거 가능.
    if cfg.base_layer_thickness > 0
        c.add('multiptc','d', cfg.base_layer_thickness, 'nh', 1, ...
              'rect', {'n', cfg.nc, 'x', cfg.P/2, 'y', cfg.P/2, ...
                       'xspan', cfg.P, 'yspan', cfg.P, 'nvm', true});
    end

    switch job.shape
        case 'cyl'
            r = job.p1;
            if r > 0
                c.add('multiptc','d', cfg.h, 'nh', 1, ...
                      'circ', {'n', cfg.nc, 'x', cfg.P/2, 'y', cfg.P/2, ...
                               'radius', r, 'nvm', true});
            end
        case 'sq'
            s = job.p1;
            if s > 0
                c.add('multiptc','d', cfg.h, 'nh', 1, ...
                      'rect', {'n', cfg.nc, 'x', cfg.P/2, 'y', cfg.P/2, ...
                               'xspan', s, 'yspan', s, 'nvm', true});
            end
        case 'cross'
            L = job.p1; W = job.p2;
            if L > 0 && W > 0
                c.add('multiptc','d', cfg.h, 'nh', 1, ...
                      'rect', {'n', cfg.nc, 'x', cfg.P/2, 'y', cfg.P/2, ...
                               'xspan', L, 'yspan', W, 'nvm', true}, ...
                      'rect', {'n', cfg.nc, 'x', cfg.P/2, 'y', cfg.P/2, ...
                               'xspan', W, 'yspan', L, 'nvm', true});
            end
        case 'xcross'
            % 45°돌아간 십자가 (×) — 두 rect 를 각각 45° 회전
            % fmm 이 rect 에 대해 'phi'(deg) 회전 파라미터를 지원한다고 가정.
            % 만약 fmm 의 인자 이름이 'angle' 또는 'theta' 이면 아래 'phi' 를 교체.
            L = job.p1; W = job.p2;
            if L > 0 && W > 0
                c.add('multiptc','d', cfg.h, 'nh', 1, ...
                      'rect', {'n', cfg.nc, 'x', cfg.P/2, 'y', cfg.P/2, ...
                               'xspan', L, 'yspan', W, 'phi', 45, 'nvm', true}, ...
                      'rect', {'n', cfg.nc, 'x', cfg.P/2, 'y', cfg.P/2, ...
                               'xspan', W, 'yspan', L, 'phi', 45, 'nvm', true});
            end
        otherwise
            error('Unknown shape: %s', job.shape);
    end

    c.compute;
    t = c.out.tl(cfg.nx+1, cfg.ny+1);
    T = c.out.Tl(cfg.nx+1, cfg.ny+1);
end

% =====================================================================
function cfg = set_defaults(cfg)
    nm = 1e-9;

    cfg = set_if_empty(cfg, 'lam0',   266*nm);
    cfg = set_if_empty(cfg, 'P',      260*nm);
    cfg = set_if_empty(cfg, 'h',      500*nm);
    cfg = set_if_empty(cfg, 'nSiO2',  1.46);
    cfg = set_if_empty(cfg, 'n_air',  1.0);
    cfg = set_if_empty(cfg, 'addpath_fmm', 'G:\내 드라이브\Reserach_Source\A. Simulation\fmm');
    cfg = set_if_empty(cfg, 'xlsx_path', '');
    cfg = set_if_empty(cfg, 'nx', 6);
    cfg = set_if_empty(cfg, 'ny', 6);
    cfg = set_if_empty(cfg, 'eta', 45);
    cfg = set_if_empty(cfg, 'base_layer_thickness', 100*nm);

    if ~isfield(cfg,'cyl') || isempty(cfg.cyl)
        cfg.cyl = struct();
    end
    if ~isfield(cfg.cyl, 'r') || isempty(cfg.cyl.r)
        cfg.cyl.r = linspace(25*nm, (cfg.P/2 - 10*nm), 40);
    end

    if ~isfield(cfg,'sq') || isempty(cfg.sq)
        cfg.sq = struct();
    end
    if ~isfield(cfg.sq, 's') || isempty(cfg.sq.s)
        cfg.sq.s = linspace(50*nm, (cfg.P - 20*nm), 40);
    end

    if ~isfield(cfg,'cross') || isempty(cfg.cross)
        cfg.cross = struct();
    end
    if ~isfield(cfg.cross, 'L') || isempty(cfg.cross.L)
        cfg.cross.L = linspace(80*nm, (cfg.P - 20*nm), 12);
    end
    if ~isfield(cfg.cross, 'W') || isempty(cfg.cross.W)
        cfg.cross.W = linspace(40*nm, 120*nm, 8);
    end

    if ~isfield(cfg,'xcross') || isempty(cfg.xcross)
        cfg.xcross = struct();
    end
    if ~isfield(cfg.xcross, 'L') || isempty(cfg.xcross.L)
        cfg.xcross.L = cfg.cross.L;      % 기본적으로 cross 와 동일 sweep
    end
    if ~isfield(cfg.xcross, 'W') || isempty(cfg.xcross.W)
        cfg.xcross.W = cfg.cross.W;
    end

    if ~isfield(cfg, 'save_mat') || isempty(cfg.save_mat)
        cfg.save_mat = sprintf('lib_%dnm.mat', round(cfg.lam0*1e9));
    end
end

function cfg = set_if_empty(cfg, f, v)
    if ~isfield(cfg, f) || isempty(cfg.(f))
        cfg.(f) = v;
    end
end
