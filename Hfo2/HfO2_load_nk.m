function [n, k, nc] = HfO2_load_nk(lam_m, xlsx_path)
% HfO2_load_nk  HfO2 레진 n,k 값을 파장(m)에 맞게 보간해서 반환
%
%   [n, k, nc] = HfO2_load_nk(lam_m)
%   [n, k, nc] = HfO2_load_nk(lam_m, xlsx_path)
%
%   INPUT
%     lam_m     : 파장 [m] (예: 266e-9)
%     xlsx_path : HfO2 n,k 엑셀 파일 경로 (기본값: 260317_HfO2레진n,k.xlsx)
%
%   OUTPUT
%     n  : 실수 굴절률
%     k  : 흡수 계수
%     nc : 복소 굴절률 n + 1i*k
%
%   엑셀 포맷: 4행부터 데이터, col1 = lambda(nm), col2 = n, col3 = k

    if nargin < 2 || isempty(xlsx_path)
        here = fileparts(mfilename('fullpath'));
        xlsx_path = fullfile(here, '260317_HfO2레진n,k.xlsx');
    end

    raw = readmatrix(xlsx_path);
    raw = raw(all(~isnan(raw), 2), :);
    lam_nm = raw(:, 1);
    n_vec  = raw(:, 2);
    k_vec  = raw(:, 3);

    lam_query_nm = lam_m * 1e9;
    if lam_query_nm < min(lam_nm) || lam_query_nm > max(lam_nm)
        warning('HfO2_load_nk:OutOfRange', ...
            '요청 파장 %.1f nm 이 데이터 범위 [%.1f, %.1f] nm 를 벗어남 - 외삽', ...
            lam_query_nm, min(lam_nm), max(lam_nm));
    end

    n = interp1(lam_nm, n_vec, lam_query_nm, 'pchip', 'extrap');
    k = interp1(lam_nm, k_vec, lam_query_nm, 'pchip', 'extrap');
    nc = n + 1i * k;
end
