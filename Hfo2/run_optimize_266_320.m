%% run_optimize_266_320.m
%   266 nm, 320 nm 두 파장에 대해 형상(원기둥/사각기둥/십자가) 라이브러리를
%   sweep 하고, (h, P) 외부 grid search 까지 수행하여 best library 를 저장.
%
%   필요사항: fmm 객체 + addpath_fmm 경로
%   결과: ./opt_266nm/, ./opt_320nm/ 폴더에 lib_*.mat, lib_*_opt.txt 저장

clc; close all; clear;
nm = 1e-9;

%% 공통 설정
common = struct();
common.nSiO2 = 1.46;
common.n_air = 1.0;
common.nx = 6; common.ny = 6; common.eta = 45;
common.addpath_fmm = 'G:\내 드라이브\Reserach_Source\A. Simulation\fmm';
common.xlsx_path   = '';   % 같은 폴더의 260317_HfO2레진n,k.xlsx 자동 사용
common.base_layer_thickness = 100*nm;
common.verbose = true;

% (h, P) 후보 — 두 파장 공통
%   P = 190 nm 단일 고정 (min feature 80 nm, min gap 60 nm 공정 규칙)
%   FOV ±40° (quadratic) Nyquist :  P <= lam0 / (2*(NA + sin40°))
%     266nm P=190 : NA ≤ 0.057  (NA=0 극한 FOV 최대 ±44°)
%     320nm P=190 : NA ≤ 0.199
common.h_list  = [400 500 600 700 800]*nm;
common.P_list  = [190]*nm;

% 형상 sweep 해상도
common.nR = 30; common.nS = 30; common.nL = 10; common.nW = 10;

% 위상 bin 최적화 옵션
common.opt = struct( ...
    'N_bin',  64, ...
    'T_min',  0.5, ...
    'prefer', {{'cyl','sq','cross','xcross'}}, ...
    'allow_wrap', true);

% 점수(참고용 로그): best 선정 자체는 (coverage, mean_T) lexicographic 으로
% A_Rect_search 안에서 2단계 정렬됩니다. score_fun 은 summary/log 에만 쓰임.
common.score_fun = @(o) o.coverage * 10 + o.mean_T;

%% (1) 266 nm
cfg266 = common;
cfg266.lam0 = 266*nm;
cfg266.save_dir = 'opt_266nm';
res266 = A_Rect_search(cfg266);

%% (2) 320 nm
cfg320 = common;
cfg320.lam0 = 320*nm;
cfg320.save_dir = 'opt_320nm';
res320 = A_Rect_search(cfg320);

%% 요약
fprintf('\n\n##################  SUMMARY  ##################\n');
fprintf('  266 nm best : h=%dnm  P=%dnm  cov=%.1f%%  T=%.3f  -> %s\n', ...
    round(res266.best.h*1e9), round(res266.best.P*1e9), ...
    res266.best.opt.coverage*100, res266.best.opt.mean_T, res266.best.txt_path);
fprintf('  320 nm best : h=%dnm  P=%dnm  cov=%.1f%%  T=%.3f  -> %s\n', ...
    round(res320.best.h*1e9), round(res320.best.P*1e9), ...
    res320.best.opt.coverage*100, res320.best.opt.mean_T, res320.best.txt_path);
fprintf('################################################\n');

% 두 파장 결과를 한 .mat 에 묶어 저장
save('search_results_266_320.mat', 'res266', 'res320');

%% (3) 266/320 nm 동시 성능 비교 (같은 h,P 에서 best 찾기)
fprintf('\n>>> 두 파장 결합 분석 자동 실행...\n');
cmb = compare_266_320(res266, res320);
