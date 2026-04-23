%% Anti-aliased Lens 
clc ; close all ; clear ; 
% It is connected to the "metadimer_optimization.m" 

%% Parameter 
% Path setting 
% basic paramter 
nm = 1e-9 ; um = 1e-6 ; mm = 1e-3;
lam0 = 532*nm ; k0 = 2*pi/lam0 ; 
ni = 1;
% lens paramter 
P = 380*nm ; 
ratio = lam0/P ; 
D = 114*um ; NA = 0.6 ; kt_max = k0*NA ; %d/2=odd여야 빵꾸안남
m=D/P;
theta = asin(NA/ni) ; 
f = D/2/tan(theta);
DOF = lam0 / (1 - cos(asin(NA)))/um; 

if P > lam0/2/NA
    warning('Period should be smaller than %.2f nm', lam0/2/NA/nm); % Nyquist sampling criterion
end
fprintf('Diameter = %.2f mm\n',D/mm);
fprintf('Focal length = %.5f mm\n',f/mm);
fprintf('NA = %.2f \n',NA);
fprintf('maximum NA = %.2f \n', lam0/2/P);
fprintf('Depth of focus = %.2f um\n',DOF);
%% 

%%%%%%%%%%%%%%%%%%%%%%%%% Grid and Phase map %%%%%%%%%%%%%%%%%%%%
% Hexagonal grid 
ps = 0;
% x_range = -D/2:P:D/2 ; 
x_range = -D/2+P/2:P:D/2-P/2 ;
% y_range = -D/2:P:D/2 ; 
y_range = -D/2+P/2:P:D/2-P/2 ;
[X_hexa,Y_hexa] = meshgrid(x_range,y_range) ; 
R_hexa = sqrt(X_hexa.^2 + Y_hexa.^2) ; 
% phase_map_hexa = mod(ps + ni*k0*(f-sqrt(R_hexa.^2 + f.^2)),2*pi) - pi ; 
phase_map_hexa = mod(ps - ni*k0*(R_hexa.^2 / (2*f)), 2*pi) - pi;
phase_map_hexa(R_hexa>D/2) = 0 ; 

figure('position',[100 580 500 400]);
imagesc(x_range,y_range,phase_map_hexa);
axis xy;
shading interp;
custom_colormap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1);  % 파란색-흰색
                   ones(128,1) linspace(1,0,128)' linspace(1,0,128)']; % 흰색-빨간색

colormap(custom_colormap);

colorbar;

% %%
% %%%%%errormap%%%%%
% linesize = ceil(D/P+1);
linesize = ceil(D/P);
totalsize = linesize.^2;
data = load('260223_KU_nanoconnect.txt') ; 
h = 630*nm ; 
factor = 0*nm ; % minus 10 nm of structure 

i = 1 ;
radius = zeros(totalsize,1);
atomphase = zeros(totalsize,1);
error = zeros(totalsize,1);
complexlens = zeros(totalsize,1);

for jr = 1:totalsize
    phase = phase_map_hexa(jr) ; 
    [Min,ind] = min(abs(phase - (mod(data(:,2)+pi,2*pi)-pi))) ;
    radius(i) = data(ind,1) - factor/2 ;
    atomphase(i) = data(ind,2);
    error(i) = (phase - atomphase(i));
    complexlens(i) = data(ind,3)*exp(1i*atomphase(i));
    i = i+1 ;
end
% reshape
error_2D_abs = reshape(abs(error), [linesize, linesize]);
error_2D_signed = reshape(error, [linesize, linesize]);

% pupil 함수를 사용한 원형 마스크
circ = pupil(linesize, linesize, linesize/2);  % or D/(2*P) depending on your scale

% 마스크 적용
error_2D_abs_masked = error_2D_abs .* circ;
error_2D_signed_masked = error_2D_signed .* circ;

% error_2D = reshape(abs(error), [linesize, linesize]);
figure;
% imagesc(error_2D);
imagesc(error_2D_abs_masked);
colorbar;
colormap(jet);
title('Phase Error Map');
caxis([0 1]);
xlabel('X Index');
ylabel('Y Index');
axis image;

% error_2D = reshape(error, [linesize, linesize]);
figure;
% imagesc(error_2D);
imagesc(error_2D_signed_masked);
colorbar;
colormap(custom_colormap);
title('Phase Error Map');
caxis([-3.1416 3.1416]);
xlabel('X Index');
ylabel('Y Index');
axis image;
%%%

%% Fresnel Diffraction
circ = pupil(linesize, linesize, D/2);

complex_metalens = reshape(complexlens, [linesize, linesize]);
inc_angle=0;

circcomplexmeta = circ.*complex_metalens;


u0 = circcomplexmeta.*exp(1i*k0*sind(inc_angle));
view = 45;
u2=RS_2D_Diffraction_FW(D, P, lam0, f, u0); %rs\
N= length (u0);
% u2=ASM_2D_FW(D, lam0, f, u0); %asm
center_idx = round(length(u2)/2);  % Ensure the center index is an integer

start_idx = max(1, round(center_idx - view/2));  % Ensure at least 1
end_idx = min(length(u2), round(center_idx + view/2));  % Ensure no index exceeds array length
I_o = abs(u2(start_idx:end_idx, start_idx:end_idx)).^2;

% Extract the corresponding coordinates from xx and yy for plotting
x_view = x_range(start_idx:end_idx);
y_view = y_range(start_idx:end_idx);

% Plotting the result
figure('position',[100 80 1000 800]);
imagesc(x_view/um, y_view/um, I_o);  % Use x_view and y_view for axis scaling
axis xy;
colormap('hot');
colorbar;
title('Intensity');
% pbaspect([1 f/D 1]);
%% ================= 기존 시뮬레이션 2D & 1D 코드 (상단 동일) =================
figure('position',[100 80 1000 800]);

% --- (1) 시뮬레이션 2D Intensity (Zoomed) ---
subplot(2,1,1);
imagesc(x_view/um, y_view/um, I_o);
axis image; axis xy;
colormap(hot);
colorbar;
title('Intensity (Zoomed)','FontWeight','normal');
xlabel('x (\mum)'); ylabel('y (\mum)');

% --- (2) 시뮬레이션 1D 프로파일 (가로) ---
[maxVal, maxIdx] = max(I_o(:));
[row_c, col_c]   = ind2sub(size(I_o), maxIdx);

profile_calc = I_o(row_c, :);
profile_calc = profile_calc / max(profile_calc);
x_axis_um    = x_view / um;

subplot(2,1,2); hold on;
hCalc = plot(x_axis_um, profile_calc, 'k-', 'LineWidth', 1.8, ...
    'DisplayName','Calculated');

% ---- 시뮬레이션 FWHM ----
above_half_c = profile_calc >= 0.5;
edges_c      = diff(above_half_c);
rise_c       = find(edges_c == 1, 1, 'first');
fall_c       = find(edges_c == -1,1,'last');
if ~isempty(rise_c) && ~isempty(fall_c)
    x1c = interp1(profile_calc(rise_c:rise_c+1), x_axis_um(rise_c:rise_c+1), 0.5);
    x2c = interp1(profile_calc(fall_c:fall_c+1), x_axis_um(fall_c:fall_c+1), 0.5);
    FWHM_calc = abs(x2c - x1c);
    % (얇은 회색 선)
    plot([x1c x1c],[0 1],'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1);
    plot([x2c x2c],[0 1],'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1);
    % FWHM 라벨(옵션: 지저분하면 주석)
    text(mean([x1c x2c]), 0.82, sprintf('FWHM_{calc}=%.2f \\mum',FWHM_calc), ...
        'Color',[0.25 0.25 0.25],'HorizontalAlignment','center','FontSize',8);
else
    FWHM_calc = NaN;
end
%% --- (추가) 3×FWHM 원 내부 적분 ----------------------------------------
% ────────────────────────────────────────────────────────────────────────
% 준비:  row_c, col_c          ->> 최대강도 위치 (이미 구함)
%        FWHM_calc (µm)        ->> 계산된 FWHM (이미 구함)
%        I_o                   ->> intensity^2 (이미 구함)
%        x_view, y_view (m)    ->> I_o 가 써먹는 실제 좌표 (이미 구함)
% ────────────────────────────────────────────────────────────────────────

% 1) 실측 좌표(µm) 그리드 만들기
[XX_um, YY_um] = meshgrid(x_view/um, y_view/um);   % µm 단위

% 2) 중심 좌표(µm)
x0 = XX_um(row_c, col_c);
y0 = YY_um(row_c, col_c);

% 3) 각 픽셀까지의 반경 r (µm)
R_um = sqrt( (XX_um - x0).^2 + (YY_um - y0).^2 );

% 4) 원 마스크 (r ≤ 3×FWHM)
Rmax_um   = 3 * FWHM_calc;        % 원하는 반경 [µm]
mask_circ = R_um <= Rmax_um;      % 논리형 마스크

% 5) (선택) 픽셀 면적 dA [µm²]  ───────────────
dx_um = abs( x_view(2) - x_view(1) ) / um;   % x 픽셀 간격
dy_um = abs( y_view(2) - y_view(1) ) / um;   % y 픽셀 간격
dA_um2 = dx_um * dy_um;                      % 한 픽셀 면적

% 6) 적분 결과
E_tot      = sum( I_o(mask_circ), 'all' ) * dA_um2;  % ∫∫ I dA  [정규화 X]
area_um2   = pi * Rmax_um^2;                         % π r²    [이론적 면적]

% ------------------------------------------------------------
% (앞선 3×FWHM 적분 코드를 그대로 두고, E_tot·dA_um2까지 나온 직후에 붙이세요)
% ------------------------------------------------------------

% 0) 창(view) 전체 에너지 ------------- (참고용, 이미 계산했다면 생략 가능)
E_full = sum(I_o(:)) * dA_um2;    

% 1) 메타렌즈 ‘입사’ 파워 (렌즈 개구 전체) -------------------------------
dx_inc_um = P / um;                 % 한 픽셀 폭 = 주기 P = 0.22 µm
dA_inc_um2 = dx_inc_um^2;           % 픽셀 면적    [µm²]

% ① 평면파 진폭 = 1 (phase·transmission 무시) 라면:
E_incident = sum(circ(:)) * dA_inc_um2;          % circ = 개구(0/1)

% ② 메타원소 투과계수(data(ind,3))까지 포함하고 싶다면:
% E_incident = sum(abs(circcomplexmeta(:)).^2) * dA_inc_um2;

% 2) 비율 계산 -----------------------------------------------------------
frac_view     = E_tot   / E_full;      % (선택) 시뮬레이션 창 대비
frac_incident = E_tot   / E_incident;  % 렌즈 전체 입사 대비

% 3) 출력 ---------------------------------------------------------------
fprintf('\n[3×FWHM 원 적분]\n');
fprintf(' - 반경 R           = %.3f µm\n', Rmax_um);
fprintf(' - 창(view) 대비    = %.2f %%\n', frac_view*100);
fprintf(' - 렌즈 입사 대비   = %.2f %%\n\n', frac_incident*100);


%% ================= 실측 TIF 추가 & 스타일 =================
% ---- 사용자 입력 ----
tif_name  = '2_1_sample_1_1_exposure9_0049_14.55487_1.01.tif';
channel   = 1;
Nwin      = 41;                    % 홀수 권장
pixelsize = (27.84/117)*um;        % m/px

% ---- TIF 읽기 & 중심 찾기 ----
B = double(imread(tif_name));
if ndims(B)==3, B = B(:,:,channel); end
[maxB,~] = max(B(:));
[rp_all, cp_all] = find(B==maxB);
col_exp = round(mean(cp_all)); row_exp = round(mean(rp_all));

% ---- Crop ----
if mod(Nwin,2)==0, Nwin=Nwin+1; end
halfW = (Nwin-1)/2;
rmin = max(row_exp-halfW,1); rmax = min(row_exp+halfW,size(B,1));
cmin = max(col_exp-halfW,1); cmax = min(col_exp+halfW,size(B,2));
Crop = B(rmin:rmax, cmin:cmax);

% ---- 정규화 (피크) ----
CropN = Crop / max(Crop(:));

% ---- 1D (중심 행) ----
row_mid = floor(size(CropN,1)/2)+1;
profile_meas = CropN(row_mid,:);
Nexp = length(profile_meas);
x_exp_um = ((1:Nexp) - (Nexp+1)/2) * (pixelsize/um);

% ---- (필요시 x축 매칭) : 보간하여 공통 그리드 사용 (권장) ----
% 공통 x축을 선택 (시뮬레이션 범위 안 or 둘을 합친 범위)
x_common_um = linspace( max(min(x_axis_um), min(x_exp_um)), ...
                        min(max(x_axis_um), max(x_exp_um)), 1000);

prof_calc_interp = interp1(x_axis_um, profile_calc, x_common_um, 'linear','extrap');
prof_meas_interp = interp1(x_exp_um, profile_meas,  x_common_um, 'linear','extrap');

% ---- Measured 플롯 (파란 점선) ----
hMeas = plot(x_common_um, prof_meas_interp / max(prof_meas_interp), ...
    '--','Color',[0 0.38 0.85],'LineWidth',1.8,'DisplayName','Measured');

% ---- 실측 FWHM (원래 좌표 기준) ----
prof_meas_norm = profile_meas / max(profile_meas);
above_half_e = prof_meas_norm >= 0.5;
edges_e      = diff(above_half_e);
rise_e       = find(edges_e==1,1,'first');
fall_e       = find(edges_e==-1,1,'last');

if ~isempty(rise_e) && ~isempty(fall_e)
    x1e = interp1(prof_meas_norm(rise_e:rise_e+1), x_exp_um(rise_e:rise_e+1), 0.5);
    x2e = interp1(prof_meas_norm(fall_e:fall_e+1), x_exp_um(fall_e:fall_e+1), 0.5);
    FWHM_meas = abs(x2e - x1e);
    % FWHM 선 (파란 점선 얇게)
    plot([x1e x1e],[0 1],'Color',[0 0.38 0.85],'LineStyle','--','LineWidth',1);
    plot([x2e x2e],[0 1],'Color',[0 0.38 0.85],'LineStyle','--','LineWidth',1);
    text(mean([x1e x2e]), 0.60, sprintf('FWHM_{meas}=%.2f \\mum',FWHM_meas), ...
        'Color',[0 0.38 0.85],'HorizontalAlignment','center','FontSize',8);
else
    FWHM_meas = NaN;
end

%% ================= 축/레이블/스타일 공통 마감 =================
xlabel('x (\mum)');
ylabel('Normalized intensity');
title('Calculated vs. Measured Line Profile','FontWeight','normal');

% 축 범위: 두 데이터 모두 커버
xmin_all = min([x_axis_um(1), x_exp_um(1)]);
xmax_all = max([x_axis_um(end), x_exp_um(end)]);
xlim([xmin_all xmax_all]);

ylim([0 1.05]);
set(gca,'FontName','Arial','FontSize',9,'LineWidth',1,...
        'TickDir','out','Box','off');

% 범례 정리
legend([hCalc hMeas],'Location','northeast','Box','off');

% FWHM 비교 콘솔
if ~isnan(FWHM_calc) && ~isnan(FWHM_meas)
    fprintf('FWHM(calc)=%.3f µm, FWHM(meas)=%.3f µm, ratio=%.3f\n', ...
        FWHM_calc, FWHM_meas, FWHM_meas/FWHM_calc);
end

% Figure 배경/저장
set(gcf,'Color','w');

% 벡터 저장 (필요시)
% exportgraphics(gcf,'Profile_Comparison.pdf','ContentType','vector');
% print(gcf,'Profile_Comparison','-depsc','-painters');

% %% ===== 2D Intensity + 1D Profile (추가) =====
% % I_o : 이미 start_idx:end_idx 영역에서 잘라낸 intensity^2
% % x_view, y_view : 해당 영역의 좌표 (단위: m)
% 
% figure('position',[100 80 1000 800]);
% 
% % (1) 2D 이미지
% subplot(2,1,1);
% imagesc(x_view/um, y_view/um, I_o);
% axis xy; axis image;
% colormap(hot);
% colorbar;
% title('Intensity (Zoomed)');
% xlabel('x [\mum]');
% ylabel('y [\mum]');
% 
% % (2) 중심 또는 최대 강도 위치 선택
% % 옵션 A: 단순히 가운데
% % row_c = round(size(I_o,1)/2);
% % col_c = round(size(I_o,2)/2);
% 
% % 옵션 B: 실제 최대값 위치 (권장)
% [maxVal, maxIdx] = max(I_o(:)); %#ok<NASGU>
% [row_c, col_c] = ind2sub(size(I_o), maxIdx);
% 
% % 가로 1D 프로파일 (x-방향)
% profile_x = I_o(row_c, :);
% profile_x = profile_x / max(profile_x);   % 정규화
% 
% x_axis_um = x_view/um;   % x 좌표 (µm)
% 
% subplot(2,1,2);
% plot(x_axis_um, profile_x, 'LineWidth', 1.6);
% grid on; box on;
% xlabel('x [\mum]');
% ylabel('Normalized Intensity');
% title('Central Horizontal Line Profile');
% xlim([x_axis_um(1) x_axis_um(end)]);
% 
% % ---------- FWHM 계산 ----------
% above_half = profile_x >= 0.5;
% edges = diff(above_half);
% rise_idx  = find(edges == 1);
% fall_idx  = find(edges == -1);
% 
% if ~isempty(rise_idx) && ~isempty(fall_idx)
%     i1 = rise_idx(1);
%     i2 = fall_idx(end);
% 
%     % 선형 보간 (0.5 지점)
%     x1 = interp1(profile_x(i1:i1+1), x_axis_um(i1:i1+1), 0.5);
%     x2 = interp1(profile_x(i2:i2+1), x_axis_um(i2:i2+1), 0.5);
%     FWHM_um = abs(x2 - x1);
% 
%     hold on;
%     yline(0.5,'r--');
%     plot([x1 x1],[0 1],'r--');
%     plot([x2 x2],[0 1],'r--');
%     plot([x1 x2],[0.5 0.5],'ro-','MarkerFaceColor','r');
%     text(mean([x1 x2]), 0.65, sprintf('FWHM = %.2f \\mum', FWHM_um), ...
%         'HorizontalAlignment','center','Color','r');
% 
%     fprintf('FWHM (x-profile) = %.3f µm\n', FWHM_um);
% else
%     fprintf('FWHM 계산 실패: 프로파일이 0.5를 두 번 교차하지 않습니다.\n');
% end
% 
% %% ===== (추가) 실측 TIF PSF 프로파일 오버레이 =====
% % ---- 사용자 입력 ----
% tif_name  = '2_1_sample_1_1_exposure9_0049_14.55487_1.01.tif';
% channel   = 1;          % 컬러 TIF이면 특정 채널 선택
% Nwin      = 41;         % 한 변 픽셀 수 (홀수 권장, 짝수면 자동 +1)
% pixelsize = (27.84/111)*um;   % (m/px) = 0.251... µm (당신 코드)
% % pixelsize = 0.22*um;  % 다른 값 시험하려면 주석 해제
% 
% % ---- TIF 읽기 ----
% B = double(imread(tif_name));
% if ndims(B) == 3
%     B = B(:,:,channel);
% end
% 
% % ---- 최대값 위치(복수) → 평균 좌표 ----
% [maxB, ~] = max(B(:));
% [rows_peak, cols_peak] = find(B == maxB);
% col_c = round(mean(cols_peak));
% row_c = round(mean(rows_peak));
% 
% % ---- Crop (홀수 보정) ----
% if mod(Nwin,2)==0, Nwin = Nwin + 1; end
% halfW = (Nwin - 1)/2;
% 
% r_min = max(row_c - halfW, 1);
% r_max = min(row_c + halfW, size(B,1));
% c_min = max(col_c - halfW, 1);
% c_max = min(col_c + halfW, size(B,2));
% 
% % 만약 이미지 가장자리라서 크기가 줄어들면 패딩 없이 있는 영역만 사용
% Crop = B(r_min:r_max, c_min:c_max);
% 
% % ---- PSF 정규화 (피크) ----
% CropN = Crop / max(Crop(:));
% 
% % ---- 1D 가로 프로파일 (중심 행) ----
% row_mid = round(size(CropN,1)/2);
% profile_exp = CropN(row_mid, :);
% 
% % ---- 좌표 (µm) ----
% Nexp = length(profile_exp);
% x_exp_um = ((1:Nexp) - (Nexp+1)/2) * (pixelsize/um);  % 중심 0 기준
% 
% % ---- (기존 subplot 유지) ----
% subplot(2,1,2);
% hold on;
% 
% % ---- 실측 프로파일 플로팅 ----
% plot(x_exp_um, profile_exp, 'Color', [0 0.45 0.85], 'LineWidth', 1.5, ...
%     'DisplayName', 'Measured (TIF)');
% 
% % ---- 실측 FWHM 계산 ----
% profE = profile_exp / max(profile_exp);        % (보강) 최대 1
% above_half_e = profE >= 0.5;
% edges_e = diff(above_half_e);
% rise_e  = find(edges_e == 1);
% fall_e  = find(edges_e == -1);
% 
% FWHM_exp_um = NaN;
% if ~isempty(rise_e) && ~isempty(fall_e)
%     j1 = rise_e(1);
%     j2 = fall_e(end);
%     x1e = interp1(profE(j1:j1+1), x_exp_um(j1:j1+1), 0.5);
%     x2e = interp1(profE(j2:j2+1), x_exp_um(j2:j2+1), 0.5);
%     FWHM_exp_um = abs(x2e - x1e);
% 
%     % 표시
%     plot([x1e x1e],[0 1],'b--');
%     plot([x2e x2e],[0 1],'b--');
%     plot([x1e x2e],[0.5 0.5],'bo-','MarkerFaceColor','b');
%     text(mean([x1e x2e]), 0.35, sprintf('Exp FWHM = %.2f \\mum', FWHM_exp_um), ...
%         'HorizontalAlignment','center','Color','b','FontSize',9);
%     fprintf('Experimental FWHM = %.3f µm\n', FWHM_exp_um);
% else
%     fprintf('실측 FWHM 계산 실패 (0.5 교차 부족)\n');
% end
% 
% % ---- 범례 / 축 범위 조정 ----
% % 시뮬레이션 x축 변수: x_axis_um (이미 존재)
% xmin_all = min([x_axis_um(1), x_exp_um(1)]);
% xmax_all = max([x_axis_um(end), x_exp_um(end)]);
% xlim([xmin_all xmax_all]);
% 
% legend({'Simulation','Measured'},'Location','best');
% 
% % ---- (선택) 정규화 방식 통일 팁 ----
% % 만약 두 프로파일의 면적(∑)을 같게 맞추고 싶다면:
% % profile_x = profile_x / trapz(x_axis_um, profile_x);
% % profile_exp = profile_exp / trapz(x_exp_um, profile_exp);
% % 그런 뒤 다시 플롯 (현재는 peak 정규화)
% 
% % ---- (선택) 두 FWHM 비교 콘솔 출력 ----
% if exist('FWHM_um','var') && ~isnan(FWHM_um) && ~isnan(FWHM_exp_um)
%     fprintf('FWHM ratio (Exp / Sim) = %.3f\n', FWHM_exp_um / FWHM_um);
% end

% %% 1D Rayleigh-Sommerfeld Diffraction
% u0_1D = u0(round(length(u0)/2), :);
% % u0_1D = u0_1D(2:end-1); 
% 
% resolv_z = 456;
% z_d = linspace(0,f*2,resolv_z);
% U_1D = [];
% 
% for i = 1:length(z_d)
%     z_ = z_d(i);
%     temp = RS_1D_Diffraction(D, P, lam0, z_, u0_1D);
%     U_1D = [U_1D; temp];
% end
% 
% I_1D = abs(U_1D).^2;

% Plot 1D -----------------------------------------------------------------
% view_x = floor(m/40); % scaling
% view_z = 15;
% X_ = x_range(length(x_range)/2+1-view_x:length(x_range)/2+view_x);
% z_ = z_d(length(z_d)/2+1-view_z:length(z_d)/2+view_z);
% I_1D_z = I_1D(length(I_1D(:,1))/2+1-view_z:length(I_1D(:,1))/2+view_z,length(I_1D(1,:))/2+1-view_x:length(I_1D(1,:))/2+view_x);
% 
% figure;
% % ('position',[600 80 500 400]);
% imagesc(X_,z_, I_1D_z/max(max(I_1D_z)));
% axis xy;
% colormap('hot');
% colorbar;
% caxis([0 1]);
% title('Intensity')
% pbaspect([1 f/D 1]);

%% Metaatom GDS (Hexa)
% path setting

% clear temp xy X Y; clc ; 
name = '260223_KU_nanoconnect_gds';
nc = 10 ; theta = linspace(0,2*pi,nc+1) ; theta = theta(1:nc) ; % number of corners - 1; 
n = nnz(phase_map_hexa) ; % number of metaatoms 
data = load('260223_KU_nanoconnect.txt') ; 
h = 600*nm ; factor = 0*nm ; % minus 10 nm of structure 
idx = find(phase_map_hexa~=0) ; % indices of metaatom

% Find the radius 
i = 1 ;

radius = zeros(n,1);
for jr = 1:n 
    phase = phase_map_hexa(idx(jr)) ; 
    [Min,ind] = min(abs(phase - (mod(data(:,2)+pi,2*pi)-pi))) ;
    radius(i) = data(ind,1) - factor/2 ;
    i = i+1 ;
end

temp = zeros(n,nc,2);
for i = 1:n 
    x_cent = X_hexa(idx(i)) ; 
    y_cent = Y_hexa(idx(i)) ; 
    [X,Y] = pol2cart(theta,radius(i)) ;
    temp(i,:,:) = [X ; Y]' + repmat([x_cent y_cent],[nc,1]); 
end
X = [temp(:,:,1) temp(:,1,1)] ; % ?? ?? ?? 
Y = [temp(:,:,2) temp(:,1,2)] ; 
xy = zeros(nc*2,n);

% xy = zeros(2(nc+1),N); % xy = [corners(x,y),total # of building block]
for i=1:n
    for j=1:nc
        xy(2*(j-1)+1,i) = X(i,j);
        xy(2*j,i) = Y(i,j);
    end
end

xy = xy*1e+9 ; 
xy = round(xy) ;

fid = fopen([pwd ,'\',name,'.txt'],'w');

fprintf(fid,'HEADER 3;\r\n');
fprintf(fid,'BGNLIB;\r\n');
fprintf(fid,['LIBNAME ',name,'.txt;\r\n']);
fprintf(fid,'UNITS 1.000000e+000 1.000000e-009;\r\n');
fprintf(fid,'BGNSTR;\r\n');
fprintf(fid,['STRNAME ',name,';\r\n']);
%%
for i=1:n
    if radius(i) ~= 0
        fprintf(fid,'BOUNDARY\r\n');
        fprintf(fid,'LAYER 46;\r\n');
        fprintf(fid,'DATATYPE 46;\r\n');
        fprintf(fid,'XY\r\n');
        fprintf(fid,'%d\t:\t%d\r\n',xy(:,i)');
        fprintf(fid,'ENDEL\r\n');

    end
end


fprintf(fid,'ENDSTR\r\n');
fprintf(fid,'ENDLIB\r\n');

fclose(fid);
