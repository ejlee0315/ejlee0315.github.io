function result = RS_1D_Diffraction(L, P, lambda, z_, u0)
% L = Diameter
% P = Sampling period
% z_d = Size of Window in the z-direction
% u0 = Input 1D complex amplitude
% lambda = wavelength
% U_z = Output 1D complex amplitude
% h = impulse response function
% H = Transfer function

%% Initial Value ----------------------------------------
k0 = 2*pi/lambda; % wavenumber
m = L/P; % number of pixel
dx = P; % period of x

xx = linspace(-L/2+P/2,L/2-P/2,m); % Coordinates x

X = xx(1)-xx(end:-1:2);
XX = xx(1:end) - xx(1);

X = [X XX];
R = sqrt(X.^2 + z_^2);
N = length(u0);

%% 1D Rayleigh-Sommerfeld Diffraction -------------------
u = zeros(1, 2*N-1);
u(1,1:N) = u0; 

hk = besselh(1,1,k0*R); % 1D Hankel function
h = 1j/2*k0*z_./R .* hk; % impulse response function in 1D RS diffraction

% U = fftshift(fft(fftshift(u)));  % Does not work when uses fftshift and ifftshift
% H = fftshift(fft(fftshift(h)));
% S = ifftshift(ifft(ifftshift(U.*H)))*dx;

U = fft(u);
H = fft(h);
S = ifft(U.*H)*dx;

U_z = S(N:end); % 1D complex amplitude at the z_d

result = U_z;
end

% Reference
% Shen, F. & Wang, A. Fast-Fourier-transform based numerical integration method for the Rayleigh-Sommerfeld diffraction formula. Appl. Opt. 45, 1102–1110 (2006).
% made by Seong-Won Moon
