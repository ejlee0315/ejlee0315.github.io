function result = RS_2D_Diffraction_FW(L, P, lambda, z_, u1)
% L = Diameter
% P = Sampling period
% z_ = Size of Window in the z-direction
% u1 = Input 1D complex amplitude
% lamda = wavelength
% U_z = Output 1D complex amplitude
% h = impulse response function
% H = Transfer function

%% Coordinate ------------------------------------------
k0 = 2*pi/lambda; % wavenumber
N = length(u1); % number of pixel
dx = P; dy = P; % period of x and y direction

xx1 = linspace(-L/2+P/2,L/2-P/2,N); % Coordinates x
yy1 = fliplr(linspace(-L/2+P/2,L/2-P/2,N)); % Coordinates y

xx2 = linspace(-L/2+P/2,L/2-P/2,N); % Coordinates x
yy2 = fliplr(linspace(-L/2+P/2,L/2-P/2,N)); % Coordinates y


X = xx1(1)-xx2(end:-1:2);
XX = xx1(1:end) - xx2(1);
X = [X XX];

Y = yy1(1)-yy2(end:-1:2);
YY = yy1(1:end) - yy2(1);
Y = [Y YY];

[X_, Y_] = meshgrid(X,Y);
X_ = X_';
Y_ = fliplr(Y_');

R = sqrt(X_.^2 + Y_.^2 + z_^2);

%% 2D Rayleigh-Sommerfeld Diffraction -------------------
u = zeros(2*N-1);
u(1:N,1:N) = u1;

h = exp(1j*k0*R).*z_./(2*pi*R.^3).*(1-1j*k0*R); % impulse response function in 2D RS diffraction

U = fft2(u);  % Does not work when uses fftshift and ifftshift
H = fft2(h);

S = ifft2(U.*H)*dx*dy;

u2 = S(N:end,N:end);  % 2D complex amplitude at the z_d

result = u2;
end


% Reference
% Shen, F. & Wang, A. Fast-Fourier-transform based numerical integration method for the Rayleigh-Sommerfeld diffraction formula. Appl. Opt. 45, 1102–1110 (2006).
% made by Seong-Won Moon
