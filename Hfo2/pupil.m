function circ = pupil(M, N, radius)
%make aperture function
% M = size of matrix
% N = size of aperture matrix
% radius = size of aperture

r = radius;
dx = 2*radius/N;
total_size = M*dx;

xx = linspace(-total_size/2+dx/2, total_size/2-dx/2, M);
yy = linspace(-total_size/2+dx/2, total_size/2-dx/2, M);
[X, Y] = meshgrid(xx, yy);

circ = sqrt(X.^2 + Y.^2);

circ(circ>r) = 0;
circ(circ>0) = 1;
end

