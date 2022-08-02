function [U,X,Y] = fourier2real(UF)
% FOURIER2REAL Convert a grid function in Fourier space into Real space.
%	FOURIER2REAL(UF) returns a time-dependent grid function in Real space.
%
%	Input arguments:
%		UF (3D array): time-dependent grid function in Fourier space.
%			The time is stored in the third coordinate.
%
%	Returns:
%		U (3D array): time-dependent grid function in Real space.
%			At the k-th time, the code computes ifft2(UF(:,:,k)).
%		X, Y (2D array): X, Y coordinates.
%
%	See also: RK4FFT2 ETDRK4FFT2 EIGENSOLVER
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/13 (v1.0)

N = size(UF,1); % Get grid points in each direction
x = 2*pi*(-N/2:N/2-1)'/N;
[X,Y] = meshgrid(x,x);
U = zeros(size(UF));

L = size(UF,3);
for k = 1:L
    U(:,:,k) = ifft2(UF(:,:,k)); 	% IFFT2 at each time
end
