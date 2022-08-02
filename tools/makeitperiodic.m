function [UU,XX,YY] = makeitperiodic(U,X,Y)
% MAKEITPERIODIC Patch a grid function to make it visibly periodic.
%	MAKEITPERIODIC(U,X,Y) adds an extra column and row to U,X,Y,
%		to make U truly periodic. If U is in Fourier space, then this 
%		function shifts the (0,0) wavenumberto the center of the spectrum. 
%		This is useful, in particular, when plotting in Real or Fourier
%		space.
%
%	Input arguments:
%		U (3D array): time-dependent grid function. The time is stored in
%			the third coordinate.
%		X, Y (2D array): If U is in Real space, then these are the X, Y
%			coordinates. If U is in Fourier space, then these are the KX, KY
%			2D wavenumbers in the default order by Matlab.
%		
%	Returns:
%		UU (3D array): If U is of size N-by-N-by-M, then UU is of size
%			(N+1)-by-(N+1)-by-M containing the missing entries that make
%			U a visibly periodic function. If U is in Fourier space, then UU is
%			shifted so that the wavenumber (0,0) appears in the center of the
%			spectrum.
%		XX, YY (2D array): If X, Y are of size N-by-N, then XX, YY are of size
%			(N+1)-by-(N+1). If X, Y are coordinates in real space, then XX, YY
%			is a meshgrid-type pair that contains, in addition to the components
%			of X and Y, the coordinate x=pi and y=pi. If X, Y are integer
%			wavenumbers then XX, YY contains also the wavenumber kx = N/2 and
%			ky = N/2.
%			
%	See also RK4FFT2 ETDRK4FFT2 EIGENSOLVER
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/12 (v1.0)

N = size(X,1);

if X(1,1) == 0  % Input is in Fourier space
    X = fftshift(X);    Y = fftshift(Y);
    for l = 1:size(U,3)
        U(:,:,l) = fftshift(U(:,:,l));
    end
end

c = -X(1,1);
XX = [X, c*ones(N,1); X(1,:), c];
YY = [Y, Y(:,1); c*ones(1,N), c];
if length(size(U)) == 3
    UU = zeros(size(U)+[1 1 0]);
    for l = 1:size(U,3)
        UU(:,:,l) = [U(:,:,l), U(:,1,l); U(1,:,l), U(1,1,l)];
    end
else
    UU = [U(:,:), U(:,1); U(1,:), U(1,1)];
end
