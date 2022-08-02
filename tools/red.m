function [R,G] = red(s,UF,KX,KY)
% RED s-level Radial Energy Density (RED)
%	RED(s,UF,KX,KY) returns the decay of the shell-summed Fourier coefficients.
%	
%	Input arguments:
%		s (float): level of the RED
%		UF (3D array): time-dependent grid function at times given in the 
%			third coordinate.
%		KX, KY (2D array): 2D wavenumbers in the default order by Matlab.
%
%	Returns:
%		G (2D array): Radial energy density at each time. If UF is of size
%			N-by-N-by-M, then G is an M-by-(N/4) array, with each row
%			correspoding to a different time, and each column to a different R.
%		R (1D array): Radius for the RED between R and N/2. The annulus has width
%			equal to 2, so R is a 1-by-N/4 array.
%	
%	See also RK4FFT2 ETDRK4FFT2
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/13 (v1.1)


for l =1:size(UF,3)
    UF(:,:,l) = fftshift(UF(:,:,l)); 	% Shift everything so wavenumber (0,0) 
    									% is at the center of the spectrum.
end
KX = fftshift(KX); KY = fftshift(KY);

N = size(UF,1);
KX = KX(:); KY = KY(:);
m = sqrt(KX.^2 + KY.^2);
R = 2:2:N/2;
G = zeros(size(UF,3),N/4);

for l = 1:length(R)
    jin = find( (R(l)-2 <= m) & (m < R(l)) ); % Find all values R-2<= m <R
    pos = [KX(jin),KY(jin)]+N/2+1;  % Find associated wavenumbers, 
    								% shift to coincide with Matlab indexing
    elem = sub2ind([N N],pos(:,2),pos(:,1));
    for row = 1:size(UF,3)
        Uaux = UF(:,:,row);
        G(row,l) = sum( (1+KX(jin).^2+KY(jin).^2).^s.*abs(Uaux(elem)).^2 )/N^2;
    end
end

