function snorm = hsnorm(UF,KX,KY,s)
% HSNORM H^s Sobolev norm.
%	HSNORM(UF,KX,KY,s) returns the H^s discrete Sobolev norm
%	of UF.
%
%	Input arguments:
%		UF (2D array): Grid function in Fourier space.
%		KX, KY (2D array): 2D wavenumbers in the default order by Matlab.
%		s (float): s power in the definition of the norm.
%
%	Returns:
%		snorm (float): H^s norm of UF.
%
%	See also RED NORM
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/13 (v1.0)

N = size(UF,2);
UF = UF(:); KX = KX(:); KY = KY(:);

snorm = (1+KX.^2+KY.^2).^s.*abs(UF).^2;
snorm = (4*pi^2/N^4)*sum(snorm);
snorm = sqrt(snorm);
