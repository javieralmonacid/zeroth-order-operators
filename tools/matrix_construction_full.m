function H = matrix_construction_full(PM,BM,r)
% MATRIX_CONSTRUCTION_FULL Matrix for eigenvalue problem.
%	MATRIX_CONSTRUCTION_FULL(PM,BM,r) returns the discretization
%		matrix corresponding to the operator
%				P(D) - r*beta(x)
%		where \Delta is the standard Laplacian in the 2-torus.
%		
%	Input arguments:
%		PM (2D array): Evaluation of the symbol p(\xi) in the
%			2D wavenumbers, e.g., PM = p(KX,KY).
%		BM (2D array): Evaluation of the function beta in the
%			2D X, Y coordinates, i.e., BM = beta(X,Y).
%		r (float): roughness parameter.
%
%	Returns
%		H (2D matrix, full): If PM and BM are of size N-by-N, then
%			H is N^2-by-N^2 full matrix.
%
%	See also: eigensolver
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/18 (v1.0)
%
N = size(PM,2);

% Diagonal parts
PM = PM(:);
BM = BM(:);

% Final matrix
H = fft(eye(N));
H = kron(H,H);
H = diag(PM)-r*H*diag(BM)*H'/N^2;
