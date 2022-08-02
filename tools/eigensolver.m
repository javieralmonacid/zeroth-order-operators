function [lambda,VF,KX,KY] = eigensolver(r,beta,nu,N,eigvals)
% EIGENSOLVER Eigenvalue solver
%   EIGENSOLVER(r,beta,nu,N,eigvals) computes the first 'eigvals'
%   of the operator 
%
%           <D>^{-1} D_y + 1i*nu*\Delta - r*beta 
%
%   using a pseudo-spectral method by means of a DFT with N modes per 
%   direction. If eigvals << N^2, then the algorithm the built-in function
%   "eigs". If eigvals is not provided, or eigvals = N^2, then it uses "eig".
%
%   Input arguments:
%       r (float): roughness parameter,
%       beta (function handle): beta function in the defintion of 
%           the operator.
%       nu (float, nonegative): viscosity. 
%       N (int): discretization points on the periodic interval [-pi,pi].
%       eigvals (int, <= N^2): amount of eigenvalues to compute.
%
%   Returns:
%       lambda (eigvals-by-1 array): eigenvalues of the operator.
%           They are ordered in ascending order of magnitude, with
%           eigenvalues with positive real part appearing first,
%           then those with negative real part. If nu = 0, the eigenvalues
%           are real, so they are ordered in ascending order.
%       VF (N-by-N-by-eigvals array): eigenfunctions in Fourier space
%           (unshifted). VF(:,:,k) is the eigenfunction corresponding to the
%           k-th eigenvalue.
%       KX, KY (N-by-N array): 2D wavenumbers in the default order by Matlab.
%
%   See also ETDRK4FFT2 RK4FFT2 MAKEITPERIODIC EIGS
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/18 (v1.1)

% Matrix construction
P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2) - 1i*nu*(kx.^2+ky.^2);
x = 2*pi*(-N/2:N/2-1)'/N; y = x; 
[X,Y] = meshgrid(x,y);
BM = beta(X,Y);

ky = [0:N/2-1 -N/2 -N/2+1:-1]'; kx = ky; 
[KX,KY] = meshgrid(kx,ky);
PM = P(KX,KY);

F2 = matrix_construction_full(PM,BM,r);

% Eigenvalue solver
if nargin > 4
    F2 = sparse(F2); % To be able to use 'eigs' instead of 'eig'
    try
        [VF,lambda] = eigs(F2,[],eigvals,0,'Display',1); 
    catch ME
        if strcmp(ME.identifier,'MATLAB:eigs:SingularA')
            % 0 is a eigenvalue, replacing initial shift by nearby value
            warning('Matrix is singular. Replacing shift by 1e-10i')
            [VF,lambda] = eigs(F2,[],eigvals,1e-10i); 
        end 
    end
else
    % If eigvals is not provided, compute all using 'eig'
    warning(['Computing ALL ',num2str(N^2),' eigenvalues for nu = ',num2str(nu)])
    [VF,lambda] = eig(F2);
end
lambda = diag(lambda);

% Sort eigenfunctions
if nu ~= 0
    % OLD ORDER
    % Sorted in increasing order of magnitude then phase
    %[lambda,ii] = sort(lambda);
    % NEW ORDER
    % Positive real part first, then negative real part
    pos = find(real(lambda)>=0);
    neg = setdiff(1:length(lambda),pos)';
    ii = [pos; neg];
    lambda = lambda(ii);
else
    % Spectrum appears to be continuous with no imaginary part.
    [lambda,ii] = sort(lambda,'ComparisonMethod','real'); 
end
VF = VF(:,ii);

% Reshape into 2D arrays and IFFT
VF = reshape(VF,N,N,length(lambda));