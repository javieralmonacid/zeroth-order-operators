function [UF,KX,KY,tsteps] = etdrk4fft2(P,beta,r,omega0,N,dt,T,fun,u0,snaps)
% ETDRK4FFT2 ETDRK4-FFT2 based transient solver
%   ETDRK4FFT2(P,beta,r,omega0,N,dt,T,fun,u0,snaps) returns the solution
%   to the problem in the 2-torus
%
%       i u_t - (P(x,y,D) + r*beta(x,y))* u = fun(x,y)*exp(-i*omega0*t),
%       u = u0 at initial time.
%
%   The time discretization is based on the Exponential Time-Differencing
%   Runge-Kutta 4-th order (ETDRK4) from Kassam & Trefethen (2005). In turn, 
%   the space pseudo-spectral discretization uses the built-in function FFT2.
%
%   Input arguments:
%       P (function handle): symbol of the pseudo-differential operator P,
%           for example, P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2).
%       beta (function handle): boundary-like smooth periodic function.
%       r (float, nonnegative): parameter related to the roughness of the
%           solution. The higher the value, the more singular is.
%       omega0 (float): forcing frequency.
%       N (int): number of gridpoints in the periodic interval [-pi,pi]
%       dt (float, positive): time step.
%       T (float, positive): final time of simulation.
%       fun (function handle): complex forcing function in real space.
%       u0 (function handle / N-by-N matrix): Initial data.
%       snaps (int, nonzero): snapshots of the evolution between time 0
%           and time T, equally spaced if (T/snaps)/dt is an integer value.
%
%   Returns:
%       UF (N-by-N-by-(snaps+1) array): Solution in Fourier Space to the IVP
%           at times 0,...,T. UF(:,:,k) returns a 2D array corresponding to the
%           solution at the k-th snapshot. The solution at initial and final
%           times is always included.
%       KX, KY (N-by-N array): 2D wavenumbers in the default order by Matlab.
%       tsteps (1-by-N array): recorded times.
%   
%   Added to v1.3: u0 can also be a N-by-N matrix.
%
%   See also FFT2 IFFT2 RK4FFT2 MAKEITPERIODIC
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/06 (v1.4)
    
ky = [0:N/2-1 -N/2 -N/2+1:-1]'; kx = ky;
[KX,KY] = meshgrid(kx,ky);
LIN = -1i*P(KX,KY);
clear kx ky;

% ETDRK4 quantities
E = exp(dt*LIN); E2 = exp(dt*LIN/2);
M = 64; % Number of discretization points in the circle
LIN = reshape(LIN,N^2,1);
rts = exp(2*(1:M)*pi*1i/M); % Roots of unity
LR = dt*LIN(:,ones(M,1)) + rts(ones(N^2,1),:); % MOST EXPENSIVE PART
clear LIN rts;
Q     = dt*mean( (exp(LR/2)-1)./LR ,2);  % Trapezoidal rule
alpha = dt*mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2);
breta = dt*mean( (2+LR+exp(LR).*(-2+LR))./LR.^3 ,2);
gamma = dt*mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2);
clear LR;
Q = reshape(Q,N,N); alpha = reshape(alpha,N,N); 
breta = reshape(breta,N,N); gamma = reshape(gamma,N,N);

x = 2*pi*(-N/2:N/2-1)'/N; y = x;
[X,Y] = meshgrid(x,y); 
F = fft2(fun(X,Y)); 
CX = beta(X,Y);
if r == 0
    NON = @(vv,tt) -1i*F.*exp(-1i*omega0*tt);
else
    NON = @(vv,tt) 1i*r*fft2(CX.*ifft2(vv)) - 1i*F.*exp(-1i*omega0*tt);
end

if isa(u0,'function_handle')
    U = fft2(u0(X,Y));
else
    U = fft2(u0);
end

% Time-stepping loop
nmax = round(T/dt); nplt = floor((T/snaps)/dt);

UF = U; tsteps = 0;

textprogressbar('Time-stepping progress:    ');
for n = 1:nmax
    textprogressbar(ceil(100*n/nmax))
    t = (n-1)*dt;
    NU = NON(U,t);
    A = E2.*U + Q.*NU;
    NA = NON(A,t+dt/2);
    B = E2.*U + Q.*NA;
    NB = NON(B,t+dt/2);
    C = E2.*A + Q.*(2*NB-NU);
    NC = NON(C,t+dt);
    U = E.*U + alpha.*NU + 2*breta.*(NA+NB) + gamma.*NC;
    if mod(n,nplt)==0
        UF = cat(3,UF,U);
        tsteps = cat(1,tsteps,n*dt);
    end
end
textprogressbar(' Finished!');

% This is in case the endtime was not recorded.
if mod(nmax,nplt) ~= 0
    UF = cat(3,UF,U);
    tsteps = cat(1,tsteps,t);
end

tsteps = tsteps';
