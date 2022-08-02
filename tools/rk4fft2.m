function [UF,KX,KY,tsteps] = rk4fft2(P,beta,r,omega0,N,dt,T,fun,u0,snaps)
% RK4FFT2 RK4-FFT2 based transient solver
%   RK4FFT2(P,beta,r,omega0,N,dt,T,fun,u0,snaps) returns the solution
%   to the problem in the 2-torus
%
%       i u_t - (P(x,y,D) + r*beta(x,y))* u = fun(x,y)*exp(-i*omega0*t),
%       u = u0 at initial time.
%
%   The time discretization is based on the classical fourth-order Runge-Kutta
%   method. In turn, the space pseudo-spectral discretization uses the 
%   built-in function FFT2.
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
%   See also FFT2 IFFT2 ETDRK4FFT2 MAKEITPERIODIC
%
%   Author: Javier Almonacid
%           Department of Mathematics
%           Simon Fraser University
%   Date:   2020/05/06 (v1.4)

x = 2*pi*(-N/2:N/2-1)'/N; y = x;
ky = [0:N/2-1 -N/2 -N/2+1:-1]'; kx = ky;
[X,Y] = meshgrid(x,y); [KX,KY] = meshgrid(kx,ky);
F = fft2(fun(X,Y));
LIN = -1i*P(KX,KY);
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

nmax = round(T/dt); nplt = floor((T/snaps)/dt);

UF = U; tsteps = 0;

textprogressbar('Time-stepping progress:    ');
for n = 1:nmax
    textprogressbar(ceil(100*n/nmax))
    t = (n-1)*dt;
    k1 = dt*(LIN.*U + NON(U,t));
    k2 = dt*(LIN.*(U+k1/2) + NON(U+k1/2,t+dt/2));
    k3 = dt*(LIN.*(U+k2/2) + NON(U+k2/2,t+dt/2));
    k4 = dt*(LIN.*(U+k3) + NON(U+k3,t+dt));
    U = U + (k1 + 2*k2 + 2*k3 + k4)/6;
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
