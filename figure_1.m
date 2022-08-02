% Figure 1. Development of attractors for the pseudo-differential operator 
% given in (2.5). That is, we solve
%                  iu_t - Pu = f
% where P is the zeroth-order operator with symbol
%         p(x,k) = k_2 <k>^{-1} - 2 \cos(x_1)
% and 
%         f(x,t) = -5(-|x|^2 + i(2 x_1 + x_2)).
%
% Initially, u=0. Then, we evolve this until T=300 using a time step 
% dt = 0.5 and a 256x256 mesh. We show the results at times t = 50, 100,
% 150, 200, 250, and 300.
close all; clear; clc;

%% Compute solution and save it to disk.
% We save the data in a .mat file with the same name as this script.
% If the script finds this mat file in the folder, it will read it
% instead of computing the solution again. Delete this file if you
% want to recompute the solution.

auxst = dbstack; filename = auxst.name;
filename = [filename,'.mat'];

if ~isfile(filename)
    % Set the data and the symbol for the pseudo-differential operator.
    omega0 = 0; r = 2;
    fun = @(x,y) -5*exp(-((x+0.0).^2 + (y+0.0).^2)*1+1i*2*x + 1i*y);
    u0 = @(x,y) 0*x+0*y;
    P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2);
    beta = @(x,y) cos(x) + 0*y;
    
    % Set the discretization parameters and end time.
    N = 2^8; dt = 0.5; T = 300;
    
    % Solve the PDE using an spectral method with an ETDRK4 time
    % discretization. Since dt = 0.5 and T = 300, we can output up to 600
    % snapshots. We only choose to retrieve 6 of them. The routine will
    % output the solution as an N-by-N-by-(snapshots+1) matrix. The first
    % N-by-N layer corresponds to the initial solution. As we will not be
    % plotting this solution in Fourier space, we can omit the second and
    % third arguments in the output, as this corresponds to the spatial
    % mesh in Fourier space.
    [UF,~,~,t] = etdrk4fft2(P,beta,r,omega0,N,dt,T,fun,u0,6);

    % We first transform the solution to real space using FFT2 and then
    % make it periodic so that we can plot it.
    [U,X,Y] = fourier2real(UF);
    [U,X,Y] = makeitperiodic(U,X,Y);

    % Save MAT file
    save(filename)
else
    disp(['File exists. Loading ',filename])
    load(filename);
end

%% Main Plot
% As the solution is complex valued, we will plot only its magnitude.
figure(1)

% Get the upper limit for the color bar
umax = max(max(abs(U(:,:,end))));
umin = 0;

% As the first recorded time step is the initial solution u = 0, we plot
% the solution from the second time step onwards.
for k = 2:7
    subplot(2,3,k-1)
    surf(X,Y,abs(U(:,:,k))),
    shading interp, colormap(brewermap([],'*RdBu'));
	xlabel('x_1'); ylabel('x_2'); zlabel('|u|')
    colorbar; caxis([0 umax])
    axis([-pi pi -pi pi umin umax]),
    view([-17.4692 21.2159])
    set(0,'DefaultAxesTitleFontWeight','normal');
    title(['t = ',num2str(t(k))])
    
end
set(gcf,'Position',[358 467 1116 420])