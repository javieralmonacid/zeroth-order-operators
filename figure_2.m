% Figure 2. Convergence studies. Left: Fourth-order convergence in time.
% Right: Local spectral accuracy in space.
%
% Parameters for the test:
% Symbol: p(x,k) = k_2 <k>^{-1} - 2 \cos(x_1),
% Forcing function: f = -5*exp(-((x_1+0.9).^2 + (x_2+0.8).^2)*3 + i*2*x_1 + i*x_2),
% End time: T = 10.
%
% Expected computational usage (as measured using a 6-core Intel Core i5-9600K 
% 3.70GHz, 16GB of RAM):
% - Time convergence study: ~20 min, ~3 GB of RAM
% - Space convergence study: ~2 min, ~7 GB of RAM
close all; clear; clc;
addpath tools;

%% Convergence in time for the ETDRK4 and RK4 methods
% We consider a forcing frequency omega_0 = 0.1 and a 64-by-64 spatial
% mesh. The exact solution will be taken as the one computed using a
% significantly small time step size (dt = 2^(-10)*1e-02 \approx 1e-05).
% We then compare this solution with those obtained using time step sizes
% dt = 1, 2^(-1), 2^(-2), ... , 2^(-8).

% Get this filename and set the MAT file as filename_time.mat
auxst = dbstack; filename = auxst.name;
filename = [filename,'_time.mat'];

if ~isfile(filename)
    omega0 = 0.1; T = 10; r = 2;
    N = 2^6; dt = 1./2.^(0:8); dtex = 2^(-10)*1e-2;

    fun = @(x,y) -5*exp(-((x+0.9).^2 + (y+0.8).^2)*3+1i*2*x + 1i*y);
    u0 = @(x,y) 0*x+0*y;
    P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2);
    beta = @(x,y) cos(x)+0*y;

    % Solver: RK4. 
    % Compute the exact solution first. We only need the last snapshot of
    % the solution, so we set the last parameter 1. In this way, the output
    % only contains the initial and final solutions.
    tic
    Uex = rk4fft2(P,beta,r,omega0,N,dtex,T,fun,u0,1);
    toc
    disp('Exact solution computed')

    % Set vectors to store errors and rates of convergence.
    error1 = zeros(size(dt));
    rate1 = zeros(size(dt));

    for k = 1:length(dt)
        tic
        % For each time step size, solve the PDE using an RK4 time
        % discretization.
        U = rk4fft2(P,beta,r,omega0,N,dt(k),T,fun,u0,1);

        % Compute errors and rates of convergence. The error is measured in
        % the whole domain (this will not be the case for the convergence
        % study in space) using a discrete L2-norm. Recall that the
        % solution is returned in Fourier space.
        error1(k) = norm(U(:,:,end)-Uex(:,:,end),'fro')*2*pi/N^2; % L2-norm
        if k > 1
            rate1(k) = log(error1(k-1)/error1(k))/log(dt(k-1)/dt(k));
        end

        disp(['Error ',num2str(k),'/',num2str(length(dt)),' computed.'])
        toc
    end
    
    % Repeat the steps above but using this time the ETDRK4 time
    % discretization. This method is just slightly faster than RK4.
    % However, ETDRK4 tends to require more RAM than RK4.
    tic
    Uex = etdrk4fft2(P,beta,r,omega0,N,dtex,T,fun,u0,1);
    toc
    disp('Exact solution computed')

    error2 = zeros(size(dt));
    rate2 = zeros(size(dt));

    for k = 1:length(dt)
        tic
        U = etdrk4fft2(P,beta,r,omega0,N,dt(k),T,fun,u0,1);
        error2(k) = norm(U(:,:,end)-Uex(:,:,end),'fro')*2*pi/N^2; % L2-norm
        if k > 1
            rate2(k) = log(error2(k-1)/error2(k))/log(dt(k-1)/dt(k));
        end
        disp(['Error ',num2str(k),'/',num2str(length(dt)),' computed.'])
        toc
    end

    % Save as filename_time.mat
    save(filename)
else
    disp(['File exists. Loading ',filename])
    load(filename);
end

disp(' ')
disp('L2-norm of the error using RK4'), disp(error1)
disp('L2-norm of the error using ETDRK4'), disp(error2)
disp('Rate of convergence using RK4'), disp(rate1)
disp('Rate of convergence using ETDRK4'), disp(rate2)

% For plotting, we compute the errors relative to the size of the exact
% solution.
error1 = error1/(norm(Uex(:,:,end),'fro')*2*pi/N^2); 
error2 = error2/(norm(Uex(:,:,end),'fro')*2*pi/N^2); 

figure(21)
loglog(dt,error1,'-d',dt,error2,'-s',dt,1e-2*dt.^4,'k-.','LineWidth',2)
axis([dt(end) dt(1) 1e-12 1]), grid on
xlabel('{\Delta}t','FontSize',12)
ylabel('L^2-Relative Error in {[-\pi,\pi]^2}','FontSize',12)

legend('RK4','ETDRK4','O({\Delta}t^4)','FontSize',12,'Location','NorthWest')
set(gca,'FontSize',12)
set(gcf,'Position',[669   530   520   391])

%% Convergence in space using the ETDRK4 method.
% We consider a forcing frequency omega_0 = 0.5 and a time step size of
% dt = 0.01. The exact solution will be taken as the one computed using a
% significantly refined mesh (N = 2^10). We then compare this solution 
% with those obtained using N = 2^3, 2^4, ... , 2^7. We will also use
% different forcing functions:
%  f_1(x) = sin(x_1) cos(2x_2),
%  f_2(x) = sin(x_1) cos(2x_2) + sin(5x_1)cos(2x_2) + i*sin(5x_1)cos(4x_2),
%  f_3(x) = 0.5*exp(-2|x|^2),
%  f_4(x) = -5*exp(-((x_1+0.9)^2 + (x_2+0.8).^2)*3 + i*2*x_1 + i*x_2).
% Differently from the time convergence study, the singularities will
% force us to compute the error in a region where the solution is smooth.
% For these particular examples, such region will be taken as the square
% [-pi/4, pi/4]^2.
clear;

% Get this filename and set the MAT file as filename_space.mat
auxst = dbstack; filename = auxst.name;
filename = [filename,'_space.mat'];

if ~isfile(filename)
    % Set the parameters and forcing functions. We will use an ETDRK4 time
    % discretization.
    solver = 'etdrk4fft2'; dt = 1e-01;
    omega0 = 0.5; T = 10; r = 2; Nex = 2^10;
    nvalues = 2.^(3:7);
    funcs = {@(x,y) sin(x).*cos(2*y), ...
             @(x,y) sin(x).*cos(2*y) + sin(5*x).*cos(2*y) + 1i*sin(5*x).*cos(4*y), ...
             @(x,y) 0.5*exp(-2*(x.^2 + y.^2)), ...
             @(x,y) -5*exp(-((x+0.9).^2 + (y+0.8).^2)*3+1i*2*x + 1i*y)};
    u0 = @(x,y) 0*x+0*y;
    P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2);
    beta = @(x,y) cos(x)+0*y;
    
    % Set matrices to store L^2 and H^(-1) errors and associated rates of
    % convergence.
    error2  = zeros(length(funcs),length(nvalues));
    errorm1 = zeros(length(funcs),length(nvalues));
    rate2   = zeros(length(funcs),length(nvalues));
    ratem1  = zeros(length(funcs),length(nvalues));

    for k = 1:length(funcs)
        disp(['Forcing term ', num2str(k), '/', num2str(length(funcs)), ':'])

        % Compute the exact solution.
        tic
        UexF = feval(solver,P,beta,r,omega0,Nex,dt,T,funcs{k},u0,1);
        disp('Exact solution computed')
        UexF = UexF(:,:,end);       % Discard initial time step
        Uex = fourier2real(UexF);   % Convert into real space
        toc

        % Select only those values in the region [-pi/4,pi/4]^2
        Uex0 = Uex(3*Nex/8+1:5*Nex/8+1,3*Nex/8+1:5*Nex/8+1);

        % For each number of degrees of freedom N, compute the solution
        for l = 1:length(nvalues)
            N = nvalues(l);
            [UF,KX,KY] = feval(solver,P,beta,r,omega0,N,dt,T,funcs{k},u0,1);
            disp(['Solution ', num2str(l), '/', num2str(length(nvalues)),' computed'])

            UF = UF(:,:,end);     % Discard initial time step
            U = fourier2real(UF); % Convert into real space

            % Compute H^{-1} error in the whole domain.
            ind = (Nex/N)*(0:N-1)+1;
            errorm1(k,l) = hsnorm(UexF(ind,ind)-UF,KX,KY,-1);

            % Reduce the solution to the region [-pi/4,pi/4]^2
            U0 = U(3*N/8+1:5*N/8+1,3*N/8+1:5*N/8+1);
            x = 2*pi*(-N/8:N/8)/N;

            % Note that the solution is not periodic in the considered 
            % subregion, so we cannot use Parseval's identity to compute
            % the error in Fourier space. Hence, we compute the relative 
            % error in real space using the trapezoidal rule.
            skip = (length(Uex0)-1)/(length(U0)-1);
            aux = sqrt( trapz(x,trapz(x,abs(Uex0(1:skip:end,1:skip:end)).^2,2)) );
            error2(k,l) = sqrt( trapz(x,trapz(x,abs(Uex0(1:skip:end,1:skip:end)-U0).^2,2)) )/aux; % Relative error
            
            % Rates of convergence
            if l > 1
                ratem1(k,l) = log(errorm1(k,l-1)/errorm1(k,l))/log(nvalues(l)/nvalues(l-1));
                rate2(k,l) = log(error2(k,l-1)/error2(k,l))/log(nvalues(l)/nvalues(l-1));
            end
        end
    end

    % Save as filename_space.mat
    save(filename)
else
    disp(['File exists. Loading ',filename])
    load(filename)
    %load('/home/javier/Dropbox/sfu/MSc Thesis/nigam-zworski-01/internal-waves/old_mat_spaceconvergence.mat');
end

figure(22)
loglog(nvalues,error2(1,:),'-*',nvalues,error2(2,:),'-s',nvalues,error2(3,:),'-d',nvalues,error2(4,:),'-o','LineWidth',2)
grid on
xlabel('N','FontSize',12); 
ylabel('L^2-Relative Error in {[-\pi/4,\pi/4]^2}','FontSize',12)
legend('f_1','f_2','f_3','f_4','location','SouthWest','FontSize',12)
set(gca,'FontSize',12)
set(gcf,'Position',[669   530   520   391])
xticks(nvalues)
