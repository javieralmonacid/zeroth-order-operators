% Figure 3. Spectral convergence of the first 12 eigenvalues for nu = 0.01
% and nu = 0.001.
% 
% Parameters for the tests:
% Symbol: p(x,k) = k_2 <k>^{-1} + i\nu|k|^2 - 0.5 (\cos(x_1)+\sin(x_2)),
% nu = 0.01 and nu = 0.001,
% N = 12, 16, 20, ... , 64 in each case.
% Nex = 80.
%
% Expected computational usage (as measured using a 6-core Intel Core i5-9600K 
% 3.70GHz, 16GB of RAM):
% Time: ~4 minutes
% Memory: ~6 GB 
close all; clear; clc;
addpath tools;

%% Convergence study of the eigensolver for nu = 0.01

% Get this filename and set the MAT file as filename_1.mat
auxst = dbstack; filename = auxst.name;
filename = [filename,'_1.mat'];

if ~isfile(filename)
    % Set the parameters to solve the viscous eigenvalue problem.
    r = 0.5; beta = @(x,y) cos(x)+sin(y); nu = 1e-2; 
    eigvals = 12; % First 12 eigenvalues
    N = (12:4:64)'; Nex = 80;
    
    % Define vectors to store eigenvalues and errors.
    lambda = zeros(eigvals,length(N));
    errorl = zeros(eigvals,length(N));
    erroref = zeros(eigvals,length(N));
    
    % Compute "exact" eigenvalues and eigenfunctions
    [lambdaex,VFex,KXex,KYex] = eigensolver(r,beta,nu,Nex,eigvals);
    disp('Exact eigenvalues computed')
    % We will need to fftshift the KX and KY variables so that we can
    % project the eigenfunctions onto coarser meshes.
    KXexsh = fftshift(KXex);
    KYexsh = fftshift(KYex);
    
    % For each value of N, compute eigenvalues and eigenfunctions.
    for k = 1:length(N)
        disp(['Computing step ',num2str(k),'/',num2str(length(N))])
        [L,VF,KX,KY] = eigensolver(r,beta,nu,N(k),eigvals);
        lambda(:,k) = L;
        
        % Magnitude of the error between exact and computed eigenvalues.
        errorl(:,k) = abs(lambdaex-L);

        % L2-norm of the error between exact and computed eigenfunctions.
        KXsh = fftshift(KX);
        KYsh = fftshift(KY);
        for l = 1:eigvals
            VFsh = fftshift(VFex(:,:,l));
            VFexint = interp2(KXexsh,KYexsh,VFsh,KXsh,KYsh,'spline');
            erroref(l,k) = norm(abs(VFexint)-abs(fftshift(VF(:,:,l))),'fro');
        end
    end
    
    % Save the data as figure_3_1.mat
    save(filename);
else
    disp(['File exists. Loading ',filename])
    load(filename)
end

%% Plot the error for nu = 0.01
figure(31)
loglog(N,errorl,'-s','LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
grid on
set(gca,'FontSize',12)
xlabel('N','FontSize',12)
ylabel('Eigenvalue Error Magnitude','FontSize',12)
set(gcf,'Position',[784   410   528   402])

%% Convergence study of the eigensolver for nu = 0.001
clear;

% Get this filename and set the MAT file as filename_2.mat
auxst = dbstack; filename = auxst.name;
filename = [filename,'_2.mat'];

if ~isfile(filename)
    % Repeat the steps above but this time with nu = 1e-3.
    r = 0.5; beta = @(x,y) cos(x)+sin(y); nu = 1e-3; 
    eigvals = 12; % First 12 eigenvalues
    N = (12:4:64)'; Nex = 80;

    lambda = zeros(eigvals,length(N));
    errorl = zeros(eigvals,length(N));
    erroref = zeros(eigvals,length(N));
    
    % Compute "exact" eigenvalues
    [lambdaex,VFex,KXex,KYex] = eigensolver(r,beta,nu,Nex,eigvals);
    disp('Exact eigenvalues computed')
    KXexsh = fftshift(KXex);
    KYexsh = fftshift(KYex);
    
    % Compute eigenvalues for the rest of the refinements
    for k = 1:length(N)
        disp(['Computing step ',num2str(k),'/',num2str(length(N))])
        [L,VF,KX,KY] = eigensolver(r,beta,nu,N(k),eigvals);
        lambda(:,k) = L;
        
        % Eigenvalue error
        errorl(:,k) = abs(lambdaex-L);

        % Eigenfunction error
        KXsh = fftshift(KX);
        KYsh = fftshift(KY);
        for l = 1:eigvals
            VFsh = fftshift(VFex(:,:,l));
            VFexint = interp2(KXexsh,KYexsh,VFsh,KXsh,KYsh,'spline');
            erroref(l,k) = norm(abs(VFexint)-abs(fftshift(VF(:,:,l))),'fro');
        end
    end
    
    % Save the data as figure_3_2.mat
    save(filename);
else
    disp(['File exists. Loading ',filename])
    load(filename)
end

%% Plot the error for nu = 0.001
figure(32)
loglog(N,errorl,'-s','LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
grid on
set(gca,'FontSize',12)
xlabel('N','FontSize',12)
ylabel('Eigenvalue Error Magnitude','FontSize',12)
set(gcf,'Position',[784   410   528   402])