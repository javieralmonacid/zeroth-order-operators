% Eigenvalues and eigenfunctions
% r = 0.5, beta(x) = cos(x_1)

close all; clear; clc; 

auxst = dbstack; filename = auxst.name;
filename = [filename,'.mat'];

if ~isfile(filename)
    r = 0.5;
    beta = @(x,y) cos(x) + 0*y;
    nu = (1:100)*1e-4;
    nu = flip(nu);
    % Select only a few for visualization
    nu = nu(8:10:78);
    N = 64;
    eigvals = 8;
    
    L = zeros(eigvals,length(nu));
    VF = cell(length(nu),1);
    VR = cell(length(nu),1);
    E0 = zeros(eigvals,N/4,length(nu));
    
    for k = 1:length(nu)
        tic
        disp(['Step ',num2str(k),'/',num2str(length(nu))])
        [lambda,V,KX,KY] = eigensolver(r,beta,nu(k),N,eigvals);
        L(:,k) = lambda;
        VF{k} = V;
        [VR{k},X,Y] = fourier2real(V);
        [R,E0(:,:,k)] = red(0,V,KX,KY);
        toc
    end
    
    save(filename);
else
    disp(['File exists. Loading ',filename])
    load(filename)
end

%% EIGENVALUES AS NU \TO 0

%L = L(:,8:10:78);

figure(9), clf
set(gcf,'Position',[784   410   1056   402])

subplot(2,1,1)
for k = 1:size(L,2)
    scatter(real(L(:,k)),imag(L(:,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    text(real(L(:,k))+0.5e-2,imag(L(:,k)),strsplit(num2str(1:size(L,1))),'FontSize',12,'Clipping','on')
    drawnow, hold on
end
grid on
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
ylabel('Im(\lambda)','FontSize',12)

subplot(2,1,2)
for k = 1:size(L,2)
    scatter(real(L(:,k)),imag(L(:,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    text(real(L(:,k))+0.3e-2,imag(L(:,k)),strsplit(num2str(1:size(L,1))),'FontSize',12,'Clipping','on')
    axis([-0.2 0.2 -0.202 -0.196])
    drawnow, hold on
end
grid on
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
ylabel('Im(\lambda)','FontSize',12)
hold off

%% RADIAL ENERGY DENSITY

figure(121), clf
set(gcf,'Position',[784   410   1056   201])
nus = 1:8;

subplot(1,4,1)
% REDO MATRIX
E = zeros(length(nus),N/4);
for k = 1:length(nus)
    E(k,:) = E0(1,:,k);
end
loglog(R,E,'LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
grid on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_0(\phi_1)}','FontSize',12)
axis([R(1) N/2 1e-20 1])
xticks([10 30])

subplot(1,4,2)
% REDO MATRIX
E = zeros(length(nus),N/4);
for k = 1:length(nus)
    E(k,:) = E0(2,:,k);
end
loglog(R,E,'LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
grid on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_0(\phi_2)}','FontSize',12)
axis([R(1) N/2 1e-20 1])
xticks([10 30])

subplot(1,4,3)
% REDO MATRIX
E = zeros(length(nus),N/4);
for k = 1:length(nus)
    E(k,:) = E0(5,:,k);
end
loglog(R,E,'LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
grid on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_0(\phi_5)}','FontSize',12)
axis([R(1) N/2 1e-20 1])
xticks([10 30])

subplot(1,4,4)
% REDO MATRIX
E = zeros(length(nus),N/4);
for k = 1:length(nus)
    E(k,:) = E0(6,:,k);
end
loglog(R,E,'LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
grid on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_0(\phi_6)}','FontSize',12)
axis([R(1) N/2 1e-20 1])
xticks([10 30])

%% EIGENFUNCTIONS IN REAL SPACE

VR = VR{end}; VF = VF{end}; % Select smallest viscosity
disp(['Eigenfunctions to be plot correspond to nu = ',num2str(nu(end))])

% Prepare the functions for plotting
[VR,X,Y]   = makeitperiodic(VR,X,Y);
[VF,KX,KY] = makeitperiodic(VF,KX,KY);

figure(131), clf
set(gcf,'Position',[784   410   1056   201])
plt = 1;
for k = [1 2 5 6]
    subplot(1,4,plt)
    surf(X,Y,abs(VR(:,:,k))),
    shading interp, colormap(brewermap([],'*RdBu'));
	xlabel('x_1'); ylabel('x_2');
    colorbar; axis tight
    title(['$|\phi_',num2str(k),'|$'],'Interpreter', 'latex','FontSize',12)
    view([0 90])
    plt = plt + 1;
end

%% EIGENFUNCTIONS IN FOURIER SPACE AND LONG-TERM EVOLUTION

figure(14), clf
set(gcf,'Position',[784   245   1056   366])
ef = [1 2 5 6];
plot_position = [1 2 5 6];
for k = 1:length(ef)
    p(k) = subplot(2,4,plot_position(k));
    surf(KX,KY,abs(VF(:,:,ef(k)))),
    shading interp, colormap(brewermap([],'*RdBu'));
	xlabel('$k_1$','Interpreter','latex'); ylabel('$k_2$','Interpreter','latex');
    colorbar; axis tight
    title(['$|\widehat{\phi_',num2str(ef(k)),'}|$'],'Interpreter', 'latex','FontSize',12)
    view([0 90])
end

load('figure_7.mat')
p(5) = subplot(2,4,[3 4 7 8]);
surf(KX1,KY1,abs(UF1(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('$k_1$','Interpreter','latex'); ylabel('$k_2$','Interpreter','latex'); view([0 90])
colorbar; axis tight; set(gca,'ColorScale','log')
title('$|\widehat{u_1}|$ ($t=2000$)', 'Interpreter', 'latex', 'FontSize', 12)
