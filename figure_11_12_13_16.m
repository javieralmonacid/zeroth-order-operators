% Figure 11. Evolution of the first seven eigenvalues for Test 3 
% ($r=0.55$, $\beta(x) = \cos(x_1-2x_2)+\sin(2x_2)$). Left: first 
% three eigenvalues when the viscosity decreases from $10^{-2}$ to $4.4 
% \cdot 10^{-3}$ (below this viscosity, these eigenvalues become the 
% first, second and fifth eigenvalues, respectively).  Right: all seven 
% eigenvalues as $\nu$ decreases from $10^{-2}$ to $3 \cdot 10^{-4}$.
%
% Figure 12. Radial energy density $E_0$ (log-log scale) of different 
% viscous eigenfunctions $\phi_j$. Each curve represents a fixed value of 
% $\nu$ within the chosen range. In general, as the viscosity decreases, 
% the curves move to the right, which shows that the viscous 
% eigenfunctions become less regular.
%
% Figure 13. Magnitude of eigenmodes in real space $\phi_j$ for Tests 1, 
% 2 and 3. The shape of some of the eigenmodes resembles that of the 
% attractors in Figure 7.
%
% Figure 16. Left half: Magnitude of some eigenmodes in frequency space 
% $\widehat{\phi_j}$ of $\ip{D}^{-1} D_{x_2} + i\nu\Delta 
% - 0.55\left(\cos(x_1-2x_2) + \sin(2x_2)\right)$ 
% ($\nu = 3 \cdot 10^{-4}$). Right half: Long-term evolution in frequency 
% space (this is the Fourier transform of Figure \ref{fig:lte}-right).

close all; clear; clc; 
addpath tools;

auxst = dbstack; filename = auxst.name;
filename = [filename,'.mat'];

if ~isfile(filename)
    r = 0.55;
    beta = @(x,y) cos(x-2*y) + sin(2*y);
    nu = (1:100)*1e-4;
    nu = flip(nu);
    % Select a few viscosities for visualization
    nu = nu([1:8:57, 58:4:98]);
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

%% REORDERING EIGENVALUES

% SORT BY MAGNITUDE
for k = 1:size(L,2)
    [lambda,ind] = sort(L(:,k));
    L(:,k) = L(ind,k);
    VF{k} = VF{k}(:,:,ind);
    VR{k} = VR{k}(:,:,ind);
end
E0 = E0(ind,:,:);
L = L(1:7,:);  % WE'LL WORK WITH ONLY THE FIRST THREE ONES.

% SECOND SORTING
for k = 1:size(L,2)
    ind1 = find(abs(real(L(:,k))) <= 0.5e-3);
    ind2 = find(real(L(:,k)) > 0.5e-3);
    ind3 = find(real(L(:,k)) < -0.5e-3);
    ind = [ind1;ind2;ind3];
    L(:,k) = L(ind,k);
    VF{k} = VF{k}(:,:,ind);
    VR{k} = VR{k}(:,:,ind);
end
E0 = E0(ind,:,:);
%% EIGENVALUES (7)

%{
figure(1), clf
set(gcf,'Position',[784   410   1056   402])

subplot(1,4,1)
for k = 1:8:57
    scatter(real(L(1:3,k)),imag(L(1:3,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    text(real(L(1,k))+0.2e-3,imag(L(1,k)),'1','FontSize',12,'Clipping','on')
    text(real(L(2,k))+0.2e-3,imag(L(2,k)),'2','FontSize',12,'Clipping','on')
    text(real(L(3,k))-0.3e-3,imag(L(3,k)),'3','FontSize',12,'Clipping','on')
    axis([-1e-3 1e-3 -0.07 0])
    drawnow, hold on
end
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
ylabel('Im(\lambda)','FontSize',12)
grid on

subplot(1,4,[2 3 4])
for k = 1:8:57
    scatter(real(L(1:3,k)),imag(L(1:3,k)),36,'o','filled','MarkerFaceColor',[0.741, 0.741, 0.741])
    hold on
end
for k = 58:4:98
    scatter(real(L(1:3,k)),imag(L(1:3,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    text(real(L(1,k))+0.4e-3,imag(L(1,k)),'1','FontSize',12,'Clipping','on')
    text(real(L(2,k))+0.4e-3,imag(L(2,k)),'2','FontSize',12,'Clipping','on')
    text(real(L(3,k))-0.5e-3,imag(L(3,k)),'3','FontSize',12,'Clipping','on')
    axis([-8e-3 8e-3 -0.07 0])
    drawnow, hold on
end
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
%ylabel('Im(\lambda)','FontSize',12)
grid on
%}

figure(11), clf
set(gcf,'Position',[784   410   1056   402])

subplot(1,4,1)
for k = 1:8
    scatter(real(L(1:3,k)),imag(L(1:3,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    text(real(L(1,k))+0.2e-6,imag(L(1,k)),'1','FontSize',12,'Clipping','on')
    text(real(L(2,k))+0.2e-6,imag(L(2,k)),'2','FontSize',12,'Clipping','on')
    text(real(L(3,k))-0.3e-6,imag(L(3,k)),'3','FontSize',12,'Clipping','on')
    axis([-1e-6 1e-6 -0.07 0])
    drawnow, hold on
end
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
ylabel('Im(\lambda)','FontSize',12)
grid on

subplot(1,4,[2,3,4])
for k = 1:8
    scatter(real(L(1:3,k)),imag(L(1:3,k)),36,'o','filled','MarkerFaceColor',[0.741, 0.741, 0.741])
    hold on
    scatter(real(L(4:7,k)),imag(L(4:7,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    text(real(L(4:5,k))+1.5e-3,imag(L(4:5,k)),strsplit(num2str(4:5)),'FontSize',12,'Clipping','on')
    text(real(L(6:7,k))-2e-3,imag(L(6:7,k)),strsplit(num2str(6:7)),'FontSize',12,'Clipping','on')
end
counter = 1;
for k = 9:length(nu)
    scatter(real(L(:,k)),imag(L(:,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    if mod(counter,2) ~= 0
        text(real(L(1:4,k))+1.5e-3,imag(L(1:4,k)),strsplit(num2str(1:4)),'FontSize',12,'Clipping','on')
        text(real(L(5:7,k))-2.5e-3,imag(L(5:7,k)),strsplit(num2str(5:7)),'FontSize',12,'Clipping','on')
    end
    axis([-0.04 0.04 -0.2 0])
    drawnow, hold on
    counter = counter+1;
end
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
%ylabel('Im(\lambda)','FontSize',12)
grid on

%{
for k = 1:size(L,2)
    scatter(real(L(1,k)),imag(L(1,k)),36,'o','filled','MarkerFaceColor',[0.741, 0.741, 0.741])
    scatter(real(L(2:7,k)),imag(L(2:7,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    if mod(k,2) == 0
        text(real(L(1:4,k))+1.5e-3,imag(L(1:4,k)),strsplit(num2str(1:4)),'FontSize',12,'Clipping','on')
        text(real(L(5:7,k))-2.5e-3,imag(L(5:7,k)),strsplit(num2str(5:7)),'FontSize',12,'Clipping','on')
    end
    drawnow, hold on
end
grid on, axis([-0.04 0.04 -0.15 0])
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
ylabel('Im(\lambda)','FontSize',12)
%}

%% RADIAL ENERGY DENSITY
figure(123), clf
set(gcf,'Position',[784   410   1056   201])

plt = 1;
for k = [1 2 5]
    subplot(1,3,plt)
    %loglog(R,E0(k,:,9),'k--','LineWidth',2)
    %hold on
    for l = 9:length(nu)
        loglog(R,E0(k,:,l),'LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
        hold on
        %pause
    end
    %loglog(R,E0(k,:,end),'--','LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
    
    %loglog(R,7e1*R.^(-5),'k','LineWidth',1)
    
    
    grid on
    set(gca,'FontSize',12)
    xlabel('R','FontSize',12)
    ylabel(['RED {E_0(\phi_',num2str(k),')}'],'FontSize',12)
    axis([R(1) N/2 1e-20 1])
    plt = plt+1;
    xticks([10 20 30])
end

%% EIGENFUNCTIONS IN REAL SPACE

[VVR,XX,YY] = makeitperiodic(VR{end},X,Y);
[VVF,KKX,KKY] = makeitperiodic(VF{end},KX,KY);

figure(133), clf
set(gcf,'Position',[784   410   1056   201])
plt = 1;
for k = [1 2 5]
    subplot(1,3,plt)
    surf(XX,YY,abs(VVR(:,:,k))),
    shading interp, colormap(brewermap([],'*RdBu'));
	xlabel('x_1'); ylabel('x_2');
    colorbar; axis tight;
    view([0 90])
    title(['$|\phi_',num2str(k),'|$'],'Interpreter', 'latex','FontSize',12)
    plt = plt + 1;
end

%% EIGENFUNCTIONS IN FOURIER SPACE AND LONG-TERM EVOLUTION

figure(16), clf
set(gcf,'Position',[784   245   1056   366])
ef = [1 2 5];
plot_position = [1 2 6];
for k = 1:3
    p(k) = subplot(2,4,plot_position(k));
    surf(KKX,KKY,abs(VVF(:,:,ef(k)))),
    shading interp, colormap(brewermap([],'*RdBu'));
	ylabel('$k_2$','Interpreter','latex');
    if k == 3 || k == 1
        xlabel('$k_1$','Interpreter','latex');
    end
    colorbar
    axis tight
    view([0 90])
    title(['$|\widehat{\phi_',num2str(ef(k)),'}|$'],'Interpreter', 'latex','FontSize',12)
end
hold off

load('figure_7.mat')
p(4) = subplot(2,4,[3 4 7 8]);
surf(KX3,KY3,abs(UF3(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('$k_1$','Interpreter','latex'); ylabel('$k_2$','Interpreter','latex'); view([0 90])
colorbar; axis tight; set(gca,'ColorScale','log')
title('$|\widehat{u_3}|$ ($t=1000$)', 'Interpreter', 'latex', 'FontSize', 12)

