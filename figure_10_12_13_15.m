% Eigenvalues and eigenfunctions

close all; clear; clc; 

auxst = dbstack; filename = auxst.name;
filename = [filename,'.mat'];

if ~isfile(filename)
    r = 0.45;
    beta = @(x,y) cos(x-2*y) + sin(2*y);
    nu = (1:100)*1e-4;
    nu = flip(nu);
    % Select a few for visualization
    nu = nu([1:4:81, 82:2:98]);
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
L = L(1:7,:);  % WE'LL WORK WITH ONLY THE FIRST SEVEN ONES.

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

%% EIGENVALUES AS NU \TO 0

figure(10), clf
set(gcf,'Position',[784   410   1056   402])

subplot(1,4,1)
for k = 1:size(L,2)
    scatter(real(L(1,k)),imag(L(1,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    if mod(k,2) == 0
        text(real(L(1,k))+0.2e-6,imag(L(1,k)),'1','FontSize',12,'Clipping','on')
    end
    drawnow, hold on
    axis([-1e-6 1e-6 -1.2e-2 0])
end
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
ylabel('Im(\lambda)','FontSize',12)
grid on

subplot(1,4,[2 3 4])
for k = 1:size(L,2)
    scatter(real(L(1,k)),imag(L(1,k)),36,'o','filled','MarkerFaceColor',[0.741, 0.741, 0.741])
    scatter(real(L(2:7,k)),imag(L(2:7,k)),36,'o','filled','MarkerFaceColor',[0.6350 0.0780 0.1840])
    if mod(k,2) ~= 0
        text(real(L(1:4,k))+1.5e-3,imag(L(1:4,k)),strsplit(num2str(1:4)),'FontSize',12,'Clipping','on')
        text(real(L(5:7,k))-2e-3,imag(L(5:7,k)),strsplit(num2str(5:7)),'FontSize',12,'Clipping','on')
    end
    drawnow, hold on
end
grid on, axis([-0.04 0.04 -0.15 0])
set(gca,'FontSize',12)
xlabel('Re(\lambda)','FontSize',12)
%ylabel('Im(\lambda)','FontSize',12)

%% RADIAL ENERGY DENSITY

figure(122), clf
set(gcf,'Position',[784   410   1056   201])

plt = 1;
for k = [1 2 3 4]
    subplot(1,4,plt)
    for l = length(nu)-10:length(nu)
        loglog(R,E0(k,:,l),'LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
        hold on
    end
    grid on
    set(gca,'FontSize',12)
    xlabel('R','FontSize',12)
    ylabel(['RED {E_0(\phi_',num2str(k),')}'],'FontSize',12)
    axis([R(1) N/2 1e-20 1])
    plt = plt+1;
    xticks([10 30])
end

%% EIGENFUNCTIONS IN REAL SPACE

[VVR,XX,YY] = makeitperiodic(VR{end},X,Y);
[VVF,KKX,KKY] = makeitperiodic(VF{end},KX,KY);

figure(132), clf
set(gcf,'Position',[784   410   1056   201])
plt = 1;
for k = [1 2 3 4]
    subplot(1,4,plt)
    surf(XX,YY,abs(VVR(:,:,k))),
    shading interp, colormap(brewermap([],'*RdBu'));
	xlabel('x_1'); ylabel('x_2');
    colorbar; axis tight
    view([0 90])
    title(['$|\phi_',num2str(k),'|$'],'Interpreter', 'latex','FontSize',12)
    plt = plt + 1;
end

%% EIGENFUNCTIONS IN FOURIER SPACE AND LONG-TERM EVOLUTION

figure(15), clf
set(gcf,'Position',[784   245   1056   366])
plot_position = [1 2 5 6];
for k = [1 2 3 4]
    p(k) = subplot(2,4,plot_position(k));
    surf(KKX,KKY,abs(VVF(:,:,k)))
    shading interp, colormap(brewermap([],'*RdBu'));
	ylabel('$k_2$','Interpreter','latex');
    if k == 3 || k == 4
        xlabel('$k_1$','Interpreter','latex');
    end
    colorbar; axis tight;
    title(['$|\widehat{\phi_',num2str(k),'}|$'],'Interpreter', 'latex','FontSize',12)
    view([0 90])
end
hold off

load('figure_7.mat')
p(5) = subplot(2,4,[3 4 7 8]);
surf(KX2,KY2,abs(UF2(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('$k_1$','Interpreter','latex'); ylabel('$k_2$','Interpreter','latex'); view([0 90])
colorbar; axis tight; set(gca,'ColorScale','log')
title('$|\widehat{u_2}|$ ($t=1000$)', 'Interpreter', 'latex', 'FontSize', 12)
