% Figure 5. Radial energy density $E_0$ at several times ($r = 0.25$, 
% $\beta(x) = \cos(x_1)+\sin(x_2)$).
%
% Figure 6. Radial energy density for different values of $s$. This 
% suggests that $E_s(u_N)(R)$ decays faster than $R^{2s}$ for any 
% $s \leq 0$.

close all; clear; clc;
addpath tools;

auxst = dbstack; filename = auxst.name;
filename = [filename,'.mat'];

if ~isfile(filename)
    beta = @(x,y) cos(x)+sin(y);

    fun = @(x,y) -5*exp(-((x+0.0).^2 + (y+0.0).^2)*3+1i*2*x + 1i*y);
    u0 = @(x,y) 0*x+0*y;
    P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2);

    N = 2^10; dt = 0.05; T = 2000; snaps = 20; r = 0.25; omega0 = 0;
    
    [UF,KX,KY,t] = etdrk4fft2(P,beta,r,omega0,N,dt,T,fun,u0,snaps);
    save(filename)
else
    disp(['File exists. Loading ',filename])
    load(filename)
    t = t(3:2:end)';
    UF = UF(:,:,3:2:end);
end

[R,G] = red(0,UF,KX,KY);

%% RED E0 (Figure 5)

figure(5)
semilogy(R,G,'LineWidth',1.2,'Color',[0.6350 0.0780 0.1840])
axis([R(1) N/2 1e-25 1e+10])
text(75,1e-10,'t = 200','BackgroundColor','white')
text(140,1e-10,'t = 400','BackgroundColor','white')
text(200,1e-10,'t = 600','BackgroundColor','white')
text(265,1e-10,'t = 800','BackgroundColor','white')
text(320,1e-10,'t = 1000','BackgroundColor','white')
text(380,1e-10,'t = 1200','BackgroundColor','white')
text(440,1e-10,'t = 1400','BackgroundColor','white')
text(470,1e-03,'t = 1600','BackgroundColor','white')
text(450,1e+04,'t = 1800','BackgroundColor','white')
text(470,1e+08,'t = 2000')

grid on
set(gcf,'Position',[784   410   1056   402])
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_0}','FontSize',12)

%saveas(gcf,'/home/javier/Dropbox/sfu/research/nigam-zworski-01/tex/paper/figs/fig27_red0','epsc')

%% RED E_s, s \in \{-1/2, -1, -3/2\} (Figure 6)

figure(6), clf
set(gcf,'Position',[779 498 1066 302])

subplot(1,3,1)
[R,G] = red(-1/2,UF(:,:,[1 4 8 10]),KX,KY);
loglog(R,1e+08*R.^(-1),'--k','LineWidth',1.2)
legend('O(R^{-1})','FontSize',10,'Location','SouthWest')
hold on
loglog(R,G,'LineWidth',1.2,'HandleVisibility','off','Color',[0.6350 0.0780 0.1840])
axis([R(1) 1e+03 1e-20 1e+10])
grid on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_{-1/2}}','FontSize',12)
text(90,1e-18,'t = 200','BackgroundColor','white','Rotation',90)
text(265,1e-18,'t = 800','BackgroundColor','white','Rotation',90)
text(550,1e-17,'t = 1600','BackgroundColor','white','Rotation',90)
text(700,1e-04,'t = 2000','Rotation',90)
xticks([1e1 1e2])

subplot(1,3,2)
[R,G] = red(-1,UF(:,:,[1 4 8 10]),KX,KY);
loglog(R,1e+08*R.^(-2),'--k','LineWidth',1.2)
legend('O(R^{-2})','FontSize',10,'Location','NorthEast')
hold on
loglog(R,G,'LineWidth',1.2,'HandleVisibility','off','Color',[0.6350 0.0780 0.1840])
axis([R(1) 1e+03 1e-20 1e+10])
grid on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_{-1}}','FontSize',12)
text(90,1e-18,'t = 200','BackgroundColor','white','Rotation',90)
text(265,1e-18,'t = 800','BackgroundColor','white','Rotation',90)
text(550,1e-18,'t = 1600','BackgroundColor','white','Rotation',90)
text(700,1e-07,'t = 2000','Rotation',90)
xticks([1e1 1e2])

subplot(1,3,3)
[R,G] = red(-3/2,UF(:,:,[1 4 8 10]),KX,KY);
loglog(R,1e+08*R.^(-3),'--k','LineWidth',1.2)
legend('O(R^{-3})','FontSize',10,'Location','NorthEast')
hold on
loglog(R,G,'LineWidth',1.2,'HandleVisibility','off','Color',[0.6350 0.0780 0.1840])
axis([R(1) 1e+03 1e-20 1e+10])
grid on
set(gca,'FontSize',12)
xlabel('R','FontSize',12)
ylabel('RED {E_{-3/2}}','FontSize',12)
text(90,1e-18,'t = 200','BackgroundColor','white','Rotation',90)
text(265,1e-18,'t = 800','BackgroundColor','white','Rotation',90)
text(550,1e-18,'t = 1600','BackgroundColor','white','Rotation',90)
text(700,1e-10,'t = 2000','Rotation',90)
xticks([1e1 1e2])
