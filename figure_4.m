% Figure 4. Evolution of the squared $\rm L^2$-norm for several values of 
% $r$ and $\beta(x) = \cos(x_1)+\sin(x_2)$.

close all; clear; clc;
addpath tools;

auxst = dbstack; filename = auxst.name;
filename = [filename,'.mat'];

if ~isfile(filename)
    beta = @(x,y) cos(x)+sin(y);

    fun = @(x,y) -5*exp(-((x+0.0).^2 + (y+0.0).^2)*3+1i*2*x + 1i*y);
    u0 = @(x,y) 0*x+0*y;
    P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2);

    N = 2^8; dt = 0.05; T = 2000; omega0 = 0; snaps = 10; 
    r = [0.2 0.4 0.6 0.8 1.0 1.4 1.6 2.0];
    
    norms = zeros(length(r),snaps+1);
    
    for j = 1:length(r)
        [UF,KX,KY,t] = etdrk4fft2(P,beta,r(j),omega0,N,dt,T,fun,u0,snaps);
        disp(['r = ',num2str(r(j)),' computed'])
        for l = 1:snaps+1
            norms(j,l) = 4*pi^2*norm(UF(:,:,l),'fro')^2/N^4;
        end
    end
    
    % Save as figure_4.mat
    save(filename)
else
    disp(['File exists. Loading ',filename])
    load(filename)
end

%% Plot
figure(4)
plot(t,norms(1,:),'-.o','LineWidth',1.5)
hold on
plot(t,norms(2,:),'-.s','LineWidth',1.5)
plot(t,norms(3,:),'-.d','LineWidth',1.5)
plot(t,norms(4,:),'-.v','LineWidth',1.5)
plot(t,norms(5,:),'--o','LineWidth',1.5)
plot(t,norms(6,:),'--s','LineWidth',1.5)
plot(t,norms(7,:),'--d','LineWidth',1.5)
plot(t,norms(8,:),'--vk','LineWidth',1.5)

legend('r = 0.2','r = 0.4','r = 0.6', 'r = 0.8', 'r = 1.0', 'r = 1.4', 'r = 1.6', 'r = 2.0','Location','NorthWest','FontSize',12,'NumColumns',2)
xlabel('T','FontSize',12)
ylabel('Squared {L^2-Norm}','FontSize',12)
set(gca,'FontSize',12)
grid on
set(gcf,'Position',[784   410   528   402])