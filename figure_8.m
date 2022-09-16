% Figure 8. Energy manifolds $\Sigma_j$, $j \in \{1,2,3\}$ viewed from 
% the side (top row), and from the top (bottom row). This last view 
% reveals holes in $\Sigma_3$, but not in $\Sigma_1$ and $\Sigma_2$.

clear; close all;
addpath tools;

N = 128;
x = 2*pi*(-N/2:N/2)/N;
[X,Y,T] = meshgrid(x);  % Mesh for the 3-torus

% Plot all values such that S = 0.
% Numerically we plot all values such that
% -riv <= S <= riv, wirth riv (range of isovalues)
% a small number close to 0.
riv = 1e-14;    % Range of isovalues
number_iv = 11; % Number of isovalues

figure(8)
set(gcf,'Position',[353 257 1299 739])

subplot(2,3,1)
S = 0.5*cos(X)-sin(T);
for iv = linspace(-riv,riv,number_iv)    
    p = patch(isosurface(X,Y,T,S,iv));
    p.EdgeColor='none';
    p.FaceColor = [5/255,48/255,97/255];
    hold on
end
disp('Plot 1/6 Generated')
axis([-pi pi -pi pi -pi pi])
xlabel('x_1','FontSize',14)
ylabel('x_2','FontSize',14)
zlabel('\eta','FontSize',14)
camlight('right')
lighting gouraud
grid on
view([22.9427,15.9842])
title('$\Sigma_1$','Interpreter','latex','FontSize',16)

subplot(2,3,2)
S = 0.45*(cos(X-2*Y) + sin(2*Y)) - sin(T);
for iv = linspace(-riv,riv,number_iv)
    p = patch(isosurface(X,Y,T,S,iv));
    p.EdgeColor='none';
    p.FaceColor = [5/255,48/255,97/255];
    hold on
end
disp('Plot 2/6 Generated')
axis([-pi pi -pi pi -pi pi])
xlabel('x_1','FontSize',14)
ylabel('x_2','FontSize',14)
zlabel('\eta','FontSize',14)
lightangle(45,45)
lighting gouraud
grid on
view([22.9427,15.9842])
title('$\Sigma_2$','Interpreter','latex','FontSize',16)

subplot(2,3,3)
S = 0.55*(cos(X-2*Y) + sin(2*Y)) - sin(T);
for iv = linspace(-riv,riv,number_iv)
    p = patch(isosurface(X,Y,T,S,iv));
    p.EdgeColor='none';
    p.FaceColor = [5/255,48/255,97/255];
    hold on
end
disp('Plot 3/6 Generated')
axis([-pi pi -pi pi -pi pi])
xlabel('x_1','FontSize',14)
ylabel('x_2','FontSize',14)
zlabel('\eta','FontSize',14)
lightangle(45,45)
grid on
view([22.9427,15.9842])
title('$\Sigma_3$','Interpreter','latex','FontSize',16)

subplot(2,3,4)
S = 0.5*cos(X)-sin(T);
for iv = linspace(-riv,riv,number_iv)
    p = patch(isosurface(X,Y,T,S,iv));
    p.EdgeColor='none';
    p.FaceColor = [5/255,48/255,97/255];
    hold on
end
disp('Plot 4/6 Generated')
axis([-pi pi -pi pi -pi pi])
xlabel('x_1','FontSize',14)
ylabel('x_2','FontSize',14)
zlabel('\eta','FontSize',14)
camlight('right')
lighting gouraud
grid on
view([0,90])
title('$\Sigma_1$','Interpreter','latex','FontSize',16)

subplot(2,3,5)
S = 0.45*(cos(X-2*Y) + sin(2*Y)) - sin(T);
for iv = linspace(-riv,riv,number_iv)
    p = patch(isosurface(X,Y,T,S,iv));
    p.EdgeColor='none';
    p.FaceColor = [5/255,48/255,97/255];
    hold on
end
disp('Plot 5/6 Generated')
axis([-pi pi -pi pi -pi pi])
xlabel('x_1','FontSize',14)
ylabel('x_2','FontSize',14)
zlabel('\eta','FontSize',14)
camlight
lighting gouraud
grid on
view([0,90])
title('$\Sigma_2$','Interpreter','latex','FontSize',16)

subplot(2,3,6)
S = 0.55*(cos(X-2*Y) + sin(2*Y)) - sin(T);
for iv = linspace(-riv,riv,number_iv)
    p = patch(isosurface(X,Y,T,S,iv));
    p.EdgeColor='none';
    p.FaceColor = [5/255,48/255,97/255];
    hold on
end
disp('Plot 6/6 Generated')
axis([-pi pi -pi pi -pi pi])
xlabel('x_1','FontSize',14)
ylabel('x_2','FontSize',14)
zlabel('\eta','FontSize',14)
camlight
lighting gouraud
grid on
view([0,90])
title('$\Sigma_3$','Interpreter','latex','FontSize',16)
