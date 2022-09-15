% Long-term evolution. Plots in real and Fourier space

close all; clear; clc; 

auxst = dbstack; filename = auxst.name;
filename = [filename,'.mat'];

if ~isfile(filename)
    fun = @(x,y) -5*exp(-((x+0.0).^2 + (y+0.0).^2)*1+1i*2*x + 1i*y); % Centred Gaussian
    omega0 = 0;
    u0 = @(x,y) 0*x+0*y;
    P = @(kx,ky) ky./sqrt(1+kx.^2+ky.^2);

    % Test 1
    r = 0.5; beta = @(x,y) cos(x) + 0*y;    
    N = 512; dt = 0.5; T = 2000; snaps = 1;
    [UF1,KX1,KY1,tsteps1] = etdrk4fft2(P,beta,r,omega0,N,dt,T,fun,u0,snaps);
    [U1,X1,Y1] = fourier2real(UF1);
    [UF1,KX1,KY1] = makeitperiodic(UF1,KX1,KY1);
    [U1,X1,Y1] = makeitperiodic(U1,X1,Y1);
    disp('Test 1 completed.')

    % Test 2
    r = 0.45; beta = @(x,y) cos(x-2*y) + sin(2*y);
    N = 1024; dt = 0.5; T = 1000; snaps = 1;
    [UF2,KX2,KY2,tsteps2] = etdrk4fft2(P,beta,r,omega0,N,dt,T,fun,u0,snaps);
    [U2,X2,Y2] = fourier2real(UF2);
    [UF2,KX2,KY2] = makeitperiodic(UF2,KX2,KY2);
    [U2,X2,Y2] = makeitperiodic(U2,X2,Y2);
    disp('Test 2 completed.')

    % Test 3
    r = 0.55; beta = @(x,y) cos(x-2*y) + sin(2*y);
    N = 1024; dt = 0.5; T = 1000; snaps = 1;
    [UF3,KX3,KY3,tsteps3] = etdrk4fft2(P,beta,r,omega0,N,dt,T,fun,u0,snaps);
    [U3,X3,Y3] = fourier2real(UF3);
    [UF3,KX3,KY3] = makeitperiodic(UF3,KX3,KY3);
    [U3,X3,Y3] = makeitperiodic(U3,X3,Y3);
    disp('Test 3 completed.')
    
    save(filename);
else
    disp(['File exists. Loading ',filename])
    load(filename);
end

%% Plot in real space
figure(7)
set(gcf, 'Position', [665 549 1201 270])

subplot(1,3,1)
surf(X1,Y1,abs(U1(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('x_1'); ylabel('x_2'); view([0 90])
colorbar; axis tight;
title('$|u_1|$ ($t=2000$)', 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,2)
surf(X2,Y2,abs(U2(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('x_1'); ylabel('x_2'); view([0 90])
colorbar; axis tight; 
title('$|u_2|$ ($t=1000$)', 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,3)
surf(X3,Y3,abs(U3(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('x_1'); ylabel('x_2'); view([0 90])
colorbar; axis tight;
title('$|u_3|$ ($t=1000$)', 'Interpreter', 'latex', 'FontSize', 12)

%% Plot in Fourier space
%{
figure(17)
set(gcf, 'Position', [665 549 1201 270])

subplot(1,3,1)
surf(KX1,KY1,abs(UF1(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('k_1'); ylabel('k_2'); view([0 90])
colorbar; axis tight; set(gca,'ColorScale','log')
title('$|\widehat{u_1}|$ ($t=2000$)', 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,2)
surf(KX2,KY2,abs(UF2(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('k_1'); ylabel('k_2'); view([0 90])
colorbar; axis tight; set(gca,'ColorScale','log')
title('$|\widehat{u_2}|$ ($t=1000$)', 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,3)
surf(KX3,KY3,abs(UF3(:,:,2))),
shading interp, colormap(brewermap([],'*RdBu'));
xlabel('k_1'); ylabel('k_2'); view([0 90])
colorbar; axis tight; set(gca,'ColorScale','log')
title('$|\widehat{u_3}|$ ($t=1000$)', 'Interpreter', 'latex', 'FontSize', 12)
%}

