function [v,err_proj] = l2projection(uf,eff,efr)

% uf and ef are in Fourier space, ef has to be fftshifted
uf = fftshift(uf);
for l = 1:size(eff,3)
    eff(:,:,l) = fftshift(eff(:,:,l));
end

m = log2(size(eff,2)); n = log2(size(uf,1));
if m < n
    % Reduce u
    ind = 2^(n-m)*(0:2^m-1)+1;
    uf = uf(ind,ind);
elseif m > n
    % Reduce ef
    ind = 2^(m-n)*(0:2^n-1)+1;
    eff = eff(ind,ind,:);
end

l = size(eff,3);
b = zeros(l,1);
A = zeros(l,l);
for j = 1:l
    b(j) = sum(uf.*conj(eff(:,:,j)),'all');
    for k = 1:l
        A(k,j) = sum(eff(:,:,j).*conj(eff(:,:,k)),'all');
    end
end

disp('System to solve:')
disp(abs(A));
disp(b);

coef = A\b
v = zeros(size(efr,1),size(efr,2));
for k = 1:l
    v = v + coef(k)*efr(:,:,k);
end

vf = zeros(size(eff,1),size(eff,2));
for k = 1:l
    vf = vf + coef(k)*eff(:,:,k);
end
% PROJECTION ERROR
N = size(uf,2); 

err_proj = abs(uf-vf).^2;
err_proj = (4*pi^2/N^4)*sum(err_proj,'all');
err_proj = sqrt(err_proj);