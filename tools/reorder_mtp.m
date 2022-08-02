function [L,ind] = reorder_mtp(eigenvalues,cutoff)

% Sort by magnitude (i.e. undo previous ordering)
[L,~] = sort(eigenvalues);

% Second sorting: magnitude-then-phase
ind1 = find(abs(real(L)) <= cutoff);
ind2 = find(real(L) > cutoff);
ind3 = find(real(L) < cutoff);
ind = [ind1; ind2; ind3];
L = L(ind);


