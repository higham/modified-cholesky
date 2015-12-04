%test_modchol_ldlt  Quick test of modchol_ldlt.

A = [1 1 0 1
     1 1 1 0
     0 1 1 1
     1 0 1 1]

% A = randn(4); A = A + A';  % Or try random symemtrix A.

A_eigs = eig(A) % Check definitness of A.

[L, D, P, D0] = modchol_ldlt(A) 

residual = norm(P*A*P' - L*D0*L',1)/norm(A,1) % Should be order rho*eps.

A_pert = P'*L*D*L'*P      % Modified matrix: symmetric pos def.
A_pert_eigs = eig(A_pert) % Should all be >= sqrt(eps)*norm(A,'fro').


[L1, D1, P1, D01, rho] = modchol_ldlt_m(A);
rel_diffs = [
norm(L-L1,1)/norm(L,1)
norm(D-D1,1)/norm(D,1)
norm(P-P1,1)/norm(P,1)
norm(D0-D01,1)/norm(D0,1)];

fprintf('Max relative difference between matrices computed by\n')
fprintf('modchol_ldlt and modchol_ldlt_m is %g\n', max(rel_diffs))
