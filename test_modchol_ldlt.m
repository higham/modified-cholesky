%test_modchol_ldlt  Quick test of modchol_ldlt.

A = [1 1 0 1
     1 1 1 0
     0 1 1 1
     1 0 1 1]

% A = randn(4); A = A + A';  % Or try random symemtrix A.

A_eigs = eig(A) % Check definitness of A.

[L, D, P, D0, rho] = modchol_ldlt(A) 

residual = norm(P*A*P' - L*D0*L',1)/norm(A,1) % Should be order rho*eps.

A_pert = P'*L*D*L'*P      % Modified matrix: symmetric pos def.
A_pert_eigs = eig(A_pert) % Should all be >= sqrt(eps)*norm(A,'fro').
