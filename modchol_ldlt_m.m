function [L, DMC, P, D, rho] = modchol_ldlt_m(A,delta)
%modchol_ldlt_m  Modified Cholesky algorithm based on LDL' factorization.
%   [L D,P,D0,rho] = modchol_ldlt_m(A,delta) computes a modified
%   Cholesky factorization P*(A + E)*P' = L*D*L', where 
%   P is a permutation matrix, L is unit lower triangular,
%   and D is block diagonal and positive definite with 1-by-1 and 2-by-2 
%   diagonal blocks.  Thus A+E is symmetric positive definite, but E is
%   not explicitly computed.  Also returned is a block diagonal D0 such
%   that P*A*P' = L*D0*L'.  If A is sufficiently positive definite then 
%   E = 0 and D = D0.  Rho is the growth factor for the factorization.
%   The algorithm sets the smallest eigenvalue of D to the tolerance
%   delta, which defaults to sqrt(eps)*norm(A,'fro').
%   The LDL' factorization is compute using a symmetric form of rook 
%   pivoting proposed by Ashcraft, Grimes and Lewis.
%
%   This routine does not exploit symmetry and is not designed to be
%   efficient.  The function modchol_ldlt should generally be used in
%   preference to this one.

%   Reference:
%   S. H. Cheng and N. J. Higham. A modified Cholesky algorithm based
%   on a symmetric indefinite factorization. SIAM J. Matrix Anal. Appl.,
%   19(4):1097-1110, 1998. doi:10.1137/S0895479896302898,

%   Authors: Bobby Cheng and Nick Higham, 1996; revised 2015.

if ~ishermitian(A), error('Must supply symmetric matrix.'), end
if nargin < 2, delta = sqrt(eps)*norm(A,'fro'); end

n = max(size(A));
k = 1;
DMC = eye(n);
D = eye(n);
L = eye(n);
pp = 1:n;
normA = norm(A(:),inf);
rho = normA;

alpha = (1 + sqrt(17))/8;

while k < n
      [lambda, r] = max( abs(A(k+1:n,k)) );
      r = r(1) + k;

      if lambda > 0
          swap = 0;
          if abs(A(k,k)) >= alpha*lambda
              s = 1;
          else
              j = k;
              pivot = 0;
              lambda_j = lambda;
              while ~pivot
                    [temp, r] = max( abs(A(k:n,j)) );
                    r = r(1) + k - 1;
                    temp = A(k:n,r); temp(r-k+1) = 0;
                    lambda_r = max( abs(temp) );
                    if alpha*lambda_r <= abs(A(r,r))
                       pivot = 1;
                       s = 1;
                       A( [k, r], : ) = A( [r, k], : );
                       L( [k, r], : ) = L( [r, k], : );
                       A( :, [k, r] ) = A( :, [r, k] );
                       L( :, [k, r] ) = L( :, [r, k] );
                       pp( [k, r] ) = pp( [r, k] );
                    elseif lambda_j == lambda_r
                       pivot = 1;
                       s = 2;
                       swap = 2;
                       A( [k, j], : ) = A( [j, k], : );
                       L( [k, j], : ) = L( [j, k], : );
                       A( :, [k, j] ) = A( :, [j, k] );
                       L( :, [k, j] ) = L( :, [j, k] );
                       pp( [k, j] ) = pp( [j, k] );
                       k1 = k+1;
                       A( [k1, r], : ) = A( [r, k1], : );
                       L( [k1, r], : ) = L( [r, k1], : );
                       A( :, [k1, r] ) = A( :, [r, k1] );
                       L( :, [k1, r] ) = L( :, [r, k1] );
                       pp( [k1, r] ) = pp( [r, k1] );
                    else
                       j = r;
                       lambda_j = lambda_r;               
                    end
              end
          end
          
          if s == 1
 
             D(k,k) = A(k,k);
             A(k+1:n,k) = A(k+1:n,k)/A(k,k);
             L(k+1:n,k) = A(k+1:n,k);
             i = k+1:n;
             A(i,i) = A(i,i) - A(i,k) * A(k,i);
             A(i,i) = 0.5 * (A(i,i) + A(i,i)');

          elseif s == 2

             E = A(k:k+1,k:k+1);
             D(k:k+1,k:k+1) = E;
             C = A(k+2:n,k:k+1);
             temp = C/E;
             L(k+2:n,k:k+1) = temp;
             A(k+2:n,k+2:n) = A(k+2:n,k+2:n) - temp*C';
             A(k+2:n,k+2:n) = 0.5 * (A(k+2:n,k+2:n) + A(k+2:n,k+2:n)');

          end

          if k+s <= n
             rho = max(rho, max(max(abs(A(k+s:n,k+s:n)))) );
          end
 
      else  % Nothing to do.
     
         s = 1;
         D(k,k) = A(k,k);

      end

      % Modified Cholesky perturbations.
      if s == 1

         if D(k,k) <= delta
            DMC(k,k) = delta;
         else
            DMC(k,k) = D(k,k);
         end
      
      elseif s == 2 
          
         E = D(k:k+1,k:k+1);
         [U,T] = eig(E);
         for ii = 1:2
             if T(ii,ii) <= delta
                T(ii,ii) = delta;
             end
         end
         temp = U*T*U';
         DMC(k:k+1,k:k+1) = (temp + temp')/2;  % Ensure symmetric.

      end

      k = k + s;

      if k == n
         D(n,n) = A(n,n);
         if D(k,k) <= delta
            DMC(k,k) = delta;
         else
            DMC(k,k) = D(k,k);
         end
         break
      end

end

if nargout >= 3, P = eye(n); P = P(pp,:); end
rho = rho/normA;
