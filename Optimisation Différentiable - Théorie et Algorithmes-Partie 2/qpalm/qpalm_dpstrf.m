%%
% function [A,piv,rank,info] = qpalm_dpstrf (A,uplo,tol)
%
% Qpalm_dpstrf is a gateway function to the Lapack routine dpstrf, which
% computes the (possibly singular) Cholesky factorization with complete
% pivoting of a real symmetric positive semidefinite matrix A.
%
% The factorization has the form
%    P'*A*P = U'*U,  if uplo = 'U',
%    P'*A*P = L*L',  if uplo = 'L',
% where U is an upper triangular matrix and L is lower triangular, and P
% is a permutation matrix stored in the vector piv. This algorithm does
% not attempt to check whether A is positive semidefinite.
%
% On entry
% ^^^^^^^^
% A: a positive semidefinite matrix whose Cholesky factor is required;
% uplo: the character 'u' or 'l', which specifies whether the upper or
%   lower part of A must be used and must contain the Cholesky factor on
%   return;
% tol: a tolerance used to terminate the algorithm and the rank of the
%   matrix; if tol < 0, then n*eps*max(diag(A)) will be used, where
%   n=size(A,1).
%
% On return
% ^^^^^^^^^
% A: if info == 0, its upper (if uplo == 'u') or lower (if uplo == 'l')
%   part contains the Cholesky factor of A on entry and the other part
%   is zeroed;
% piv: pivoting vector of length n, which contains a permutation of
%   [1:n] and is such that the nonzero entries of P are P(piv(i),i) = 1,
%   with i in [1:n];
% info: information flag
%   < 0: the (-info)-th argument had an illegal value;
%   = 0: algorithm completed successfully;
%   > 0: the matrix A is either rank deficient with computed rank as
%        returned in 'rank', or is indefinite.
