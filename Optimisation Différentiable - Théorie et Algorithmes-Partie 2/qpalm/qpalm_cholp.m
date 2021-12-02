function [L,piv,r] = qpalm_cholp (A,tol)

%%
% [L,piv,r] = qpalm_cholp (A,tol)
%
% Computes the Cholesky factorization of the matrix A with pivotting (to
% ensure stability): P*A*P' = L*L', where P is the permutation matrix
% and L is lower triangular. Only the lower triangular part of A is used
% and the lower triangular part of L is filled in. The matrix A must be
% positive semidefinite, otherwise the factorization fails, obviously,
% and this is highlighted by a negative value of r (the value of the
% pivot).
%
% The permutation matrix P is orthogonal (P*P' = identity) and has the
% same rows and columns as the identity matrix, but in another order.
% This order is specified by the vector piv, which is defined by
% . piv = P*[1:n]'.
% We deduce from this that:
% . P'*piv = [1:n]',
% . the nonzero entries of P are P(i,piv(i))=1, for i = 1:n,
% . the i-th row of P is the piv(i)-th row of the identity matrix,
% . for an n-dimentional vector v, P*v = v(piv) and P'*v(piv) = v,
% . to get w = P*v, write w = v(piv),
% . to get v = P'*w, write v(piv) = w,
% . to get N = P*M*P', write N = M(piv,piv),
% . to get N = P'*M*P, write N(piv,piv) = M.
% If the permutation matrix P1 is represented by the vector piv1 and the
% permutation matrix P2 is represented by the vector piv2, then
% . P1*P2 is represented by piv2(piv1).
%
% The factorization also tries to guess the rank of A and returns it in
% the scalar r. The returned number is highly sensitive to rounding
% errors and monitored by the positive precision threshold tol, which
% aims at measuring a too small number (the value n*eps is recommended
% by Hammarling, Higham et Lucas (2007)). If the first pivot p1 is
% nonpositive or if current pivot is less than tol*p1, the factorization
% is interrupted and the rank r is determined. In that case, L(1:r,1:r)
% is lower triangular and L(:,r+1:n) is zero.
%
% Qpalm_cholp is a straightforward implementation of the algorithm
% described in section 10.3 of the book by Higham (2002).
%
% References
% - Hammarling, Higham, and Lucas (2007). "LAPACK-style codes for
%   pivoted Cholesky and QR updating". In B. Kågström, editor. Applied
%   Parallel Computing - State of the Art in Scientific Computing.
%   Lecture Notes in Computer Science 4699, pages 137–146. Springer,
%   Berlin, Heidelberg.
% - N.J. Higham (2002). "Accuracy and Stability of Numerical
%   Algorithms", second edition. SIAM Publication, Philadelphia.

%-----------------------------------------------------------------------
%
% Author: Jean Charles Gilbert, INRIA.
% Copyright 2014, INRIA.
%
% QPALM is distributed under the terms of the Q Public License version
% 1.0.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
% License version 1.0 for more details.
%
% You should have received a copy of the Q Public License version 1.0
% along with this program.  If not, see
% <http://doc.trolltech.com/3.0/license.html>.
%
%-----------------------------------------------------------------------

% Initialize

  n     = size(A,1);	% matrix order
  L     = zeros(n,n);	% Cholesky factor
  piv   = [1:n]';	% vertical pivot vector
  r     = n;		% the rank is a priori n
  range = 1:n;		% range of the remaining matrix to factor

% Pivoted algorithm described by Higham (2002)

  % Loop

  for k = 1:n-1

    km1 = k-1;

    % Determine the pivot and adapt piv

    [p,I] = max(diag(A(range,range)));	% pivot
    l = km1+min(I); 			% pivot index
    pivk   = piv(k);
    piv(k) = piv(l);
    piv(l) = pivk;

    % Test for exit

    if k == 1
      p1 = p;		% first pivot
      if p <= 0
        r = p;
        return
      end
    elseif p < -tol*p1	% nonpositive definiteness
      r = p;
      return
    elseif p <= tol*p1	% stopping criterion
      r = k-1;
      return
    end

    % Swap

    if k ~= l

      % Swap the rows k/l of L

      if (k > 1) & (km1 > 0)
        range1 = 1:km1;
        auxkm1      = L(k,range1);
        L(k,range1) = L(l,range1);
        L(l,range1) = auxkm1;
      end

      % Swap the rows/columns k/l of A

      auxnmkp1   = A(k,range);
      A(k,range) = A(l,range);
      A(l,range) = auxnmkp1;

      auxnmkp1   = A(range,k);
      A(range,k) = A(range,l);
      A(range,l) = auxnmkp1;

    end

    % Range of the next remaining matrix to factor

    range = k+1:n;

    % Get a new column of L

    p = sqrt(p);	% square root of the pivot
    L(k,k) = p;
    L(range,k) = A(range,k)/p;
    A(range,range) = A(range,range) - L(range,k)*L(range,k)';

  end

  % Consider the last element

  p = A(n,n);
  if n == 1
    if p <= 0
      r = p;
    else
      L(1,1) = sqrt(p);
      r = 1;
    end
  else
    if p < -tol*p1
      r = p;
    elseif p <= tol*p1
      r = n-1;
    else
      L(n,n) = sqrt(p);
      r = n;
    end
  end

% Return

  return
