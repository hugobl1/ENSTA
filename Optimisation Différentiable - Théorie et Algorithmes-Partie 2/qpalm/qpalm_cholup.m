function [L,r,piv] = qpalm_cholup (L,r,v,tol)

%%
% [L,r,piv] = qpalm_cholup (L,r,v,tol)
%
% Computes the Cholesky factor of the rank one correction of L*L', where
% L is a lower triangular square matrix of rank r (the last size(L,2)-r
% columns are supposed to be zero). On return L is the lower triangular
% Cholesky factor of P*(L*L'+v*v')*P', where v is the vector on entry
% and P is a permutation matrix, and r gives the rank of the returned L.
%
% It is assumed (without verifying it) that the r first diagonal
% elements of L (on entry) are nonzero (r is assumed to be nonnegative).
% Hence the rank of the resulting L cannot decrease and can only
% increase by at most one.
%
% If r ~= 0, the new L is computed thanks to r Givens rotations and piv
% is empty.
%
% If r == 0 and v ~= 0, the returned L is zero except its first column
% that contains v or -v. It is only in this latter case that P may
% differ from the identity matrix (i.e. piv is empty): a row pivoting
% may occur to put in L(1,1) the largest |v(i)|. 
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
% On entry
% ^^^^^^^^
% L: lower triangular matrix of rank r (value on entry); it is assumed
%   that the last size(L,2)-r columns of L vanish;
% r: declared rank of L;
% v: vector of length size(L,1);
% tol: positive tolerance for detecting a vanishing column of the
%   resulting L, hence its rank.
%
% On return
% ^^^^^^^^^
% L: Cholesky factor of L*L'+v*v';
% r: rank of the returned matrix L;
% piv: pivoting vector (empty iff there is no pivoting).

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

% A priori, no pivoting

  piv = [];

% Case when v == 0: just return

  if all(abs(v)<=tol); return; end

% Initialization

  n   = size(L,1);
  piv = [];		% a priori no pivoting

% Case when v ~= 0 and r == 0: pivot row 1 and row k where |v(k)| is maximal, next put v or -v in L(:,1)

  if r == 0
    [p,I] = max(abs(v));	% locate pivot
    k = min(I); 		% pivot index
    if k ~= 1			% then pivoting occors
      piv = [1:n]';		% vertical pivot vector
      piv1   = piv(1);
      piv(1) = piv(k);
      piv(k) = piv1;
      v1     = v(1);
      v(1)   = v(k);
      v(k)   = v1;
    end
    if v(1) < 0; v = -v; end
    L(:,1) = v;
    r = 1;
    return
  end

% Case when v ~= 0 and r ~= 0

  L = [v(:) L];

  for i = 1:r
    x = L(i,i);
    y = L(i,i+1);
    aux = norm([x y]);	% nonzero by the assumtion
    c = x/aux;
    s = y/aux;
    L(i:n,i:i+1) = L(i:n,i:i+1)*[c -s; s c];
  end

% See whether the rank has increased

  L(:,n+1) = [];
  if (r<n) & any(abs(L(:,r+1))>tol)
    r = r+1;
  end

% Exit

  return
