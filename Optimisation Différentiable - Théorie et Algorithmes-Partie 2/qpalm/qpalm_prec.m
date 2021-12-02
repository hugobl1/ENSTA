function [pv,vpv,flag] = qpalm_prec (prec,v)

%%
% [pv,vpv,flag] = qpalm_prec (prec,v)
%
% When flag == 0, qpalm_prec returns P*v in 'pv', where P is a
% preconditioned matrix represented by the vector or structure 'prec'
% and 'v' is a vector. When flag == 1, a zero curvature direction has
% been found and qpalm_prec returns that direction in 'pv'.
%
% On entry
% ^^^^^^^^
% prec: structure representing the preconditioning matrix P, defined by
%   the function qpalm_prmk, and depending on the value of
%   'options.cgprec'
% v: vector (of appropriate dimension) to precondition
%
% On return
% ^^^^^^^^^
% pv: P*v (if flag == 0) or a direction of zero curvature satisfying
%   v'*pv >= 0 (if flag == 1);
% vpv: v'*P*v (if flag == 0) or [] (if flag == 1);
% flag: tells whether a zero curvature direction has been detected
%   == 0: no zero curvature direction detected
%   == 1: a zero curvature direction has been detected.

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

  global QPALM_options QPALM_values

%-----------------------------------
% Initialization
%-----------------------------------

  nv = length(v);

%-----------------------------------
% Diagonal preconditioner
%-----------------------------------

  if QPALM_options.cgprec == QPALM_values.cg_prec_diag

    pmin = min(prec.diag);		% detect the minimum curvature directions
    imin = find(prec.diag <= pmin);	% does not work with '[pmin,imin] = min(prec.diag)'
    if pmin > QPALM_options.ubtol	% then precondition all the components
      pv  = v./prec.diag;		% pv = preconditioned direction
      vpv = v'*pv;
      flag = 0;
    else				% e^imin are directions of zero curvature
      pv = zeros(nv,1);
      if all(abs(v(imin)) <= QPALM_options.ubtol)	% then ignore the components in 'imin' in the preconditioning
        pv(imin) = 0;					% impose to work in the complementary space
        cimin = setdiff([1:nv],imin);			% complement of imin
        pv(cimin)  = v(cimin)./prec.diag(cimin);	% precondition only in vect{e^cimin}
        vpv = v(cimin)'*pv(cimin);
        flag = 0;
      else						% then return the first zero curvature direction with nonzero slope
        I = imin(find(abs(v(imin)) > QPALM_options.ubtol));	% is nonempty
        imin1 = I(1);
        if v(imin1) > 0
          pv(imin1) = 1;
        else
          pv(imin1) = -1;
        end
        vpv = [];
        flag = 1;
      end
    end

%-----------------------------------
% Cholesky preconditioner
%-----------------------------------

  elseif QPALM_options.cgprec == QPALM_values.cg_prec_chol

    if isempty(prec.piv)	% the Cholesky factorization uses no pivoting matrix

      % Compute pv = M\v using the Cholesky factorization of M = L*L'; L is in prec.L

      opts.LT = true;			% L is lower triangular
      pv = linsolve(prec.L,v,opts);	% pv = L\v
      opts.TRANSA = true;		% use L'
      pv = linsolve(prec.L,pv,opts); 	% pv = L'\(L\v)
      vpv = v'*pv;
      flag = 0;

    else			% the Cholesky factorization uses the pivoting vector prec.piv

      r = prec.rank;

      if r == nv		% the preconditioner is nonsingular

        % Compute pv = M\v using the Cholesky factorization of P*M*P' = L*L'

        pv = v(prec.piv);		% pv = P*v
        pv = pv(:);			% make sure that pv is a column vector
        opts.LT = true;			% use the lower triangular part of L in linsolve
        pv = linsolve(prec.L,pv,opts);	% pv = L\(P*v)
        opts.TRANSA = true;		% use L' instead of L in linsolve
        pv = linsolve(prec.L,pv,opts);	% pv = L'\(L\(P*v))
        pv(prec.piv) = pv; 		% pv = P'*L'\(L\(P*v))) = -M\v
        vpv = v'*pv;
        flag = 0;

      else			% the preconditioner is singular

        % get the two nonzero blocks L11 and L21 of L

        L11 = prec.L(1:r,1:r);
        L21 = prec.L(r+1:nv,1:r);

        % compute Pv = P*v and prepare for checking whether the LS has a solution

        Pv = v(prec.piv);
        opts.LT = true;
        opts.TRANSA = false;
        w1 = linsolve(L11,Pv(1:r),opts);	% w1 = L11\(P*v)(1:r)
        z = L21*w1-Pv(r+1:nv);

        % proceed differently depending on whether the LS has a solution or not

        if norm(z) <= QPALM_options.ubtol	% the LS has a solution, hence compute one solution
          opts.TRANSA = true;
          pv = [linsolve(L11,w1,opts);zeros(nv-r,1)];	% pv = [L11'\(L11\(P*v)(1:r)); 0]
          pv(prec.piv) = pv; 				% pv = P'*[L11'\(L11\(P*v)(1:r)); 0]
          vpv = v'*pv;
          flag = 0;
        else					% the LS has no solution, hence compute a candidate for an unbounded direction (its opposite actually)
          opts.TRANSA = true;
          pv = [linsolve(L11,L21'*z,opts);-z];	% pv = [L11'\(L21'*z);-z]
          pv(prec.piv) = pv;			% pv = P'*[L11'\(L21'*z);-z]
          vpv = [];
          flag = 1;
        end

%       % OLD STUFF
%       % Compute a cheap direction of unboundedness P' * [d1; d2]

%       L11 = prec.L(1:r,1:r);
%       L21 = prec.L(r+1:nv,1:r);
%       opts.LT = true;
%       opts.TRANSA = false;
%       d1 = linsolve(L11,v(1:r),opts);		% d1 = L11\v(1:r)
%       d2 = L21*d1+v(r+1:nv);			% d2 = L21*(L11\v(1:r)) + v(r+1:nv)
%       opts.TRANSA = true;
%       d1 = linsolve(L11,L21'*d2,opts);	% d1 = L11'\(L21'*d2)
%       pv = [d1; d2];				% permuted direction of unboundedness [d1; d2]
%       pv = pv/norm(pv);			% normalized permuted direction of unboundedness [d1; d2]
%       pv(prec.piv) = pv;			% return P'*[d1; d2]
%       vpv = v'*pv;
%       flag = 1;

      end

    end

  end

% Exit

  return
