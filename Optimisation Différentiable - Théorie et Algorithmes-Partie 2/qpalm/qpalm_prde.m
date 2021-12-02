function [dv] = qpalm_prde (prec,v)

%%
% [dv] = qpalm_prde (prec,v)
%
% Qpalm_prec makes the reverse operation than the one in qpalm_prec: it
% de-preconditions or deconditions the vector v, using the
% preconditioner described in 'prec'. This operation is only used for
% realizing some expensive printings.
%
% On entry
% ^^^^^^^^
% prec: structure representing the preconditioning matrix P, defined by
%   the function qpalm_prmk, and depending on the value of
%   'options.cgprec'
% v: vector (of appropriate dimension) to decondition
%
% On return
% ^^^^^^^^^
% dv: P\v, where P is the preconditioner, hence a matrix that would like
%   to resemble the inverse of the matrix of the linear system to solve.

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
% Diagonal preconditioner
%-----------------------------------

  if QPALM_options.cgprec == QPALM_values.cg_prec_diag

    dv  = prec.diag.*v;		% dv = deconditioned direction

%-----------------------------------
% Perfect Cholesky preconditioner
%-----------------------------------

  elseif QPALM_options.cgprec == QPALM_values.cg_prec_chol


    if isempty(prec.piv)	% the Cholesky factorization uses no pivoting matrix

      dv = prec.L*(prec.L'*v);	% dv = L*L'*v

    else			% the Cholesky factorization uses the pivoting vector prec.piv; whatever is the rank of M, one
				% computes dv = M*v using the Cholesky factor L such that M = P'*L*L'*P

        dv = v(prec.piv);		% dv = P*v
        dv = prec.L*(prec.L'*dv);	% dv = L*L'*P*v
        dv(prec.piv) = dv; 		% dv = P'*L*L'*P*v

    end

  end

% Exit

  return
