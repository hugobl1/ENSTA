function prec = qpalm_prmk (VB,WI)

%%
% prec = qpalm_prmk (VB,WI)
%
% Qpalm_prmk is a preconditioner maker for the CG iterations. It returns
% in 'prec' a vector or a structure, depending on the type of
% preconditioner specified by options.cgprec.
%
% On entry
% ^^^^^^^^
% VB: indices of the free variables x;
% WI: indices of the free variables y;
%
% On return
% ^^^^^^^^^
% prec: variable or matrix that depends on the type of preconditioner
%   specified by options.cgprec. The following matrix is used:
%
%     M := H(VB,VB)+r*AE(:,VB)'*AE(:,VB)+r*AI(WI,VB)'*AI(WI,VB);
%
%   It is assumed, however, that, in any case, prec.rank < means that
%   the matrix M is not positive semidefinite.
%
%   if options.cgprec == 'diag': prec is a structure with the following
%     fields:
%     - prec.diag is the diagonal of M;
%     - prec.rank is the estimated rank of M; if <0, the matrix M is not
%       positive semidefinite.
%   if options.cgprec == 'chol': prec is a structure with the following
%     fields:
%     - prec.L is the (possibly singular) lower triangular matrix L such
%       that L*L' is
%       . M when M is positive definite and
%       . P*M*P' with pivoting matrix P specified by the vector prec.piv
%         otherwise;
%     - prec.piv is empty when M is positive definite and is a pivoting
%       vector representing the pivoting matrix P otherwise;
%     - prec.rank is the estimated rank of M, if <0, the matrix M is not
%       positive semidefinite.

% Global variables used
% ^^^^^^^^^^^^^^^^^^^^^
% QPALM_info: structure giving various information; those that are used in
%   this function are the following
%
%   QPALM_info.rlag: value of the augmentation parameter;

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

  global QPALM_mi QPALM_me
  global QPALM_H QPALM_AI QPALM_AE
  global QPALM_info QPALM_options QPALM_values

%-----------------------------------
% Initializations
%-----------------------------------

  fout = QPALM_options.fout;		% output channel
  verb = QPALM_options.verb;		% verbosity level

  r    = QPALM_info.rlag;

  prec.diag = [];

%-----------------------------------
% Diagonal preconditioner
%-----------------------------------

  if QPALM_options.cgprec == QPALM_values.cg_prec_diag

    if verb >= 3; fprintf(fout,'\n  > Recomputing the diagonal preconditionner'); end

    % make prec.diag

    prec.diag = diag(QPALM_H(VB,VB));
    if QPALM_me
      prec.diag = prec.diag + r*sum(QPALM_AE(:,VB).*QPALM_AE(:,VB),1)';
    end
    if QPALM_mi
      prec.diag = prec.diag + r*sum(QPALM_AI(WI,VB).*QPALM_AI(WI,VB),1)';
    end

    % detect singularity or indefiniteness

    if any(prec.diag < 0)	% detect indefiniteness
      prec.rank = min(prec.diag);
    else
      prec.rank = sum(prec.diag >= QPALM_options.ubtol);	% estimate the rank
    end

%-----------------------------------
% Perfect Cholesky preconditioner
%-----------------------------------

  elseif QPALM_options.cgprec == QPALM_values.cg_prec_chol

    if verb >= 3; fprintf(fout,'\n  > Cholesky factorization'); end

    % matrix M of the LS to solve

    M = QPALM_H(VB,VB);
    if QPALM_me
      M = M + r*(QPALM_AE(:,VB)'*QPALM_AE(:,VB));
    end
    if QPALM_mi
      M = M + r*(QPALM_AI(WI,VB)'*QPALM_AI(WI,VB));
    end

    % make prec

    [prec.L,prec.piv,prec.rank,info_dpstrf] = qpalm_dpstrf (M,'l',QPALM_options.ubtol);

%> OLD STUFF
%>    [prec.L,p] = chol(M,'lower');
%>
%>    if ~p	% M is positive definite
%>
%>      prec.piv  = [];
%>      prec.rank = size(M,1);
%>
%>    else	% M is probably not positive definite
%>
%>      % Try using 'qpalm_cholp', which is a personnal Matlab implementation of a Cholesky factorization with pivoting (such a
%>      % factorization does not exist in Matlab, 2012). Deduce the solution if it finds that the rank of M is n; or compute a
%>      % direction of unboundedness, otherwise.
%>
%>      [prec.L,prec.piv,prec.rank] = qpalm_cholp (M,QPALM_options.ubtol);	% Cholesky factorization of P*M*P' = L*L' with pivoting matrix P
%>
%>    end

  end

% Exit

  % To check that the factorization is correct, execute the following statements. The error matrix M-MM should be close to zero

%>   MM = prec.L*prec.L';		% MM = L*L'
%>   if ~isempty(prec.piv); MM(prec.piv,prec.piv) = MM; end	% MM = P'*L*L'*P
%>   norm(M-MM)
%>   keyboard

  return
