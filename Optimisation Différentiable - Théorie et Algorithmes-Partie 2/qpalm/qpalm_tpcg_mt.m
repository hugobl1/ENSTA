function [prec] = qpalm_tpcg_mt (prec)

%%
% [prec] = qpalm_tpcg_mt (prec)
%
% Qpalm_tpcg_mt is a truncated preconditioned conjugate gradient solver
% for the convex piecewise quadratic optimization problem
%
%   minimize (in x)  phi(x)
%   subject to       lbx <= x <= ubx,
%
% where phi(x) is the minimum value of augmented Lagrangian
% lr(x,y,lmi,lme) for y in [lbi,ubi]. This function phi is piecewise
% quadratic and differentiable; its value and gradient are computed by
% qpalm_phi.
%
% Qpalm_tpcg_mt uses preconditioned conjugate gradient iterations to
% explore a face of [lbx,ubx]. It stops as soon as a solution has been
% found or a bound has been hit in x or a boundary of a quadratic piece
% of phi is encountered (or equivalently, the updated y hits a boundary
% of [lbi,ubi]). In the two latter case, the corresponding components of
% x and/or y are fixed to the boundary values (this is indicated by an
% update of the sets QPALM_IB and QPALM_II).
%
% On entry
% ^^^^^^^^

% Global variables used
% ^^^^^^^^^^^^^^^^^^^^^
% QPALM_info: structure giving various information; those that are used
%   in this function are the following
%
%   QPALM_info.cgit = total number of conjugate-gradient iterations so
%     far, a value that is updated during the run;
%   QPALM_info.phi = value of phi at the final point;
%   QPALM_info.flag (int8) specifies the status of the TCG solver on
%     return:
%     == QPALM_values.stop_on_nonconvexity: the Hessian H probably has
%          negative eigenvalues,
%     == QPALM_values.a_bound_is_hit: the CG iterations encountered a
%          bound before having reached optimality, in which case the
%          returned iterate is the one at the bound (the index of the
%          bound is not returnedm since thater may be more than one such
%          limiting bounds),
%     == QPALM_values.stop_on_max_cgit: max CG iterations reached,
%     == QPALM_values.stop_on_unboundedness: the problem looks unbounded
%          (see QPALM_info.ubdd),
%     == QPALM_values.relative_minimum: a minimum in the current active
%          face has been found up to the required accuracy;
%   QPALM_info.rq: Raylieght quotient in a direction of negative
%     curvature (valid in case QPALM_info.flag ==
%     QPALM_values.stop_on_nonconvexity |
%     QPALM_values.stop_on_unboundedness);
%   QPALM_info.rqmin, QPALM_info.rqmax: minimum and maximum Rayleigh
%     quotients, which are updated during the run (this is usefull in
%     case the CG solve is done in several phases).
%   QPALM_info.ubdd = (meaningful if QPALM_info.flag ==
%     QPALM_values.stop_on_unboundedness) this is a direction in x along
%     which the QP goes to go to -Inf; the direction is feasible for the
%     constraints and normalized, i.e., |QPALM_info.ubdd| = 1; to
%     improve the feasibility of the direction, decrease the value of
%     QPALM_options.ubtol.

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

  global QPALM_n QPALM_nb QPALM_mi QPALM_me
  global QPALM_g QPALM_H QPALM_AI QPALM_AE
  global QPALM_x QPALM_y QPALM_lmi QPALM_lme
  global QPALM_info QPALM_options QPALM_values
  global QPALM_IB QPALM_II

% Just return if there is no free variables

  if QPALM_nb & all(QPALM_IB ~= 0); return; end

% Parameters

  mt1 = 1.e-2;	% parameter intervening in the sufficient decresed condition

% Initializations

  QPALM_info.flag = QPALM_values.relative_minimum;	% will be used for diagonosis inside, so reset it

  fout   = QPALM_options.fout;				% output channel
  verb   = QPALM_options.verb;				% verbosity level
  rntols = QPALM_options.gtol^2;			% residual norm tolerance squared

  if QPALM_mi+QPALM_me; r = QPALM_info.rlag; end

  iter = 0;						% local iteration counter
  iter_max = QPALM_options.cgit - QPALM_info.cgit;	% max iterations allowed

  d = zeros(QPALM_n,1);					% displacement in x
  p = zeros(QPALM_n,1);					% displacement in the residual
  if QPALM_mi
    q = zeros(QPALM_mi,1);				% displacement in y
  else
    q = [];
  end

  alphamax = inf;

  decmax = 0;				% maximal positive decrease obtained so far

  % unboundedness_direction is updated by qpalm_prmk when this one is called; in any case, ~unboundedness_direction
  % is used to decide whether making CG iterations is appropriate

% unboundedness_direction = false;

% Initial bound activity

  if QPALM_nb

    Il  = find(QPALM_IB == -1);		% indices of the variables at the lower bound
    Iu  = find(QPALM_IB == +1);		% indices of the variables at the upper bound
%   Il  = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables at the lower bound
%   Iu  = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables at the upper bound
    WB  = union(Il,Iu)';		% indices of the bound variables
    VB  = setdiff([1:QPALM_n]',WB);	% indices of the free variables (vertical vector)

    QPALM_info.wb = length(WB);		% nb of bound variables

    if verb >= 3; qpalm_print_activb (fout,verb); end

  else

    VB  = [1:QPALM_n]';
    QPALM_info.wb = 0;			% nb of bound variables

  end

  if QPALM_mi
    WI = find(QPALM_II ~= 0); WI = WI(:);
    VI = setdiff([1:QPALM_mi]',WI); VI = VI(:);
    if verb >= 3; qpalm_print_activi(fout,verb); end
  else
    WI = [];
    VI = [];
  end

% Make the preconditioner. Detect indefiniteness or a direction of unboundedness if appropriate

  if (QPALM_options.cgprec ~= QPALM_values.cg_prec_none) & (QPALM_info.ccgp == 1)
    prec = qpalm_prmk (VB,WI);
%VB
%L=prec.L
%rank=prec.rank
%piv=prec.piv
%M = QPALM_H(VB,VB);
%if QPALM_mi; M = M + r*(QPALM_AI(WI,VB)'*(QPALM_AI(WI,VB))); end
%if QPALM_me; M = M + r*(QPALM_AE(:,VB)'*(QPALM_AE(:,VB))); end
%MF = prec.L*prec.L';
%MF(prec.piv,prec.piv) = MF;
%if norm(M-MF)>1.e-10; keyboard; end
    if prec.rank < 0
      QPALM_info.flag = QPALM_values.stop_on_nonconvexity;
      if verb > 0, fprintf(fout,'\n\n### qpalm_tpcg_mt: the QP is nonconvex\n'); end
      QPALM_info.rq = prec.rank;
      return
    end
  end

% Initial phi (= lr since y = yp), the residual, and its squared norm 'rns'

  [phi,gphi] = qpalm_al ();
  rns = gphi(VB)'*gphi(VB);
  rni = log10(rns);

% Initial printings

  if verb >= 3
    fprintf(fout,'\n  > Variables %i',QPALM_n-QPALM_info.wb);
    if QPALM_options.cgprec == QPALM_values.cg_prec_chol; fprintf(fout,', matrix rank %i',prec.rank); end
    fprintf(fout,', max iterations %i',iter_max);
    fprintf(fout,', residual tol %9.3e',QPALM_options.gtol);
    fprintf(fout,'\n  > Iter    |d|    Rayleigh  stepsize          AL            |res|');
    fprintf(fout,'\n  >    0                               %18.11e  %9.3e',phi,sqrt(rns));
  end

%------------------
% CG iteration loop
%------------------

  while 1

    % Stop on residual

    if rns <= rntols
      QPALM_info.flag = QPALM_values.relative_minimum;
      if verb >= 3, fprintf(fout,'\n  > stop on convergence'); end
      break
    end

    % Increase iteration counter

    if iter >= iter_max
      QPALM_info.flag = QPALM_values.stop_on_max_cgit;
      if verb >= 3, fprintf(fout,'\n  > too many CG iterations'); end
      break
    end
    iter = iter+1;
    if verb >= 3, fprintf(fout,'\n  > %4i',iter); end

    % Displacement d(VB) in x

    if iter == 1
      if QPALM_options.cgprec == QPALM_values.cg_prec_none
        d(VB) = -gphi(VB);
      else
        [pg,rns,zc_direction] = qpalm_prec (prec,gphi(VB));
        d(VB) = -pg;
      end
    else
      if QPALM_options.cgprec == QPALM_values.cg_prec_none
        d(VB) = -gphi(VB) + (rns/rns0)*d(VB);
      else
        [pg,rns,zc_direction] = qpalm_prec (prec,gphi(VB));
        if zc_direction
          d(VB) = -pg;
        else
          d(VB) = -pg + (rns/rns0)*d(VB);
        end
      end
    end

    % Displacement p(VB) in the residual

    if (QPALM_options.cgprec ~= QPALM_values.cg_prec_none) & zc_direction
      p(VB) = 0;	% the computation in the else-part should give the same result
    else
      p(VB) = QPALM_H(VB,VB)*d(VB);
      if QPALM_mi; p(VB) = p(VB) + r*(QPALM_AI(WI,VB)'*(QPALM_AI(WI,VB)*d(VB))); end
      if QPALM_me; p(VB) = p(VB) + r*(QPALM_AE(:,VB)'*(QPALM_AE(:,VB)*d(VB))); end
    end

    % Displacement q in y (when x changes, only the components VI of y change)

    if QPALM_mi
      q(VI) = QPALM_AI(VI,VB)*d(VB);
    end

%   % Max stepsize due to the bounds in x and y. Put in (IIl, IIu, IBl, IBu) the indices of the hit bounds

%   if QPALM_nb+QPALM_mi
%     [alphamax,IBl,IBu,IIl,IIu] = qpalm_maxstep (d,q,VB,VI);
%   end

    % Rayleigh quotient and detection of nonconvexity

    dns = d(VB)'*d(VB);	% d norm squared
    pd  = p(VB)'*d(VB);	% product saved since used three times
    rq  = pd/dns;	% Rayleigh quotient
    if verb >= 3, fprintf(fout,'  %7.1e  %8.1e',sqrt(dns),rq); end
    if rq < 0
      QPALM_info.flag = QPALM_values.stop_on_nonconvexity;
      if verb > 0, fprintf(fout,'\n\n### qpalm_tpcg_mt: the QP is nonconvex\n'); end
      QPALM_info.rq = rq;
      break
    end
    if (QPALM_options.cgprec == QPALM_values.cg_prec_none) | ~zc_direction
      QPALM_info.rqmin = min(QPALM_info.rqmin,rq);	% min Rayleigh quotient encountered
      QPALM_info.rqmax = max(QPALM_info.rqmax,rq);	% max Rayleigh quotient encountered
      if (verb >= 5) & (QPALM_options.cgprec ~= QPALM_values.cg_prec_none)	% if zc_direction, dt(below) = 0
        dt = qpalm_prde (prec,d(VB));
        rqp = pd/(dt'*d(VB));
        QPALM_info.rqpmin = min(QPALM_info.rqpmin,rqp);	% min Rayleigh quotient of the preconditioned system encountered
        QPALM_info.rqpmax = max(QPALM_info.rqpmax,rqp);	% max Rayleigh quotient of the preconditioned system encountered
      end
    end

    % Unboundedness

    if (alphamax == inf) & (rq <= QPALM_options.ubtol)
      QPALM_info.flag = QPALM_values.stop_on_unboundedness;
      if verb >= 3, fprintf(fout,'\n  > Unbounded problem'); end
      QPALM_info.ubdd = d/sqrt(dns);
      QPALM_info.rq = rq;
      break
    end

    % Stepsize

    if (QPALM_options.cgprec ~= QPALM_values.cg_prec_none) & zc_direction

      alpha = alphamax;
      QPALM_info.flag = QPALM_values.a_bound_is_hit;

    else

      alpha = rns/pd;

      if alpha > alphamax		% a bound is hit
        alpha = alphamax;
        QPALM_info.flag = QPALM_values.a_bound_is_hit;
      end

    end

    % New phi

    dec = alpha*(gphi'*d) + (alpha^2/2)*pd;
    phi = phi + dec;

    % New x, y, and residual

    QPALM_x(VB) = QPALM_x(VB) + alpha*d(VB);
    if QPALM_mi; QPALM_y(VI) = QPALM_y(VI) + alpha*q(VI); end
    gphi(VB) = gphi(VB) + alpha*p(VB);
    rns0 = rns;
    rns = gphi(VB)'*gphi(VB);

    % Printings

    if verb >= 3
      fprintf(fout,'  %8.2e  %18.11e  %9.3e',alpha,phi,sqrt(rns));
      if (QPALM_options.cgprec ~= QPALM_values.cg_prec_none) & zc_direction; fprintf(fout,'  zero curvature direction'); end
    end

    % Test whether the CG should be interrupted

    dec    = -dec;			% positive decrease
    test1  = (dec <= mt1*decmax);
    decmax = max(decmax,dec);		% update decmax
stop
    % If a bound has been hit, update the preconditioner, WI, VI, and break

    if QPALM_info.flag == QPALM_values.a_bound_is_hit

      Bhit = union(IBl,IBu);
      Ihit = union(IIl,IIu);

      % update the preconditioner if appropriate

      if QPALM_options.cgprec ~= QPALM_values.cg_prec_none
        prec = qpalm_prup (prec,VB,Bhit,WI,Ihit);
      end

      % update WB, VB, and QPALM_IB if appropriate (TODO: complete Ibl and IBu considering dxmin)

      if QPALM_nb & ~isempty(Bhit)
        if verb >= 3; fprintf(fout,'\n  > variable bounds hit:'); end
        if ~isempty(IBl)
          if verb >= 3; fprintf(fout,' %0i%',-IBl); end
          WB = union(WB,IBl);
          VB = setdiff(VB,IBl);
          QPALM_IB(IBl) = -1;
        end
        if ~isempty(IBu)
          if verb >= 3; fprintf(fout,' %0i%',IBu); end
          WB = union(WB,IBu);
          VB = setdiff(VB,IBu);
          QPALM_IB(IBu) = +1;
        end
      end

      % update WI, VI, and QPALM_II if appropriate

      if QPALM_mi & ~isempty(Ihit)
        if verb >= 3; fprintf(fout,'\n  > Inequality bounds hit:'); end
        if ~isempty(IIl)
          if verb >= 3; fprintf(fout,' %0i%',-IIl); end
          WI = union(WI,IIl);
          VI = setdiff(VI,IIl);
          QPALM_II(IIl) = -1;
        end
        if ~isempty(IIu)
          if verb >= 3; fprintf(fout,' %0i%',IIu); end
          WI = union(WI,IIu);
          VI = setdiff(VI,IIu);
          QPALM_II(IIu) = +1;
        end
      end

      % break

      break;

    end

  end

% Before returning

  QPALM_info.cgit = QPALM_info.cgit + iter;
  QPALM_info.phi  = phi;
  QPALM_info.glan = sqrt(rns);

  % compute the exact residual norm and estimate the precision measure

  if QPALM_info.flag == QPALM_values.relative_minimum
    [phi,gphi] = qpalm_al ();
    rne = gphi(VB)'*gphi(VB);
    QPALM_info.gopt = sqrt(rne);	% useful in case there is no GP phase next
    if verb >= 3, fprintf(fout,'\n  > |exact residual| = %8.2e',QPALM_info.gopt); end
    if (rne == 0) & (rns == 0)
      diff = 0;
    else
      rne = log10(rne);
      rn = log10(rns);
      diff = rn-rni;
    end
    if diff
      accu = (rne-rni)/diff;
      QPALM_info.accu = max(QPALM_info.accu,accu);	% catch the best accuracy obtained when a solution has been found
    else
      QPALM_info.accu = 1;
    end
    if verb >= 3, fprintf(fout,', accuracy = %9.2e',QPALM_info.accu); end
  end

  if (verb >= 3) & (iter > 0) & ((QPALM_options.cgprec == QPALM_values.cg_prec_none) | ~zc_direction)
    fprintf(fout,'\n  > Estimated spectrum and condition number: [%8.2e,%8.2e] and %8.2e',QPALM_info.rqmin,QPALM_info.rqmax,QPALM_info.rqmax/QPALM_info.rqmin);
    if (verb >= 5) & (QPALM_options.cgprec ~= QPALM_values.cg_prec_none)
      fprintf(fout,'\n  > ......... for the preconditioned system: [%8.2e,%8.2e] and %8.2e',QPALM_info.rqpmin,QPALM_info.rqpmax,QPALM_info.rqpmax/QPALM_info.rqpmin);
    end
  end

% Exit

  return
