function [] = qpalm_gp ()

%%
% [] = qpalm_gp ()
%
% Qpalm_gp checks optimality of the current augmented Lagrangien
% subproblem, using the projected gradient and, in case optimality is
% not reached, it realizes one step of gradient projection (a GP step)
% for the convex quadratic optimization problem
%
%      minimize (in x)  phi(x)
%      subject to       lbx <= x <= ubx,
%
% where phi(x) is the minimum value of augmented Lagrangian
% lr(x,y,lmi,lme) for y in [lbi,ubi]. This function phi is
% differentiable and its value and gradient are computed by qpalm_phi.
% The GP step consists in a linesearch that approximately minimizes the
% piecewise quadratic nonconvex function h : (alpha>0) ->
% phi(P(x-alpha*gphi)), where x is the entry point, gphi is the gradient
% of phi, and P is the orthogonal projector on [lbx,ubx].
%
% The GP phase also recomputes the bound index set in QPALM_IB and the
% inequality index set in QPALM_II. The latter may change from the one
% obtained at the end of a previous minimization phase, since the y
% computed by the minimization phase may differ from the projetion of
% (AI*x+lmi/r) on [lbi,ubi].
%
% On entry
% ^^^^^^^^
%
% On return
% ^^^^^^^^^
% yp: the projection of (AI*x+lmi/r) on [lbi,ubi];
%
% Global
% ^^^^^^
% QPALM_x: updated primal variable.
% QPALM_y: the projetion of (AI*x+lmi/r) on [lbi,ubi];

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
  global QPALM_g QPALM_H QPALM_lbx QPALM_ubx QPALM_AI QPALM_lbi QPALM_ubi QPALM_AE
  global QPALM_x QPALM_lmi
  global QPALM_info QPALM_options QPALM_values
  global QPALM_IB QPALM_II

% Parameters

  si = 1.e-2;	% interpolation safeguard, must be in ]0,0.5[
  se = 1.e+1;	% extrapolation safeguard, must be > 1

% Initializations

  mie = QPALM_mi+QPALM_me;

  if mie; r = QPALM_info.rlag; end

  fout = QPALM_options.fout;		% output channel
  verb = QPALM_options.verb;		% verbosity level

  QPALM_info.wb = 0;

  IB0 = QPALM_IB;				% remember initial activity to see whether change occurs
  II0 = QPALM_II;				% remember initial activity to see whether change occurs
  QPALM_info.no_change_in_activity = false;	% a priori the activity will change

% Initial activity on x

  if QPALM_nb

    if verb >= 3; qpalm_print_activb (fout,verb); end

    Ilb0 = find(QPALM_IB < 0);
    Iub0 = find(QPALM_IB > 0);

%   Ilb0 = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables at the lower bound
%   Iub0 = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables at the upper bound
%   QPALM_IB       = zeros(QPALM_n,1);
%   QPALM_IB(Ilb0) = -1;
%   QPALM_IB(Iub0) = +1;

  end

% Compute the initial value 'phi0' of phi, its gradient 'gphi' at x, and the intial activity in x

  if QPALM_mi & (verb >= 3); qpalm_print_activi (fout,verb); end

  [phi0,gphi] = qpalm_phi (true);

%-------------------------------------------------------------------------------------------------------------------------------
% Detect optimality at the current x.
%
% It is only here that optimality of the current AL subproblem is detected. Since the method generates no multiplier, a primal
% test is used, namely, it is checked whether the projected gradient 'pgphi' vanishes or, equivalently, whether the gradient
% 'gphi' is in the dual of the tangent cone to the feasible set at x. Recall that x is always feasible for the bound constraints
% and that the projected gradient 'pgphi' of phi at x is the projection of the gradient of phi on the opposite of the tangent
% cone to [lbx,ubx] at x. It is given by
%
%               { min(0,gphi(i))  if lbx(i) = x(i)
%    pgphi(i) = { gphi(i)         if lbx(i) < x(i) < ubx(i)
%               { max(0,gphi(i))  if          x(i) = ubx(i).
%
%-------------------------------------------------------------------------------------------------------------------------------

  % One computes d = -pgphi instead, since it is that direction that is useful below

  d = -gphi;

  if QPALM_nb
    d(Ilb0) = max(0,d(Ilb0));
    d(Iub0) = min(0,d(Iub0));
  end

  % Test convergence

  QPALM_info.gopt = norm(d,inf);		% inf-norm of the projected gradient of phi
  if verb >= 3; fprintf(fout,'\n  > |Projected gradient| = %12.5e',QPALM_info.gopt); end
  if QPALM_info.gopt <= QPALM_options.gtol
    QPALM_info.flag = QPALM_values.success;	% a solution of the AL subproblem has been found
    QPALM_info.phi = phi0;			% return the value of phi
    if QPALM_nb; QPALM_info.wb = length(union(Ilb0,Iub0)); end	% return the nb of active x-bounds
    if verb >= 3; fprintf(fout,' (optimality reached)'); end
    return
  end

%-------------------------------------------------------------------------------------------------------------------------------
% Detect nonconvexity, which occurs when d'*H*d < 0.
%
% This test is only done along the direction d (it implements a sufficient condition, which may not be necessarily satisfied
% when the quadratic function is nonconvex), hence it must be repeated each time qpalm_gp is called.
%-------------------------------------------------------------------------------------------------------------------------------

  curv_init = d'*QPALM_H*d;
  QPALM_info.rq = curv_init/(d'*d);	% Rayleigh quotient

  if curv_init < 0
    QPALM_info.flag = QPALM_values.stop_on_nonconvexity;
    if verb > 0, fprintf(fout,'\n\n### qpalm_gp: the QP is nonconvex\n'); end
    return
  end

%-------------------------------------------------------------------------------------------------------------------------------
% Detect unboundedness, which occurs when
% . H*d = 0,
% . d is in the asymptotic cone of [lbx,ubx],
% . AI*d is in the asymptotic cone of [lbi,ubi],
% . AE*d = 0.
%
% This test is only done along the direction d, hence it must be repeated each time qpalm_gp is called.
%-------------------------------------------------------------------------------------------------------------------------------

  while 1	% not a loop, just used to make exiting easily from the unboundedness detection

    % Initialization

    curv_unbd_tol = QPALM_options.ubtol*(d'*d);
    if QPALM_mi	% also useful for curv_init below
      AId = QPALM_AI*d;
      AIx = QPALM_AI*QPALM_x+QPALM_lmi/r;
    end
    if QPALM_me	% also useful for curv_init below
      AEd = QPALM_AE*d;
      rAEds = r*(AEd'*AEd);
    end

    % Contribution of the cost function to 'curv_unbd': pursue detection if H*d vanishes

    curv_unbd = curv_init;
    if (curv_unbd > curv_unbd_tol); break; end

    % Contribution of the inequality constraints to 'curv_unbd': pursue detection if AId is in the asymptotic cone of [lbi,ubi]

    if QPALM_mi
      v = AId(find(((AId<0)&(QPALM_lbi>-inf))|((AId>0)&(QPALM_ubi<inf))));
      if ~isempty(v)
        curv_unbd = curv_unbd + r*(v'*v);
      end
      if (curv_unbd > curv_unbd_tol); break; end
    end

    % Contribution of the equality constraints to 'curv_unbd': pursue detection if AEd vanishes

    if QPALM_me
      curv_unbd = curv_unbd + rAEds;
      if (curv_unbd > curv_unbd_tol); break; end
    end

    % Here, 'curv_unbd' is small. Hence if d is in the asymptotic cone of [lbx,ubx], then phi is unbounded below (along d
    % actually)

    if all((d<=0)|(QPALM_ubx==inf)) & all((d>=0)|(QPALM_lbx==-inf))
      QPALM_info.flag = QPALM_values.stop_on_unboundedness;
      if verb >= 3, fprintf(fout,'\n  > unbounded problem'); end
      QPALM_info.ubdd = d/norm(d);
      return
    end

    break

  end

%-------------------------------------------------------------------------------------------------------------------------------
% Piecewise linesearch on h : (alpha>0) -> phi(P(x-alpha*gphi)).
%
% We have assumed that this is better than computing the exact minimizer, when the dimension is large, despite h is piecewise
% quadratic. The reason is that there may be many pieces when the dimension is large, while having the exact minimizer has
% litte interest; what matters is to guarantee convergence rapidly.
%-------------------------------------------------------------------------------------------------------------------------------

  % resolution in alpha, deduced from the resolution in x (given by options.dxmin)

  damin = QPALM_options.dxmin/norm(d,inf);

  % Compute the initial stepsize trial, starting by computing the one obtained by the initial curvature.
  % Put in 'curv_init' the initial curvature of phi along d

  if QPALM_mi
    v = AId(find((AIx<QPALM_lbi)|((AIx==QPALM_lbi)&(AId<0))|(AIx>QPALM_ubi)|((AIx==QPALM_ubi)&(AId>0))));
    if ~isempty(v)
      curv_init = curv_init + r*(v'*v);
    end
  end

  if QPALM_me
    curv_init = curv_init + rAEds;
  end

  % Initial stepsize. Here, the curvature 'curv_init' is >= 0 and 'curv_unbd' is > 0. If 'curv_init > 0', an initial stepsize
  % can be computed, which is the one that gives the minimum of the quadratics Q interpolating the conditions Q'(0)=slope0 and
  % Q''(0)=curv_init: alpha = -slope0/curv_init. If 'curv_init == 0', do the same thing using 'curv_unbd' instead of
  % 'curv_init'.

  slope0 = -d'*d;		% slope of h at 0 = -pgphi'*pgphi = -d'*d
  if curv_init > 0
    alpha = -slope0/curv_init;
    alpha_type = 'curv-init';	% to memorize that alpha has been computed by using curv_init
  else
    alpha = -slope0/curv_unbd;
    alpha_type = 'curv_unbd';	% to memorize that alpha has been computed by using curv_unbd
  end

%>   % Setting the initial stepsize trial to the last-hit stepsize is not a good idea (too many stepsize trials without improvement
%>   % in the difficult cases)
%>   alpha_hit = zeros(QPALM_n,1);
%>   I = find((d>0)&(QPALM_ubx<inf));
%>   if ~isempty(I); alpha = max((QPALM_ubx(I)-QPALM_x(I))./d(I)); end
%>   I = find((d<0)&(QPALM_lbx>-inf));
%>   if ~isempty(I); alpha = max(alpha,max((QPALM_lbx(I)-QPALM_x(I))./d(I))); end
%>   alpha_type = 'last_hit';		% to memorize that alpha has been computed by using curv_unbd
%>   if alpha == 0; stop; end

%>  % Initial stepsize. Here, the curvature 'curv_init' is >= 0. If 'curv_init > 0', an initial stepsize can be computed, which is
%>  % the one that gives the minimum of the quadratics Q interpolating the conditions Q'(0)=slope0 and Q''(0)=curv_init: alpha =
%>  % -slope0/curv_init. If 'curv_init == 0', take the largest alpha such that x*alpha*d hits a new bound.
%>
%>  slope0 = -d'*d;		% slope of h at 0 = -pgphi'*pgphi = -d'*d
%>  if curv_init > 0
%>    alpha = -slope0/curv_init;
%>    alpha_type = 'min';		% to memorize that alpha has been computed by quadratic interpolation
%>  else
%>    Iub = find(d>0);
%>    if isempty(Iub)
%>      alphau = 0;
%>    else
%>      alphau = max((QPALM_ubx(Iub)-QPALM_x(Iub))./d(Iub));
%>    end
%>    Ilb = find(d<0);
%>    if isempty(Ilb)
%>      alphal = 0;
%>    else
%>      alphal = max((QPALM_lbx(Ilb)-QPALM_x(Ilb))./d(Ilb));
%>    end
%>    alpha = max(alphal,alphau);
%>    alpha_type = 'bound';	% to memorize that alpha has been set to the last hit on the bounds
%>  end

  if QPALM_options.gpin == QPALM_values.inc_bkpt;	% take for initial stepsize the first break point that increases phi

    Ixl = find((d<0)&(QPALM_lbx>-inf)&(QPALM_x>QPALM_lbx+QPALM_options.dxmin));
    alpha_xl = (QPALM_lbx(Ixl)-QPALM_x(Ixl))./d(Ixl);
    Ixu = find((d>0)&(QPALM_ubx<inf)&(QPALM_x<QPALM_ubx-QPALM_options.dxmin));
    alpha_xu = (QPALM_ubx(Ixu)-QPALM_x(Ixu))./d(Ixu);
    Ixlu = [Ixl;Ixu];
    alpha_xlu = [alpha_xl;alpha_xu];
    [alpha_x,Ix] = sort(alpha_xlu);		% the sorted break points: alpha_x = alpha_xlu(Ix)
    phix_ = phi0;
    xx = QPALM_x;
    alpha_ = 0;
    dd = d;
    nhit = length(Ixlu);

    if nhit

      alpha_small = alpha;			% this the stepsize computed by the inital curvature

      for i = 1:nhit
        alpha = alpha_x(i);
        xx = xx+(alpha-alpha_)*dd;
        dd(Ixlu(i)) = 0;
        phixx = qpalm_phi (false,xx);
        if phixx > phix_; break; end
	alpha_ = alpha;
      end

      alpha_type = sprintf('break point %i',i);	% to memorize that alpha has been computed by using curv_unbd

      if i > 1
        alpha_small = min(alpha_small,alpha_x(i-1));	% include the previous break point as a rescue stepsize
      end
      alpha_small = max(alpha_small,damin);
      si = max(1.e-5,min(si,alpha_small/alpha));	% allow alpha_small to be considered by the interpolation safeguard

    end

  end


  % Initialization of the linesearch

  x0 = QPALM_x;
  itls = 0;				% linesearch iteration counter

  % Piecewise lineasearch
  
  if QPALM_options.plsn == QPALM_values.armijo

    % Armijo's rule

%>    if verb >= 3
%>      fprintf(fout,'\n  > Armijo search');
%>      fprintf(fout,'\n  > . damin = %10.4e, slope0 = %12.5e, initial stepsize -> %s',damin,slope0,alpha_type);
%>      fprintf(fout,'\n  > . iter   stepsize      decrease      FD slope');
%>    end
%>  
%>    [x,alpha,phi,info] = qpalm_armijo (x,lmi,lme,phi0,alpha,d,gphi,slope0);
%>
%>    if QPALM_info.flag == QPALM_values.stepsize_failure;
%>      alphamax = alpha;
%>      alphamin = 0;
%>    end

    w1 = 0.01;	% Armijo's constant

    if verb >= 3
      fprintf(fout,'\n  > Armijo search');
      fprintf(fout,'\n  > . damin = %10.4e, slope0 = %12.5e, initial stepsize -> %s, si = %7.1e',damin,slope0,alpha_type,si);
      fprintf(fout,'\n  > . iter   stepsize      decrease      FD slope');
    end

    while 1

      itls = itls+1;
      QPALM_x = max(QPALM_lbx,min(QPALM_ubx,x0-alpha*gphi));	% projection of x0-alpha*gphi on [lbx,ubx]
      phi = qpalm_phi (false);
      pdiff = phi-phi0;

      if pdiff <= w1*(gphi'*(QPALM_x-x0));	% Armijo's rule is satisfied, hence the stepsize is fine
        if verb >= 3
          fprintf(fout,'\n  > . %4i  %10.4e   %12.5e  %12.5e',itls,alpha,pdiff,pdiff/alpha);
          fprintf(fout,'\n  > ..................... initial AL = %18.11e',phi0);
          fprintf(fout,'\n  > ....................... final AL = %18.11e',phi);
        end
        break
      end

      if verb >= 3; fprintf(fout,'\n  > . %4i  %10.4e]  %12.5e  %12.5e',itls,alpha,pdiff,pdiff/alpha); end

      alpha = interpol2(phi0,slope0,alpha,phi,si*alpha,(1-si)*alpha);	% find a new trial stepsize by quadratic interpolation

      if alpha < damin
        QPALM_info.flag = QPALM_values.stepsize_failure;
        alphamax = alpha;
        alphamin = 0;
        break
      end

    end

    QPALM_info.gp_dec = pdiff;

  else

    % Goldstein's rule

    g1 = 0.01;	% first Goldstein's constant
    g2 = 0.99;	% second Goldstein's constant, must be in ]g1,1[

    if verb >= 3
      fprintf(fout,'\n  > Goldstein search');
      fprintf(fout,'\n  > . damin = %10.4e, slope0 = %12.5e, initial stepsize -> %s',damin,slope0,alpha_type);
      fprintf(fout,'\n  > . iter   stepsize      decrease      FD slope');
    end

    alphamin = 0;	% left bound on the security interval
    alphamax = inf;	% right bound on the security interval

    while 1

      itls = itls+1;
      QPALM_x = max(QPALM_lbx,min(QPALM_ubx,x0-alpha*gphi));	% projection of x0-alpha*gphi on [lbx,ubx]
      phi = qpalm_phi (false);
      pdiff = phi-phi0;
      slope = gphi'*(QPALM_x-x0);

      if pdiff > g1*slope;	% First Goldstein's rule is NOT satisfied
        if verb >= 3; fprintf(fout,'\n  > . %4i  %10.4e]  %12.5e  %12.5e',itls,alpha,pdiff,pdiff/alpha); end
        alphamax = alpha;
        if alphamax-alphamin < damin
% pdiff
% g1
% g2
% slope
% g1*slope
% g2*slope
% np=500
% for i=1:np+1
%   a(i) = ((i-1)/np)*alphamax;
%   xx = max(QPALM_lbx,min(QPALM_ubx,x0-a(i)*gphi));
%   gg = gphi'*(xx-x0);
%   gg1(i) = g1*gg;
%   gg2(i) = g2*gg;
%   ph(i) = qpalm_phi (false,xx)-phi0;
% end
% a;
% ph;
% clf
% hold on
% %axis([0 alphamax -2.e-28 2.e-28]);
% plot(a,ph);
% plot(a,gg1,'r');
% plot(a,gg2,'r');
% figure(gcf)
          QPALM_info.flag = QPALM_values.stepsize_failure;
          break
        end
        % find a new trial stepsize by quadratic interpolation
        alpha = interpol2 (phi0,slope0,alpha,phi,(1-si)*alphamin+si*alphamax,si*alphamin+(1-si)*alphamax);
      else			% First Goldstein's rule is satisfied
        if pdiff >= g2*slope;		% Second Goldstein's rule is satisfied
          if verb >= 3
            fprintf(fout,'\n  > . %4i  %10.4e   %12.5e  %12.5e',itls,alpha,pdiff,pdiff/alpha);
            fprintf(fout,'\n  > ..................... initial AL = %18.11e',phi0);
            fprintf(fout,'\n  > ....................... final AL = %18.11e',phi);
          end
          break
        else
          if verb >= 3; fprintf(fout,'\n  > . %4i [%10.4e   %12.5e  %12.5e',itls,alpha,pdiff,pdiff/alpha); end
          alphamin = alpha;
          if alphamax == inf
            alpha = extrapol2 (phi0,slope0,alpha,phi,se*alpha);
          else
            if alphamax-alphamin < damin
              QPALM_info.flag = QPALM_values.stepsize_failure;
              break
            end
            alpha = interpol2 (phi0,slope0,alpha,phi,(1-si)*alphamin+si*alphamax,si*alphamin+(1-si)*alphamax);
          end
        end
      end

    end

    QPALM_info.gp_dec = pdiff;

  end

  % Failure of the piecewise linesearch on dxmin

  if QPALM_info.flag == QPALM_values.stepsize_failure

    if verb >= 3; fprintf(fout,'\n  > too small trial stepsize interval (%11.5e)',alphamax-alphamin); end

    % Check whether, by chance, the optimality conditions are satisfied at the new point

    [phi,gphi] = qpalm_phi (true);

    if QPALM_nb
      Ilb = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables at the lower bound
      Iub = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables at the upper bound
      QPALM_IB      = zeros(QPALM_n,1);
      QPALM_IB(Ilb) = -1;
      QPALM_IB(Iub) = +1;
      gphi(Ilb) = min(0,gphi(Ilb));
      gphi(Iub) = max(0,gphi(Iub));	% projected gradient at the x
    end

    QPALM_info.gopt = norm(gphi,inf);
    if verb >= 3; fprintf(fout,'\n  > |Projected gradient| = %12.5e',QPALM_info.gopt); end
    if QPALM_info.gopt <= QPALM_options.gtol
      QPALM_info.flag = QPALM_values.success;	% a solution has been found
      QPALM_info.phi = phi;
      if verb >= 3; fprintf(fout,' (optimality reached)'); end
    else
      QPALM_x = x0;				% return the initial point and phi if failure
      QPALM_info.phi = phi0;
      if verb >= 3; fprintf(fout,' (stop on dxmin during linesearch)'); end
    end

    return

  end

%-------------------------------------------------------------------------------------------------------------------------------
% Detect too restrictive required precision.
%
% It is at this point that QPALM checks whether the precision required on the solution by 'options.gtol' can be reached. At this
% point, x on entry does not satisfy the required precision and a new point x has been computed by the GP phase. It is then
% decided that if
%  - the preceeding minimisation phase found a minimum in the previous activated face (this is detected by QPALM_info.flag ==
%    QPALM_values.relative_minimum) and
%  - the present GP phase did not change the activated phase (this is detected by a comparison with the activated phase on
%    entry, stored in Ilb0 and Iub0),
% then returning to the minimization phase is useless. In that case, a new test of optimality is done (at the new x) and if this
% one is not satified (this is likely since the GP phase is not aimed at making good progress on optimality), failure of the
% minimization process is declared.
%-------------------------------------------------------------------------------------------------------------------------------

  % Recompute activity

  if QPALM_nb
    Ilb = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables at the lower bound
    Iub = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables at the upper bound
    QPALM_IB      = zeros(QPALM_n,1);
    QPALM_IB(Ilb) = -1;
    QPALM_IB(Iub) = +1;
    if (verb >= 3); qpalm_print_activb (fout,verb); end
  end

  if QPALM_mi
    qpalm_yp ();
    if (verb >= 3); qpalm_print_activi (fout,verb); end
  end

  % GP phase is successful if the previous minimization phase did not stop at a minimizer of the previous activated face (this
  % is only the case the first time qpalm_gp is called, next it is always called after a minimisation phase that has reached a
  % minimum on the activated face).

  if QPALM_info.flag ~= QPALM_values.relative_minimum;
    QPALM_info.flag = QPALM_values.gp_successful;
    if ((QPALM_nb==0) | all(QPALM_IB == IB0)) & ((QPALM_mi==0) | all(QPALM_II == II0))
      QPALM_info.no_change_in_activity = true;
    end
    return
  end

  % Check whether the activaty in (x,y) is the same as the one on entry. If not, the GP phase is successful.

  if QPALM_nb & any(QPALM_IB ~= IB0)
    QPALM_info.flag = QPALM_values.gp_successful;
    return
  end
  if QPALM_mi & any(QPALM_II ~= II0)
    QPALM_info.flag = QPALM_values.gp_successful;
    return
  end
  QPALM_info.no_change_in_activity = true;	% no change in activity

  % Here, the required precision is too stringent, the next minimization phase is unlikely to make significant progress. Hence
  % check whether, by chance, optimality has been reached at the new x by computing the projected gradient 'pgphi'
  %
  %               { min(0,gphi(i))  if lbx(i) = x(i)
  %    pgphi(i) = { gphi(i)         if lbx(i) < x(i) < ubx(i)
  %               { max(0,gphi(i))  if          x(i) = ubx(i).

  [phi,gphi] = qpalm_phi (true);
  QPALM_info.phi = phi;

  if QPALM_nb
    Ilb = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables at the lower bound
    Iub = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables at the upper bound
    gphi(Ilb) = min(0,gphi(Ilb));
    gphi(Iub) = max(0,gphi(Iub));		% projected gradient at the new x
  end

  % Test convergence

  QPALM_info.gopt = norm(gphi,inf);
  if verb >= 3; fprintf(fout,'\n  > |Projected gradient| = %12.5e',QPALM_info.gopt); end
  if QPALM_info.gopt <= QPALM_options.gtol
    QPALM_info.flag = QPALM_values.success;		% a solution has been found
    if verb >= 3; fprintf(fout,' (optimality reached)'); end
  else
    QPALM_info.flag = QPALM_values.gtol_unreachable;	% a solution has not been found
    if verb >= 3; fprintf(fout,' (too stringent optimality tolerance)'); end
  end

%fprintf('\n');
%fprintf('|x-x0|         = %12.5e\n',norm(QPALM_x-x0,inf));
%grad0 = QPALM_H*x0;
%grad = QPALM_H*QPALM_x;
%fprintf('|QPALM_H*x-QPALM_H*x0|     = %12.5e\n',norm(grad-grad0,inf));
%grad0 = QPALM_g+QPALM_H*x0;
%grad = QPALM_g+QPALM_H*QPALM_x;
%aix0 = QPALM_AI*x0;
%aix = QPALM_AI*QPALM_x;
%fprintf('|QPALM_AI*x-QPALM_AI*x0|   = %12.5e\n',norm(aix-aix0,inf));
%r = QPALM_info.rlag;
%yp0 = max(QPALM_lbi,min(QPALM_ubi,aix0+QPALM_lmi/r));
%yp = max(QPALM_lbi,min(QPALM_ubi,aix+QPALM_lmi/r));
%fprintf('|yp-yp0|       = %12.5e\n',norm(yp-yp0,inf));
%ri0 = aix0-yp0;
%ri = aix-yp;
%fprintf('|ri-ri0|       = %12.5e\n',norm(ri-ri0,inf));
%grad0 = grad0 + QPALM_AI'*(QPALM_lmi+r*ri0);
%grad = grad + QPALM_AI'*(QPALM_lmi+r*ri);
%fprintf('|grad-grad0|   = %12.5e\n',norm(grad-grad0,inf));
%re0 = QPALM_AE*x0-QPALM_be;
%re = QPALM_AE*QPALM_x-QPALM_be;
%fprintf('|re-re0|       = %12.5e\n',norm(re-re0,inf));
%grad0 = grad0 + QPALM_AE'*(QPALM_lme+r*re0);
%grad = grad + QPALM_AE'*(QPALM_lme+r*re);
%fprintf('|grad-grad0|   = %12.5e\n',norm(grad-grad0,inf));
%Ilb = find(x0<=QPALM_lbx+QPALM_options.dxmin);
%Iub = find(QPALM_ubx-QPALM_options.dxmin<=x0);
%grad0(Ilb) = min(0,grad0(Ilb));
%grad0(Iub) = max(0,grad0(Iub));
%Ilb = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);
%Iub = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);
%grad(Ilb) = min(0,grad(Ilb));
%grad(Iub) = max(0,grad(Iub));
%fprintf('|pg-pg|        = %12.5e\n',norm(grad-grad0,inf));
%keyboard
%stop

%-------------------------------------------------------------------------------------------------------------------------------
% Exit
%-------------------------------------------------------------------------------------------------------------------------------

  return

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t] = interpol2 (f0,fp0,t1,f1,tmin,tmax)

%%
% [t] = interpol2 (f0,fp0,t1,f1,tmin,tmax)
%
% Interpol2 computes the minimizer t of a quadratic function f
% interpolating the values f(0) = f0, f'(0) = fp0, and f(t1) = f1. This
% minimizer is then projected on [tmin,tmax]. It is supposed that
% 1) t1 > 0 (otherwise t = tmin)
% 2) tmin < tmax (otherwise t = tmin).

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

% Special cases

  if (t1 == 0) | (tmin >= tmax)
    t = tmin;
    return
  end

% Find the minimizer t

  a = (f1-f0-fp0*t1)/(t1*t1);	% 2*curvature of the quadratics
  if a <= 0
    if fp0*tmin + a*tmin*tmin <= fp0*tmax + a*tmax*tmax
      t = tmin;
    else
      t = tmax;
    end
    return
  else
    t = -fp0/(2*a);
  end

% Projection

  if t < tmin
    t = tmin;
  elseif t > tmax
    t = tmax;
  end

% exit

  return

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t] = extrapol2 (f0,fp0,t1,f1,tmin)

%%
% [t] = interpol2 (f0,fp0,t1,f1,tmin)
%
% Extrapol2 computes the minimizer t of a quadratic function f
% interpolating the values f(0) = f0, f'(0) = fp0, and f(t1) = f1. This
% minimizer is then projected on [tmin,Inf[. It is supposed that t1 ~= 0
% (otherwise t = tmin).

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

% Special cases

  if (t1 <= 0)
    t = tmin;
    return
  end

% Find the minimizer t

  a = (f1-f0-fp0*t1)/(t1*t1);	% 2*curvature of the quadratics
  if a <= 0
    t = tmin;
    return
  else
    t = -fp0/(2*a);
  end

t
% Projection

  if t < tmin
    t = tmin;
  end

% exit

  return

  end

