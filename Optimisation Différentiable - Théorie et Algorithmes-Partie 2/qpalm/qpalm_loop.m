function [] = qpalm_loop ()

%%
% [] = qpalm_loop ()
%
% This function implements the main optimization loop of QPALM. There is
% one loop per augmented Lagrangien iteration (hence a single loop if
% there is no linear equality or inequality constraint).

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
  global QPALM_g QPALM_lbx QPALM_ubx QPALM_AI QPALM_lbi QPALM_ubi QPALM_AE QPALM_be
  global QPALM_x QPALM_y QPALM_lmi QPALM_lme
  global QPALM_info QPALM_options QPALM_values
  global QPALM_IB QPALM_II

%-----------------------------------------------------------------------
% Initializations
%-----------------------------------------------------------------------

  fout = QPALM_options.fout;
  verb = QPALM_options.verb;

  QPALM_info.alit = 0;	% counter of AL iterations
  QPALM_info.cgit = 0;	% counter of CG iterations
  QPALM_info.gpph = 0;	% counter of GP phases
  QPALM_info.cgph = 0;	% counter of CG phases
  cgph = 0;		% value of the CG phase counter before calling qpalm_alsp
  gpph = 0;		% value of the GP phase counter before calling qpalm_alsp

  mie = QPALM_mi+QPALM_me;

  A = [QPALM_AI; QPALM_AE];

  stopping_test_on_c = true;	% initially the stopping test is on the constraint value
  cmark  = '*';
  scmark = ' ';

%-----------------------------------------------------------------------
% Prepare the AL loop
%-----------------------------------------------------------------------

  if mie	% no AL loop otherwise
    alit      = 0;
    QPALM_info.rlag = QPALM_options.rlag;	% initial r
    rlag_prev = 0;				% no previous r, no rlag reduction
    max_rlag_red = 2;				% max number of rlag reductions
    L         = 0;				% Lipschitz constant estimate
  end

  QPALM_IB = zeros(QPALM_n,1);
  if QPALM_nb
    Ilb0 = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables at the lower bound
    Iub0 = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables at the upper bound
    QPALM_IB(Ilb0) = -1;
    QPALM_IB(Iub0) = +1;
  end

  % TODO initialize correctly

  if QPALM_mi; QPALM_II = zeros(QPALM_mi,1); end	% II(i) = -1 (resp. = +1) if lower (resp. upper) bound active on y, II(i) = 0 otherwise

%-----------------------------------------------------------------------
% Dry printings in case of AL iterations
%-----------------------------------------------------------------------

  if (verb == 2) & mie
    fprintf(fout,'\n%s\nALit     r   ',QPALM_values.dline);
    if QPALM_nb+QPALM_mi; fprintf(fout,'  GPph  MINph'); end
    fprintf(fout,'  CGit');
    cgit = QPALM_info.cgit;	% used to count the number of CG iterations used to solve this AL subproblem (QPALM_info.cgit is a cumulated counter)
    if QPALM_nb; fprintf(fout,'  |WB|'); end
    if QPALM_mi; fprintf(fout,'  |WI|'); end
    fprintf(fout,'     |c_QP|   ');
    if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'    |c_CFQP|  '); end
    fprintf(fout,'    rho ');
    if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'     rho'''); end
    fprintf(fout,'      L');
  end

%-----------------------------------------------------------------------
% Augmented Lagrangian iteration loop
%-----------------------------------------------------------------------

  while 1

    %-------------------------------------------------------------------
    % Increment AL iteration counter (alit = 1 at the first iteration)
    %-------------------------------------------------------------------

    if mie	% no AL loop otherwise
      if alit >= QPALM_options.alit
        QPALM_info.flag = QPALM_values.stop_on_max_alit;
        break;
      end
      alit = alit+1;
    end

    %-------------------------------------------------------------------
    % Printings
    %-------------------------------------------------------------------

    if verb >= 3
      if mie
        fprintf(fout,'\n%s\nAL iter = %0i, r = %7.1e\n%s',QPALM_values.dline,alit,QPALM_info.rlag,QPALM_values.dline);
      else
        fprintf(fout,'\n%s',QPALM_values.dline);
      end
    end

    %-------------------------------------------------------------------
    % Minimise an AL subproblem in x subject to bounds
    %-------------------------------------------------------------------

%QPALM_lmi
%QPALM_lme
    qpalm_alsp ();

    if (verb == 2) & mie
      fprintf(fout,'\n%4i  %7.1e',alit,QPALM_info.rlag);
      if QPALM_nb+QPALM_mi
        fprintf(fout,'  %4i  %5i',QPALM_info.gpph-gpph,QPALM_info.cgph-cgph);
	gpph = QPALM_info.gpph;	% previous counter
	cgph = QPALM_info.cgph;	% previous counter
      end
      fprintf(fout,'  %4i',QPALM_info.cgit-cgit);
      cgit = QPALM_info.cgit;	% reset the counter
      if QPALM_nb; fprintf(fout,'  %4i',QPALM_info.wb); end
      if QPALM_mi; fprintf(fout,'  %4i',QPALM_info.wi); end
    end

    switch QPALM_info.flag
      case QPALM_values.stop_on_max_gpph
        break;
      case QPALM_values.stop_on_unboundedness
        break;
      case QPALM_values.stop_on_max_cgit;
        break;
      case QPALM_values.stop_on_max_cgph;
        break;
      case QPALM_values.stop_on_max_cput;
        break;
      case QPALM_values.stop_on_nonconvexity;
        break;
    end

    % Terminate if mie == 0

    if ~mie; break; end

    %-------------------------------------------------------------------
    % Evaluate the constraints (c) and constraint changes (cc)
    %-------------------------------------------------------------------

    % remember previous values if not the first iteration

    if alit > 1
      c_ = c;
      cn_  = cn;
      if alit > 2; cc_ = cc; end
    end

    % constraints

    c = [];
    if QPALM_mi
      AIx = QPALM_AI*QPALM_x;
      c = AIx-QPALM_y;
    end
    if QPALM_me
      c = [c; QPALM_AE*QPALM_x-QPALM_be];
    end
    cn = norm(c);

    % constraint changes

    if alit > 1; cc = c-c_; end

    %-------------------------------------------------------------------
    % New multipliers
    %-------------------------------------------------------------------

    if QPALM_mi; QPALM_lmi = QPALM_lmi + QPALM_info.rlag*c(1:QPALM_mi); end
    if QPALM_me; QPALM_lme = QPALM_lme + QPALM_info.rlag*c(QPALM_mi+1:QPALM_mi+QPALM_me); end

    %-------------------------------------------------------------------
    % Evaluate feasibility of the closest feasible problem
    %-------------------------------------------------------------------

    if QPALM_options.rctl ~= QPALM_values.r_constraint
      if alit > 1; scn_ = scn; end			% remember previous value if not the first iteration
      sc1 = [];
      sc2 = [];
      if mie
        sc1 = A'*c;						% projection of A'*c on the opposite of the tangent cone to [lb,ub]
        I = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables x at the lower bound
        sc1(I) = min(0,sc1(I));
        I = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables x at the upper bound
        sc1(I) = max(0,sc1(I));
      end
      if QPALM_mi
        sc2 = AIx;					% projection of AI*x on [li,bi] - y
        I = find(AIx<=QPALM_lbi+QPALM_options.dymin);	% indices of the variables AI*x at the lower bound
        sc2(I) = QPALM_lbi(I);
        I = find(AIx>=QPALM_ubi-QPALM_options.dymin);	% indices of the variables AI*x at the upper bound
        sc2(I) = QPALM_ubi(I);
        sc2 = sc2 - QPALM_y;
      end
      scn = norm([sc1;sc2]);
    end

    %-------------------------------------------------------------------
    % Check whether the constraint norm decreases
    %-------------------------------------------------------------------

    if alit > 1;
      rho_c = cn/cn_;
      if QPALM_options.rctl ~= QPALM_values.r_constraint; rho_sc = scn/scn_; end
      if stopping_test_on_c | (alit == 2)
        rho = rho_c;
      else
        rho = rho_sc;
      end
      L = max(L,rho*QPALM_info.rlag);
    end

    % printings

    if verb == 2
      if alit == 1;
        if QPALM_info.flag == QPALM_values.success
          fprintf(fout,'  %11.5e%s ',cn,cmark);
          if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,' %11.5e ',scn); end
        else
          fprintf(fout,' (%11.5e)',cn);
          if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'(%11.5e)',scn); end
        end
        fprintf(fout,'                             ');
      else
%     elseif alit == 2
        if QPALM_info.flag == QPALM_values.success
          fprintf(fout,'  %11.5e%s',cn,cmark);
          if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'  %11.5e%s',scn,scmark); end
          fprintf(fout,'  %7.1e',rho_c);
          if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'  %7.1e',rho_sc); end
          fprintf(fout,'  %7.1e',L);
        else
          fprintf(fout,' (%11.5e%s)',cn,cmark);
          if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'(%11.5e%s)',scn,scmark); end
          fprintf(fout,'(%7.1e)',rho_c);
          if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'(%7.1e)',rho_sc); end
          fprintf(fout,'(%7.1e)',L);
        end
%     else
%       if QPALM_info.flag == QPALM_values.success
%         fprintf(fout,'  %11.5e%s  %11.5e%s  %10.4e  %11.5e  %7.1e',cn,cmark,scn,scmark,rho_c,rho_sc,L);
%       else
%         fprintf(fout,' (%11.5e%s)(%11.5e%s)(%10.4e)(%11.5e)(%7.1e)',cn,cmark,scn,scmark,rho_c,rho_sc,L);
%       end
      end
    elseif verb >= 3;
      if QPALM_info.flag == QPALM_values.success
        es = '= ';
      else
        es = '~ ';
      end
      fprintf(fout,'\n\n  |AL''| %s%9.3e,  |c| %s%9.3e,  |shifted c| %s',es,QPALM_info.gopt,es,cn,es);
      if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'%9.3e',scn); end
      if alit > 1
        if success_ & (QPALM_info.flag == QPALM_values.success)
          es = '= ';
        else
          es = '~ ';
        end
        fprintf(fout,'\n  rho_c %s%9.3e',es,rho_c);
        if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,',  rho_sc %s%9.3e',es,rho_sc); end
        fprintf(fout,',  estimated L %s%11.5e',es,L);
      end
%     fprintf(fout,'\n');
    end
    success_ = (QPALM_info.flag == QPALM_values.success);

%> %<<<
%> %===
%> % Compute rate of convergence to sbar (to remove)
%> 
%>   % smallest feasible shift for problem randqp10k
%> 
%>   sbar = [...
%>     -0.000000000000000e+00
%>      9.244193382542121e-04
%>      7.710271643781375e-02
%>     -0.000000000000000e+00
%>     -0.000000000000000e+00
%>     -7.417528378605331e-02
%>     -1.245679683464694e-01
%>      4.787046367159019e-01
%>      8.975633732312741e-02
%>     -4.370441461438220e-01
%>      6.831319277423180e-01
%>      1.543578321738044e-01
%>      1.966239179524698e-02
%>     -5.183744246437756e-01
%>     -2.250428030751985e-01];
%> 
%>   % shifted constraints
%> 
%>   shc = c+sbar;
%>   if alit > 1; shcn0 = shcn; end
%>   shcn = norm(shc);
%>   if alit > 1; fprintf(fout,'  rho_shcn %s%9.3e\n',es,shcn/shcn0); end
%> 
%> %>>>

    %-------------------------------------------------------------------
    % Stopping tests (both can be used at the iteration at which one
    % switches from stopping_test_on_c to ~stopping_test_on_c, so that
    % it is not an if-then-else construct)
    %-------------------------------------------------------------------

    if stopping_test_on_c
      if cn <= QPALM_options.feas;
        QPALM_info.flag = QPALM_values.success;
        break;
      end
      if alit > 1;
        log_rho_c = log10(rho_c);
        if QPALM_options.rctl ~= QPALM_values.r_constraint; log_rho_sc = log10(rho_sc); end
        if alit == 2
          min_log_rho_c = log_rho_c;
          if QPALM_options.rctl ~= QPALM_values.r_constraint; min_log_rho_sc = log_rho_sc; end
        else
          if log_rho_c >= QPALM_options.dcrf*min_log_rho_c
            if QPALM_options.rctl == QPALM_values.r_constraint
              if verb >= 3
                fprintf(fout,'\n\n  Too stringent feasibility threshold or infeasible QP');
                fprintf(fout,'\n  In the latter case, set options.rctl to ''constraint_change'' to solve the CFQP');
              end
              QPALM_info.flag = QPALM_values.ctol_unreachable;
              break;
            else
              stopping_test_on_c = false;
              cmark  = ' ';
              scmark = '*';
              if verb >= 3; fprintf(fout,'\n\n  The stopping test is now on the constraints of the CFQP'); end
            end
          end
          min_log_rho_c = min(min_log_rho_c, log_rho_c);
          if QPALM_options.rctl ~= QPALM_values.r_constraint; min_log_rho_sc = min(min_log_rho_sc,log_rho_sc); end
        end
      end
    end

    if ~stopping_test_on_c
      if scn < QPALM_options.feass;
        QPALM_info.flag = QPALM_values.success_shifted;
        QPALM_info.s = -c;
        QPALM_info.feass = scn;
        break;
      end
      if (alit > 2) & (QPALM_options.rctl ~= QPALM_values.r_constraint)
        log_rho_sc = log10(rho_sc);
        if log_rho_sc >= QPALM_options.dcrf*min_log_rho_sc
          if verb >= 3; fprintf(fout,'\n  Too stringent feasibility threshold'); end
          QPALM_info.flag = QPALM_values.ctol_unreachable;
          break;
        end
        min_log_rho_sc = min(min_log_rho_sc,log_rho_sc);
      end
    end

    %-------------------------------------------------------------------
    % Update r
    %-------------------------------------------------------------------

%   if QPALM_options.rctl == QPALM_values.r_constraint_change;
%     fprintf(fout,'\n  ''QPALM_options.rctl == constraint_change'' is not implemented');
%     return
%   end

    if QPALM_options.rctl ~= QPALM_values.r_fixed
      if QPALM_info.accu < QPALM_options.accu	% the CG solver lacks of precision
        if rlag_prev > 0		% there is a previous r
          QPALM_info.rlag = rlag_prev;
          rlag_prev = -1;		% there is no previous smaller r any more, and r has been reduced once
          if verb >= 3; fprintf(fout,'\n  not enough accuracy => previous smaller r is recovered and no longer increased'); end
        elseif rlag_prev > -max_rlag_red;	% r has been reduces -rlag_prev times but not yet max_rlag_red times, hence can still be reduced
          QPALM_info.rlag = 1.e-1*QPALM_info.rlag;
          rlag_prev = rlag_prev-1;	% r has been reduced once more
          if verb >= 3; fprintf(fout,'\n  not enough accuracy => r is reduced'); end
        else
          QPALM_info.flag = QPALM_values.fail_on_accuracy;
          if verb >= 3; fprintf(fout,'\n  not enough accuracy, r too many times reduced => failure is declared'); end
          break;
        end
      elseif (alit > 1) ...		% no way to update r if alit == 1
           & (rho > QPALM_options.dcr) ...	% an increase of r is desirable since the convergence rate is not sufficient
           & (rlag_prev >= 0)		% r has never been decreased
        rlag_prev = QPALM_info.rlag;
        QPALM_info.rlag = (rho/QPALM_options.dcr)*QPALM_info.rlag;
      end
    end

  end

%-----------------------------------------------------------------------
% Put information
%-----------------------------------------------------------------------

  if mie
    QPALM_info.alit = alit;
    QPALM_info.lips = L;
  end

%-----------------------------------------------------------------------
% Exit
%-----------------------------------------------------------------------

  return
