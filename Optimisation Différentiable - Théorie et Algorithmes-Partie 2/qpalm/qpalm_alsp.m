function [] = qpalm_alsp ()

%%
% [] = qpalm_alsp ()
%
% This function aims at solving an augmented Lagrangian subproblem in
% the solver QPALM.

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

  global QPALM_nb QPALM_mi QPALM_me
  global QPALM_g QPALM_AI QPALM_AE
  global QPALM_x QPALM_y
  global QPALM_info QPALM_options QPALM_values
  global QPALM_IB QPALM_II

% Parameters

  % There are two parameters that monitor the Moré-Toraldo strategy:
  % - the parameter 'mt1' is used to see whether the decrease QPALM_info.gp_dec (negative) obtained in the last GP phase is large
  %   enough, i.e., suffciently negative (this is eta_2 in the 1991 Moré-Toraldo paper, which is set to 0.25 in their numerical
  %   experiments); if not, a CG phase is done next; this parameter must be sufficiently large to avoid zigzaging (0.25 is not
  %   large enough);
  % - the parameter 'mt2' is used to see whether the decrease in the successive GP phases is not too large with respect to the
  %   one obtained in the previous CG phase, in which case the problem is suspected to be unbounded and a CG phase is done next
  %   to detect that possibility (GP phases may indeed not be able to detect unboundedness).

  mt1 = 1.0;	% must be >0 (if larger, there is less consecutive GP phases)
  mt2 = 0.1;	% must be >0 (if larger, there is more consecutive GP phases)

% Initializations

  fout = QPALM_options.fout;		% output channel
  verb = QPALM_options.verb;		% verbosity level

  mie = QPALM_mi+QPALM_me;

  bounds = QPALM_nb+QPALM_mi;	% nb of bounds on x and y
  np = 0;			% phase counter
  QPALM_y = [];			% it is the GP phase that first determines y

  prec = [];			% no preconditioner has been made yet

  % reset the values of the min and max Rayleigh quotients of H, for the next QP solve (same thing for the preconditioned H, in
  % case verb >= 5 and preconditioning is required)

  QPALM_info.rqmin =  Inf;
  QPALM_info.rqmax = -Inf;
  if (verb >= 5) & (QPALM_options.cgprec ~= QPALM_values.cg_prec_none)
    QPALM_info.rqpmin =  Inf;
    QPALM_info.rqpmax = -Inf;
  end

% Decide which phase type to start first

  if bounds
    phase_type = QPALM_values.gp_phase;
  else
    phase_type = QPALM_values.cg_phase;
  end

% Initialization

  QPALM_info.flag   = QPALM_values.success;	% be sure that QPALM_info.flag ~= QPALM_values.relative_minimum
  QPALM_info.ccgp   = 0;			% counter of consecutve CG phases
  QPALM_info.accu   = 1;			% accuracy of the CG linear solver
  gp_dec_min        = 0;			% max negative phi decrease obtained in the current series of successive GP phases
  gp_dec_tot        = 0;			% total negative phi decrease obtained in the current series of successive GP phases
  QPALM_info.cg_dec = 0;			% negative phi decrease obtained in the las CG phase

%-------------------------------------------------------------------------------------------------------------------------------
% Loop on the phases
%-------------------------------------------------------------------------------------------------------------------------------

  while 1

    %-------------------------------------------------------------------
    % test cput
    %-------------------------------------------------------------------

    if (QPALM_options.cput < inf) & (etime(clock,QPALM_info.time0) > QPALM_options.cput);
      QPALM_info.flag = QPALM_values.stop_on_max_cput;
      break
    end

    %-------------------------------------------------------------------
    % increase phase counter
    %-------------------------------------------------------------------

    np = np+1;
    if (np > QPALM_options.phase)
      QPALM_info.flag = QPALM_values.stop_on_max_phase;
      break
    end

    switch phase_type
    case QPALM_values.gp_phase

    %-------------------------------------------------------------------
    % Gradient projection phase.
    %-------------------------------------------------------------------

      % Increase the counter of GP phases

      if QPALM_info.gpph >= QPALM_options.gpph
        QPALM_info.flag = QPALM_values.stop_on_max_gpph;
        break
      end
      QPALM_info.gpph = QPALM_info.gpph+1;

      if verb >= 3
        if np > 1; fprintf(fout,'\n'); end
        fprintf(fout,'\n  (%i: gradient projection phase)\n',np);
      end

      % GP phase

      qpalm_gp ();

      % Stop solving the current AL subproblem if the GP phase is not successful

      if QPALM_info.flag ~= QPALM_values.gp_successful; break; end

      % Decide which phase to do next

      if any(QPALM_IB == 0)				% otherwise, there is no variable to change in the CG phase
        gp_dec_tot = gp_dec_tot + QPALM_info.gp_dec;
        if QPALM_options.gpsucc ~= QPALM_values.gp_mt	% 'hit-and-fix' strategy: a GP phase is always followed by a CG phase
          phase_type = QPALM_values.cg_phase;
        elseif QPALM_info.no_change_in_activity		% no change in constraint activity -> CG phase next
          if verb >= 3; fprintf(fout,'\n  > No change in constraint activity -> CG phase'); end
          phase_type = QPALM_values.cg_phase;
        elseif QPALM_info.gp_dec > mt1*gp_dec_min	% not enough decrease in the last GP phase -> CG phase next
          if verb >= 3; fprintf(fout,'\n  > Too small decrease (%12.5e > %12.5e*%8.2e) -> CG phase',QPALM_info.gp_dec,gp_dec_min,mt1); end
          phase_type = QPALM_values.cg_phase;	
        elseif gp_dec_tot < mt2*QPALM_info.cg_dec	% too large total decrease in the successive GP phases -> CG phase next
          if verb >= 3
            fprintf(fout,'\n  > Too large total decrease (%12.5e < %12.5e*%8.2e) -> CG phase',gp_dec_tot,QPALM_info.cg_dec,mt2);
          end
          phase_type = QPALM_values.cg_phase;	
        end
      end

      if phase_type == QPALM_values.cg_phase
        gp_dec_min = 0;
      else
        gp_dec_min = min(gp_dec_min,QPALM_info.gp_dec);	% max phi decrease obtained in the current series of successive GP phases
        QPALM_info.ccgp = 0;				% reset the counter of consecutve CG phases (still useful?)
      end

    case QPALM_values.cg_phase

    %-------------------------------------------------------------------
    % Conjugate gradient phase
    %-------------------------------------------------------------------

      if QPALM_info.cgph >= QPALM_options.cgph
        QPALM_info.flag = QPALM_values.stop_on_max_cgph;
        break
      end
      QPALM_info.cgph = QPALM_info.cgph+1;	% increase the counter of CG phases
      QPALM_info.ccgp = QPALM_info.ccgp+1;	% increase the counter of consecutve CG phases

      if verb >= 3
        if bounds
          if (np > 1), fprintf(fout,'\n'); end
          fprintf(fout,'\n  (%i: conjugate gradient phase)\n',np);
        else
          fprintf(fout,'\n  > Conjugate gradient iterations',np);
        end
      end

      % CG phase

      if QPALM_options.cgas == QPALM_values.cg_hf;	% 'hit-and-fix' strategy
        [prec] = qpalm_tpcg (prec);
      else						% 'more-toraldo' strategy
        [prec] = qpalm_tpcg_mt (prec);
stop
      end

      % Decide which phase to do next

      switch QPALM_info.flag
      case QPALM_values.relative_minimum
        if bounds
          phase_type = QPALM_values.gp_phase;	% do a gradient projection phase next
          gp_dec_tot = 0;
        else
          QPALM_info.flag = QPALM_values.success;
          break					% no need to pursue if there is no bound on x
        end
      case QPALM_values.a_bound_is_hit
        if all(QPALM_IB ~= 0)				% do a GP phase next if there is no more free variable
          phase_type = QPALM_values.gp_phase;
          gp_dec_tot = 0;
        end
      otherwise
        break;
      end

    end

  end

% Exit

  QPALM_info.wi = length(find(QPALM_II~=0));

% Exit

  return
