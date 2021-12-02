function [] = qpalm_finish (lm)

%%
% [] = qpalm_finish (lm)
%
% QPALM termination procedure. It prints the output status

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
  global QPALM_g QPALM_H QPALM_lbx QPALM_ubx QPALM_AI QPALM_lbi QPALM_ubi QPALM_AE
  global QPALM_info QPALM_options QPALM_values
  global QPALM_x

% Initializations

  fout = QPALM_options.fout;
  verb = QPALM_options.verb;

  mie = QPALM_mi+QPALM_me;

% Cost

  QPALM_info.f = QPALM_x'*(QPALM_g+0.5*(QPALM_H*QPALM_x));

% Print the results

  fprintf(fout,'\n%s\nQPALM optimization solver (exit point)',QPALM_values.dline);
  fprintf(fout,'\n  Status %i ',QPALM_info.flag);

  switch QPALM_info.flag
  case QPALM_values.success

    fprintf(fout,'(a solution to the QP has been found)');
    fprintf(fout,'\n  Cost                                         %12.5e',QPALM_info.f);
    if QPALM_nb+QPALM_mi
      fprintf(fout,'\n  Constraint activity');
      if QPALM_nb; fprintf(fout,'\n  . nb of active bounds       %19i',QPALM_info.wb); end
      if QPALM_mi; fprintf(fout,'\n  . nb of active inequalities %19i',QPALM_info.wi); end
    end
    fprintf(fout,'\n  Optimality conditions (inf norm)');
    fprintf(fout,'\n  . projected gradient of the Lagrangian        ');
    if QPALM_info.pgln < inf
      fprintf(fout,'%7.2e',QPALM_info.pgln);
    else
      fprintf(fout,'unavailable');
    end
    if QPALM_nb+mie
      fprintf(fout,'\n  . feasibility                                 %7.2e',QPALM_info.feas);
    end

  case QPALM_values.success_shifted

    fprintf(fout,'(a solution to the closest feasible QP has been found)');
    fprintf(fout,'\n  Cost                                         %12.5e',QPALM_info.f);
    if QPALM_nb+QPALM_mi
      fprintf(fout,'\n  Constraint activity');
      if QPALM_nb; fprintf(fout,'\n  . nb of active bounds       %19i',QPALM_info.wb); end
      if QPALM_mi; fprintf(fout,'\n  . nb of active inequalities %19i',QPALM_info.wi); end
    end
    fprintf(fout,'\n  Optimality conditions (inf norm)');
    fprintf(fout,'\n  . projected gradient of the Lagrangian        ');
    if QPALM_info.pgln < inf
      fprintf(fout,'%7.2e',QPALM_info.pgln);
    else
      fprintf(fout,'unavailable');
    end
    fprintf(fout,'\n  . feasibility                                 %7.2e',QPALM_info.feas);
    fprintf(fout,'\n  . closest feasibility                         %7.2e',QPALM_info.feass);
    fprintf(fout,'\n  . feasibility shift                           %7.2e',norm(QPALM_info.s,inf));

  case QPALM_values.fail_on_argument

    fprintf(fout,'(erroneous argument)');

  case QPALM_values.stop_on_max_alit

    fprintf(fout,'(max nb of augmented Lagrangian iterations reached)');

  case QPALM_values.stop_on_max_phase

    fprintf(fout,'(max nb of phases reached)');

  case QPALM_values.stop_on_max_gpph

    fprintf(fout,'(max nb of gradient projection phases reached)');

  case QPALM_values.stop_on_max_cgph

    fprintf(fout,'(max nb of conjugate gradient phases reached)');

  case QPALM_values.stop_on_max_cgit

    fprintf(fout,'(max nb of conjugate gradient iterations reached)');

  case QPALM_values.stop_on_max_cput

    fprintf(fout,'(max CPU time reached)');

  case QPALM_values.stop_on_unboundedness

    if mie
      fprintf(fout,'(the closest feasible problem is unbounded)');
    else
      fprintf(fout,'(unbounded problem)');
    end
    fprintf(fout,'\n  Found an l2-normalized feasible direction d=info.ubdd satisfying');
    fprintf(fout,'\n  . g''*d                  = %12.5e (must be negative)',QPALM_g'*QPALM_info.ubdd);
    fprintf(fout,'\n  . d''*H*d                = %12.5e (must vanish)', ...
      QPALM_info.ubdd'*(QPALM_H*QPALM_info.ubdd)/(QPALM_info.ubdd'*QPALM_info.ubdd));
    if QPALM_nb
      fprintf(fout,'\n  . |P_{[l,u]^inf}(d)|    = %12.5e (must vanish)', ...
        norm([min(0,QPALM_info.ubdd(find(QPALM_lbx>-inf)));max(0,QPALM_info.ubdd(find(QPALM_ubx<inf)))]));
    end
    if QPALM_mi
      aid = QPALM_AI*QPALM_info.ubdd;
      fprintf(fout,'\n  . |P_{[l,u]^inf}(AI*d)| = %12.5e (must vanish)', ...
        norm([min(0,aid(find(QPALM_lbi>-inf)));max(0,aid(find(QPALM_ubi<inf)))]));
    end
    if QPALM_me
      fprintf(fout,'\n  . |AE*d|                = %12.5e (must vanish)',norm(QPALM_AE*QPALM_info.ubdd));
    end

  case QPALM_values.stop_on_nonconvexity

    fprintf(fout,'(H is not positive semidefinite)');

  case QPALM_values.stepsize_failure

    fprintf(fout,'(stop on dxmin during linesearch)');
    fprintf(fout,'\n  Cost                                   %12.5e',QPALM_info.f);
    fprintf(fout,'\n  Optimality conditions (inf norm)');
    fprintf(fout,'\n  . gradient of the cost function         %11.5e',QPALM_info.gopt);

  case QPALM_values.gtol_unreachable

    fprintf(fout,'(too stringent optimality tolerance)');
    fprintf(fout,'\n  Cost                                         %12.5e',QPALM_info.f);
    fprintf(fout,'\n  Optimality conditions (inf norm)');
    fprintf(fout,'\n  . projected gradient of the Lagrangian        ');
    if QPALM_info.pgln < inf
      fprintf(fout,'%7.2e',QPALM_info.pgln);
    else
      fprintf(fout,'unavailable');
    end
    if QPALM_nb+mie
      fprintf(fout,'\n  . feasibility                                 %7.2e',QPALM_info.feas);
    end

  case QPALM_values.ctol_unreachable

    fprintf(fout,'(too stringent feasibility tolerance');
    if QPALM_options.rctl == QPALM_values.r_constraint; fprintf(fout,' or infeasible QP');end
    fprintf(fout,')');
    fprintf(fout,'\n  Cost                                         %12.5e',QPALM_info.f);
    fprintf(fout,'\n  Optimality conditions (inf norm)');
    fprintf(fout,'\n  . projected gradient of the Lagrangian        ');
    if QPALM_info.pgln < inf
      fprintf(fout,'%7.2e',QPALM_info.pgln);
    else
      fprintf(fout,'unavailable');
    end
    if QPALM_nb+mie
      fprintf(fout,'\n  . feasibility                                 %7.2e',QPALM_info.feas);
    end

  case QPALM_values.fail_on_accuracy

    fprintf(fout,'(too much rounding error in CG solve - measured by options.accu)');

  end

  fprintf(fout,'\n  Counters');
  if QPALM_nb+QPALM_mi; fprintf(fout,'\n  . nb of gradient-projection phases  %11i',QPALM_info.gpph); end;
  fprintf(fout,'\n  . nb of conjugate-gradient phases  %12i',QPALM_info.cgph);
  if verb >= 3
    fprintf(fout,'\n  . nb of conjugate-gradient iterations   %7i',QPALM_info.cgit);
    if mie
      fprintf(fout,'\n  Augmented Lagrangian (AL) algorithm');
      fprintf(fout,'\n  . nb of AL iterations                   %7i',QPALM_info.alit);
      fprintf(fout,'\n  . augmentation parameter r                   %12.5e',QPALM_info.rlag);
      fprintf(fout,'\n  . Lipschitz modulus L                        %12.5e',QPALM_info.lips);
    end
  end
  fprintf(fout,'\n  CPU time (before final diagnosis)             %5.3f sec',QPALM_info.cput);

  fprintf(fout,'\n%s\n',QPALM_values.eline);

% Warning on (those set off in qpalm_prelim)

  warning ('on','MATLAB:nearlySingularMatrix');		% 'Matrix is close to singular or badly scaled' from linsolve

% Exit

  return
