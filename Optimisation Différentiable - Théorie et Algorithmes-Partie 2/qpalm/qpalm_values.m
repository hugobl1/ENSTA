function [] = qpalm_values ()

%%
% [] = qpalm_values ()
%
% Set the constant values used in Qpalm.

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

  global QPALM_values

% Diagnosis values

  QPALM_values.success               = int8(  0);	% solution found
  QPALM_values.success_shifted       = int8(  1);	% solution found (nonzero smallest feasible shift)
  QPALM_values.stop_on_unboundedness = int8(  2);	% unbounded QP
  QPALM_values.ctol_unreachable      = int8(  3);	% the feasibility tolerance options.feas cannot be reached
  QPALM_values.gtol_unreachable      = int8(  4);	% the optimality tolerance options.gtol cannot be reached

  QPALM_values.fail_on_argument      = int8(  5);	% an argument is wrong
  QPALM_values.fail_on_mex           = int8(  6);	% cannot generate a MEX-gile
  QPALM_values.stop_on_max_alit      = int8(  7);	% max augmented Lagrangian iterations
  QPALM_values.stop_on_max_phase     = int8(  8);	% max phases
  QPALM_values.stop_on_max_gpph      = int8(  9);	% max gradient projection phases
  QPALM_values.stop_on_max_cgph      = int8( 10);	% max conjugate gradient phases
  QPALM_values.stop_on_max_cgit      = int8( 11);	% max conjugate gradient iterations
  QPALM_values.stop_on_max_cput      = int8( 12);	% max CPU time reached
  QPALM_values.stop_on_nonconvexity  = int8( 13);	% H is not positive semidefinite
  QPALM_values.fail_on_accuracy      = int8( 14);	% too many rounding error in CG solve - measured by options.accu

% Running values

  QPALM_values.gp_phase              = int8(100);	% characterizes a gradient projection phase
  QPALM_values.cg_phase              = int8(101);	% characterizes a conjugate gradient phase
  QPALM_values.cg_prec_none          = int8(102);	% CG iterations are not preconditioned
  QPALM_values.cg_prec_diag          = int8(103);	% CG iterations are preconditioned by a diagonal preconditioner
  QPALM_values.cg_prec_chol          = int8(104);	% CG iterations are preconditioned by a Cholesky preconditioner
  QPALM_values.a_bound_is_hit        = int8(105);	% a linear solver has hit a bound
  QPALM_values.relative_minimum      = int8(106);	% a linear solver has found a minimum on the active face
  QPALM_values.gp_successful         = int8(107);	% the gradient projection phase has found a new iterate
  QPALM_values.stepsize_failure      = int8(108);	% a stepsize cannot be determined in linesearch
  QPALM_values.r_fixed               = int8(109);	% maintain r fixed
  QPALM_values.r_constraint          = int8(110);	% update r to have the a convergence rate 'options.dcr' on the constraints
  QPALM_values.r_constraint_change   = int8(111);	% update r to have the a convergence rate 'options.dcr' on the constraint changes
  QPALM_values.curv_init             = int8(112);	% the initial stepsize trial in the GP is computed using the initial curvature
  QPALM_values.inc_bkpt              = int8(113);	% the initial stepsize trial in the GP is computed using the minimum breakpoint
  QPALM_values.armijo                = int8(114);	% Armijo's piecewise lineasearch
  QPALM_values.goldstein             = int8(115);	% Goldstein's piecewise lineasearch
  QPALM_values.gp_one                = int8(116);	% Do several GP phases until Moré-Toraldo test is satisfied
  QPALM_values.gp_mt                 = int8(117);	% Do several GP phases until Moré-Toraldo test is satisfied
  QPALM_values.cg_hf                 = int8(128);	% CG active set strategy: hit-and-fix
  QPALM_values.cg_mt                 = int8(119);	% CG active set strategy: Moré-Toraldo
  QPALM_values.cg_mt_out             = int8(120);	% CG active set strategy: Moré-Toraldo - stop CG since also out of B bounds
  QPALM_values.cg_mt_badface         = int8(121);	% CG active set strategy: Moré-Toraldo - stop CG since also wrong activated face

% Other values

  QPALM_values.eline = '================================================================================';
  QPALM_values.dline = '--------------------------------------------------------------------------------';

% Return

return
