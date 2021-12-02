function [x,lm,info] = qpalm (g,H,lbx,ubx,AI,lbi,ubi,AE,be,x0,lm0,options)

%%
% [x,lm,info] = qpalm (g,H,lbx,ubx,AI,lbi,ubi,AE,be,x0,lm0,options)
%
% QPALM solves a convex quadratic optimization problem written as
% follows (the notation is explained below)
%
%      minimize (in x)  g'*x + 0.5*x'*H*x
%      subject to       lbx <=    x <= ubx
%                       lbi <= AI*x <= ubi
%                       AE*x = be
%
% QPALM prints a few warning messages using the 'warning' function and a
% message identifier starting by 'qpalm:'. These messages can therefore
% be disabled or enabled by using "warning('off','qpalm')" or
% "warning('on','qpalm')".
%
% Input arguments
% ^^^^^^^^^^^^^^^
% g: n vector, specifying the gradient at x=0 of the objective;
% H: n x n symmetric positive semidefinite matrix, specifying the
%   Hessian of the objective function;
% lbx, ubx: n vectors used to specify the lower (lbx) and upper (ubx)
%   bounds on x (see the model above); a component of lbx (resp. ubx)
%   having the infinite value -inf (resp. inf) indicates that the
%   corresponding bound is not present; there must hold ubx-lbx >=
%   2*options.dxmin, hence variables cannot be maintained fixed (remove
%   variable i from the set of variables instead of imposing lbx(i) =
%   ubx(i), this is more efficient);
% AI: mi x n matrix used to defined the linear inequality constraints
%   (see the model above);
% lbi, ubi: mi vectors used to specify the lower (lbi) and upper (ubi)
%   bounds on AI*x (see the model above); a component of lbi (resp. ubi)
%   having the infinite value -inf (resp. inf) indicates that the
%   corresponding bound is not present; there must hold ubi-lbi >=
%   2*options.dymin (transform inequality constraints lbi(i) <=
%   AI(i,:)*x <= ubi(i) into an equality constraint AI(i,:)*x = lbi(i)
%   instead of imposing lbi(i) = ubi(i), this is more efficient); the
%   two bounds lbi(i) and ubi(i) cannot be infinite simultaneously
%   (remove inequality constraint i instead, this is more efficient);
% AE: me x n matrix used to defined the linear equality constraints (see
%   the model above);
% be: me vector used to defined the right hand side of the linear
%   equality constraints (see the model above);
% x0: n vector, specifying an initial guess of the primal solution (this
%   may be set to empty);
% lm0 is an n+mi+me vector, giving an initial guess of the dual solution
%   or Lagrange multiplier associated with the bound constraints
%   (components 1:n), the linear inequality contraints (components
%   n+1:n+mi) and the equality contraints (components n+mi+1:n+mi+me);
% options: rich structure of values, which can be used for tuning the
%   behavior of the solver; a default value is provided if an option is
%   not set; keep the default value of an option whose meaning is not
%   clear to you (many options are provided for the specialists in
%   numerical optimization);
%
%   options.accu: a number in [0,1], specifying the required accuracy of
%     the computed residual at the final iteration, compared with the
%     exact residual, when the linear systems are solved by CG
%     iterations (the closer to 1 the more accurate the residual must
%     be, default value is 1.e-1); this is an attempt to control the
%     rounding errors occuring during the CG iterations, due to the
%     ill-conditioning of the linear systems to solve; this
%     ill-conditioning is sometimes introduced by the approach itself,
%     when the augmentation parameter r is getting too large to try to
%     reach the convergence rate specified by options.dcr (it is then
%     decided to prevent increasing r or even decreasing it);
%   options.alit (only used when mi+me > 0): maximum number of augmented
%     Lagrangian iterations (a positive number, default 100);
%   options.cgit: maximum number of conjugate-gradient iterations (a
%     positive number, default inf);
%   options.cg_as: specifies the active-set technique used to defined
%     the CG phases;
%     = 'hit-and-fix': each time a bound in (x,y) is hit, fix the
%       corresponding variables and restart the CG iterations on the new
%       activated face; pursue up to convergence on the last activated
%       face; then do a gradient-projection step.
%     = 'more-toraldo': use a technique similar to the one proposed by
%       Moré-Toraldo in SIAM Journal on Optimization 1 (1991) 93-113; it
%       consists in minimizing the quadratic function on the affine hull
%       of a particular face up to the satisfaction of a criterion, then
%       projecting on the bounds.
%   options.cgprec: specifies which proconditioner must be used for the
%     CG iterations if any;
%     = 'none' (default): CG iterations are not preconditioned;
%     = 'diag': CG iterations are preconditioned using the diagonal of
%       the matrix of linear systems;
%     = 'chol': CG iterations are preconditioned, using the Cholesky
%       factorization of the matrix of linear systems, resulting in
%       solving them in a single iteration (in exact arithmetic);
%     = 'xxxx' (not implemented): CG iterations are preconditioned,
%       using the Cholesky factorization of the Hessian of the augmented
%       Lagrangian Hr = H+r*A'*A; the factorization is updated when r/r-
%       is not in the interval [1/options.rfct,options.rfct], where r-
%       is the augmentation parameter of the last AL iteration at which
%       a factorization was formed;
%     Qpalm can deal with the cases where the preconditioners are
%     singular.
%   options.cput: limit in second on the CPU time allowed to solve the
%     problem (must be > 0, default inf); if exceeded, Qpalm stops
%     without having found the solution; the timer is called each time a
%     new phase is started;
%   options.ubtol: small nonnegative number used to detect a direction
%     of unboundedness (default n*eps); a recession/asymptotic direction
%     of the feasible set with a Rayleigh quotient of the Hessein of the
%     augmented Lagrangian less than options.ubtol is declared to be a
%     direction of unboundedness;
%   options.dcr (only used when mi+me > 0): desired convergence rate (a
%     number in the interval ]0,1[, default 1.e-1); QPALM updates the
%     augmentation parameter r such that the constraints (if
%     options.rctl == 'constraint') or the constraints shifted by the
%     smallest feasible shift (if options.rctl == 'constraint-change')
%     converges to zero at a linear convergence rate equal to
%     options.dcr; a too mall value induces a large augmentation
%     parameter, hence ill-conditioning, that may prevent convergence;
%   options.dcrf (only used when mi+me > 0): fraction of the desired
%     convergence rate (a number in the interval [0,1[, default 1.e-3),
%     which is used to appreciate whether feasibility is impossible and
%     whether rounding error prevails; a small number forces the solver
%     to get feasibility, but this generally induces a large
%     augmentation parameter, hence ill-conditioning, that may prevent
%     convergence;
%   options.dxmin: resolution in x; two points x and x' such that
%     norm(x-x',inf) <= options.dxmin are considered identical (a
%     positive number, default 1.e-15); this threshold is used in
%     linesearch and for detecting the activity of the bounds on x;
%   options.dymin (only used when mi > 0): resolution in y = AI*x; two
%     points y and y' such that norm(y-y',inf) <= options.dymin are
%     considered identical (a positive number, default 1.e-15); this
%     threshold is used for detecting the activity of the bounds on
%     AI*x;
%   options.feas (only used when mi+me > 0): required accuracy on the
%     constraint feasibility (a positive number, default 1.e-8),
%   options.feass (only used when mi+me > 0): required accuracy on the
%     feasibility of the shifted constraints (a positive number, default
%     options.feas),
%   options.fout: output channel; set 1 for the standard ouput (screen,
%     default); for a specific file, use "options.fout = fopen(...)"
%     before calling QPALM;
%   options.gpph (only used when mi+me > 0): maximum number of gradient
%     projection phases (a positive number, default inf);
%   options.gp_init: specifies how to choose the first stepsize trial in
%     the piecewise linesearch of the GP phose:
%     = 'curv-init': use the first curvature if this one is >0 to
%       determine the initial stepsize trial (the default),
%     = 'inc-bkpt': use the first break point at which phi increases;
%   options.gp_pls: type of piecewise linesearch to use in the gradient
%     projection phase:
%     = 'armijo': the Armijo linesearch, which is more robust than the
%       Goldstein ruel (see below) when the rounding errors prevail (the
%       default),
%     = 'goldstein': the Goldstein linesearch;
%   options.gp_succ: specifies whether one or more GP phases must be done
%     successively:
%     = 'one': after each GP phase, not yielding optimality, launch a CG
%       phase (default);
%     = 'mt': do several GP phases successively until optimality is
%       reached or until the Moré-Toraldo test is realized or until a
%       decrease larger than a fraction of the one of the previous
%       minimization phase has been obtained; if optimality is not
%       reached, do a CG phase next;
%   options.gtol: required smallness on the infinite norm of the
%     projected gradient of phi, used to detect optimality when solving
%     the AL subproblems (a positive number, default 1.e-8);
%   options.phase: maximum number of phases (a positive number, default
%     inf); a phase can be a "minimisation phase" on an active face of
%     the feasible set or a "gradient projection phase", which selects
%     another face to explore; the algorithm alternates between these
%     phases;
%   options.rctl (only used when mi+me > 0): type of control for the
%     augmented parameter r; this option has also an impact on the
%     stopping criteria used by the solver;
%     = 'fixed': r is maintained fixed to its initial value options.rlag
%       during all the run; this option is not recommended but may be
%       useful for an algorithmician wanting to test the solver or the
%       AL algorithm; the stopping criterion focuses on the constraint
%       values if possible and, otherwise, on the shifted constraint
%       values;
%     = 'constraint': r is updated to have the desired convergence rate
%       options.dcr of the constraint values, by examining the ratio
%       |c|/|c_| of the Euclidean norm of the constraint values at the
%       current ieration (i.e., |c|) and the previous iteration (i.e.,
%       |c_|); if the problem is infeasible r will blow up so that this
%       option is not recommended unless the user knows that the problem
%       is feasible, in which case this is the best control type; the
%       stopping criterion focuses on the constraint values;
%     = 'constraint-change' (default): r is updated to have the desired
%       convergence rate options.dcr of the change in the constraint
%       values, by examining |c-c_|/|c_-c__|; this is the recommended
%       control type if the user does not know whether the QP is
%       feasible; the stopping criterion focuses on the constraint
%       values if possible and, otherwise, on the shifted constraint
%       values;
%   options.rfct: when options.cgprec is set to 'chol', the Cholesky
%     factorization of the Hessian of the augmented Lagrangian Hr =
%     H+r*A'*A is updated when r/r- is not in the interval
%     [1/options.rfct,options.rfct], where r- is the augmentation
%     parameter of the last AL iteration at which a factorization was
%     formed;
%   options.rlag (only used when mi+me > 0): initial value of the
%     augmentation parameter (a positive number, default 1);
%   options.verb: verbosity of the solver (a nonnegative number,
%     default 0):
%     == 0: works silently (default);
%     >= 1: error messages;
%     >= 2: + initial setting, final status, and one line per AL
%           iteration (if any);
%     >= 3: + information of the AL subproblem iterations;
%     >= 4: + more details (with possibly long lists);
%     >= 5: + more details (with possibly expensive computations).
%
% Output arguments
% ^^^^^^^^^^^^^^^^
% x: n vector, giving the computed solution (if info.flag == 0);
% lm is an (n+mi+me) vector, giving the computed dual solution or
%   Lagrange multiplier associated with the bounds-on-x in lm(1:n),
%   inequality constraints in lm(n+1:n+mi), and equality contraints in
%   lm(n+mi+1:n+mi+me) (if info.flag = 0); a negative multiplier with
%   component in [1:n+mi] is associated with a lower bound and a
%   positive one is associated with an upper bound;
% info: structure giving various information (only those that are
%   meaningful outside the solver are listed below):
%
%   info.alit = nb of AL iterations realized;
%   info.cgit = total number of conjugate-gradient iterations;
%   info.cgph = total number of conjugate-gradient phases;
%   info.cput = CPU time spent in QPALM (in sec);
%   info.f = (meaningful if info.flag == 0 or 1) value of the cost
%     function at the final point;
%   info.flag (int8) specifies the status of the QP solver on return:
%     ==  0: a solution of the QP is found up to the required accuracy,
%     ==  1: a solution of the closest feasible QP is found (with a
%            nonzero smallest feasible shift) up to the required
%            accuracy,
%     ==  2: the (closest feasible) problem is likely to be unbounded,
%     ==  3: the feasibility tolerance options.feas cannot be reached,
%     ==  4: the optimality tolerance options.gtol cannot be reached,
%     ==  5: erroneous argument,
%     ==  6: failure in the generation of a MEX-file,
%     ==  7: options.alit exceeded without reaching optimality,
%     ==  8: options.phase exceeded without reaching optimality,
%     ==  9: options.gpph exceeded without reaching optimality,
%     == 10: options.cgph exceeded without reaching optimality,
%     == 11: options.cgit exceeded without reaching optimality,
%     == 12: options.cput exceeded without reaching optimality,
%     == 13: H is not positive semidefinite,
%     == 14: CG computations lack of accuracy (measured by options.accu);
%   info.gopt = inf-norm of the projected gradient of the objective;
%   info.gpph = total number of gradient-projection phases;
%   info.rq = (meaningful if info.flag == 6) Rayleigh quotient
%     d'*H*d/(d'*d), where d is a direction of unboundedness;
%   info.s = (meaningful if info.flag == 1) shift vector of length mi;
%     the constraint l(n+1:n+mi) <= AI*x + info.s <= u(n+1:n+mi) is
%     feasible in x (up to machine precision), while the constraint
%     l(n+1:n+mi) <= AI*x <= u(n+1:n+mi) is probably infeasible in x;
%   info.ubdd = (meaningful if info.flag == 5) this is a direction in x
%     along which the closest feasible QP goes to -inf; the direction is
%     normalized, i.e., |info.ubdd| = 1; to improve the feasibility of
%     the direction, decrease the value of options.ubtol.
%
% References
% ^^^^^^^^^^
% [1] F. Delbos, J.Ch. Gilbert (2005). Global linear convergence of an
%     augmented Lagrangian algorithm for solving convex quadratic
%     optimization problems. Journal of Convex Analysis, 12, 45–69.
% [2] A. Chiche, J.Ch. Gilbert (2014). How the augmented Lagrangian
%     algorithm can deal with an infeasible convex quadratic
%     optimization problem. INRIA Research Report RR-8583.
% [3] J.Ch. Gilbert, É. Joannopoulos. OQLA/QPALM - Convex quadratic
%     optimization solvers using the augmented Lagrangian approach, able
%     to deal with infeasibility. INRIA Research Report (to appear).

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

  global QPALM_n QPALM_mi QPALM_me								% dimensions
  global QPALM_g QPALM_H QPALM_lbx QPALM_ubx QPALM_AI QPALM_lbi QPALM_ubi QPALM_AE QPALM_be	% data
  global QPALM_x QPALM_lmx QPALM_lmi QPALM_lme							% optimal variables
  global QPALM_info QPALM_options QPALM_values							% info-options-values

% Initialisation of the timer

  QPALM_info.time0 = clock;

% Set the constant values used in Qpalm

  qpalm_values();
  QPALM_info.flag = QPALM_values.success;

% Set absent input arguments to empty, build globals with frequently used variables, and impose the entry vectors to be column
% vectors

  if nargin < 12, options = struct();
    if nargin < 11, lm0 = [];
      if nargin < 10, x0 = [];
        if nargin < 9, be = [];
          if nargin < 8, AE = [];
            if nargin < 7, ubi = [];
              if nargin < 6, lbi = [];
                if nargin < 5, AI = [];
                  if nargin < 4, ubx = [];
                    if nargin < 3, lbx = [];
                      if nargin < 2
                        fprintf('\n### qpalm: the first 2 arguments are required\n\n');
                        x = x0;
                        lm = lm0;
                        info.flag = QPALM_values.fail_on_argument;
                        return
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end

% Create globals for frequently used variables and impose the entry vectors to be column vectors

  QPALM_g       = g(:);
  QPALM_H       = H;
  QPALM_lbx     = lbx(:);
  QPALM_ubx     = ubx(:);
  QPALM_AI      = AI;
  QPALM_lbi     = lbi(:);
  QPALM_ubi     = ubi(:);
  QPALM_AE      = AE;
  QPALM_be      = be(:);
  QPALM_options = options;

  x0  = x0(:);
  lm0 = lm0(:);

% Preliminaries

  qpalm_prelim (x0,lm0);

  if QPALM_info.flag
    if QPALM_options.verb; fprintf(QPALM_options.fout,'\n%s\n',QPALM_values.eline); end
    x    = QPALM_x;
    lm   = [QPALM_lmx;QPALM_lmi;QPALM_lme];
    info = QPALM_info;
    return
  end

% Solve the QP

  qpalm_loop ();

% Compute the optimality conditions

  qpalm_optim ();
  lm = [QPALM_lmx;QPALM_lmi;QPALM_lme];

% Conclude

  QPALM_info.cput = etime (clock,QPALM_info.time0);
  if QPALM_options.verb >= 2
    qpalm_finish (lm);
  end
  x    = QPALM_x;
  info = QPALM_info;

% Exit

  return
