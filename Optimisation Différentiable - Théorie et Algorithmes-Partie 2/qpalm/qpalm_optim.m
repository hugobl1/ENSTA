function [] = qpalm_optim ();

%%
% [] = qpalm_optim ();
%
% Global variables
% ^^^^^^^^^^^^^^^^
% QPALM_lmx (return): multipliers associated with the bound constraints.
% QPALM_info.feas (return):
% QPALM_info.flag (return):
% QPALM_info.pgln (return):

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
  global QPALM_g QPALM_H QPALM_lbx QPALM_ubx QPALM_AI QPALM_lbi QPALM_ubi QPALM_AE QPALM_be
  global QPALM_info QPALM_options QPALM_values
  global QPALM_x QPALM_lmx QPALM_lmi QPALM_lme

% Initialization

  pgln = [];

% Feasibility

  feas = 0;
  if QPALM_nb
    feas = max(feas,max(QPALM_lbx-QPALM_x));
    feas = max(feas,max(QPALM_x-QPALM_ubx));
  end
  if QPALM_mi
    AIx = QPALM_AI*QPALM_x;
    feas = max(feas,max(QPALM_lbi-AIx));
    feas = max(feas,max(AIx-QPALM_ubi));
  end
  if QPALM_me
    feas = max(feas,norm(QPALM_AE*QPALM_x-QPALM_be,inf));
  end
  QPALM_info.feas = feas;

% Gradient of the Lagrangian

  pglx = QPALM_g+QPALM_H*QPALM_x;
  if QPALM_mi; pglx = pglx+QPALM_AI'*QPALM_lmi; end
  if QPALM_me; pglx = pglx+QPALM_AE'*QPALM_lme; end

% Multiplier associated with the bounds

  if (QPALM_info.flag == QPALM_values.success) | (QPALM_info.flag == QPALM_values.success_shifted) | (QPALM_info.flag == QPALM_values.ctol_unreachable)
    QPALM_lmx = -pglx;
    I = find((QPALM_lbx+QPALM_options.dxmin < QPALM_x) & (QPALM_x < QPALM_ubx-QPALM_options.dxmin));	% indices of the free variables x
    QPALM_lmx(I) = 0;
  else
    QPALM_lmx = NaN*ones(QPALM_n,1);
  end

% Projected gradient of the Lagrangian in x

  % Norm of the projected gradient of the Lagrangian in x
  %              { min(0,glx(i))  if lbx(i) = x(i)
  %    pglx(i) = { gphi(i)        if lbx(i) < x(i) < ubx(i)
  %              { max(0,glx(i))  if          x(i) = ubx(i).

  if QPALM_nb
    Bl = find(QPALM_x<=QPALM_lbx+QPALM_options.dxmin);	% indices of the variables x at the lower bound
    pglx(Bl) = min(0,pglx(Bl));
    Bu = find(QPALM_ubx-QPALM_options.dxmin<=QPALM_x);	% indices of the variables x at the upper bound
    pglx(Bu) = max(0,pglx(Bu));
  end

  pgln = norm(pglx,inf);

  % Add the norm of the projected gradient of the Lagrangian in y

  if QPALM_mi
    Il = find(AIx<=QPALM_lbi+QPALM_options.dymin);	% indices of the variables y at the lower bound
    pgln = max(pgln,norm(max(0,QPALM_lmi(Il)),inf));
    Iu = find(QPALM_ubi-QPALM_options.dymin<=AIx);	% indices of the variables y at the upper bound
    pgln = max(pgln,norm(min(0,QPALM_lmi(Iu)),inf));
  end

  QPALM_info.pgln = pgln;

% Take another chance of declaring success, when this new way of computing optimality is successful

  if (QPALM_info.flag ~= QPALM_values.success) & (pgln <= QPALM_options.gtol) & (feas <= QPALM_options.feas)
    QPALM_info.flag = QPALM_values.success;
  end

% Exit

return
