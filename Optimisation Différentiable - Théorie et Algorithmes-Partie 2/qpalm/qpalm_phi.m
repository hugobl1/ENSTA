function [val,grad] = qpalm_phi (withii,x)

%%
% [val,grad]    = qpalm_phi (withii,x)
%
% The input argument 'x' and the output argument 'grad' are optional.
%
% Qpalm_phi computes the value 'val' and the gradient 'grad' (if the
% argument is present) of the differentiable function
%
%   phi : x -> inf {lr(x,y,lm): y in [lbi,ubi]},
%
% where 'lr' is the augmented Lagrangian with parameter 'r', defined by
%
%   lr(x,y,lm) = g'*x + 0.5*x'*H*x + lmi'*(AI*x-y)  + 0.5*r*norm(AI*x-y)^2
%                                  + lme'*(AE*x-be) + 0.5*r*norm(AE*x-be)^2.
%
% Qpalm_phi also computes the global QPALM_y as the projection of
% (AI*x+lmi/r) on [lbi,ubi].
%
% If the input argument 'x' is present, the evaluations are made at the
% given 'x', otherwise the global 'QPALM_x' is used.
%
% The argument 'withii' specifies whether the inequality activity (in
% the global 'QPALM_II') must be recomputed (this is the case if withii
% == true).
%
% On entry
% ^^^^^^^^
% withii: logical that specifies whether the inequality activity (in
%   QPALM_II) must be recomputed (withii == true);
% QPALM_lmi (global): multipliers of the inequality constraints;
% QPALM_lme (global): multipliers of the equality constraints.
%
% On return
% ^^^^^^^^^
% QPALM_y (global): the projetion of (AI*x+lmi/r) on [lbi,ubi].
% QPALM_II (global): activity in the inequalities (recomputed when
%   withii == true)

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
  global QPALM_g QPALM_H QPALM_AI QPALM_lbi QPALM_ubi QPALM_AE QPALM_be
  global QPALM_x QPALM_y QPALM_lmi QPALM_lme
  global QPALM_info QPALM_options
  global QPALM_II

% Initialization

  if QPALM_mi+QPALM_me; r = QPALM_info.rlag; end

% Unconstrained case

  if nargin == 1
    grad = QPALM_g+QPALM_H*QPALM_x;
    val = 0.5*QPALM_x'*(QPALM_g+grad);
  else
    grad = QPALM_g+QPALM_H*x;
    val = 0.5*x'*(QPALM_g+grad);
  end

% Take into account the inequality constraints

  if QPALM_mi
    if nargin == 1
      aix = QPALM_AI*QPALM_x;					% AI*x
    else
      aix = QPALM_AI*x;						% AI*x
    end
    if withii							% computation of yp and II
      QPALM_y = aix+QPALM_lmi/r;				% AI*x + lmi/r
      Il = find(QPALM_y<=QPALM_lbi+QPALM_options.dymin);	% indices of the variables yp at the lower bound
      Iu = find(QPALM_y>=QPALM_ubi-QPALM_options.dymin);	% indices of the variables yp at the lower bound
      QPALM_y(Il) = QPALM_lbi(Il); QPALM_y(Iu) = QPALM_ubi(Iu);	% yp = projetion of (AI*x+lmi/r) on [lbi,ubi]
      QPALM_II = zeros(QPALM_mi,1);
      QPALM_II(Il) = -1;
      QPALM_II(Iu) = +1;
    else							% computation of yp only
      QPALM_y = max(QPALM_lbi,min(QPALM_ubi,aix+QPALM_lmi/r));	% yp = projetion of (AI*x+lmi/r) on [lbi,ubi]
    end
    ri = aix-QPALM_y;				% inequality constraint residual
    val = val + QPALM_lmi'*ri + 0.5*r*(ri'*ri);	% phi value update
    if nargout > 1
      grad = grad + QPALM_AI'*(QPALM_lmi+r*ri);	% phi gradient update
    end
  end

% Take into account the equality constraints

  if QPALM_me
    if nargin == 1
      re = QPALM_AE*QPALM_x-QPALM_be;		% equality constraint residual
    else
      re = QPALM_AE*x-QPALM_be;			% equality constraint residual
    end
    val = val + QPALM_lme'*re + 0.5*r*(re'*re);	% phi value update
    if nargout > 1
      grad = grad + QPALM_AE'*(QPALM_lme+r*re);	% phi gradient update
    end
  end

% Return

  return

end
