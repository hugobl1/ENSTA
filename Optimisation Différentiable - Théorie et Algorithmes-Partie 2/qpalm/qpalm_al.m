function [val,grad] = qpalm_al ()

%%
% [val]      = qpalm_al ()
% [val,grad] = qpalm_al ()
%
% Qpalm_al computes the value 'val' of the augmaented Lagrangian 'lr'
% with parameter 'r', defined by
%
%   lr(x,y,lm) = g'*x + 0.5*x'*H*x + lmi'*(AI*x-y)  + 0.5*r*norm(AI*x-y)^2
%                                  + lme'*(AE*x-be) + 0.5*r*norm(AE*x-be)^2.
%
% When required it also computes the gradient 'grad' of lr wrt x.
% 
%
% On entry
% ^^^^^^^^
% QPALM_lmi (global): multipliers of the inequality constraints
% QPALM_lme (global): multipliers of the equality constraints

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
  global QPALM_g QPALM_H QPALM_AI QPALM_AE QPALM_be
  global QPALM_info
  global QPALM_x QPALM_y QPALM_lmi QPALM_lme

% Initialization

  grad = [];

  mie = QPALM_mi+QPALM_me; 
  if mie; r = QPALM_info.rlag; end

% Unconstrained case

  grad = QPALM_g+QPALM_H*QPALM_x;
  val = 0.5*QPALM_x'*(QPALM_g+grad);

% Take into account the inequality constraints

  if QPALM_mi
    aix = QPALM_AI*QPALM_x;			% AI*x
    ri = aix-QPALM_y;				% inequality constraint residual
    val = val + QPALM_lmi'*ri + 0.5*r*(ri'*ri);	% phi value update
    if nargout > 1
      grad = grad + QPALM_AI'*(QPALM_lmi+r*ri);	% phi gradient update
    end
  end

% Take into account the equality constraints

  if QPALM_me
    re = QPALM_AE*QPALM_x-QPALM_be;		% equality constraint residual
    val = val + QPALM_lme'*re + 0.5*r*(re'*re);	% phi value update
    if nargout > 1
      grad = grad + QPALM_AE'*(QPALM_lme+r*re);	% phi gradient update
    end
  end

% Return

  return

end
