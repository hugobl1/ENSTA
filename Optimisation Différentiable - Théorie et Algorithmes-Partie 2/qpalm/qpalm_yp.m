function [] = qpalm_yp ()

%%
% [] = qpalm_yp ()
%
% Qpalm_yp computes
%
%   yp = projetion of (AI*x+lmi/r) on [lbi,ubi],
%
% It is assumed that size(AI,1) is > 0.
%
% On entry
% ^^^^^^^^
% QPALM_lmi (global): multipliers of the inequality constraints
%
% On return
% ^^^^^^^^^
% QPALM_y (global): the projetion of (AI*x+lmi/r) on [lbi,ubi];

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

  global QPALM_mi					% dimensions
  global QPALM_AI QPALM_lbi QPALM_ubi			% data
  global QPALM_x QPALM_y QPALM_lmi			% variables
  global QPALM_info QPALM_options			% info-options-values
  global QPALM_II

% Take into account the inequality constraints

  aix = QPALM_AI*QPALM_x;					% AI*x
  QPALM_y = aix+QPALM_lmi/QPALM_info.rlag;			% AI*x + lmi/r
  Il = find(QPALM_y<=QPALM_lbi+QPALM_options.dymin);		% indices of the variables yp at the lower bound
  Iu = find(QPALM_y>=QPALM_ubi-QPALM_options.dymin);		% indices of the variables yp at the lower bound
  QPALM_y(Il) = QPALM_lbi(Il); QPALM_y(Iu) = QPALM_ubi(Iu);	% yp = projetion of (AI*x+lmi/r) on [lbi,ubi]
  QPALM_II     = zeros(QPALM_mi,1);
  QPALM_II(Il) = -1;
  QPALM_II(Iu) = +1;

% Return

  return

end
