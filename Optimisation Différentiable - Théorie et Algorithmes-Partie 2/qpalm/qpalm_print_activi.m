function [] = qpalm_print_activi (fout,verb)

%%
% [] = qpalm_print_activi (fout,verb)
%
% Print active inequality bounds.

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

  global QPALM_II

% Print activity of the inequaliti constraints

  WIl = find(QPALM_II<0);
  WIln = length(WIl);
  fprintf(fout,'\n  > Active inequalities\n  > . %i lower bounds',WIln);
  if WIln & (verb >= 4)
    fprintf(fout,':');
    fprintf(fout,' %0i',WIl);
  end
  WIu = find(QPALM_II>0);
  WIun = length(WIu);
  fprintf(fout,'\n  > . %i upper bounds',WIun);
  if WIun & (verb >= 4)
    fprintf(fout,':');
    fprintf(fout,' %0i',WIu);
  end

% exit

  return
