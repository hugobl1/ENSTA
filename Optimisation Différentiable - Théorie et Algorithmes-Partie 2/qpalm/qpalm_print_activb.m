function [] = qpalm_print_activb (fout,verb)

%%
% [] = qpalm_print_activb (fout,verb)
%
% Print active bounds on the variables.

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

  global QPALM_IB

% Print activity of the bound constraints

  WBl = find(QPALM_IB<0);
  WBln = length(WBl);
  fprintf(fout,'\n  > Active bounds\n  > . %i lower bounds',WBln);
  if WBln & (verb >= 4)
    fprintf(fout,':');
    fprintf(fout,' %0i',WBl);
  end
  WBu = find(QPALM_IB>0);
  WBun = length(WBu);
  fprintf(fout,'\n  > . %i upper bounds',WBun);
  if WBun & (verb >= 4)
    fprintf(fout,':');
    fprintf(fout,' %0i',WBu);
  end

% exit

  return
