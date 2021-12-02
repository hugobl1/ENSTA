function [alphamax,IBl,IBu,IIl,IIu] = qpalm_maxstep (d,q,VB,VI)

%%
% [alphamax,IBl,IBu,IIl,IIu] = qpalm_maxstep (d,q,VB,VI)
%
% Computes the largest stepsize 'alphamax' from 'x' in the direction 'd'
% such that x+alphamax*d remains in [lbx,ubx] and y+AI*d remains in
% [lbi,ubi]. If a bound is not encountered, alphamax = inf.
%
% On return
% ^^^^^^^^^
% IIl list of the indices i in VI such that (y+alphamax*AI*d)(i) is at
%   the lower bound lbi;
% IIu list of the indices i in VI such that (y+alphamax*AI*d)(i) is at
%   the upper bound ubi;
% IBl list of the indices i in VB such that (x+alphamax*d)(i) is at the
%   lower bound lbx;
% IBu list of the indices i in VB such that (x+alphamax*d)(i) is at the
%   upper bound ubx.

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

  global QPALM_nb QPALM_mi
  global QPALM_lbx QPALM_ubx QPALM_lbi QPALM_ubi
  global QPALM_info
  global QPALM_x QPALM_y

% initialization

  alphamax = inf;	% a priori infinite stepsize is allowed
  IIl = [];		% lower y-bound (lbi) hit at alphamax
  IIu = [];		% upper y-bound (ubi) hit at alphamax
  IBl = [];		% lower x-bound (lbx) hit at alphamax
  IBu = [];		% upper x-bound (ubx) hit at alphamax

% compute alphamax

  if QPALM_nb & ~isempty(VB)
%fprintf('\n\n i      l(i)          x(i)          u(i)          d(i)        alpha max    (x+am*d)(i)\n');
    for i = VB(:)'		% make sure that the list is an horizontal vector
%fprintf('%2i  %12.5e  %12.5e  %12.5e  %12.5e',i,QPALM_lbx(i),QPALM_x(i),QPALM_ubx(i),d(i));
      if ((QPALM_lbx(i) > -inf) & (d(i) < 0))
        alpham = (QPALM_lbx(i)-QPALM_x(i))/d(i);
%fprintf('  %12.5e  %12.5e',alpham,QPALM_x(i)+alpham*d(i));
        if alpham < alphamax
          IBl = i;
          IBu = [];
          alphamax = alpham;
        elseif alpham == alphamax
          IBl = [IBl,i];
        end
      elseif ((QPALM_ubx(i) < inf) & (d(i) > 0))
        alpham = (QPALM_ubx(i)-QPALM_x(i))/d(i);
%fprintf('  %12.5e  %12.5e',alpham,QPALM_x(i)+alpham*d(i));
        if alpham < alphamax
          IBl = [];
          IBu = i;
          alphamax = alpham;
        elseif alpham == alphamax
          IBu = [IBu,i];
        end
      end
%fprintf('\n');
    end
  end

  if QPALM_mi & ~isempty(VI)
%fprintf('\n\n i      l(i)          y(i)          u(i)          q(i)        alpha max    (y+am*q)(i)\n');
    for i = VI'		% the 'for' loop requires an horizontal vector
%fprintf('%2i  %12.5e  %12.5e  %12.5e  %12.5e',i,QPALM_lbi(i),QPALM_y(i),QPALM_ubi(i),q(i));
      if ((QPALM_lbi(i) > -inf) & (q(i) < 0))
        alpham = (QPALM_lbi(i)-QPALM_y(i))/q(i);
%fprintf('  %12.5e  %12.5e',alpham,QPALM_y(i)+alpham*q(i));
        if alpham < alphamax
          IBl = [];
          IBu = [];
          IIl = i;
          IIu = [];
          alphamax = alpham;
        elseif alpham == alphamax
          IIl = [IIl,i];
        end
      elseif ((QPALM_ubi(i) < inf) & (q(i) > 0))
        alpham = (QPALM_ubi(i)-QPALM_y(i))/q(i);
%fprintf('  %12.5e  %12.5e',alpham,QPALM_y(i)+alpham*q(i));
        if alpham < alphamax
          IBl = [];
          IBu = [];
          IIl = [];
          IIu = i;
          alphamax = alpham;
        elseif alpham == alphamax
          IIu = [IIu,i];
        end
      end
%fprintf('\n');
    end
  end
%alphamax

% Exit

  return
