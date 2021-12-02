function [x,info] = qpalm_chol (g,H,info,options,values)

%%
% [x,info] = qpalm_chol (g,H,info,options,values)
%
% Qpalm_chol uses a Cholesky factorization of the positive semi-definite
% matrix H, if possible, to solve a convex quadratic optimization
% problem written as follows
%
%   min (in x)  g'*x + 0.5*x'*H*x
%
% If the problem has a solution, this one is returned in x. Otherwise,
% the problem is unbounded and qpalm_chol returns x=zeros(n,1) and a
% direction info.ubdd along which the cost goes to -Inf.

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

  global QPALM_n

% Set the output variables

  x = [];

% Initial printings

  if options.verb >= 3
    fprintf(options.fout,'\n  > %i variables',QPALM_n);
  end

% Compute the Cholesky factorization of H = L*L'

  % It is better to try 'chol' first, since it is much faster than a Matlab implementation (below); 'chol' assumes the positive
  % definiteness of the matrix however and works without pivoting

  [L,p] = chol(H,'lower');

  if ~p	% H is positive definite

    % Compute x = -H\g using the Cholesky factorization of H = L*L'
  
    if options.verb >= 3
      fprintf(options.fout,'\n  > LS solved by a Cholesky factorization',QPALM_n);
    end

    opts.LT = true;		% L is lower triangular
    x = linsolve(L,-g,opts);	% x = -L\g
    opts.TRANSA = true;		% use L'
    x = linsolve(L,x,opts); 	% x = -L'\(L\g)

  else	% H is probably not positive definite

    % Try using 'qpalm_cholp', which is a personnal Matlab implementation of a Cholesky factorization with pivoting (such a
    % factorization does not exist in Matlab, 2012). Deduce the solution if it finds that the rank of H is n; or compute a
    % direction of unboundedness, otherwise.

    if options.verb >= 3
      fprintf ('\n  > H is not positive definite to working precision');
      fprintf ('\n  > ... trying to recover by pivoting (Matlab coding, slow for large problems)');
    end

    [L,piv,rank] = qpalm_cholp (H,options.ubtol);	% Cholesky factorization of P*H*P' = L*L' with pivoting matrix P

    if options.verb >= 3
      fprintf(options.fout,'\n  > LS solved by a Cholesky factorization, with pivoting\n  > estimated rank deficiency: ',QPALM_n);
    end

    if rank == QPALM_n

      if options.verb >= 3
        fprintf(options.fout,'0');
      end

      % Compute x = -H\g using the Cholesky factorization of P*H*P' = L*L'

      x = -g(piv);		% x = -P*g
      x = x(:);			% make sure that x is a column vector
      opts.LT = true;		% use the lower triangular part of L in linsolve
      x = linsolve(L,x,opts);	% x = -L\(P*g)
      opts.TRANSA = true;	% use L' instead of L in linsolve
      x = linsolve(L,x,opts);	% x = -L'\(L\(P*g))
      x(piv) = x; 		% x = -P'*L'\(L\(P*g))) = -H\g

    else

      if options.verb >= 3
        fprintf(options.fout,'%0i',QPALM_n-rank);
      end

      % Compute a cheap direction of unboundedness P' * [d1; d2]
%fprintf(options.fout,'\nP =');
%fprintf(options.fout,'\n  %12.5e',piv);
      L11 = L(1:rank,1:rank);
      L21 = L(rank+1:QPALM_n,1:rank);
      opts.LT = true;
      opts.TRANSA = false;
      d1 = linsolve(L11,-g(1:rank),opts);	% d1 = - L11\g(1:rank)
      d2 = L21*d1-g(rank+1:QPALM_n);		% d2 = - L21*(L11\g(1:rank)) - g(rank+1:n)
      opts.TRANSA = true;
      d1 = linsolve(L11,L21'*d2,opts);		% d1 = L11'\(L21'*d2)
      info.ubdd = [d1; d2];			% permuted direction of unboundedness [d1; d2]
      info.ubdd = info.ubdd/norm(info.ubdd);	% normalized permuted direction of unboundedness [d1; d2]
      info.ubdd(piv) = info.ubdd;		% return P'*[d1; d2]

      % The "solution" is zero by convention

      x = zeros(QPALM_n,1);
      info.f = 0;
      info.gopt = norm(g,inf);
      info.flag = values.stop_on_unboundedness;

      return

    end

  end

% Before returning with a solution x

  r = g+H*x;  			% residual

  info.f = 0.5*x'*(g+r);
  info.gopt = norm(r,inf);
  info.flag = 0;

% Return

  return
