function [L, d, flag] = cholmod (M, small, big)

%
% [L, d] = cholmod (M, small, big)
% [L, d, flag] = cholmod (M, small, big)
%
% Computes a modified Cholesky factorization of the symmetric matrix M,
% such that M + E = L*D*L', where L is unit lower triangular and D =
% diag(d), and E is some nonnegative diagonal matrix. The matrix E is
% computed as small as possible such that d(i) >= small (for all i) and
% abs(L(i,j)*sqrt(d(j))) <= big (for all i>j). It is required to have
% sqrt(small) <= big. The algorithm uses no permutation.
%
% If present, the output argument flag describes the computation:
%   0 (computation is fine),
%   1 (M is not square or empty),
%   2 (small<=0 or big<=0 or sqrt(small)>big).
%
% Author: J.Ch. Gilbert (Inria), following the description in the book
% by J. Nocedal and S. Wright (1999), "Numerical Optimization", pg
% 145-149.

% Initialize output arguments

  L = [];
  d = [];
  if (nargout > 2) flag = 0; end

% Check input arguments

  n = size(M,1);
  if (n == 0) | (size(M,2) ~= n)
    if (nargout > 2) flag = 1; end
    return
  end
  if (small <= 0) | (big <= 0) | (sqrt(small) > big)
    if (nargout > 2) flag = 2; end
    return
  end

% The easy case

  if n == 1
    L = 1;
    d = small;
    return
  end

% Compute L and d column by column

  L = eye(n);

  for j=1:n

    % dd = d(j)

    dd = M(j,j);
    for k = 1:j-1
      dd = dd - L(j,k)^2*d(k);
    end

    % L(i,j) and theta (depends on j)

    theta = 0;
    for i = j+1:n
      aux = M(i,j);
      for k = 1:j-1
        aux = aux - L(i,k)*L(j,k)*d(k);
      end
      L(i,j) = aux;
      theta = max(theta, abs(aux));
    end

    % Reset dd and set d(j)

    dd = max([abs(dd), (theta/big)^2, small]);
    d(j) = dd;

    % Reset L(i,j)

    for i = j+1:n
      L(i,j) = L(i,j)/dd;
    end

  end

return
