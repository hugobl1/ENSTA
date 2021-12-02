function [prec] = qpalm_prup (prec,VB,Bhit,WI,Ihit)

%%
% [prec] = qpalm_prup (prec,VB,Bhit,WI,Ihit)
%
% Qpalm_prup updates the preconditioner 'prec' for the CG iterations,
% depending on type of preconditioner specified by options.cgprec.
%
% On entry
% ^^^^^^^^
% VB: indices of the free variables x;
% Bhit: new hitten bounds in x;
% WI: indices of the free variables yp;
% Ihit: new hitten bounds in y;
%
% On return
% ^^^^^^^^^
% prec: represents the new preconditioner in a variable or structure
%   that depends on the type of preconditioner specified by
%   options.cgprec. The following matrix is used:
%
%     M(VB,WI) := H(VB,VB)+r*AE(:,VB)'*AE(:,VB)+r*AI(WI,VB)'*AI(WI,VB),
%
%   where VB = setminus(VB_on_entry,Bhit), WI = union(WI_on_entry,Ihit);
%   if options.cgprec = 'diag': prec.diag is the diagonal of M(VB,WI).

% Global variables used
% ^^^^^^^^^^^^^^^^^^^^^
% QPALM_info: structure giving various information; those that are used
%   in this function are the following
%
%   QPALM_info.rlag: value of the augmentation parameter

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

  global QPALM_AI QPALM_AE
  global QPALM_info QPALM_options QPALM_values

% Initializations

  mi = size(QPALM_AI,1);
  r  = QPALM_info.rlag;

% Diagonal preconditioner

  if QPALM_options.cgprec == QPALM_values.cg_prec_diag

    for i = Bhit(:)'		% make sure the list is horizontal
      k = find(VB==i);
      prec.diag(k) = [];
      VB(k) = [];
    end
    for i = Ihit(:)'		% make sure the list is horizontal
      prec.diag = prec.diag + r*(QPALM_AI(i,VB).*QPALM_AI(i,VB))';
    end

% Cholesky preconditioner

  elseif QPALM_options.cgprec == QPALM_values.cg_prec_chol

    if isempty(prec.piv)	% the Cholesky factorization uses no pivoting matrix

      nvb = size(prec.L,1);	% order of L

      % consider the case when x-bounds have been hit

      for i = Bhit(:)'		% make sure the list is horizontal
        k = find(VB==i);	% index of the rwo/column to remove from M = L*L'
        prec.L = [prec.L(1:k-1  ,1:k-1) , zeros(k-1,nvb-k);
                  prec.L(k+1:nvb,1:k-1) , cholupdate(prec.L(k+1:nvb,k+1:nvb)',prec.L(k+1:nvb,k))'];
        VB(k) = [];
        nvb = nvb-1;
      end
      prec.rank = prec.rank-length(Bhit);

      % consider the case when y-bounds have been hit

      for i = Ihit(:)'		% make sure the list is horizontal
        prec.L = cholupdate(prec.L',sqrt(r)*QPALM_AI(i,VB)')';
      end

    else			% the Cholesky factorization uses the pivoting vector prec.piv and may be singular

      rk = prec.rank;
      nv = length(prec.piv);
      small = QPALM_options.ubtol*prec.L(1,1);	% threshold for detecting a small scalar

      % consider the case when x-bounds have been hit

%Bhit
      for i = Bhit(:)'		% make sure the list is horizontal
        k = find(VB==i);	% index of the rwo/column to remove from M = P'*L*L'*P (P is the permutation matrix)
        l = find(prec.piv==k);	% index of the rwo/column to remove from L*L'
        if l < rk
          v = prec.L(l+1:nv,l);
          if norm(v) <= small*sqrt(nv-l)	% if v == 0, remove row/column l from L, rank decreased by 1
            prec.L = [prec.L(1:l-1 ,1:l-1) , zeros(l-1 ,nv-l);
                      prec.L(l+1:nv,1:l-1) , prec.L(l+1:nv,l+1:nv)];
            rk = rk-1;
          else					% remove row/column l from L, update L22, pivot 2 rows if appropriate, rank may change
            [L22,r22] = qpalm_cholup (prec.L(l+1:nv,l+1:nv),rk-l,v,small);
            prec.L = [prec.L(1:l-1 ,1:l-1) , zeros(l-1,nv-l);
                      prec.L(l+1:nv,1:l-1) , L22];
            rk = l-1+r22;
          end
        elseif l == rk
          if l == nv		% simple case: remove the last row/column from L, rank decreased by 1
            prec.L = prec.L(1:nv-1,1:nv-1);
            rk = rk-1;
          else
            v = prec.L(l+1:nv,l);
            if norm(v) <= small*sqrt(nv-l)	% if v == 0, remove row/column l from L, rank decreased by 1
              prec.L = [prec.L(1:l-1 ,1:l-1) , zeros(l-1 ,nv-l);
                        prec.L(l+1:nv,1:l-1) , zeros(nv-l,nv-l)];
              rk = rk-1;
            else				% remove row l and the zero column l+1 from L, pivot 2 rows if appropriate, rank unchanged
              prec.L = [prec.L(1:l-1 ,1:l) , zeros(l-1 ,nv-l-1);
                        prec.L(l+1:nv,1:l) , zeros(nv-l,nv-l-1)];
              [vmax,I] = max(v);
              j = I(1);
              if j ~= 1					% pivot row l and l+j-1 of L
		pivlp1   = prec.piv(l+1);		% mind that piv has not been reduced in size yet (done below)
                prec.piv(l+1) = prec.piv(l+j);
                prec.piv(l+j) = pivlp1;
                range = 1:nv-1;				% mind that L has been reduced in size already
                rowl                = prec.L(l,range);
                prec.L(l,range)     = prec.L(l+j-1,range);
                prec.L(l+j-1,range) = rowl;
              end
            end
          end
        else			% here l > rk, simple case: remove row/column l from L, rank unchanged
          prec.L = [prec.L(1:l-1 ,1:l-1) , zeros(l-1 ,nv-l);
                    prec.L(l+1:nv,1:l-1) , zeros(nv-l,nv-l)];
        end
        prec.piv(l) = [];		% remove row l from P
        K = find(prec.piv>k);
        prec.piv(K) = prec.piv(K)-1;	% remove column k from P
        VB(k) = [];			% remove variable i from VB
        nv = nv-1;
      end
      prec.rank = rk;

      % consider the case when y-bounds have been hit

      for i = Ihit(:)'		% make sure the list is horizontal
	v = sqrt(r)*QPALM_AI(i,VB)';
	v = v(prec.piv);	% P*v
        [prec.L,rk,piv] = qpalm_cholup (prec.L,rk,v,small);
	if ~isempty(piv)
check7	% this section has still not been verified on an example, send a mail with your data to Jean-Charles.Gilbert@inria.fr who will do it
          prec.piv = prec.piv(piv);
	end
      end
      prec.rank = rk;

    end

  end

% Exit

  return
