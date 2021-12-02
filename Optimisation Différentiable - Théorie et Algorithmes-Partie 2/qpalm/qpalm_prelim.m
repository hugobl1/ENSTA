function [] = qpalm_prelim (x0,lm0);

%%
% [] = qpalm_prelim (x0,lm0);
%
% This function realizes the following preliminary tasks:
% - set QPALM_options (default values if absent, numeric values for
%   lexical options)
% - check the given QPALM_options
% - get the dimensions
% - initial printings

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

% Warning off (reset 'on' in qpalm_finish)

  warning ('off','MATLAB:nearlySingularMatrix');	% 'Matrix is close to singular or badly scaled' from linsolve

% Check validity of the options only if some of them have been set by the user; here 'options' is most likely to be a
% struture (with possibily no field) or an empty variable

  checkoptions = 1;
  if isempty(QPALM_options) | isempty(fieldnames(QPALM_options))
    checkoptions = 0;
  end

% Set fout and verb

  if exist('QPALM_options')
    if ~isstruct(QPALM_options)
      warning('qpalm:OptionsNotStruct','argument ''options'' is expected to be a structure\n');
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
    if isfield(QPALM_options,'fout')
      if isempty(QPALM_options.fout) | (QPALM_options.fout <= 0) | (fix(QPALM_options.fout)-QPALM_options.fout)
        if QPALM_options.verb, fprintf('\n\n### qpalm_prelim: options.fout = %i not recongnized\n',QPALM_options.fout); end
        QPALM_info.flag = QPALM_values.fail_on_argument;
        return
      end
    end
  end
  if ~isfield(QPALM_options,'fout')
    QPALM_options.fout = 1;			% default output channel
  end
  fout = QPALM_options.fout;

  if isfield(QPALM_options,'verb')
    if isempty(QPALM_options.verb) | (QPALM_options.verb < 0)
      if QPALM_options.verb, fprintf(fout,'\n\n### qpalm_prelim: options.verb must be nonnegative\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.verb = 1;			% default verbosity level
  end
  verb = QPALM_options.verb;			% to reduce time access

% Check options fields

  possible_options = {'accu', 'alit', 'cg_as', 'cgit', 'cgph', 'cgprec', 'cput', 'ubtol', 'dcr', 'dcrf', 'dxmin', 'dymin' 'feas', ...
    'feass', 'fout', 'gpph', 'gp_init', 'gp_pls', 'gp_succ', 'gtol', 'phase', 'rctl', 'rfct', 'rlag', 'verb', 'xxx'};
  actual_options = fieldnames(QPALM_options);
  actual_options_char = char(actual_options);
  for i=1:length(actual_options)
    if ~ismember(actual_options(i,:),possible_options)
      if verb
        fprintf(fout,'\n\n### qpalm_prelim: field ''%s'' of structure ''options'' is not recongnized\n', ...
	  actual_options_char(i,:));
      end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  end

% Decode entry arguments and check consistency

  QPALM_n = length(QPALM_g);
  if verb >= 2
    fprintf(fout,'\n%s',QPALM_values.eline);
    fprintf(fout,'\nQPALM solver (version 0.5, November 2014)');
    fprintf(fout,'\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^');

    fprintf(fout,'\nDescription of the problem');
    fprintf(fout,'\n  Nb of variables %14i',QPALM_n);
  end

  if (ndims(QPALM_g) > 2) | (size(QPALM_g,2) ~= 1) | (QPALM_n <= 0)
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent length for ''g'' (should be > 0)\n'); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end
  if (ndims(QPALM_H) > 2) | (size(QPALM_H,1) ~= QPALM_n) | (size(QPALM_H,2) ~= QPALM_n);
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent dimensions for ''H'' (should be %0i x %0i)\n',QPALM_n,QPALM_n); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end
  hnorm = norm(QPALM_H,'fro');
  if (hnorm > 0) & (norm(QPALM_H-QPALM_H','fro')/hnorm > 100*eps);
    warning('qpalm:NonSymmetricHessian','H doesn''t look symmetric, resetting H = (H+H'')/2\n\n');
    QPALM_H = (QPALM_H+QPALM_H')*0.5;
  end

  if isempty(QPALM_lbx)
    QPALM_lbx = -inf*ones(QPALM_n,1);
  elseif (ndims(QPALM_lbx) > 2) | (size(QPALM_lbx,1) ~= QPALM_n) | (size(QPALM_lbx,2) ~= 1)
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent length for ''lbx'' (should be %0i)\n',QPALM_n); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end
  if isempty(QPALM_ubx)
    QPALM_ubx = inf*ones(QPALM_n,1);
  elseif (ndims(QPALM_ubx) > 2) | (size(QPALM_ubx,1) ~= QPALM_n) | (size(QPALM_ubx,2) ~= 1)
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent length for ''ubx'' (should be %0i)\n',QPALM_n); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end
  QPALM_nb = sum(isfinite(QPALM_lbx)|isfinite(QPALM_ubx));	% nb of variables subject to a bound

  if verb >= 2
    if QPALM_nb
      nxl  = 0;
      nxu  = 0;
      nxlu = 0;
      for i=1:QPALM_n
        if QPALM_lbx(i) > -inf
          if QPALM_ubx(i) < inf
            nxlu = nxlu+1;
          else
            nxl = nxl+1;
          end
        else
          if QPALM_ubx(i) < inf
            nxu = nxu+1;
          end
        end
      end
      fprintf(fout,'\n  Bounds on the variables %6i lower, %0i lower+upper, %0i upper',nxl,nxlu,nxu);
    end
  end

  QPALM_mi = size(QPALM_AI,1);
  if (QPALM_mi ~= 0) & (size(QPALM_AI,2) ~= QPALM_n);
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent dimensions for ''AI'' (should be %0i x %0i)\n',QPALM_mi,QPALM_n); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end
  if isempty(QPALM_lbi)
    QPALM_lbi = -inf*ones(QPALM_mi,1); 
  elseif (ndims(QPALM_lbi) > 2) | (size(QPALM_lbi,1) ~= QPALM_mi) | (size(QPALM_lbi,2) ~= 1)
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent length for ''lbi'' (should be %0i)\n',QPALM_mi); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end
  if isempty(QPALM_ubi)
    QPALM_ubi = inf*ones(QPALM_mi,1); 
  elseif (ndims(QPALM_ubi) > 2) | (size(QPALM_ubi,1) ~= QPALM_mi) | (size(QPALM_ubi,2) ~= 1)
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent length for ''ubi'' (should be %0i)\n',QPALM_mi); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end
  if any ((QPALM_lbi == -inf) & (QPALM_ubi == inf))
    if verb, fprintf(fout,'\n\n### qpalm_prelim: some lower and upper bounds on AI*x are both infinite, discard the constraint\n'); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end

  if verb >= 2
    if QPALM_mi
      nil  = 0;
      niu  = 0;
      nilu = 0;
      for i=1:QPALM_mi
        if QPALM_lbi(i) > -inf
          if QPALM_ubi(i) < inf
            nilu = nilu+1;
          else
            nil = nil+1;
          end
        else
          if QPALM_ubi(i) < inf
            niu = niu+1;
          end
        end
      end
      fprintf(fout,'\n  Inequality constraints %7i lower, %0i lower+upper, %0i upper',nil,nilu,niu);
    end
  end

  QPALM_me = size(QPALM_AE,1);
  mie      = QPALM_mi+QPALM_me;
  nmie     = QPALM_n+mie;
  if QPALM_me
    if (ndims(QPALM_AE) > 2) | (size(QPALM_AE,2) ~= QPALM_n)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent dimensions for ''AE'' (should be %0i x %0i)\n',QPALM_me,QPALM_n); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
    if isempty(QPALM_be)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: right-hand-side ''be'' for equality constraints is missing\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
    if (ndims(QPALM_be) > 2) | (size(QPALM_be,1) ~= QPALM_me) | (size(QPALM_be,2) ~= 1)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent length for ''be'' (should be %0i)\n',QPALM_me); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  elseif ~isempty(QPALM_be)
    if verb, fprintf(fout,'\n\n### qpalm_prelim: right-hand-side ''be'' for equality constraints is not expected\n'); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end

  if verb >= 2
    if QPALM_me
      fprintf(fout,'\n  Nb equality constraints %6i',QPALM_me);
    end
    sparsedata  = [issparse(QPALM_g); issparse(QPALM_H); issparse(QPALM_AI); issparse(QPALM_AE); issparse(QPALM_be)];
    if any(sparsedata)
      fprintf(fout,'\n  The folowing data is sparse:');
      if sparsedata(1) == true; fprintf(fout,' g'); end
      if sparsedata(2) == true; fprintf(fout,' H'); end
      if sparsedata(3) == true; fprintf(fout,' AI'); end
      if sparsedata(4) == true; fprintf(fout,' AE'); end
      if sparsedata(5) == true; fprintf(fout,' be'); end
    else
      fprintf(fout,'\n  The data is full');
    end
  end

  if ~isempty(x0) & ( (ndims(x0) > 2) | (size(x0,1) ~= QPALM_n) | (size(x0,2) ~= 1) )
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent length for ''x0'' (should be %0i)\n',QPALM_n); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end

  if ~isempty(lm0) & ( (ndims(lm0) > 2) | (size(lm0,1) ~= nmie) | (size(lm0,2) ~= 1) )
    if verb, fprintf(fout,'\n\n### qpalm_prelim: inconsistent size for ''lm0'' (should be %0i)\n',nmie); end
    QPALM_info.flag = QPALM_values.fail_on_argument;
    return
  end

% Decode options and check consistency

  if isfield(QPALM_options,'alit')
    if isempty(QPALM_options.alit) | (QPALM_options.alit <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.alit must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.alit = 100;			% default maximal number of augmented Lagrangian iterations
  end

  if isfield(QPALM_options,'cgph')
    if isempty(QPALM_options.cgph) | (QPALM_options.cgph <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.cgph must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.cgph = inf;			% default maximal number of conjugate gradient phases
  end

  if isfield(QPALM_options,'gpph')
    if isempty(QPALM_options.gpph) | (QPALM_options.gpph <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.gpph must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.gpph = inf;			% default maximal number of gradient projection phases
  end

  if isfield(QPALM_options,'cput')
    if isempty(QPALM_options.cput) | (QPALM_options.cput <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.cput must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.cput = inf;			% default maximal number of augmented Lagrangian iterations
  end

  if isfield(QPALM_options,'rctl')
    switch lower(QPALM_options.rctl)
      case 'fixed'
        QPALM_options.rctl = QPALM_values.r_fixed;
      case 'constraint'
        QPALM_options.rctl = QPALM_values.r_constraint;
      case 'constraint-change'
        QPALM_options.rctl = QPALM_values.r_constraint_change;
      otherwise
        if verb
          fprintf(fout,'\n\n### qpalm_prelim: options.rctl = ''%s'' is not recognized',QPALM_options.rctl);
          fprintf(fout,'\n### qpalm_prelim: possible options are ''fixed'', ''constraint'', and ''constraint-change''\n');
        end
        QPALM_info.flag = QPALM_values.fail_on_argument;
        return
    end
  else
    QPALM_options.rctl = QPALM_values.r_constraint_change;
  end

  if isfield(QPALM_options,'rlag')
    if isempty(QPALM_options.rlag) & (QPALM_options.rlag <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.rlag must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.rlag = 1.e+00;		% default initial augmentation parameter
  end

  if isfield(QPALM_options,'dcr')
    if isempty(QPALM_options.dcr) & (QPALM_options.dcr <= 0) & (QPALM_options.dcr >= 1)
      fprintf('\n>>> qpalm-error: options.dcrt must be chosen in ]0,1[\n\n');
      QPALM_info.flag = 1;
      return
    end
  else
    QPALM_options.dcr = 1.e-1;			% default desired convergence rate of the constraint norm
  end

  if isfield(QPALM_options,'dcrf')
    if isempty(QPALM_options.dcrf) & (QPALM_options.dcrf < 0) & (QPALM_options.dcrf >= 1)
      fprintf('\n>>> qpalm-error: options.dcrft must be chosen in [0,1[\n\n');
      QPALM_info.flag = 1;
      return
    end
  else
    QPALM_options.dcrf = 1.e-3;			% default fraction of the desired convergence rate
  end

  if isfield(QPALM_options,'dxmin')
    if isempty(QPALM_options.dxmin) | (QPALM_options.dxmin <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.dxmin must be > 0\n'); end
      QPALM_info.flag = 1;
      return
    end
  else
    QPALM_options.dxmin = 1.e-15;		% default dxmin
  end

  if any(QPALM_ubx-QPALM_lbx < 2*QPALM_options.dxmin)
    if verb
      fprintf(fout,'\n\n### qpalm_prelim: some lower (lbx) and upper (ubx) bound too close w.r.t. options.dxmin (<%9.2e)\n', ...
        2*QPALM_options.dxmin);
    end
    QPALM_info.flag = 1;
    return
  end

  if isfield(QPALM_options,'dymin')
    if isempty(QPALM_options.dymin) | (QPALM_options.dymin <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.dymin must be > 0\n'); end
      QPALM_info.flag = 1;
      return
    end
  else
    QPALM_options.dymin = 1.e-15;		% default dymin
  end

  if any(QPALM_ubi-QPALM_lbi < 2*QPALM_options.dymin)
    if verb
      fprintf(fout,'\n\n### qpalm_prelim: some lower (lbi) and upper (ubi) inequality bound too close w.r.t. options.dymin (<%9.2e)\n',2*QPALM_options.dymin);
      I = find(QPALM_ubi-QPALM_lbi<2*QPALM_options.dymin);
      fprintf(fout,' i     lbi(i)        ubi(i)\n');
      for i = I(:)'	% be sure it is an horizontal vector
        fprintf(fout,'%2i  %12.5e  %12.5e\n',i,QPALM_lbi(i),QPALM_ubi(i));
      end
    end
    QPALM_info.flag = 1;
    return
  end

  if isfield(QPALM_options,'feas')
    if isempty(QPALM_options.feas) | (QPALM_options.feas <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.feas must be > 0\n'); end
      QPALM_info.flag = 1;
      return
    end
  else
    QPALM_options.feas = 1.e-8;			% default feasibility tolerance
  end

  if isfield(QPALM_options,'feass')
    if isempty(QPALM_options.feass) | (QPALM_options.feass <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.feass must be > 0\n'); end
      QPALM_info.flag = 1;
      return
    end
  else
    QPALM_options.feass = QPALM_options.feas;	% default tolerance for the shifted constraints
  end

  if isfield(QPALM_options,'gtol')
    if isempty(QPALM_options.gtol) | (QPALM_options.gtol <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.gtol must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.gtol = 1.e-8;			% default optimality tolerance
  end

% if isfield(QPALM_options,'mvpd')
%   if isempty(QPALM_options.mvpd) | (QPALM_options.mvpd <= 0)
%     if verb, fprintf(fout,'\n\n### qpalm_prelim: options.mvpd must be > 0\n'); end
%     QPALM_info.flag = QPALM_values.fail_on_argument;
%     return
%   end
% else
%   QPALM_options.mvpd = Inf;			% default maximal number of CG iterations
% end

  if isfield(QPALM_options,'phase')
    if isempty(QPALM_options.phase) | (QPALM_options.phase <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.phase must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.phase = Inf;			% default maximal number of phases
  end

  % CG tuning

  if isfield(QPALM_options,'cg_as') & (QPALM_nb+QPALM_mi)	% Active set strategy
    switch lower(QPALM_options.cg_as)		% it is faster to recognize a number than a string
      case 'hit-and-fix'
        QPALM_options.cgas = QPALM_values.cg_hf;
      case 'more-toraldo'
        QPALM_options.cgas = QPALM_values.cg_mt;
      otherwise
        if verb
          fprintf(fout,'\n\n### qpalm_prelim: options.cg_as = ''%s'' is not recognized',QPALM_options.cg_as);
          fprintf(fout,'\n### qpalm_prelim: possible options are ''hit-and-fix'' and ''more-toraldo''\n');
        end
        QPALM_info.flag = QPALM_values.fail_on_argument;
        return
    end
  else
    QPALM_options.cgas = QPALM_values.cg_hf;	% default is the hit-and-fix strategy
  end

  if isfield(QPALM_options,'accu')
    if isempty(QPALM_options.accu) | (QPALM_options.accu < 0) | (QPALM_options.accu > 1)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.accu must be in [0,1]\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.accu = 1.e-1;			% default accuracy measure
  end

  if isfield(QPALM_options,'cgit')
    if isempty(QPALM_options.cgit) | (QPALM_options.cgit <= 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.cgit must be > 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.cgit = Inf;			% default maximal number of conjugate gradient iterations
  end

  if isfield(QPALM_options,'cgprec')
    switch lower(QPALM_options.cgprec)
      case 'none'
        QPALM_options.cgprec = QPALM_values.cg_prec_none;
      case 'diag'
        QPALM_options.cgprec = QPALM_values.cg_prec_diag;
      case 'chol'
        QPALM_options.cgprec = QPALM_values.cg_prec_chol;
      otherwise
        if verb
          fprintf(fout,'\n\n### qpalm_prelim: options.cgprec = ''%s'' is not recognized',QPALM_options.cgprec);
          fprintf(fout,'\n### qpalm_prelim: possible options are ''none'', ''diag'', and ''chol''\n');
        end
        QPALM_info.flag = QPALM_values.fail_on_argument;
        return
    end
  else
    QPALM_options.cgprec = QPALM_values.cg_prec_none;
  end

  if isfield(QPALM_options,'ubtol')
    if isempty(QPALM_options.ubtol) | (QPALM_options.ubtol < 0)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.ubtol must be >= 0\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.ubtol = QPALM_n*eps;		% default unboundedness tolerance
  end

  if isfield(QPALM_options,'rfct')
    if isempty(QPALM_options.rfct) | (QPALM_options.rfct < 1)
      if verb, fprintf(fout,'\n\n### qpalm_prelim: options.rfct must be >= 1\n'); end
      QPALM_info.flag = QPALM_values.fail_on_argument;
      return
    end
  else
    QPALM_options.rfct = 10;			% default is to update the preconditionner when r/r- is not in [0.1,10]
  end

  % GP tuning

  if isfield(QPALM_options,'gp_init')		% set the type of initial stepsize for the piecewise linesearch
    switch lower(QPALM_options.gp_init)		% it is faster to recognize a number than a string
      case 'curv-init'
        QPALM_options.gpin = QPALM_values.curv_init;
      case 'inc-bkpt'
        QPALM_options.gpin = QPALM_values.inc_bkpt;
      otherwise
        if verb
          fprintf(fout,'\n\n### qpalm_prelim: options.gp_init = ''%s'' is not recognized',QPALM_options.gp_init);
          fprintf(fout,'\n### qpalm_prelim: possible options are ''curv-init'' and ''inc-bkpt''\n');
        end
        QPALM_info.flag = QPALM_values.fail_on_argument;
        return
    end
  else
    QPALM_options.gpin = QPALM_values.curv_init;% default is to use the initial curvature
  end

  if isfield(QPALM_options,'gp_pls')		% Type of piecewise linesearch
    switch lower(QPALM_options.gp_pls)		% it is faster to recognize a number than a string
      case 'armijo'
        QPALM_options.plsn = QPALM_values.armijo;
      case 'goldstein'
        QPALM_options.plsn = QPALM_values.goldstein;
      otherwise
        if verb
          fprintf(fout,'\n\n### qpalm_prelim: options.gp_pls = ''%s'' is not recognized',QPALM_options.gp_pls);
          fprintf(fout,'\n### qpalm_prelim: possible options are ''armijo'' and ''goldstein''\n');
        end
        QPALM_info.flag = QPALM_values.fail_on_argument;
        return
    end
  else
    QPALM_options.plsn = QPALM_values.armijo;	% default is the Armijo piecewise lineaserch
  end

  if isfield(QPALM_options,'gp_succ')		% Number of successive GP phases
    switch lower(QPALM_options.gp_succ)	% it is faster to recognize a number than a string
      case 'one'
        QPALM_options.gpsucc = QPALM_values.gp_one;
      case 'mt'
        QPALM_options.gpsucc = QPALM_values.gp_mt;
      otherwise
        if verb
          fprintf(fout,'\n\n### qpalm_prelim: options.gp_succ = ''%s'' is not recognized',QPALM_options.gp_succ);
          fprintf(fout,'\n### qpalm_prelim: possible options are ''one'' and ''mt''\n');
        end
        QPALM_info.flag = QPALM_values.fail_on_argument;
        return
    end
  else
    QPALM_options.gpsucc = QPALM_values.gp_one;	% default is to do one GP phase between CG phases
  end

% Initialize x and lm

  if isempty(lm0)
    lm_set_to_zero = true;
    QPALM_lmi = zeros(QPALM_mi,1);
    QPALM_lme = zeros(QPALM_me,1);
  else
    I = isnan(lm0);
    if all(I)
      lm_set_to_zero = true;
      QPALM_lmi = zeros(QPALM_mi,1);
      QPALM_lme = zeros(QPALM_me,1);
    else
      lm_set_to_zero = false;
      lm0(I) = 0;
      QPALM_lmi = lm0(QPALM_n+1:QPALM_n+QPALM_mi);
      QPALM_lme = lm0(QPALM_n+QPALM_mi+1:nmie);
    end
  end

  if isempty(x0)
    if lm_set_to_zero
      x_set = 0;				% means "initial x set to 0 and projected on [lbx,ubx]"
      QPALM_x = zeros(QPALM_n,1);		% give priority to zero (good for the QPs generated by SQP asymptotically)
    else
      x_set = 2;				% means "initial x deduced from the multiplier"
      Bl = find((lm0(1:QPALM_n)<0)&(QPALM_lbx>-inf));
      Bu = find((lm0(1:QPALM_n)>0)&(QPALM_ubx<inf));
      QPALM_x     = zeros(QPALM_n,1);		% give priority to zero (good for the QPs generated by SQP asymptotically)
      QPALM_x(Bl) = QPALM_lbx(Bl);
      QPALM_x(Bu) = QPALM_ubx(Bu);
%>% The technique below, which iterates by projeting alternatively on the bound constraints in x and the x's that satisfy the
%>% inequality bounds, has been abandoned since its linear system to solve are too expensive and the improvement is not
%>% significantly better. If a faster algorithm can be found to satisfy these two linear systems (try the one by Spingarn (1985)
%>% for instance), the conclusion couls be different
%>      if QPALM_mi & QPALM_options.x0it
%>        B  = union(Bl,Bu);		% indices of blocked variables
%>        F  = setdiff([1:QPALM_n],B);		% indices of free variables
%>        % initial x(F)
%>        QPALM_x    = zeros(QPALM_n,1);
%>        QPALM_x(B) = QPALM_x(B);
%>        I    = find((QPALM_lbx>-inf)&(QPALM_ubx<inf));
%>        QPALM_x(I) = 0.5*(QPALM_lbx(I)+QPALM_ubx(I));
%>        I    = find((QPALM_lbx>-inf)&(QPALM_ubx==inf));
%>        QPALM_x(I) = QPALM_lbx(I);
%>        I    = find((QPALM_lbx==-inf)&(QPALM_ubx<inf));
%>        QPALM_x(I) = QPALM_ubx(I);
%>        I    = find((QPALM_lbx==-inf)&(QPALM_ubx==inf));
%>        QPALM_x(I) = 0;
%>	xx = QPALM_x;
%>        % rhs
%>        Il = find((lm0(QPALM_n+1:QPALM_n+QPALM_mi)<0)&(QPALM_lbi>-inf));
%>        Iu = find((lm0(QPALM_n+1:QPALM_n+QPALM_mi)>0)&(QPALM_ubi<inf));
%>        I  = union(Il,Iu);
%>        b = [QPALM_lbi(Il);QPALM_ubi(Iu)]-QPALM_AI(I,B)*QPALM_x(B);
%>        % start iterate
%>        itermax = QPALM_options.x0it;
%>        iter = 0;
%>        while 1
%>          iter = iter+1;
%>          if (iter > itermax)
%>            QPALM_info.x0it = itermax;
%>            break
%>          end
%>%iter
%>          QPALM_x(F) = QPALM_x(F)+QPALM_AI(I,F)'*((QPALM_AI(I,F)*QPALM_AI(I,F)')\(b-QPALM_AI(I,F)*QPALM_x(F)));	% project x on {y: AI(I,F)*y(F) = b}
%>%QPALM_AI(I,:)*QPALM_x
%>%QPALM_x(F)-max(QPALM_lbx(F),min(QPALM_ubx(F),QPALM_x(F)))
%>          QPALM_x(F) = max(QPALM_lbx(F),min(QPALM_ubx(F),QPALM_x(F)));					% project x on [lbx,ubx]
%>%keyboard
%>	  if norm(QPALM_x-xx,inf) < 1.e-10*max(norm(QPALM_x,inf),norm(xx,inf))
%>            QPALM_info.x0it = iter;
%>            break
%>          end
%>          xx = QPALM_x;
%>        end
%>%       Il = find((lm0(QPALM_n+1:QPALM_n+QPALM_mi)<0)&(QPALM_lbi>-inf));
%>%       Iu = find((lm0(QPALM_n+1:QPALM_n+QPALM_mi)>0)&(QPALM_ubi<inf));
%>%       I  = union(Il,Iu);
%>%       QPALM_x(F) = QPALM_AI(I,F)\([QPALM_lbi(Il);QPALM_ubi(Iu)]-QPALM_AI(I,B)*QPALM_x(B));	% least squares solution
%>      end
%keyboard
    end
  else
    x_set = 1;						% means "initial x taken on entry and projected on [lbx,ubx]"
    QPALM_x = x0;
  end
  QPALM_x = max(QPALM_lbx,min(QPALM_ubx,QPALM_x));	% project x on [lbx,ubx]

% Print description of the problem

  if verb >= 2
    fprintf(fout,'\nRequired accuracy');
    fprintf(fout,'\n  On feasibility               %7.2e',QPALM_options.feas);
    if QPALM_options.rctl ~= QPALM_values.r_constraint; fprintf(fout,'\n  On feasibility of the CFQP   %7.2e',QPALM_options.feass); end
    fprintf(fout,'\n  On optimality                %7.2e\n',QPALM_options.gtol);

    fprintf(fout,'Tunings\n');
    if verb == 2
      fprintf(fout,'  Verbosity level           %4i\n',verb);
    else
      fprintf(fout,'  Verbosity level                         %4i\n',verb);
      fprintf(fout,'  Resolution in x                            %7.2e\n',QPALM_options.dxmin);
      if QPALM_mi; fprintf(fout,'  Resolution in AI*x                         %7.2e\n',QPALM_options.dymin); end
      fprintf(fout,'  Unboundedness tolerance                    %7.2e\n',QPALM_options.ubtol);
    end
    if QPALM_nb+QPALM_mi
      fprintf(fout,'  Gradient projection (GP) phases\n');
      fprintf(fout,'    Maximal number of GP phases %14i\n    ',QPALM_options.gpph);
      switch QPALM_options.plsn
        case QPALM_values.armijo;
          fprintf(fout,'Armijo''s');
        case QPALM_values.goldstein;
          fprintf(fout,'Goldstein''s');
      end
      fprintf(fout,' linesearch\n    Initial stepsize ');
      switch QPALM_options.gpin
        case QPALM_values.curv_init;
          fprintf(fout,'computed from the initial curvature\n');
        case QPALM_values.inc_bkpt;
          fprintf(fout,'is the first breakpoint at which phi increases\n');
      end
      if QPALM_options.gpsucc == QPALM_values.gp_mt
        fprintf(fout,'    Possibly more than one GP phase between CG phases\n');
      else
        fprintf(fout,'    One single GP phase between CG phases\n');
      end
    end
    fprintf(fout,'  Conjugate gradient (CG) phases\n');
    if verb >= 3
      fprintf(fout,'    Maximal number of CG phases %14i\n',QPALM_options.cgph);
      fprintf(fout,'    Maximal number of CG iterations %10i\n',QPALM_options.cgit);
      switch QPALM_options.cgprec
        case QPALM_values.cg_prec_none;
          fprintf(fout,'    No preconditioning\n');
        case QPALM_values.cg_prec_diag;
          fprintf(fout,'    Diagonal preconditioning\n');
        case QPALM_values.cg_prec_chol;
          fprintf(fout,'    Cholesky preconditioning\n');
      end
      if QPALM_nb+QPALM_mi
        fprintf(fout,'    Active-set strategy: ');
        switch QPALM_options.cgas
          case QPALM_values.cg_hf;
            fprintf(fout,'''hit-and-fix''\n');
          case QPALM_values.cg_mt;
            fprintf(fout,'''more-toraldo''\n');
        end
      end
    end

    if mie & (verb >= 3)
      fprintf(fout,'  Linear constraints are dealt with augmented Lagrangian iterations\n');
      fprintf(fout,'    Augmentation parameters ');
      if QPALM_options.rctl == QPALM_values.r_fixed
        fprintf(fout,'fixed\n');
      elseif QPALM_options.rctl == QPALM_values.r_constraint
        fprintf(fout,'controlled by the constraint values\n');
      elseif QPALM_options.rctl == QPALM_values.r_constraint_change
        fprintf(fout,'controlled by the constraint value changes\n');
      end
      fprintf(fout,'    Maximal number of AL iterations       %4i\n',QPALM_options.alit);
      fprintf(fout,'    Initial augmentation parameter           %7.2e\n',QPALM_options.rlag);
      fprintf(fout,'    Desired convergence rate                 %7.2e\n',QPALM_options.dcr);
      fprintf(fout,'    Desired convergence rate fraction        %7.2e\n',QPALM_options.dcrf);
    end
    if verb >= 3
      fprintf(fout,'  Maximal\n');
      fprintf(fout,'    number of phases %25i\n',QPALM_options.phase);
      fprintf(fout,'    CPU time (in sec)   %22g\n',QPALM_options.cput);
    end

    qpalm_optim ();
    fprintf(fout,'Initial accuracy');
    if verb == 2
      fprintf(fout,'\n  On optimality (inf norm)     ');
      if QPALM_info.pgln < inf
        fprintf(fout,'%7.2e',QPALM_info.pgln);
      else
        fprintf(fout,'unavailable');
      end
      if QPALM_nb+mie
        fprintf(fout,'\n  On feasibility (inf norm)    %7.2e',QPALM_info.feas);
      end
    else
      fprintf(fout,'\n  On optimality (inf norm)                   ');
      if QPALM_info.pgln < inf
        fprintf(fout,'%7.2e',QPALM_info.pgln);
      else
        fprintf(fout,'unavailable');
      end
      if QPALM_nb+mie
        fprintf(fout,'\n  On feasibility (inf norm)                  %7.2e',QPALM_info.feas);
      end
    end
    if verb >= 3
      if x_set == 0
        fprintf(fout,'\nInitial x set to 0 and projected on [lbx,ubx]');
      elseif x_set == 1
        fprintf(fout,'\nInitial x taken on entry and projected on [lbx,ubx]');
      else
        fprintf(fout,'\nInitial x deduced from the multiplier');
      end
      if lm_set_to_zero
        fprintf(fout,'\nInitial mulitplier set to 0');
      else
        fprintf(fout,'\nInitial mulitplier taken on entry');
      end
    end
  end

% Build the MEX-file for dpstrf if necessary

  if (QPALM_options.cgprec == QPALM_values.cg_prec_chol) & (exist('qpalm_dpstrf') ~= 3)	% the MEX-file for dpstrf has not been created

    % Introductory message

    if verb >= 3; fprintf(fout,'\n%s\nqpalm_prelim: no MEX-file for ''dpstrf''\n',QPALM_values.dline); end

    % Got to the Qpalm directory

    working_dir = pwd;
    qpalm_dir = which('qpalm');
    qpalm_dir = qpalm_dir(1:strfind(qpalm_dir,'/qpalm.m')-1);
    cd(qpalm_dir);

    % Try creating the MEX-file from qpalm_dpstrf.F and the Lapack library

    try
      if verb >= 3; fprintf(fout,'qpalm_prelim: Trying to create a MEX-file for ''dpstrf'' from ''qpalm_dpstrf.F''\n'); end
      mex qpalm_dpstrf.F -llapack;	% create the MEX-file for dpstrf
    catch

    % Try creating the MEX-file from qpalm_dpstrf_inline.F, which contains the source code of 'dpstrf'

      try
        if verb >= 3; fprintf(fout,'qpalm_prelim: Trying to create a MEX-file for ''dpstrf'' from ''qpalm_dpstrf_inline.F''\n'); end
        mex qpalm_dpstrf_inline.F;	% create the MEX-file for dpstrf
        newfile1 = ls('qpalm_dpstrf_inline.mex*');	% get the generated file
        newfile2 = strcat('qpalm_dpstrf.',newfile1(strfind(newfile1,'mex'):end));	% modifie the root name of the file
        movefile(newfile1,newfile2);	% rename the generated file
      catch
        if verb
          fprintf(fout,'\n### qpalm_prelim: failed to create the MEX-file for ''qpalm_dpstrf'';\n')
          fprintf(fout,'                  Cholesky precondtioning cannot be used; change options.cgprec\n')
        end
        QPALM_info.flag = QPALM_values.fail_on_mex;	% failure on MEX-file generation
      end
    end
    cd(working_dir);			% get back to the working directory
  end

% Exit

  return
