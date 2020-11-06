function popt = se2p_parse_params(opt)
 
% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% copy all params
popt = se3p_parse_params(opt);

% step size
popt.p_half = (mod(opt.P,2)==0)*opt.P/2+(mod(opt.P,2)~=0)*(opt.P-1)/2;

% Old delta (should cover grid Gaussian)
delta_old = popt.h*opt.P/2;
if isfield(opt, 'window') && strcmp(opt.window, 'kaiser')
  %delta = max(delta_old, sqrt(popt.beta*2)/popt.xi);
  delta = max(delta_old, sqrt(popt.beta/3)/popt.xi);
else
  % New delta (should cover remainder Gaussian)
  if popt.eta < 1
      deltaNew = sqrt((1-popt.eta)/2)*popt.m/popt.xi;
  else
      deltaNew = delta_old;
  end
  % Take max
  delta = max(deltaNew, delta_old);
end
EVEN_GRIDS = true;
if EVEN_GRIDS
    deltaM = 2*ceil(delta / popt.h);
else
    deltaM = ceil(2*delta / popt.h);
end
delta = popt.h*deltaM/2;

popt.Lz = opt.box(3) + 2*delta;
popt.Mz = opt.M(3) + deltaM;
%popt.Mz = ceil(popt.box(3)/popt.h)+popt.P;
%popt.Lz = popt.box(3)+2*popt.w;

% Even grids if even grids are in the z direction
%if(mod(popt.Mz,2)~=0)
%    popt.Mz = 2*ceil(popt.Mz/2);
%end
%popt.Lz = popt.h*popt.Mz;

% h should be the same in all directions
assert(abs(popt.Lz/popt.Mz-popt.h)<eps)

% sampling factor (oversampling)
if( isfield(opt,'s')), popt.s = opt.s; else popt.s=2; end;
if( isfield(opt,'s0')), popt.s0 = opt.s0; else popt.s0=2; end;
if( isfield(opt,'n')),
    popt.n = min(opt.n,ceil(popt.M/2));
else 
    popt.n=max(ceil(popt.M/2),1); 
end;

popt.R = popt.Lz;

% increase s and such that FFTN has integer size vectors.
popt.s  = ceil(popt.s*popt.Mz)/popt.Mz;
popt.s0 = ceil(popt.s0*popt.Mz)/popt.Mz;

% local pads
% FIXME: We assume that the same number of modes in x and y directions
% are oversampled.
n = popt.n;
if(n>1)
    n = min(floor((popt.M-1)/2),n); % half modes should be at most
                                    % half of M
end
% zero mode is the first element. But for simplicity we add it and
% overwrite it whenever needed.    
popt.local_pad = [1:n+1 popt.M-n+1:popt.M];

% collect
popt.PH = popt.P/2;
wbox = [0 popt.L; 0 popt.L; -popt.w  -popt.w+popt.Lz];
popt.a = wbox(3,1);
