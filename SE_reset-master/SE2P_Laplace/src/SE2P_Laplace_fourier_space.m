function [u,walltime] = SE2P_Laplace_fourier_space(eval_idx,x,f,opt)
% SE2P_LAPLACE_FOURIER_SPACE
% Compute Fourier-space part of Ewald sum for 2-periodic electrostatic potential,
% using the Spectral Ewald method.
%
% Periodicity is assumed in the first two directions (x and y),
% while the third direction (z) is assumed free.
%
% u = SE2P_Laplace_fourier_space(eval_idx, x, f, opt)
%   Return potential
%
% [u, walltime] = SE2P_Laplace_fourier_space(...)
%   Return potential and walltime
%
% Parameters:
% :param eval_idx: index of source locations where potential should be evaluated
% :param x: source locations (N×3)
% :param f: source charges (N×1)
% :param opt: Ewald options:
% :param opt.box:           Size of periodic cell [L1,L2,L3] (required)
% :param opt.M:             Grid size [M1,M2,M3] (required)
% :param opt.add_M3:        Extra number to add to M3 (default: 0), see below
% :param opt.xi:            Ewald parameter (required)
% :param opt.P:             Support of window, in number of grid points (required)
% :param opt.s:             Oversampling factor for local pad (required)
% :param opt.s0:            Oversampling factor for zero mode (required)
% :param opt.n:             Used to define which modes to include in the local pad
%                           (default: n=max(ceil(opt.M(1)/2),1))
% :param opt.window:        Window function (default: 'gaussian')
% :param opt.w:             Width of window (default: w=h*P/2)
% :param opt.m:             Gaussian shape function (default: m=0.94*sqrt(pi*P))
% :param opt.eval_x:        External points to evaluate potential in (Nex×3)
% :param opt.fast_gridding: Use fast gridding and spreading (default: true)
%
% The grid size in the free direction (M3) will be set to
%   M3 = opt.M(3) + opt.P + opt.add_M3,
% i.e. both opt.P and opt.add_M3 are added to the value provided
% in the opt struct.
%
% Valid window functions are 'gaussian', 'expsemicirc', 'kaiser_exact'
% and 'kaiser_poly'.
%
% Parameters specific to the expsemicirc, kaiser_exact and
% kaiser_poly windows:
% :param opt.betaP:   Shape parameter beta/P (default: betaP=2.5)
%
% Parameters specific to the kaiser_poly window:
% :param opt.polynomial_degree: Polynomial degree (default: 1)
%
% NOTE: This function may give unexpected results if there are
%       sources outside of opt.box in the free direction. It is
%       up to the caller to ensure that the sources are inside
%       the box.
%
% :returns: **u** -- Fourier-space potential
% :returns: **walltime** -- struct containing timings

opt = parse_params(opt);

% Initialize timing array
walltime = struct('grid',0,'fft',0,'scale',0,'ifft',0,'int',0,'pre',0,'prefft',0);

% Window function specifics
pre = [];
W_precomp_gridding = @(x,opt) SE2P_Laplace_gridding_precomp(x, opt);
if strcmp(opt.window, 'gaussian')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_thrd_mex_2p(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  W_fast_int = @(F,opt,S) SE_fg_int_split_mex_2p(0,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_mex_2p is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_plain_grid = @SE_fg_grid_mex_2p;
  W_plain_int = @SE_fg_int_mex_2p;
elseif strcmp(opt.window, 'expsemicirc') || strcmp(opt.window, 'kaiser_exact') ...
    || strcmp(opt.window, 'kaiser_poly')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_kaiser_mex_2p(x,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  W_fast_int = @(F,opt,S) SE_fg_int_split_kaiser_mex_2p(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_kaiser_mex_2p is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_plain_grid = @SE_fg_grid_kaiser_mex_2p;
  W_plain_int = @SE_fg_int_kaiser_mex_2p;

  if strcmp(opt.window, 'expsemicirc')
    % Precompute window Fourier transform
    prefft_t = tic;
    pre = SE2P_window_fft_precomp(opt);
    walltime.prefft = walltime.prefft + toc(prefft_t);
  end
else
  error('Unsupported window function');
end

% Functions for gridding and spreading
if opt.fast_gridding
  % Gridder
  pre_t = tic;
  S = W_precomp_gridding(x, opt);
  walltime.pre = walltime.pre + toc(pre_t);
  grid_fcn = @(F) W_fast_grid(x(S.perm,:), F(S.perm), opt, S);
  % Integrator
  iperm = @(u) u(S.iperm,:);
  int_fcn = @(F) iperm(W_fast_int(F, opt, S));
else
  error('Fast gridding should be used!');
  grid_fcn = @(F) W_plain_grid(x, F, opt);
  int_fcn = @(F) W_plain_int(x(eval_idx,:), F, opt);
end

% Gridding
grid_t = tic;
H = grid_fcn(f);
walltime.grid = walltime.grid + toc(grid_t);

% Fourier transform
fft_t = tic;
[G, Gr, G0] = SE2P_Laplace_fftnd(H, opt);
walltime.fft = walltime.fft + toc(fft_t);

% Scaling
scale_t = tic;
[G, Gr, G0] = SE2P_Laplace_scaling(G, Gr, G0, opt, pre);
walltime.scale = walltime.scale + toc(scale_t);

% Inverse Fourier transform
ifft_t = tic;
F = SE2P_Laplace_ifftnd(G, Gr, G0, opt);
walltime.ifft = walltime.ifft + toc(ifft_t);

% Spreading and integration
int_t = tic;
u = 4*pi*int_fcn(F);
walltime.int = walltime.int + toc(int_t);

if opt.fast_gridding
  u = u(eval_idx,:); % pick out the evaluation points
end

if nargout == 2
  walltime.total = sum(struct2array(walltime));
end

% ------------------------------------------------------------------------------
function opt = parse_params(opt)
% Check that mandatory options are present
assert(isfield(opt, 'box'), 'cell size box must be given in opt struct');
assert(isfield(opt, 'M'), 'grid size M must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald parameter xi must be given in opt struct');
assert(isfield(opt, 'P'), 'window support P must be given in opt struct');
assert(isfield(opt, 's'), 'oversampling factor s must be given in opt struct');
assert(isfield(opt, 's0'), 'zero-mode oversampling factor s0 must be given in opt struct');

% Verify assumptions on parameters
opt.L = opt.box(1);
opt.h = opt.L/opt.M(1); % step size (M contains number of subintervals)
if ~isfield(opt,'add_M3'), opt.add_M3 = 0; end
opt.Mz = opt.M(3) + opt.P + opt.add_M3; % increase grid size in z direction
% If Mx is even, make Mz even too
if mod(opt.M(1), 2) == 0
  opt.Mz = 2*ceil(opt.Mz/2);
end
opt.M(3) = opt.Mz; % for good measure
opt.Lz = opt.h*opt.Mz;
opt.box(3) = opt.Lz; % for good measure
% Check that h is the same in all directions
assert(abs(opt.h - opt.box(2)/opt.M(2)) < eps);
assert(abs(opt.h - opt.Lz/opt.Mz) < eps);

if ~isfield(opt,'window'), opt.window = 'gaussian'; end
if ~isfield(opt,'w'), opt.w = opt.h*opt.P/2; end
if ~isfield(opt,'m'), opt.m = 0.94*sqrt(pi*opt.P); end
opt.eta = (2*opt.w*opt.xi/opt.m)^2;
opt.c = 2*opt.xi^2/opt.eta;

wbox = [0 opt.L; 0 opt.L; -opt.w -opt.w+opt.Lz];
opt.a = wbox(3,1);

% Sampling factor (oversampling)
if ~isfield(opt,'s'), opt.s = 4; end % Are these mandatory or not?
if ~isfield(opt,'s0'), opt.s0 = 2; end % -"-
if ~isfield(opt,'n'), opt.n = max(ceil(opt.M(1)/2),1); end
% Increase s and s0 such that FFTN has integer size vectors.
opt.s = ceil(opt.s*opt.Mz)/opt.Mz;
opt.s0 = ceil(opt.s0*opt.Mz)/opt.Mz;

% Local pads
% FIXME: We assume that the same number of modes in x and y directions
% are oversampled.
n = opt.n;
if n > 1
    n = min(floor((opt.M(1)-1)/2), n); % half-modes should be at most half of M
end
% Zero mode is the first element. But for simplicity we add it and
% overwrite it whenever needed.
opt.local_pad = [1:n+1 opt.M(1)-n+1:opt.M(1)];

if ~isfield(opt,'fast_gridding'), opt.fast_gridding = true; end

% Options for ExpSemiCirc and Kaiser-Bessel windows
if strcmp(opt.window,'expsemicirc') || strcmp(opt.window,'kaiser_exact') ...
    || strcmp(opt.window,'kaiser_poly')
  if ~isfield(opt,'betaP'), opt.betaP = 2.5; end
  opt.beta = opt.betaP*opt.P;
end
if strcmp(opt.window,'kaiser_exact') || strcmp(opt.window,'kaiser_poly')
  opt.kaiser_scaling = 1/besseli(0,opt.beta);
end
if strcmp(opt.window,'kaiser_poly')
  if ~isfield(opt,'polynomial_degree'), opt.polynomial_degree = 1; end
end

% Half-support of window
if mod(opt.P,2)==0
  opt.p_half = opt.P/2;
else
  opt.p_half = (opt.P-1)/2;
end

% External evaluation points
if isfield(opt,'eval_x') && numel(opt.eval_x) > 0
  opt.eval_external = true;
  fprintf('WARNING: External evaluation points are not implemented\n');
else
  opt.eval_external = false;
  opt.eval_x = [];
end
% TODO: external evaluation points are not implemented yet!

opt.R = opt.Lz;
