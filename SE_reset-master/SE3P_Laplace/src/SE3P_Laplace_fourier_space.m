function [pot, force, walltime] = SE3P_Laplace_fourier_space(eval_idx, x, f, opt)
% SE3P_LAPLACE_FOURIER_SPACE
% Compute Fourier-space part of Ewald sum for 3-periodic electrostatic potential,
% using the Spectral Ewald method.
%
% pot = SE3P_Laplace_fourier_space(eval_idx, x, f, opt)
%   Return potential
%
% [pot, force] = SE3P_Laplace_fourier_space(...)
%   Return potential and force
%
% [pot, force, walltime] = SE3P_Laplace_fourier_space(...)
%   Return potential, force and walltime
%
% Unless opt.force is true, the force output will be empty.
%
% Parameters:
% :param eval_idx: index of source locations where potential should be evaluated
% :param x: source locations (N×3)
% :param f: source charges (N×1)
% :param opt: Ewald options:
% :param opt.potential:     Compute potential (default: true)
% :param opt.force:         Compute force (default: false)
% :param opt.box:           Size of periodic cell [L1,L2,L3] (required)
% :param opt.M:             Grid size [M1,M2,M3] (required)
% :param opt.base_factor:   Integer that the final grid size should be divisible by
%                           (default: 4); the grid size will be rounded up
% :param opt.xi:            Ewald parameter (required)
% :param opt.P:             Support of window, in number of grid points (required)
% :param opt.window:        Window function (default: 'gaussian')
% :param opt.w:             Width of window (default: w=h*P/2)
% :param opt.m:             Gaussian shape parameter (default: m=0.95*sqrt(pi*P))
% :param opt.eval_x:        External points to evaluate potential in (Nex×3)
% :param opt.fast_gridding: Use fast gridding and spreading (default: true)
% :param opt.fourier_differentiation: Compute forces using differentiation in
%                                     Fourier space and three IFFTs; has no
%                                     effect unless opt.force is true
%                                     (default: false)
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
% :returns: **pot** -- Fourier-space potential
% :returns: **force** -- Fourier-space force
% :returns: **walltime** -- struct containing timings

opt = parse_params(opt);

% Initialize timing array
walltime = struct('grid',0,'fft',0,'scale',0,'ifft',0,'int',0,'pre',0,'prefft',0);

% Window function specifics
pre = [];
W_precomp_gridding = @(x,opt) SE3P_Laplace_gridding_precomp(x, opt);
if strcmp(opt.window, 'gaussian')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_thrd_mex(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  W_fast_int_pot = @(F,opt,S) SE_fg_int_split_mex(0,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_mex is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_fast_int_force = @(F,opt,S) SE_fg_int_split_force_mex(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,...
                                S.zfx,S.zfy,S.zfz,S.idx);
  W_plain_grid = @SE_fg_grid_thrd_mex;
  W_plain_int_pot = @SE_fg_int_mex;
  W_plain_int_force = @SE_fg_int_force_mex;
elseif strcmp(opt.window, 'expsemicirc') || strcmp(opt.window, 'kaiser_exact') ...
    || strcmp(opt.window, 'kaiser_poly')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_kaiser_mex(x,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  W_fast_int_pot = @(F,opt,S) SE_fg_int_split_kaiser_mex(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_kaiser_mex is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_fast_int_force = @(F,opt,S) SE_fg_int_split_kaiser_force_mex(x,F,opt,...
                                S.zx,S.zy,S.zz,...
                                S.zfx,S.zfy,S.zfz,S.idx);
  W_plain_grid = @SE_fg_grid_kaiser_mex;
  W_plain_int_pot = @SE_fg_int_kaiser_mex;
  W_plain_int_force = @SE_fg_int_kaiser_force_mex;

  if strcmp(opt.window, 'expsemicirc')
    if opt.force
      error('Force calculations are not supported for the expsemicirc window');
    end
    % Precompute window Fourier transform
    prefft_t = tic;
    pre = SE3P_window_fft_precomp(opt);
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
  int_fcn_pot = @(F) iperm(W_fast_int_pot(F, opt, S));
  int_fcn_force = @(F) iperm(W_fast_int_force(F, opt, S));
else
  warning('SE3P:SlowGridding', 'Using slow gridding routines');
  grid_fcn = @(F) W_plain_grid(x, F, opt);
  int_fcn_pot = @(F) W_plain_int_pot(x(eval_idx,:), F, opt);
  int_fcn_force = @(F) W_plain_int_force(x(eval_idx,:), F, opt);
end

% Gridding
grid_t = tic;
H = grid_fcn(f);
walltime.grid = walltime.grid + toc(grid_t);

% Fourier transform
fft_t = tic;
H = fftn(H);
walltime.fft = walltime.fft + toc(fft_t);

% Scaling
scale_t = tic;
if opt.force && opt.fourier_differentiation
  [H, K1, K2, K3] = SE3P_Laplace_scaling(H, opt, pre);
  Hi = 1i * H;
  F1 = Hi .* K1;
  F2 = Hi .* K2;
  F3 = Hi .* K3;
else
  H = SE3P_Laplace_scaling(H, opt, pre);
end
walltime.scale = walltime.scale + toc(scale_t);

% Inverse Fourier transform
ifft_t = tic;
if opt.potential || (opt.force && ~opt.fourier_differentiation)
  H = ifftn(H); % this should be real to eps accuracy!
end
if opt.force && opt.fourier_differentiation
  F1 = ifftn(F1);
  F2 = ifftn(F2);
  F3 = ifftn(F3);
end
walltime.ifft = walltime.ifft + toc(ifft_t);

% Spreading and integration
pot = []; force = [];
if opt.potential
  int_t = tic;
  pot = 4*pi*int_fcn_pot(H);
  walltime.int = walltime.int + toc(int_t);
end
if opt.force
  if opt.fourier_differentiation
    int_t = tic;
    F1 = 4*pi*int_fcn_pot(F1);
    F2 = 4*pi*int_fcn_pot(F2);
    F3 = 4*pi*int_fcn_pot(F3);
    force = [F1, F2, F3];
    walltime.int = walltime.int + toc(int_t);
  else
    int_t = tic;
    force = 4*pi*int_fcn_force(H);
    walltime.int = walltime.int + toc(int_t);
  end
end

if opt.fast_gridding
  % Pick out the evaluation points
  if ~isempty(pot)
    pot = pot(eval_idx,:);
  end
  if ~isempty(force)
    force = force(eval_idx,:);
  end
end

% Strange corrections (FIXME: to be removed/resolved)
if ~opt.fourier_differentiation && ~isempty(force)
  if strcmp(opt.window, 'gaussian')
    force = force * 2; % FIXME TODO: some bug causes a factor 2 error!
  elseif strcmp(opt.window, 'kaiser_exact')
    force = force / opt.h; % FIXME TODO: some bug causes a factor 1/h error!
                           % Caused by a missing 1/h in the derivative
                           % computation (scaling error).
  elseif strcmp(opt.window, 'kaiser_poly')
    force = force / (opt.h*opt.P/2); % FIXME TODO: some bug causes a factor 1/(hP/2) error!
                                     % (Probably another scaling error)
  end
end


if nargout >= 3
  walltime.total = sum(struct2array(walltime));
end

% ------------------------------------------------------------------------------
function opt = parse_params(opt)
% Check that mandatory options are present
assert(isfield(opt, 'box'), 'cell size box must be given in opt struct');
assert(isfield(opt, 'M'), 'grid size M must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald parameter xi must be given in opt struct');
assert(isfield(opt, 'P'), 'window support P must be given in opt struct');

if ~isfield(opt,'potential'), opt.potential = true; end
if ~isfield(opt,'force'), opt.force = false; end
if ~isfield(opt,'fourier_differentiation'), opt.fourier_differentiation = false; end
if ~isfield(opt,'fast_gridding'), opt.fast_gridding = true; end
if ~isfield(opt,'base_factor'), opt.base_factor = 4; end

% Basic grid
opt.M(1) = opt.base_factor * ceil(opt.M(1) / opt.base_factor); % round up
opt.M(2) = opt.base_factor * ceil(opt.M(2) / opt.base_factor); % round up
opt.M(3) = opt.base_factor * ceil(opt.M(3) / opt.base_factor); % round up
opt.L = opt.box(1);
opt.h = opt.L/opt.M(1); % step size (M contains number of subintervals)
% Check that h is the same in all directions
assert(abs(opt.h - opt.box(2)/opt.M(2)) < eps);
assert(abs(opt.h - opt.box(3)/opt.M(3)) < eps);

% Window options
if ~isfield(opt,'window'), opt.window = 'gaussian'; end
if ~isfield(opt,'w'), opt.w = opt.h*opt.P/2; end
if ~isfield(opt,'m'), opt.m = 0.95*sqrt(pi*opt.P); end
opt.eta = (2*opt.w*opt.xi/opt.m)^2;
opt.c = 2*opt.xi^2/opt.eta;

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
