function [pot, force, walltime] = SE0P_Laplace_fourier_space(eval_idx, x, f, opt, pre_kernel)
% SE0P_LAPLACE_FOURIER_SPACE
% Compute Fourier-space part of Ewald sum for 0-periodic (free-space)
% electrostatic potential, using the Spectral Ewald method.
%
% pot = SE0P_Laplace_fourier_space(eval_idx, x, f, opt)
%   Return potential
%
% [pot, force] = SE0P_Laplace_fourier_space(...)
%   Return potential and force
%
% [pot, force, walltime] = SE0P_Laplace_fourier_space(...)
%   Return potential, force and walltime
%
% [...] = SE0P_Laplace_fourier_space(..., pre_kernel)
%   Provide precomputed kernel (from SE0P_Laplace_kernel_fft_precomp)
%
% Unless opt.force is true, the force output will be empty.
%
% Parameters:
% :param eval_idx: index of source locations where potential should be evaluated
% :param x: source locations (N×3)
% :param f: source charges (N×1)
% :param opt: Ewald options:
% :param opt.potential:         Compute potential (default: true)
% :param opt.force:             Compute force (default: false)
% :param opt.box:               Size of periodic cell [L1,L2,L3] (required)
% :param opt.M:                 Grid size [M1,M2,M3] (required)
% :param opt.xi:                Ewald parameter (required)
% :param opt.P:                 Support of window, in number of grid points (required)
% :param opt.s:                 Oversampling factor (required)
% :param opt.oversample_all:    Use the oversampled grid for the padded grid as well
%                               (default: false)
% :param opt.no_extra_support:  Don't take the "remainder window" into account when
%                               computing the extended grid (default: false)
% :param opt.window:            Window function (default: 'gaussian')
% :param opt.w:                 Width of window (default: w=h*P/2)
% :param opt.m:                 Gaussian shape function (default: m=0.94*sqrt(pi*P))
% :param opt.eval_x:            External points to evaluate potential in (Nex×3)
% :param opt.fast_gridding:     Use fast gridding and spreading (default: true)
%
% Four different grids are used in free space:
% - inner: domain containing all sources and targets (given by opt.box and opt.M)
% - extended: FGG grid (to cancel wrap effects)
% - padded: 2x FFT grid (for aperiodic convolution)
% - oversampled: FFT grid for truncated Green's function, using
%   oversampling factor s
% If opt.oversample_all is true, the oversampled grid is used
% also for the aperiodic convolution.
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
%       sources outside of opt.box. It is up to the caller to ensure
%       that the sources (and targets) are inside the box.
%
% :returns: **pot** -- Fourier-space potential
% :returns: **force** -- Fourier-space force
% :returns: **walltime** -- struct containing timings

opt = SE0P_parse_params(opt);

% Initialize timing array
walltime = struct('grid',0,'fft',0,'scale',0,'ifft',0,'int',0,'pre',0,'prefft',0);

% Window function specifics
pre_window = [];
W_precomp_gridding = @(x,opt) SE0P_Laplace_gridding_precomp(x, opt);
if strcmp(opt.window, 'gaussian')
  if opt.force
    warning('SE0P:Force', 'Force calculations are experimental and most likely incorrect');
  end
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
  W_plain_grid = @SE_fg_grid_mex;
  W_plain_int_pot = @SE_fg_int_mex;
  W_plain_int_force = @SE_fg_int_force_mex;
elseif strcmp(opt.window, 'expsemicirc') || strcmp(opt.window, 'kaiser_exact') ...
    || strcmp(opt.window, 'kaiser_poly')
  if opt.force
    error('Force calculations are only supported for the Gaussian window');
  end
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_kaiser_mex(x,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  W_fast_int_pot = @(F,opt,S) SE_fg_int_split_kaiser_mex(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_kaiser_mex is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_fast_int_force = @(F,opt,S) SE_fg_int_split_kaiser_mex(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx); % FIXME
  W_plain_grid = @SE_fg_grid_kaiser_mex;
  W_plain_int_pot = @SE_fg_int_kaiser_mex;
  W_plain_int_force = @SE_fg_int_kaiser_mex; % FIXME

  if strcmp(opt.window, 'expsemicirc')
    % Precompute window Fourier transform
    prefft_t = tic;
    pre_window = SE0P_window_fft_precomp(opt);
    walltime.prefft = walltime.prefft + toc(prefft_t);
  end
else
  error('Unsupported window function');
end

% Precompute Green's function
% FIXME: seems to be quite time-consuming!
if nargin < 5 || isempty(pre_kernel)
  prefft_t = tic;
  pre_kernel = SE0P_Laplace_kernel_fft_precomp(opt);
  walltime.prefft = walltime.prefft + toc(prefft_t);
end

% Functions for gridding and spreading
% Use extended box for gridding and spreading
mod_opt = struct();
mod_opt.potential = opt.potential;
mod_opt.force = opt.force;
mod_opt.box = opt.extended_box;
mod_opt.M = opt.extended_M;
mod_opt.h = mod_opt.box(1)/mod_opt.M(1);
mod_opt.window = opt.window;
mod_opt.P = opt.P;
mod_opt.w = mod_opt.h*opt.P/2;
mod_opt.eta = (2*mod_opt.w*opt.xi/opt.m)^2;
mod_opt.c = 2*opt.xi^2/mod_opt.eta;
if isfield(opt, 'beta')
  mod_opt.beta = opt.beta;
end
if isfield(opt, 'polynomial_degree')
  mod_opt.polynomial_degree = opt.polynomial_degree;
end
% Shift x for the extended box
x = bsxfun(@plus, x, opt.deltaB);
if opt.fast_gridding
  % Gridder
  pre_t = tic;
  S = W_precomp_gridding(x, mod_opt);
  walltime.pre = walltime.pre + toc(pre_t);
  grid_fcn = @(F) W_fast_grid(x(S.perm,:), F(S.perm), mod_opt, S);
  % Integrator
  iperm = @(u) u(S.iperm,:);
  int_fcn_pot = @(F) iperm(W_fast_int_pot(F, mod_opt, S));
  int_fcn_force = @(F) iperm(W_fast_int_force(F, mod_opt, S));
else
  error('Fast gridding should be used!');
  grid_fcn = @(F) W_plain_grid(x, F, mod_opt);
  int_fcn_pot = @(F) W_plain_int_pot(x(eval_idx,:), F, mod_opt);
  int_fcn_force = @(F) W_plain_int_force(x(eval_idx,:), F, mod_opt);
end

% Gridding
grid_t = tic;
H = grid_fcn(f);
walltime.grid = walltime.grid + toc(grid_t);

% Fourier transform
fft_t = tic;
H = fftshift(fftn(H, opt.padded_M));
walltime.fft = walltime.fft + toc(fft_t);

% Scaling
scale_t = tic;
H = SE0P_Laplace_scaling(H, opt, pre_kernel, pre_window);
walltime.scale = walltime.scale + toc(scale_t);

% Inverse Fourier transform
ifft_t = tic;
H = ifftn(ifftshift(H));
walltime.ifft = walltime.ifft + toc(ifft_t);

% Spreading and integration
% Truncate grid before integration
M = opt.extended_M;
pot = []; force = [];
if opt.potential
  int_t = tic;
  pot = int_fcn_pot(H(1:M(1), 1:M(2), 1:M(3)));
  walltime.int = walltime.int + toc(int_t);
end
if opt.force
  int_t = tic;
  force = int_fcn_force(H(1:M(1), 1:M(2), 1:M(3)));
  walltime.int = walltime.int + toc(int_t);
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

if nargout == 2
  walltime.total = sum(struct2array(walltime));
end
