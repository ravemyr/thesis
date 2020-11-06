function [pot, force, walltime] = SE1P_Laplace_fourier_space(eval_idx, x, f, opt)
% SE1P_LAPLACE_FOURIER_SPACE
% Compute Fourier-space part of Ewald sum for 1-periodic electrostatic potential,
% using the Spectral Ewald method.
%
% Periodicity is assumed in the first direction (x), while the
% two remaining directions (y and z) are assumed free.
%
% pot = SE1P_Laplace_fourier_space(eval_idx, x, f, opt)
%   Return potential
%
% [pot, force] = SE1P_Laplace_fourier_space(...)
%   Return potential and force
%
% [pot, force, walltime] = SE1P_Laplace_fourier_space(...)
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
% :param opt.add_M2:        Extra number to add to M2 (default: 6), see below
% :param opt.add_M3:        Extra number to add to M3 (default: 6), see below.
%                           If either add_M2 or add_M3 is set to the special
%                           value 'cover_remainder', the grid in that direction is
%                           extended to cover the "remainder Gaussian".
% :param opt.base_factor:   Integer that the final grid size should be divisible by
%                           (default: 4); the grid size will be rounded up
% :param opt.xi:            Ewald parameter (required)
% :param opt.P:             Support of window, in number of grid points (required)
% :param opt.s0:            Oversampling factor for zero mode (default: 2.5)
% :param opt.s:             Oversampling factor for local pad (default: 4)
% :param opt.sg:            Oversampling factor for the rest of the domain (default: 1)
% :param opt.n:             Used to define which modes to include in the local pad
%                           (default: ceil(M1/5))
% :param opt.window:        Window function (default: 'gaussian')
% :param opt.w:             Width of window (default: w=h*P/2)
% :param opt.m:             Gaussian shape function (default: m=0.95*sqrt(pi*P))
% :param opt.eval_x:        External points to evaluate potential in (Nex×3)
% :param opt.fast_gridding: Use fast gridding and spreading (default: true)
%
% The grid size in the free directions (M2 and M3) will be set to
%   M2 = opt.M(2) + opt.P + opt.add_M2,
%   M3 = opt.M(3) + opt.P + opt.add_M3,
% i.e. both opt.P and opt.add_M2 (or opt.add_M3) are added to the
% values provided in the opt struct.
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
%       sources outside of opt.box in the free directions. It is
%       up to the caller to ensure that the sources are inside
%       the box.
%
% :returns: **pot** -- Fourier-space potential
% :returns: **force** -- Fourier-space force
% :returns: **walltime** -- struct containing timings

opt = parse_params(opt);

% Initialize timing array
walltime = struct('grid',0,'fft',0,'scale',0,'ifft',0,'int',0,'pre',0,'prefft',0);

% Window function specifics
pre = [];
W_precomp_gridding = @(x,opt) SE1P_Laplace_gridding_precomp(x, opt);
if strcmp(opt.window, 'gaussian')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_thrd_mex_1p(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  W_fast_int_pot = @(F,opt,S) SE_fg_int_split_mex_1p(0,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_mex_1p is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_fast_int_force = @(F,opt,S) SE_fg_int_split_force_mex_1p(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,...
                                S.zfx,S.zfy,S.zfz,S.idx);
  W_plain_grid = @SE_fg_grid_mex_1p;
  W_plain_int_pot = @SE_fg_int_mex_1p;
  W_plain_int_force = @SE_fg_int_force_mex_1p;
elseif strcmp(opt.window, 'expsemicirc') || strcmp(opt.window, 'kaiser_exact') ...
    || strcmp(opt.window, 'kaiser_poly')
  if opt.force
    error('Force calculations are only supported for the Gaussian window');
  end
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_kaiser_mex_1p(x,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  W_fast_int_pot = @(F,opt,S) SE_fg_int_split_kaiser_mex_1p(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_kaiser_mex_1p is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
  W_fast_int_force = @(F,opt,S) SE_fg_int_split_kaiser_mex_1p(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx); % FIXME
  W_plain_grid = @SE_fg_grid_kaiser_mex_1p;
  W_plain_int_pot = @SE_fg_int_kaiser_mex_1p;
  W_plain_int_force = @SE_fg_int_kaiser_mex_1p; % FIXME

  if strcmp(opt.window, 'expsemicirc')
    % Precompute window Fourier transform
    prefft_t = tic;
    pre = SE1P_window_fft_precomp(opt);
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
  warning('SE1P:SlowGridding', 'Using slow gridding routines');
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
[G, Gr, G0] = SE1P_Laplace_fftnd(H, opt);
walltime.fft = walltime.fft + toc(fft_t);

% Scaling
scale_t = tic;
[G, Gr, G0] = SE1P_Laplace_scaling(G, Gr, G0, opt, pre);
walltime.scale = walltime.scale + toc(scale_t);

% Inverse Fourier transform
ifft_t = tic;
F = SE1P_Laplace_ifftnd(G, Gr, G0, opt);
walltime.ifft = walltime.ifft + toc(ifft_t);

% Spreading and integration
pot = []; force = [];
if opt.potential
  int_t = tic;
  pot = 4*pi*int_fcn_pot(F);
  walltime.int = walltime.int + toc(int_t);
end
if opt.force
  int_t = tic;
  force = 4*pi*int_fcn_force(F);
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
if ~isfield(opt,'fast_gridding'), opt.fast_gridding = true; end
if ~isfield(opt,'base_factor'), opt.base_factor = 4; end
% The following estimate for add_M2 and add_M3 works well when
% the periodic cell is a cube. Otherwise is might be completely off.
if ~isfield(opt,'add_M2'), opt.add_M2 = 6; end
if ~isfield(opt,'add_M3'), opt.add_M3 = 6; end

% Basic grid
opt.M(1) = opt.base_factor * ceil(opt.M(1) / opt.base_factor); % round up
opt.L = opt.box(1);
opt.h = opt.L / opt.M(1); % step size (M contains number of subintervals)

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

% Extend grid in free directions (y and z directions)
if (ischar(opt.add_M2) && strcmp(opt.add_M2, 'cover_remainder')) || ...
   (ischar(opt.add_M3) && strcmp(opt.add_M3, 'cover_remainder'))
  % Method used in legacy 0P code
  deltaM_grid = opt.P; % to cover the gridding window
  deltaM_rem = deltaM_grid; % to cover the "remainder Gaussian"
  if strcmp(opt.window,'expsemicirc') || strcmp(opt.window,'kaiser_exact') ...
      || strcmp(opt.window,'kaiser_poly')
    deltaM_rem = sqrt(8*opt.betaP)/opt.xi / opt.h;
  elseif opt.eta < 1
    deltaM_rem = sqrt(2*(1-opt.eta))*opt.m/opt.xi / opt.h;
  end
  deltaM = max(deltaM_grid, deltaM_rem);
  if isfield(opt, 'no_extra_support') && opt.no_extra_support == true
    deltaM = deltaM_grid;
  end
end
% y direction
if ischar(opt.add_M2) && strcmp(opt.add_M2, 'cover_remainder')
  opt.My = opt.M(2) + deltaM;
else
  % Default method
  opt.My = opt.M(2) + opt.P + opt.add_M2;
end
% z direction
if ischar(opt.add_M3) && strcmp(opt.add_M3, 'cover_remainder')
  opt.Mz = opt.M(3) + deltaM;
else
  % Default method
  opt.Mz = opt.M(3) + opt.P + opt.add_M3;
end
opt.My = opt.base_factor * ceil(opt.My / opt.base_factor); % round up
opt.Mz = opt.base_factor * ceil(opt.Mz / opt.base_factor); % round up
opt.Ly = opt.h * opt.My; % adjust box side length
opt.Lz = opt.h * opt.Mz; % adjust box side length
% Store values in given vectors as well
opt.M(2) = opt.My;
opt.M(3) = opt.Mz;
opt.box(2) = opt.Ly;
opt.box(3) = opt.Lz;
% Check that h is the same in all directions
assert(abs(opt.h - opt.Ly/opt.My) < eps);
assert(abs(opt.h - opt.Lz/opt.Mz) < eps);

wbox = [0 opt.L; -opt.w -opt.w+opt.Ly; -opt.w -opt.w+opt.Lz];
opt.free_offset = wbox(2:3,1);

% Sampling factor (oversampling)
if ~isfield(opt,'sg'), opt.sg = 1; end % global oversampling
if ~isfield(opt,'s'), opt.s = 4; end % oversampling for local pad
if ~isfield(opt,'s0'), opt.s0 = 2.5; end % oversampling for zero-mode
if ~isfield(opt,'n'), opt.n = ceil(opt.M(1)/5); end % size of local pad
% Increase s and s0 such that FFTN has integer size vectors.
opt.s = max(ceil(opt.s*opt.My)/opt.My, ceil(opt.s*opt.Mz)/opt.Mz);
opt.s0 = max(ceil(opt.s0*opt.My)/opt.My, ceil(opt.s0*opt.Mz)/opt.Mz);
opt.sl = opt.s; % FIXME: in GAUSSIAN code s is called sl (which
                % actually seems like a better name), and n is called nl

% Local pads, only if local oversampling is different from global oversampling
if opt.sg ~= opt.sl
  n = opt.n;
  % The local pad consists of k = -n:n, and the grid is -M/2:(M/2-1)
  % [if M is even] or -(M-1)/2:(M-1)/2 [if M is odd]. Thus, n can
  % at most be M/2-1 [if M is even] or (M-1)/2 [if M is odd].
  % These are both captured by floor((M-1)/2).
  n = min(floor((opt.M(1)-1)/2), n);
  % Zero mode is the first element, excluded here
  opt.local_pad = [2:n+1 opt.M(1)-n+1:opt.M(1)];
else
  opt.local_pad = 1;
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

opt.R = sqrt(opt.Ly^2 + opt.Lz^2);
