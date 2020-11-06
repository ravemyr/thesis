function opt = SE0P_parse_params(opt)
% opt = SE0P_parse_params(opt)
%
% Parse parameters and set up free-space Ewald grid sizes:
% inner = domain containing all sources and targets
% extended = FGG grid (to cancel wrap effects)
% padded = 2x FFT grid (for aperiodic convolution)
% oversampled = FFT grid for truncated Green's function
%
% For documentation, see SE0P_Laplace_fourier_space.

% Check that mandatory options are present
assert(isfield(opt, 'box'), 'cell size box must be given in opt struct');
assert(isfield(opt, 'M'), 'grid size M must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald parameter xi must be given in opt struct');
assert(isfield(opt, 'P'), 'window support P must be given in opt struct');

if ~isfield(opt,'potential'), opt.potential = true; end
if ~isfield(opt,'force'), opt.force = false; end
if ~isfield(opt,'fast_gridding'), opt.fast_gridding = true; end
if ~isfield(opt,'base_factor'), opt.base_factor = 4; end
% The following value for add_M works well when the periodic cell is a cube.
% Otherwise is might be completely off.
if ~isfield(opt,'add_M'), opt.add_M = [6 6 6]; end

% Inner grid is the given input
inner_box = opt.box;
inner_M = opt.M;
opt.L = opt.box(1);
opt.h = opt.L / opt.M(1); % step size (M contains number of subintervals)
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

% Extended grid
if ischar(opt.add_M) && strcmp(opt.add_M, 'cover_remainder')
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
else
  % Default method
  deltaM = opt.P + opt.add_M;
end
opt.extended_M = inner_M + deltaM;
opt.extended_M = opt.base_factor * ceil(opt.extended_M / opt.base_factor); % round up
opt.extended_box = opt.h * opt.extended_M; % adjust box side length

% Padded grid
opt.padded_box = opt.extended_box * 2;
opt.padded_M = opt.extended_M * 2;

% Oversampled grid
if ~isfield(opt,'s'), opt.s = 2.8; end % oversampling factor
opt.oversampling = opt.s; % FIXME: oversampling is used as an alias for s
overM = opt.base_factor * ceil(opt.oversampling * opt.extended_M / opt.base_factor); % round up
% TODO FIXME: would it make more sense to put ./ and .* in the
% two lines below?
actual_oversampling = overM / opt.extended_M; % least-squares solution?
if abs(actual_oversampling - opt.oversampling) > 10*eps(actual_oversampling)
  warning('FSE:OversamplingIncreased',...
          'Oversampling factor increased to %g achieve integer grid', ...
          actual_oversampling);
end
opt.oversampled_box = opt.extended_box * actual_oversampling;
opt.oversampled_M = overM;

if isfield(opt, 'oversample_all') && opt.oversample_all == true
  opt.padded_box = opt.oversampled_box;
  opt.padded_M = opt.oversampled_M;
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

opt.R = norm(opt.extended_box);
opt.deltaLhalf = opt.h * deltaM/2;
