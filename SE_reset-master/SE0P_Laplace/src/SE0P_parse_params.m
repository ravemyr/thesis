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
assert(isfield(opt, 's'), 'oversampling factor s must be given in opt struct');
opt.oversampling = opt.s; % FIXME: oversampling is used as an alias for s

if ~isfield(opt,'potential'), opt.potential = true; end
if ~isfield(opt,'force'), opt.force = false; end

% Inner grid is the given input
inner_box = opt.box;
inner_M = opt.M;
opt.L = opt.box(1);
opt.h = opt.L/opt.M(1); % step size (M contains number of subintervals)
% Check that h is the same in all directions
assert(abs(opt.h - opt.box(2)/opt.M(2)) < eps);
assert(abs(opt.h - opt.box(3)/opt.M(3)) < eps);

% Window options
if ~isfield(opt,'window'), opt.window = 'gaussian'; end
if ~isfield(opt,'w'), opt.w = opt.h*opt.P/2; end
if ~isfield(opt,'m'), opt.m = 0.94*sqrt(pi*opt.P); end
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
% This deltaB should cover the grid window
deltaB_grid = opt.h * opt.P / 2;
% This deltaB should cover the remainder window
deltaB_rem = deltaB_grid;
if strcmp(opt.window,'expsemicirc') || strcmp(opt.window,'kaiser_exact') ...
    || strcmp(opt.window,'kaiser_poly')
  deltaB_rem = sqrt(2*opt.betaP)/opt.xi;
elseif opt.eta < 1
  deltaB_rem = sqrt((1-opt.eta)/2)*opt.m/opt.xi;
end
% Take the largest
deltaB = max(deltaB_grid, deltaB_rem);
if isfield(opt, 'no_extra_support') && opt.no_extra_support == true
  deltaB = deltaB_grid;
end
% If inner grid is even, make extended grid even too
EVEN_GRIDS = (mod(inner_M(1), 2) == 0);
if EVEN_GRIDS
  deltaM = 2*ceil(deltaB / opt.h);
else
  deltaM = ceil(2*deltaB / opt.h);
end
deltaB = opt.h * deltaM / 2;
opt.extended_box = inner_box + 2*deltaB;
opt.extended_M = inner_M + deltaM;

% Padded grid
opt.padded_box = opt.extended_box * 2;
opt.padded_M = opt.extended_M * 2;

% Oversampled grid
% Ensure integer oversampling factor
if EVEN_GRIDS
  overM = 2 * ceil(opt.oversampling * opt.extended_M / 2);
else
  overM = ceil(opt.oversampling * opt.extended_M);
end
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

if ~isfield(opt,'fast_gridding'), opt.fast_gridding = true; end

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
opt.deltaB = deltaB;
