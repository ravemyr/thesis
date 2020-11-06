function se0p_opt = se0p_parse_params(opt)
% Setup Free-Space Ewald grid sizes
% inner = domain containing all sources and targets
% extended = FGG grid (to cancel wrap effects)
% padded = 2x FFT grid (for aperiodic convolution)
% oversampled = FFT grid for truncated Green's function

se0p_opt = opt;

inner_M = opt.M;
inner_box = opt.box;
h = opt.box(1) / opt.M(1);
popt = se3p_parse_params(opt);
P = popt.P;
m = popt.m;
xi = popt.xi;
eta = popt.eta;

se0p_opt.L = opt.box(1);
se0p_opt.h = h;
se0p_opt.p_half = se0p_opt.P/2;

if isfield(opt, 'beta')
  se0p_opt.beta = opt.beta * P;
else
  se0p_opt.beta = 2.5 * P;
end

% Old delta (should cover grid Gaussian)
delta_old = h*P/2;
if isfield(opt, 'window') && strcmp(opt.window, 'kaiser')
  %delta = delta_old;
  %delta = max(delta_old, sqrt(popt.beta*2)/xi);
  delta = max(delta_old, sqrt(popt.beta/4)/xi);
else
  % New delta (should cover remainder Gaussian)
  if popt.eta < 1
      deltaNew = sqrt((1-eta)/2)*m/xi;
      %deltaNew = delta_old;
  else
      deltaNew = delta_old;
  end
  % Take max
  delta = max(deltaNew, delta_old);
end

if isfield(opt, 'no_extra_support') && opt.no_extra_support == true;
    delta = delta_old;
end
%fprintf('P=%d, eta=%.2f, delta1=%.2f, delta2=%.2f, delta=%.2f\n', ...
%        P, eta, delta_old, deltaNew, delta);

EVEN_GRIDS = true;
if EVEN_GRIDS
    deltaM = 2*ceil(delta / h);
else
    deltaM = ceil(2*delta / h);
end
delta = h*deltaM/2;

se0p_opt.extended_box = inner_box + 2*delta;
se0p_opt.extended_M = inner_M + deltaM;

se0p_opt.padded_M = se0p_opt.extended_M * 2;
se0p_opt.padded_box = se0p_opt.extended_box * 2;

% Ensure integer oversampling rate
if EVEN_GRIDS 
    overM = ceil(opt.oversampling * se0p_opt.extended_M/2)*2;
else
    overM = ceil(opt.oversampling * se0p_opt.extended_M);
end
actual_oversampling = overM / se0p_opt.extended_M;
if abs(actual_oversampling - opt.oversampling) > 10*eps(actual_oversampling)
    warning('FSE:OversamplingIncreased',...
            'Oversampling rate increased to %g achieve integer grid', ...
            actual_oversampling);
end
se0p_opt.oversampling = opt.oversampling;
se0p_opt.oversampled_M = overM; 
se0p_opt.oversampled_box = se0p_opt.extended_box * actual_oversampling;

se0p_opt.R = norm(se0p_opt.extended_box);
se0p_opt.delta = delta;

if isfield(opt, 'oversample_all') && opt.oversample_all == true;
    se0p_opt.padded_M = se0p_opt.oversampled_M;
    se0p_opt.padded_box = se0p_opt.oversampled_box;
end

se0p_opt.m = m;
se0p_opt.eta = eta;
se0p_opt.c = popt.c;
