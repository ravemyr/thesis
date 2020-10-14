% Spectral Ewald 0P (free-space) Laplace, basic accuracy/convergence computation

clear
rng(1);

N = 10000; % number of source particles
L = 3; % box side length
box = [L L L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 20; % Set M0 to an even number, the rest is automatic

opt.M = M0*box;
opt.xi = pi*M0/12;
opt.rc = 6/opt.xi;
opt.rc = min(opt.rc, L);
opt.box = box;
opt.s = 2.8; % oversampling factor
% FIXME: an alias for opt.s is opt.oversampling

%% Direct summation
t = tic();
ref = SE0P_Laplace_direct_full(1:Neval, x, f);
tdirect = toc(t);

%% Spectral Ewald
SE0P_warnings('off');
windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
for w=1:numel(windows)
  fprintf('0P Laplace Ewald, window=%s\n', windows{w});
  opt.window = windows{w};
  if strcmp(windows{w}, 'gaussian')
    opt.P = 24;
    if isfield(opt, 'betaP'), rmfield(opt, 'betaP'), end
  else
    opt.P = 16;
    opt.betaP = 2.5;
  end
  if strcmp(windows{w}, 'kaiser_poly')
    opt.polynomial_degree = 9;
  end
  t = tic();
  pre_kernel = SE0P_Laplace_kernel_fft_precomp(SE0P_parse_params(opt));
  [uf, ~, walltime] = SE0P_Laplace_fourier_space(1:Neval, x, f, opt, pre_kernel);
  ur = SE0P_Laplace_real_space(1:Neval, x, f, opt);
  us = -f * opt.xi * 2/sqrt(pi);
  u = uf + ur + us;
  tSE = toc(t);
  rms_error = rms(u-ref) / rms(ref);
  fprintf('rms_error = %.16g  [time elapsed: %g sec]\n\n', rms_error, tSE);
end
