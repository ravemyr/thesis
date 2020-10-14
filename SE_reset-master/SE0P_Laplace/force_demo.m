% Spectral Ewald 0P (free-space) Laplace, force computation

clear
rng(1);

N = 10; % number of source particles
L = 1; % box side length
box = [L L L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 28; % Set M0 to an even number, the rest is automatic

opt.M = M0*box;
opt.xi = pi*M0/12;
opt.rc = 6/opt.xi;
opt.box = box;
opt.s = 3.2; % oversampling factor
% FIXME: an alias for opt.s is opt.oversampling

%% Direct summation
t = tic();
ref = SE0P_Laplace_direct_full_force(1:Neval, x, f);
tdirect = toc(t);

%% Spectral Ewald
SE0P_warnings('off');
%windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
windows = {'gaussian'};
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
  opt.potential = false; opt.force = true;
  t = tic();
  pre_kernel = SE0P_Laplace_kernel_fft_precomp(SE0P_parse_params(opt));
  [~, uf, walltime] = SE0P_Laplace_fourier_space(1:Neval, x, f, opt, pre_kernel);
  ur = SE0P_Laplace_real_space_force(1:Neval, x, f, opt);
  u = uf + ur;
  tSE = toc(t);
  rms_error = rms(u(:)-ref(:)) / rms(ref(:));
  fprintf('rms_error = %.16g  [time elapsed: %g sec]\n\n', rms_error, tSE);
end
