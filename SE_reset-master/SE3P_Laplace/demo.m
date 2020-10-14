% Spectral Ewald 3P Laplace, basic accuracy/convergence computation

clear
rng(1);

N = 100; % number of source particles
L = 1; % box side length
box = [L L L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 24; % Set M0 to an even number, the rest is automatic

opt.box = box;
opt.M = M0*opt.box;
opt.xi = pi*M0/12;
opt.rc = 6/opt.xi;
assert(opt.rc <= L, 'rc (%g) cannot be larger than L (%g)', opt.rc, L);

%% Direct summation, real space
Dopt.box = opt.box; Dopt.xi = opt.xi; Dopt.real_cutoff = opt.rc;
Dopt.layers = 1; % assuming rc <= L, this is enough for realspace
t = tic();
ref_real = SE3P_Laplace_direct_real_rc_mex(1:Neval, x, f, Dopt);
tDirReal = toc(t);
% NB: The only difference between SE3P_Laplace_direct_real_rc_mex
% and SE3P_Laplace_direct_real_mex is that the former removes any
% interaction where r > rc.
%% Direct summation, Fourier space
Dopt.layers = ceil((max(opt.M)-1)/2);
t = tic();
ref_fd = SE3P_Laplace_direct_fd_mex(1:Neval, x, f, Dopt);
tDirFour = toc(t);
%% Put together
ref_self = SE3P_Laplace_direct_self_mex(1:Neval, x, f, Dopt);
ref = ref_fd + ref_real + ref_self;

%% Spectral Ewald, real space
t = tic();
ur = SE3P_Laplace_real_space(1:Neval, x, f, opt);
tSEReal = toc(t);
us = -f * opt.xi * 2/sqrt(pi);
%% Spectral Ewald, Fourier space
windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
for w=1:numel(windows)
  fprintf('3P Laplace Spectral Ewald, window: %s\n', windows{w});
  opt.window = windows{w};
  if strcmp(windows{w}, 'gaussian')
    opt.P = 32;
    if isfield(opt, 'betaP'), rmfield(opt, 'betaP'), end
  else
    opt.P = 16;
    opt.betaP = 2.5;
  end
  if strcmp(windows{w}, 'kaiser_poly')
    opt.polynomial_degree = 9;
  end
  t = tic();
  [uf, ~, walltime] = SE3P_Laplace_fourier_space(1:Neval, x, f, opt);
  tSEFour = toc(t);
  u = uf + ur + us;
  rms_error = rms(u-ref) / rms(ref);
  fprintf('  RMS error: %.16g\n', rms_error);
  fprintf('  Time: %g s (Fourier), %g s (real), total %g s\n\n', ...
          tSEFour, tSEReal, tSEFour+tSEReal);
end
