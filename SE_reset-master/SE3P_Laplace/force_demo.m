% Spectral Ewald 3P Laplace, force computation

clear
rng(1);

N = 100; % number of source particles
L = 1; % box side length
box = [L L L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 28; % Set M0=M/L, the rest is automatic

opt.box = box;
opt.M = M0*opt.box;
opt.xi = pi*M0/12;
opt.rc = 6/opt.xi;
assert(opt.rc <= L, 'rc (%g) cannot be larger than L (%g)', opt.rc, L);

%% Direct summation, real space
Dopt.box = opt.box; Dopt.xi = opt.xi; Dopt.real_cutoff = opt.rc;
Dopt.layers = 1; % assuming rc <= L, this is enough for realspace
t = tic();
ref_real = SE3P_Laplace_direct_real_rc_force_mex(1:Neval, x, f, Dopt);
tDirReal = toc(t);
% NB: The only difference between SE3P_Laplace_direct_real_rc_force_mex
% and SE3P_Laplace_direct_real_force_mex is that the former removes any
% interaction where r > rc.
%% Direct summation, Fourier space
Dopt.layers = ceil((max(opt.M)-1)/2);
t = tic();
ref_fd = SE3P_Laplace_direct_fd_force_mex(1:Neval, x, f, Dopt);
tDirFour = toc(t);
%% Put together
ref = ref_fd + ref_real;

%% Spectral Ewald, real space
opt.potential = false; opt.force = true;
%opt.fourier_differentiation = true;
t = tic();
[~, ur] = SE3P_Laplace_real_space(1:Neval, x, f, opt);
tSEReal = toc(t);
%% Spectral Ewald, Fourier space
%windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
windows = {'gaussian', 'kaiser_exact', 'kaiser_poly'};
for w=1:numel(windows)
  fprintf('3P Laplace Spectral Ewald, window: %s\n', windows{w});
  opt.window = windows{w};
  if strcmp(windows{w}, 'gaussian')
    opt.P = 32;
    if isfield(opt, 'betaP'), rmfield(opt, 'betaP'); end
  else
    opt.P = 16;
    opt.betaP = 2.5;
  end
  if strcmp(windows{w}, 'kaiser_poly')
    opt.polynomial_degree = 9;
  end
  t = tic();
  [~, uf, walltime] = SE3P_Laplace_fourier_space(1:Neval, x, f, opt);
  tSEFour = toc(t);
  u = uf + ur;
  tSE = tSEReal + tSEFour;
  rms_error = rms(u(:)-ref(:)) / rms(ref(:));
  fprintf('  RMS error: %.16g\n', rms_error);
  fprintf('  Time: %g s (Fourier), %g s (real), total %g s\n\n', ...
          tSEFour, tSEReal, tSE);
end
