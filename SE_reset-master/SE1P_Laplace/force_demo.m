% Spectral Ewald 1P Laplace, force computation

clear
rng(1);

N = 100; % number of source particles
L = 1; % box side length
box = [L 1.1*L 1.2*L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 28; % Set M0=M/L, the rest is automatic

opt.box = box;
opt.M = M0*opt.box;
opt.xi = pi*M0/12;

%% Direct summation
Dopt.box = opt.box; Dopt.xi = opt.xi;
Dopt.layers = ceil((max(opt.M)-1)/2);
t = tic();
ref = SE1P_Laplace_direct_fd_force_mex(1:Neval, x, f, Dopt);
ref0 = SE1P_Laplace_direct_k0_force_mex(1:Neval, x, f, Dopt);
ref = ref+ref0;
tDirFour = toc(t);

%% Spectral Ewald
%windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
windows = {'gaussian'};
for w=1:numel(windows)
  fprintf('1P Laplace Spectral Ewald, window: %s\n', windows{w});
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
  opt.potential = false; opt.force = true;
  t = tic();
  [~, uf, walltime] = SE1P_Laplace_fourier_space(1:Neval, x, f, opt);
  tSEFour = toc(t);
  rms_error = rms(uf(:)-ref(:)) / rms(ref(:));
  fprintf('  RMS error: %.16g\n', rms_error);
  fprintf('  Time: %g s (Fourier)\n\n', tSEFour);
end
