% Spectral Ewald 2P Laplace, basic accuracy/convergence computation

clear
rng(1);

N = 20; % number of source particles
L = 2; % box side length
box = [L L L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 20; % Set M0 to an even number, the rest is automatic

opt.M = M0*box;
opt.xi = pi*M0/12;
opt.box = box;
opt.s = 4; % oversampling factor on local pad
opt.n = 8; % local pad
opt.s0 = 2; % oversampling factor of zero mode

%% Direct summation
Dopt.layers = ceil((max(opt.M)-1)/2);
Dopt.xi = opt.xi; Dopt.box = opt.box;
t = tic();
ref = SE2P_Laplace_direct_fd_mex(1:Neval, x, f, Dopt);
ref0 = SE2P_Laplace_direct_k0_mex(1:Neval, x, f, Dopt);
ref = ref+ref0;
tdirect = toc(t);

%% Spectral Ewald
windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
for w=1:numel(windows)
  fprintf('2P Laplace Ewald, window=%s\n', windows{w});
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
  [uf, walltime] = SE2P_Laplace_fourier_space(1:Neval, x, f, opt);
  tSE = toc(t);
  rms_error = rms(uf-ref) / rms(ref);
  fprintf('rms_error = %.16g  [time elapsed: %g sec]\n\n', rms_error, tSE);
end
