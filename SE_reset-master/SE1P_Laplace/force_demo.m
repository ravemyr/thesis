% Spectral Ewald 1P Laplace, force computation

clear
rng(1);

N = 20; % number of source particles
L = 1; % box side length
box = [L 1.1*L 1.2*L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 20; % Set M0 to an even number, the rest is automatic

opt.M = M0*box;
opt.xi = pi*M0/12;
opt.box = box;
%opt.sg = 1; % global oversampling factor
opt.s = 4; % oversampling factor on local pad
opt.n = 4; % local pad
opt.s0 = 2.5; % oversampling factor of zero mode
% FIXME: seems s and n are called sl and nl for GAUSSIAN 1P

%% Direct summation
Dopt.layers = ceil((max(opt.M)-1)/2);
Dopt.xi = opt.xi; Dopt.box = opt.box;
t = tic();
ref = SE1P_Laplace_direct_fd_force_mex(1:Neval, x, f, Dopt);
ref0 = SE1P_Laplace_direct_k0_force_mex(1:Neval, x, f, Dopt);
ref = ref+ref0;
tdirect = toc(t);

%% Spectral Ewald
%windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
windows = {'gaussian'};
for w=1:numel(windows)
  fprintf('1P Laplace Ewald, window=%s\n', windows{w});
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
  opt.potential = false; opt.force = true;
  t = tic();
  [~, uf, walltime] = SE1P_Laplace_fourier_space(1:Neval, x, f, opt);
  tSE = toc(t);
  rms_error = rms(uf(:)-ref(:)) / rms(ref(:));
  fprintf('rms_error = %.16g  [time elapsed: %g sec]\n\n', rms_error, tSE);
end
