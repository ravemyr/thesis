function [nG, nK] = determine_single_n(M, verbose)
% [nG, nK] = determine_single_n(M, verbose)
% Determine the size of the local pad (opt.n) for given M.
% Returns n for Gaussian window and Kaiser window.

if ~exist('M', 'var') || isempty(M)
  M = 36;
end
if ~exist('verbose', 'var') || isempty(verbose)
  verbose = true;
end

rng(1);

N = 100; % number of source particles
L = 1; % box side length
box = [L L L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

opt.M = M*box/L;
opt.xi = pi*(M/L)/12;
opt.box = box;

windows = {'gaussian', 'kaiser_poly'};
for w=1:numel(windows)
  cprintf(verbose, '\n== Window %s ==\n', windows{w});
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

  % Reference
  Ropt = opt;
  Ropt.n = max(ceil(opt.M(1)/2),1);
  ref = SE1P_Laplace_fourier_space(1:Neval, x, f, Ropt);

  n = 0;
  while true
    opt.n = n;
    uf = SE1P_Laplace_fourier_space(1:Neval, x, f, opt);
    rms_error = rms(uf-ref) / rms(ref);
    cprintf(verbose, 'n=%d; error: %.16g\n', n, rms_error);
    if rms_error < 1e-15
      break
    end
    n = n + 1;
  end
  if strcmp(windows{w}, 'gaussian')
    nG = n;
  else
    nK = n;
  end
end
