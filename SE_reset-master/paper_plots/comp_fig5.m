function varargout = comp_fig5(ref)
% Computations for Figure 5
% Show relative RMS error and runtime as a function of P

  rng(1);

  % Parameters
  N = 100; % number of source particles
  L = 1; % box side length
  bP = 2.5; % Kaiser shape parameter
  M = 28; % uniform grid size
  xi = pi*M/12; % Ewald split parameter

  n_rep = 1; % number of repetitions (for timing)
  USE_SE_REF = false;

  % Set up option struct
  per = [3 2 1 0];
  opt = cell(size(per));
  for i=1:numel(per)
    opt{i} = set_params(per(i), L, M, xi, bP);
  end

  % Random sources
  [x, f] = NEW_vector_system(N, opt{1}.box);
  Q = sum(f.^2);
  A = sqrt(Q*xi*L)/L;

  % Compute reference solutions
  if ~exist('ref', 'var') || isempty(ref)
    ref = cell(1,4);
    for i=1:numel(per)
      fprintf('Computing reference solution for %dP ...\n', per(i));
      ref{per(i)+1} = compute_reference(per(i), x, f, opt{i}, USE_SE_REF);
    end
  end

  % Compute Spectral Ewald solutions, different window functions
  windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
  Pvec = 2:2:28;
  u_se = cell(numel(per), numel(windows), numel(Pvec));
  err = NaN(numel(per), numel(windows), numel(Pvec)); % error
  time = cell(numel(per), numel(windows), numel(Pvec)); % runtime
  for i=1:numel(per)
    fprintf('\n');
    for w=1:numel(windows)
      fprintf('%dP Laplace Spectral Ewald, window=%s ', per(i), windows{w});
      for k=1:numel(Pvec)
        opti = opt{i};
        opti.window = windows{w};
        opti.P = Pvec(k);
        opti.polynomial_degree = max(min(opti.P+1, 9), 4);
        %if opti.P > 16 && ~strcmp(opti.window, 'gaussian')
        if opti.P > 16 && strcmp(opti.window, 'kaiser_poly')
          continue;
        end
        fprintf('.');
        [u, t] = compute_spectral_ewald(per(i), x, f, opti, n_rep);
        u_se{i,w,k} = u;
        time{i,w,k} = t;
        % Compute error
        err(i,w,k) = rms(u - ref{per(i)+1}) / rms(ref{per(i)+1});
      end
      fprintf('\n');
    end
  end

  save('data/fig5.mat', 'N', 'L', 'bP', 'xi', 'M', 'per', 'opt', ...
       'x', 'f', 'ref', 'windows', 'Pvec', 'u_se', 'err', 'time');
  plot_fig5;

  if nargout > 0
    varargout{1} = ref;
  end
end

function opt = set_params(periodicity, L, M, xi, bP)
  opt.box = [L L L];
  opt.M = M * opt.box;
  opt.xi = xi;
  opt.rc = 6/opt.xi;
  opt.betaP = bP;
end

function u = compute_reference(periodicity, x, f, opt, USE_SE)
  N = size(x,1);
  if USE_SE
    opt.window = 'gaussian';
    opt.P = 32;
  else
    ref_opt.box = opt.box; ref_opt.xi = opt.xi;
    ref_opt.layers = ceil((max(opt.M)-1)/2);
  end
  if periodicity == 3
    if USE_SE
      u = SE3P_Laplace_fourier_space(1:N, x, f, opt);
    else
      u = SE3P_Laplace_direct_fd_mex(1:N, x, f, ref_opt);
    end
  elseif periodicity == 2
    if USE_SE
      u = SE2P_Laplace_fourier_space(1:N, x, f, opt);
    else
      uf = SE2P_Laplace_direct_fd_mex(1:N, x, f, ref_opt);
      u0 = SE2P_Laplace_direct_k0_mex(1:N, x, f, ref_opt);
      u = uf+u0;
    end
  elseif periodicity == 1
    if USE_SE
      u = SE1P_Laplace_fourier_space(1:N, x, f, opt);
    else
      uf = SE1P_Laplace_direct_fd_mex(1:N, x, f, ref_opt);
      u0 = SE1P_Laplace_direct_k0_mex(1:N, x, f, ref_opt);
      u = uf+u0;
    end
  elseif periodicity == 0
    if USE_SE
      SE0P_warnings('off');
      u = SE0P_Laplace_fourier_space(1:N, x, f, opt);
    else
      u = SE0P_Laplace_direct_full(1:N, x, f);
      ur = SE0P_Laplace_real_space(1:N, x, f, opt);
      us = -f * opt.xi * 2/sqrt(pi);
      u = u - ur - us;
    end
  end
end

function [u, time] = compute_spectral_ewald(periodicity, x, f, opt, n_rep)
  N = size(x,1);
  if periodicity == 3
    for n=1:n_rep
      [u, ~, time] = SE3P_Laplace_fourier_space(1:N, x, f, opt);
    end
  elseif periodicity == 2
    for n=1:n_rep
      [u, ~, time] = SE2P_Laplace_fourier_space(1:N, x, f, opt);
    end
  elseif periodicity == 1
    for n=1:n_rep
      [u, ~, time] = SE1P_Laplace_fourier_space(1:N, x, f, opt);
    end
  elseif periodicity == 0
    SE0P_warnings('off');
    pre_kernel = SE0P_Laplace_kernel_fft_precomp(SE0P_parse_params(opt));
    for n=1:n_rep
      [u, ~, time] = SE0P_Laplace_fourier_space(1:N, x, f, opt, pre_kernel);
    end
  end
end
