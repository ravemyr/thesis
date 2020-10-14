function comp_figT1
% Computations for Figure T1
% Show relative RMS error and runtime as a function of P

  rng(1);

  % Parameters
  N = 100; % number of source particles
  L = 1; % box side length
  bP = 2.5; % Kaiser shape parameter
  Mvec = 2:2:40; % uniform grid size
  xi = 2*pi;

  n_rep = 1; % number of repetitions (for timing)
  USE_KAISER_REF = false;

  % Set up option struct
  per = [3 2 1 0];
  opt = cell(size(per));
  for i=1:numel(per)
    opt{i} = set_params(per(i), L, 0, xi, bP);
  end

  % Random sources
  [x, f] = NEW_vector_system(N, opt{1}.box);
  Q = sum(f.^2);
  A = sqrt(Q*xi*L)/L

  % Compute reference solutions
  ref = cell(size(per));
  for i=1:numel(per)
    fprintf('Computing reference solution for %dP ...\n', per(i));
    ref{i} = compute_reference(per(i), x, f, opt{i}, USE_KAISER_REF);
  end

  % Compute Spectral Ewald solutions, different window functions
  windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
  Pvec = [8 16 24];
  u_se = cell(numel(per), numel(windows), numel(Pvec), numel(Mvec));
  err = NaN(numel(per), numel(windows), numel(Pvec), numel(Mvec)); % error
  time = cell(numel(per), numel(windows), numel(Pvec), numel(Mvec)); % runtime
  for i=1:numel(per)
    fprintf('\n');
    for j=1:numel(Mvec)
    for w=1:numel(windows)
      fprintf('%dP Laplace Spectral Ewald, window=%s ', per(i), windows{w});
      for k=1:numel(Pvec)
        opti = opt{i};
        opti.M = Mvec(j) * opti.box
        opti.window = windows{w};
        opti.P = Pvec(k);
        opti.polynomial_degree = max(min(opti.P+1, 9), 4);
        opti = modify_M(per(i), opti);
        %if opti.P > 16 && ~strcmp(opti.window, 'gaussian')
        if opti.P > 16 && strcmp(opti.window, 'kaiser_poly')
          continue;
        end
        fprintf('.');
        [u, t] = compute_spectral_ewald(per(i), x, f, opti, n_rep);
        u_se{i,w,k,j} = u;
        time{i,w,k,j} = t;
        % Compute error
        err(i,w,k,j) = rms(u - ref{i}) / rms(ref{i});
      end
      fprintf('\n');
    end
    end
  end

  save('data/figT1.mat', 'N', 'L', 'bP', 'xi', 'Mvec', 'per', 'opt', ...
       'x', 'f', 'ref', 'windows', 'Pvec', 'u_se', 'err', 'time');
  plot_figT1;
end

function opt = set_params(periodicity, L, M, xi, bP)
  opt.box = [L L L];
  opt.M = M * opt.box;
  opt.xi = xi;
  opt.rc = 6/opt.xi;
  opt.betaP = bP;
  if periodicity == 2
    opt.s0 = 2;
    opt.s = 3.5;
    opt.n = 6;
  elseif periodicity == 1
    opt.s0 = 2.6;
    opt.s = 4;
    opt.n = 8;
  elseif periodicity == 0
    opt.s = 2.8;
  end
end

function opt = modify_M(periodicity, opt)
% if periodicity == 2
%   opt.add_M3 = 6; % FIXME: why is this needed?
% end
end

function u = compute_reference(periodicity, x, f, opt, USE_KAISER)
  N = size(x,1);
  if USE_KAISER
    opt.window = 'expsemicirc';
    opt.P = 20;
  else
    ref_opt.xi = opt.xi; ref_opt.box = opt.box;
    %ref_opt.layers = ceil((max(opt.M)-1)/2);
    ref_opt.layers = ceil((28-1)/2);
  end
  if periodicity == 3
    if USE_KAISER
      u = SE3P_Laplace_fourier_space(1:N, x, f, opt);
    else
      u = SE3P_Laplace_direct_fd_mex(1:N, x, f, ref_opt);
    end
  elseif periodicity == 2
    if USE_KAISER
      u = SE2P_Laplace_fourier_space(1:N, x, f, opt);
    else
      uf = SE2P_Laplace_direct_fd_mex(1:N, x, f, ref_opt);
      u0 = SE2P_Laplace_direct_k0_mex(1:N, x, f, ref_opt);
      u = uf+u0;
    end
  elseif periodicity == 1
    if USE_KAISER
      u = SE1P_Laplace_fourier_space(1:N, x, f, opt);
    else
      uf = SE1P_Laplace_direct_fd_mex(1:N, x, f, ref_opt);
      u0 = SE1P_Laplace_direct_k0_mex(1:N, x, f, ref_opt);
      u = uf+u0;
    end
  elseif periodicity == 0
    if USE_KAISER
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
      [u, time] = SE2P_Laplace_fourier_space(1:N, x, f, opt);
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
