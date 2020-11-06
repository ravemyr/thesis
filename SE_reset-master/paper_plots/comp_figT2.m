function comp_figT2
% Computations for Figure T2
% Show relative RMS error and runtime as a function of P

  rng(1);

  % Parameters
  N = 100; % number of source particles
  L = 1; % box side length
  bP = 2.5; % Kaiser shape parameter
  M = 12; % uniform grid size

  % Random sources
  box = [L L L];
  [x, f] = NEW_vector_system(N, box);
  Q = sum(f.^2);

  Pvec = 2:2:28;
  etol = compute_etol(Pvec) / 10;
  % Determine Ewald split parameter
  xiv = (pi/sqrt(2)) * (M/L) ./ sqrt(lambertw(4*Q./(pi^2 * L^2 * etol.^2 * M)));
  kinf = pi*M/L;
  %lBoundXi = sqrt(lambertw(1./etol * sqrt(Q/(2*L^3))));
  %Mvec = L*sqrt(3)*lBoundXi/pi .* sqrt(lambertw(4*Q^(2/3) ./ (3*L^2*(pi*lBoundXi.*etol.^2).^(2/3))));

  per = [3 2 1 0];
  n_rep = 1; % number of repetitions (for timing)
  USE_SE_REF = false;

  % Compute Spectral Ewald solutions, different window functions
  %windows = {'gaussian', 'expsemicirc', 'kaiser_exact', 'kaiser_poly'};
  windows = {'gaussian', 'kaiser_exact'};
  u_se = cell(numel(per), numel(windows), numel(Pvec));
  err = NaN(numel(per), numel(windows), numel(Pvec)); % error
  time = cell(numel(per), numel(windows), numel(Pvec)); % runtime
  for i=1:numel(per)
    fprintf('[%dP]\n', per(i));
    for k=1:numel(Pvec)
      xi = xiv(k);
      opti = set_params(per(i), box, M, xi, bP);
      fprintf('%d', Pvec(k));
      refi = compute_reference(per(i), x, f, opti, USE_KAISER_REF);
      for w=1:numel(windows)
        opti.window = windows{w};
        opti.P = Pvec(k);
        opti.polynomial_degree = max(min(opti.P+1, 9), 4);
        if opti.P > 16 && strcmp(opti.window, 'kaiser_poly')
          continue;
        end
        fprintf('.');
        [u, t] = compute_spectral_ewald(per(i), x, f, opti, n_rep);
        u_se{i,w,k} = u;
        time{i,w,k} = t;
        % Compute error
        err(i,w,k) = rms(u - refi) / rms(refi);
      end
    end
    fprintf('\n');
  end

  save('data/figT2.mat', 'N', 'L', 'bP', 'xi', 'M', 'per', ...
       'x', 'f', 'windows', 'Pvec', 'u_se', 'err', 'time');
  plot_figT2;
end

function opt = set_params(periodicity, box, M, xi, bP)
  opt.box = box;
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

function etol = compute_etol(P)
  etol = zeros(size(P));
  for k=1:numel(P)
    if P(k) == 2
      etol(k) = 5e-4;
    elseif P(k) == 4
      etol(k) = 5e-4;
    elseif P(k) == 6
      etol(k) = 5e-6;
    elseif P(k) == 8
      etol(k) = 5e-8;
    elseif P(k) == 10
      etol(k) = 1e-10;
    elseif P(k) == 12
      etol(k) = 1e-12;
    elseif P(k) == 14
      etol(k) = 1e-13;
    elseif P(k) == 16
      etol(k) = 2e-14;
    elseif P(k) >= 18
      etol(k) = 5e-16;
    end
  end
end
