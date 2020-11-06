% Timing test, Fourier space part (no error computation)

clear
rng(1);

% Parameters
N = 1e4; % number of source particles
L = 1; % box side length
bP = 2.5; % Kaiser shape parameter
%windows = {'gaussian'};
windows = {'kaiser_poly'};
M = 128; % uniform grid size
%Pvec = 4;
Pvec = 16;
per = [0];
xi = pi*M/12; % Ewald split parameter

n_rep = 5; % number of repetitions (for timing)

% Set up option struct
opt = cell(size(per));
for i=1:numel(per)
  opt{i} = set_params(per(i), L, M, xi, bP);
end

% Random sources
[x, f] = NEW_vector_system(N, opt{1}.box);

% Compute Spectral Ewald solutions, different window functions
time = cell(numel(per), numel(windows), numel(Pvec)); % runtime
for i=1:numel(per)
  for w=1:numel(windows)
    for k=1:numel(Pvec)
      fprintf('%dP, window: %s, P=%d, M=%d\n', per(i), windows{w}, Pvec(k), M);
      opti = opt{i};
      opti.window = windows{w};
      opti.P = Pvec(k);
      %opti.polynomial_degree = max(min(opti.P+1, 9), 4);
      opti.polynomial_degree = 9;
      if opti.P > 16 && strcmp(opti.window, 'kaiser_poly')
        continue;
      end
      [u, t] = compute_spectral_ewald(per(i), x, f, opti, n_rep);
      time{i,w,k} = t;
      fprintf('  Fourier time: %g sec\n', t);
    end
  end
end

function opt = set_params(periodicity, L, M, xi, bP)
  opt.box = [L L L];
  opt.M = M * opt.box;
  opt.xi = xi;
  opt.rc = 6/opt.xi;
  opt.betaP = bP;
end

function [u, time] = compute_spectral_ewald(periodicity, x, f, opt, n_rep)
  N = size(x,1);
  if periodicity == 3
    [u, ~, time] = SE3P_Laplace_fourier_space(1:N, x, f, opt);
    t = tic();
    for n=1:n_rep
      [u, ~, ~] = SE3P_Laplace_fourier_space(1:N, x, f, opt);
    end
    time = toc(t) / n_rep;
  elseif periodicity == 2
    [u, ~, time] = SE2P_Laplace_fourier_space(1:N, x, f, opt);
    t = tic();
    for n=1:n_rep
      [u, ~, ~] = SE2P_Laplace_fourier_space(1:N, x, f, opt);
    end
    time = toc(t) / n_rep;
  elseif periodicity == 1
    [u, ~, time] = SE1P_Laplace_fourier_space(1:N, x, f, opt);
    t = tic();
    for n=1:n_rep
      [u, ~, ~] = SE1P_Laplace_fourier_space(1:N, x, f, opt);
    end
    time = toc(t) / n_rep;
  elseif periodicity == 0
    SE0P_warnings('off');
    pre_kernel = SE0P_Laplace_kernel_fft_precomp(SE0P_parse_params(opt));
    [u, ~, time] = SE0P_Laplace_fourier_space(1:N, x, f, opt, pre_kernel);
    t = tic();
    for n=1:n_rep
      [u, ~, ~] = SE0P_Laplace_fourier_space(1:N, x, f, opt, pre_kernel);
    end
    time = toc(t) / n_rep;
  end
end
