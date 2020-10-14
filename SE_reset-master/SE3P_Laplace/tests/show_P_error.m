function show_P_error

rng(1);

N = 100; % number of source particles
L = 1; % box side length
box = [L L L]; % periodic box
[x, f] = NEW_vector_system(N, box); % random sources
Neval = N; % number of evaluation points

% Ewald parameters
M0 = 28;
opt.M = M0*box;
%opt.xi = pi*M0/12;
opt.xi = 6.3;
opt.box = box;
opt.betaP = 2.5;

%% Reference solution
Dopt.layers = ceil((opt.M-1)/2);
Dopt.xi = opt.xi; Dopt.box = opt.box;
t = tic();
uref = compute_direct(x, f, Dopt, Neval);
tdirect = toc(t);

n_rep = 3; % number of repetitions
exclude_pre = true; % exclude precomputation from timing

%% Iteration for Spectral Ewald
windows = {'gaussian', 'expsemicirc', 'kaiser_exact'};%, 'kaiser_poly'};
Plist = 2:2:28;
errors = zeros(numel(windows), numel(Plist));
times = zeros(numel(windows), numel(Plist));
times2 = zeros(numel(windows), numel(Plist));
for w=1:numel(windows)
  opt.window = windows{w};
  for k=1:numel(Plist)
    opt.P = Plist(k);
    [errors(w,k), times(w,k), times2(w,k)] = compute_SE(x, f, opt, Neval, uref, n_rep, exclude_pre);
  end
end

%% Plotting
style = {'o-', '*-', 's-', '^-'};
lw = 1;

sfigure(1); clf; hold on
for w=1:numel(windows)
  plot(Plist, errors(w,:), style{w}, 'DisplayName', ...
    replace(windows{w}, '_', ' '), 'LineWidth', lw)
end
set(gca, 'yscale', 'log')
xlabel('P')
ylabel('RMS error (relative)')
legend('show', 'Location', 'SouthWest')
grid on
title('Error vs P')

sfigure(2); clf; hold on
for w=1:numel(windows)
  plot(Plist, times(w,:), style{w}, 'DisplayName', ...
    replace(windows{w}, '_', ' '), 'LineWidth', lw)
end
set(gca, 'yscale', 'log')
xlabel('P')
ylabel('Time [s]')
legend('show', 'Location', 'SouthWest')
grid on
title('Total time excluding precomp vs P')

sfigure(3); clf; hold on
for w=1:numel(windows)
  plot(errors(w,:), times(w,:), style{w}, 'DisplayName', ...
    replace(windows{w}, '_', ' '), 'LineWidth', lw)
end
set(gca, 'xscale', 'log')
xlabel('RMS error (relative)')
ylabel('Time [s]')
legend('show', 'Location', 'SouthWest')
grid on
title('Total time excluding precomp vs error')

sfigure(4); clf; hold on
for w=1:numel(windows)
  plot(Plist, times2(w,:), style{w}, 'DisplayName', ...
    replace(windows{w}, '_', ' '), 'LineWidth', lw)
end
set(gca, 'yscale', 'log')
xlabel('P')
ylabel('Time [s]')
legend('show', 'Location', 'SouthWest')
grid on
title('Gridding and gathering time vs P')

sfigure(5); clf; hold on
for w=1:numel(windows)
  plot(errors(w,:), times2(w,:), style{w}, 'DisplayName', ...
    replace(windows{w}, '_', ' '), 'LineWidth', lw)
end
set(gca, 'xscale', 'log')
xlabel('RMS error (relative)')
ylabel('Time [s]')
legend('show', 'Location', 'SouthWest')
grid on
title('Gridding and gathering time vs error')

% ------------------------------------------------------------------------------
function uref = compute_direct(x, f, opt, N)
  uref = SE3P_Laplace_direct_fd_mex(1:N, x, f, opt);
  % Or use the Kaiser window as reference solution

% ------------------------------------------------------------------------------
function [err, tim, tim2] = compute_SE(x, f, opt, N, ref, n_rep, exclude_pre)
  % ignore n_rep for now...
  [u, time] = SE3P_Laplace_fourier_space(1:N, x, f, opt);
  err = rms(u-ref) / rms(ref);
  tim = time.total - exclude_pre*(time.pre + time.prefft);
  tim2 = time.grid + time.int;
