clear
rng(1);

addpath('../SE3P_Stokes/src');
addpath('../bin');
addpath('../util');
addpath('../SE_fast_gridding');
warning('off');

walltime = struct('g3p',[],'k3p',[],'g0p',[],'k0p',[]);

L = 1;
box = [L L L];
N = 10000;
[x, f] = vector_system(N, box, 3);

% Maxmimize point-to-point distance
x(1,:) = 0;
x(2,:) = box;

M0 = 32; % Set M0, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0/12;
opt.oversampling = 1 + sqrt(3);
opt.rc = 6 / opt.xi;
opt.box = box;
opt.beta = 2.3;

%% SE3P Stokes
% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = box;

% compute FD Ewald sum
%ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
opt.P = 18;
opt.window = 'kaiser';
[ref, time] = se3p_fourier_space(1:N,x,f,opt.xi,opt);

Pl = 4:2:24;
[err_3p_kaiser, err_3p_gaussian] = deal([]);
for P=Pl
  opt.P = P;
  opt.window = 'kaiser';

  [uf, time] = se3p_fourier_space(1:N,x,f,opt.xi,opt);
  walltime.k3p(end+1) = time.grid + time.fft + time.scale + time.int;

  % compute RMS error (first vector component)
  err_3p_kaiser(end+1) = rms(uf-ref)/rms(ref);

  opt.window = 'gaussian';
  [uf, time] = se3p_fourier_space(1:N,x,f,opt.xi,opt);
  walltime.g3p(end+1) = time.grid + time.fft + time.scale + time.int;

  % compute RMS error (first vector component)
  err_3p_gaussian(end+1) = rms(uf-ref)/rms(ref);
end

loglog(err_3p_kaiser(Pl<=18), walltime.k3p(Pl<=18), 'b^-', err_3p_gaussian, walltime.g3p, 'ro-')


%% SE0P Stokes
addpath('../SE0P_Stokes/src');
% Direct
%u = stokeslet_direct(x, f, box);
%ur = stokeslet_real_space(x, f, opt);
%us = -4*opt.xi/sqrt(pi)*f;
%ref = u - ur - us;
opt.P = 18;
opt.window = 'kaiser';
pre = stokeslet_precomp(opt);
[ref, time] = stokeslet_fourier_space(x, f, opt, pre);

% Precomp
Pl = 4:2:24;
[err_0p_kaiser, err_0p_gaussian] = deal([]);
for P=Pl
  opt.P = P;
  opt.window = 'kaiser';
  pre = stokeslet_precomp(opt);
  % Ewald
  [uf, time] = stokeslet_fourier_space(x, f, opt, pre);
  walltime.k0p(end+1) = time.grid + time.fft + time.scale + time.int;

  err_0p_kaiser(end+1) = rms(uf-ref) / rms(ref);

  opt.window = 'gaussian';
  pre = stokeslet_precomp(opt);
  [uf, time] = stokeslet_fourier_space(x, f, opt, pre);
  walltime.g0p(end+1) = time.grid + time.fft + time.scale + time.int;

  err_0p_gaussian(end+1) = rms(uf-ref) / rms(ref);
end

hold on
loglog(err_0p_kaiser(Pl<=18), walltime.k0p(Pl<=18), 'k.-', err_0p_gaussian, walltime.g0p, 'r*-')
xlabel('RMS error')
ylabel('Time [s]')
legend('3P Kaiser', '3P Gaussian', '0P Kaiser', '0P Gaussian')
grid on
