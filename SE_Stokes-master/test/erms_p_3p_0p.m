clear
rng(1);


addpath('../SE3P_Stokes/src');
addpath('../bin');
addpath('../util');
addpath('../SE_fast_gridding');

L = 1;
box = [L L L];
N = 10;
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
ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);

Pl = 4:2:24;
[err_3p_kaiser, err_3p_gaussian] = deal([]);
for P=Pl
  opt.P = P;
  opt.window = 'kaiser';
  uf = se3p_fourier_space(1:N,x,f,opt.xi,opt);

  % compute RMS error (first vector component)
  err_3p_kaiser(end+1) = rms(uf-ref)/rms(ref);

  opt.window = 'gaussian';
  uf = se3p_fourier_space(1:N,x,f,opt.xi,opt);

  % compute RMS error (first vector component)
  err_3p_gaussian(end+1) = rms(uf-ref)/rms(ref);
end

semilogy(Pl, err_3p_kaiser, 'b^-', Pl, err_3p_gaussian, 'ro-')


%% SE0P Stokes
addpath('../SE0P_Stokes/src');
% Direct
tic
u = stokeslet_direct(x, f, box);
toc
ur = stokeslet_real_space(x, f, opt);
us = -4*opt.xi/sqrt(pi)*f;
ref = u - ur - us;

% Precomp
Pl = 4:2:24;
[err_0p_kaiser, err_0p_gaussian] = deal([]);
for P=Pl
  opt.P = P;
  opt.window = 'kaiser';
  tic
  pre = stokeslet_precomp(opt);
  toc
  % Ewald
  tic
  uf = stokeslet_fourier_space(x, f, opt, pre);
  toc

  err_0p_kaiser(end+1) = rms(uf-ref) / rms(ref);

  opt.window = 'gaussian';
  tic
  pre = stokeslet_precomp(opt);
  toc
  tic
  uf = stokeslet_fourier_space(x, f, opt, pre);
  toc

  err_0p_gaussian(end+1) = rms(uf-ref) / rms(ref);
end

hold on
semilogy(Pl, err_0p_kaiser, 'k.-', Pl, err_0p_gaussian, 'r*-')
xlabel('P')
ylabel('RMS error')
legend('3P Kaiser', '3P Gaussian', '0P Kaiser', '0P Gaussian')
grid on
