% Demo file for stokeslet computation

clear
rng(1);

N = 100; % number of source particles
L = 1; % box side length
box = [L L L]; % periodic box

% Ewald parameters


%%
close all

% grid
opt.box = box;

%Ewald params
M0 = 32; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*opt.box;
opt.xi = M0*pi/12;

% charge-neutral system
[x, f] = SE_charged_system(N,box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = box;
%%
% Direct computation for reference
t = tic();
ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
tDirReal = toc(t);

%%
windows = {'gaussian', 'kaiser_exact', 'kaiser_poly'};
for w=1:numel(windows)
  fprintf('3P Stokes Spectral Ewald, window: %s\n', windows{w});
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
  t = tic();
  u = SE3P_Stokes(1:N, x, f, opt);
  tSEFour = toc(t);
  
  rms_error = rmse(u-ref) / rmse(ref);
  fprintf('  RMS error: %.16g\n', rms_error);
  fprintf('  Time: %g s (Fourier), %g s (real), total %g s\n\n', ...
          tSEFour);
end
