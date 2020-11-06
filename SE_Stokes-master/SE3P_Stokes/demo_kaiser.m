% Stokeslet Spectral Ewald, basic accuracy/convergence computation
% Dag Lindbo, dag@kth.se

clear all,  close all

box = [1 1 1]; % domain
N = 10;        % numbver of charged particles
xi = 4;        % ewald parameter

P = 12; % support

% grid
SE_opt.M = 32*box;
SE_opt.box = box;
SE_opt.beta = 2.5;
SE_opt.window = 'kaiser';

% charge-neutral system
[x f] = SE_charged_system(N,box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (SE_opt.M(1)-1)/2;
ED_opt.xi = xi;
ED_opt.box = box;

% compute FD Ewald sum
ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);

SE_opt.P = P;
u = se3p_fourier_space(1:N,x,f,xi,SE_opt);

% compute RMS error (first vector component)
err = rms(u-ref)/rms(ref)
