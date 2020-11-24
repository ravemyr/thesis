%param selection
% Emanuel Ravemyr 19/11 -2020

%%setup
clear
rng(1);
N = 100; % number of source particles


%% Parameter selection
L = 1; % box side length
box = [L L L]; % periodic box
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

%% Direct solution

ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
%% Estimate

F = sum(abs(f).^2);
est = @(M,xi) 2*L^2*sqrt(F)*(2*sqrt(pi)*M/2 - 3*xi*L)*exp(-(pi/(xi*L)*M/2)^2)/sqrt(pi);
A = @(F,xi,L) (F*xi/L)^(1/2);
MM = [16,24,32,40,48,56,64];
err = zero(1,len(MM));
for i = MM
    M0 = i;
    opt.M = M0*opt.box;
    opt.xi = M0*pi/12;
    u = SE3P_Stokes(1:N, x, f, opt);
    rms_error = rms(u-ref) / rms(ref);
end