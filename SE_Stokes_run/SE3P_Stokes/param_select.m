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


M0 = 128; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*opt.box;
opt.xi = M0*pi/12;
opt.P = opt.M/2;
opt.betaP = 2.5;

% charge-neutral system
[x, f] = SE_charged_system(N,box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = box;

%% Direct solution

ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
%% Estimate

F = sum(abs(f).^2);
est = @(M,xi,L,F) 2.*L^2.*sqrt(F).*(2*sqrt(pi)*M/2 + 3*xi*L).*exp(-((M/2)*pi/(xi*L)).^2)/sqrt(pi);
MM = [8,16,32];
err = [];
ee = [];
for i = MM
    M0 = i;
    opt.M = M0*opt.box;
    opt.P = i/2;
    ED_opt.layers = (opt.M(1)-1)/2;
    %u = SE3P_Stokes(1:N, x, f, opt);
    u = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
    err = [err rms(u-ref)];
    ee = [ee, est(opt.M, opt.xi,L,F)];
end
disp('error rate as M increases')
disp(err)
disp('error estimate')
disp(ee)
semilogy(MM./2,err([1,4,7]))
exportgraphics(gcf,'error_kplot.png')
%semilogy(MM./2,ee)
%exportgraphics(gcf,'error_est.png')



