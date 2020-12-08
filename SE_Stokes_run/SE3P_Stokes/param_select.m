%param selection
% Emanuel Ravemyr 19/11 -2020

%%setup
clear
rng(1);
N = 10; % number of source particles


%% Parameter selection
L = 1; % box side length
box = [L L L]; % periodic box
opt.box = box;

%Ewald params


M0=192; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*opt.box;
opt.xi = M0*pi/12;
opt.P = M0;
opt.betaP = 2.5;

% charge-neutral system
[x, f] = SE_charged_system(N,box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = box;

%% Direct solution
if(~exist('refval.mat'))
	disp('Generating reference solution')
	ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
	save('refval.mat','ref');
else
	disp('Using existing reference solution')
        reffile = load('refval.mat');
	ref = reffile.ref;
end
%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
%% Estimate

F = sum(norm(f).^2);
est = @(M,xi,L,f) 2*L^2*sqrt(f)*(2*sqrt(pi)*M/2 + 3*xi*L)*exp(-((M/2)*pi/(xi*L))^2)/sqrt(pi);
MM = [50,64,80,96,112];
err = [];
ee = [];

for i = MM
    M0 = i;
    opt.M = M0*opt.box;
    opt.P = 32;
    ED_opt.layers = (opt.M(1)-1)/2;
    %u = SE3P_Stokes(1:N, x, f, opt);
    u = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
    err = [err rmse(u-ref)/rmse(ref)];
    esti = est(opt.M(1),opt.xi,L,F)/rmse(ref);
    ee = [ee, esti];
end
disp(rmse(ref))
disp('error rate as M increases')
disp(err)
disp('estimate error comparison')
disp(ee)
semilogy(MM./2,err,'.-')
hold on
semilogy(MM./2,ee,'.--')
legend('relative error','error estimate','Location','North East')
exportgraphics(gcf,'error_kplot.png')
%semilogy(MM./2,ee)
%exportgraphics(gcf,'error_est.png','Resolution',1500)



