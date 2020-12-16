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


M0 = 192; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*[1,1,1];
opt.xi = 15;
opt.betaP = 2.5;
opt.window = 'kaiser_exact';
% charge-neutral system
[x, f] = SE_charged_system(N,box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = box;

%% Direct solution
if(exist('refval.mat'))
	refval = load('refval.mat');
	ref = refval.refv(:,1:3);
else
	ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
end
	%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
%% Compare solutions with changing P
MM = 48:8:56;
str = {};
for M = [128 MM]
	opt.M = M*[L,L,L];
	PP = 4:2:32;
	rms_err = [];
	for P = PP
    		opt.P = P;
    		u = SE3P_Stokes(1:N, x, f, opt);
    		rms_err = [rms_err rmse(u-ref)/rmse(ref)];
	end
	semilogy(PP,rms_err)
	hold on
	str = [str strcat('M = ',num2str(M))];
end
legend(str)
exportgraphics(gcf,'error_P.png')



