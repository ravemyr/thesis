%rms comparison with A
clear
rng(1);
N = 10; % number of source particles


%% Parameter selection
L = 2; % box side length
box = [L L L]; % periodic box
opt.box = box;

%Ewald params


M0 = 128; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*[1,1,1];
opt.xi = 25;
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
str = {};
M = 96; 
opt.M = M*[L,L,L];
PP = 4:2:32;
rms_err = [];
A = @(F,xi,L) F^(1/2)*(xi)^(1)*(L^(0));
F = sum(norm(f).^2);
a = A(F,opt.xi,L);
for N = [10,50,100,200]
for L = [1,2,3]
	opt.box = [L,L,L];
	opt.M = M*opt.box;
	[x,f] = SE_charged_system(N,box,'vector')
	F = sum(norm(f).^2);
	for xi = 4:2:16
	for P = PP
    		opt.P = P;
    		u = SE3P_Stokes(1:N, x, f, opt);
    		rms_err = [rms_err rmse(u-ref)/a];
	end
	semilogy(PP,rms_err,'.b')
	hold on
end
end
end
exportgraphics(gcf,'rms_a.png')
