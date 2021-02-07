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


M0 = 128; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*opt.box;
opt.xi = 30;
opt.betaP = 2.5;
%opt.window = 'kaiser_exact';
opt.window = 'kaiser_exact';
opt.polynomial_degree = 9;
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
MM = 32:4:64;
str = {};
est = @(M,xi,L,f) sqrt(f)*(xi^3*L^2/(pi^4*(M/2)^(3/2)))*exp(-(pi*(M/2)/(xi*L))^2);
e_vec = [];
F = sum(norm(f.^2));
for M = [112 MM]
	opt.M = M*[L,L,L];
	PP = [4:1:10,10:2:32];
	rms_err = [];
    
	for P = PP
    		opt.P = P;
    		u = SE3P_Stokes(1:N, x, f, opt);
    		rms_err = [rms_err rmse(u-ref)/rmse(ref)];
    end
    e = est(M,opt.xi,L,F);
    e_vec = [e_vec, e];
	semilogy(PP,rms_err)
   
	hold on
	semilogy(PP,ones(1,length(PP))*e)
	str = [str strcat('M = ',num2str(M))];
end

semilogy(PP,10*exp(-2.5.*PP),'--')
str = [str, 'estimate'];
legend(str)
e = est(48,opt.xi,L,F);
xlabel('P')
ylim([10^-14,1])
exportgraphics(gcf,'error_P.png')


