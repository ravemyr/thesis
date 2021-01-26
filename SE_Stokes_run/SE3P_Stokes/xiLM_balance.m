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


M0 = 144; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*opt.box;
opt.xi = 20;
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
%ref = SE3P_Stokes(1:N,x,f,opt);
%% Compare solutions with changing P
MM = 48:8:56;
str = {};
est = @(M,xi,L,f) sqrt(f)*(xi^3*L^2/(pi^4*(M/2)^(3/2)))*exp(-(pi*(M/2)/(xi*L))^2);
e_vec = [];

xx = [5,8,10,20];
LL = [2,1];
for xi = xx
    opt.xi = xi;
    for L = LL
	    opt.box = [L,L,L];
ED_opt.xi = xi;
ED_opt.box = [L,L,L];
	[x,f] = SE_charged_system(N,opt.box,'vector');
        ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
        F = sum(norm(f.^2));	
        for M = [96 MM]
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
    end
end
str = [str, 'estimate'];
%legend(str)
xlabel('P')
exportgraphics(gcf,'xiLM.png')



