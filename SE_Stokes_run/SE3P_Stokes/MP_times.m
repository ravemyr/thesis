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
opt.xi = 5;
opt.betaP = 2.5;
opt.c = sqrt(0.91);
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
	opt.P = 25;
	opt.M = 256*[1,1,1];
	ref = SE3P_Stokes(1:N,x,f,opt);
%	ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
end
	%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
rms_ref = rmse(ref);
%% Compare solutions with changing P
MM = 8:4:24;
str = {};
%est = @(M,xi,L,f) sqrt(f)*(xi^3*L^2/(pi^4*(M/2)^(3/2)))*exp(-(pi*(M/2)/(xi*L))^2);
est = @(M,xi,L,f) sqrt(f)*(4/(3^(1/4)*L*pi))*exp(-(pi*M/(xi*L*2))^2);
e_vec = [];
F = sum(norm(f.^2));
A = @(a,b,c) sqrt(a)*b*(b*c);
PP = [2:1:10,10:2:32];
times = zeros(length(MM),length(PP));
for M = MM
    Midx = (find(M==MM));
	opt.M = M*[1,1,1];
	rms_err = [];
	for P = PP
            Pidx = (find(P==PP));
    		opt.P = P;
            t = tic;
    		u = SE3P_Stokes(1:N, x, f, opt);
		if(Pidx>1)
            times(Midx,Pidx) = toc(t) + sum(times(Midx,1:Pidx-1));
    	    else
		    times(Midx,Pidx) = toc(t);
	    end	    
	    rms_err = [rms_err rmse(u-ref)/rms_ref];	
    end
    e = est(M,opt.xi,L,F);
    e_vec = [e_vec, e];
	semilogy(times(Midx,:),rms_err,'*')  
	hold on
	str = [str strcat('M = ',num2str(M))];
end
disp(times(1,:))
%legend(str)
xlabel('time (s)')
ylim([10^-14,1])
exportgraphics(gcf,'P-time-error.png')


