%param selection
% Emanuel Ravemyr 19/11 -2020

%%setup
clear
rng(1);
N = 10000; % number of source particles


%% Parameter selection
L = 1; % box side length
box = [L L L]; % periodic box
opt.box = box;

%Ewald params


M0 = 128; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*[1,1,1];
opt.xi = 20;
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
	opt.P = 32;
	opt.M = 200*[1,1,1];
	ref = SE3P_Stokes(1:N,x,f,opt);

	save('refmat.mat','ref');
%	ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
end
	%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
rms_ref = rmse(ref);
%% Compare solutions with changing P

str = {};
e_vec = [];
ee_vec = [];
eee_vec = [];
F = sum(norm(f.^2));
A = @(a,b,c) sqrt(a)*b*(b*c)^0;
disp(num2str(A(F,opt.xi,opt.box(1))/rms_ref))

opt.M = 112*[1,1,1];
PP = [2:2:16];
rms_err = [];
time = zeros(2,length(PP));
for i = 1:2
    k=0;
    rms_err = [];
    if(i==1)
        opt.window = 'gaussian';
    else
        opt.window = 'kaiser_poly';
    end
    for P = PP
        k = k+1;
        opt.P = P;
        t = tic;
        for it = 1:60
            u = SE3P_Stokes(1:N, x, f, opt);
        end
        tt = toc(t);
        rms_err = [rms_err rmse(u-ref)];
        time(i,k) = tt/it;
    end
    semilogy(time,rms_err)
end

legend('Gaussian','Kaiser\_poly','Location','Best')
ylim([10^-14,1])
exportgraphics(gcf,'winf_comp.png')


