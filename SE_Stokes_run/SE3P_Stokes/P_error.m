%param selection
% Emanuel Ravemyr 19/11 -2020

%%setup
clear
rng(1);
N = 1000; % number of source particles


%% Parameter selection
L = 2; % box side length
box = [L L L]; % periodic box
opt.box = box;

%Ewald params


M0 = 128; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*[1,1,1];
opt.xi = 10;
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
	opt.M = 300*[1,1,1];
	ref = SE3P_Stokes(1:N,x,f,opt);
%	ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
end
	%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
rms_ref = rmse(ref);
%% Compare solutions with changing P
MM = 4.*(8:4:20);
str = {};
est_old = @(M,xi,L,f) sqrt(f)*(xi^3*L^2/(pi^4*(M/2)^(3/2)))*exp(-(pi*(M/2)/(xi*L))^2);
%est = @(M,xi,L,f) sqrt(f)*(4/(3^(1/4)*L*pi))*exp(-(pi*M/(xi*L*2))^2);
est = @(M,xi,L,F) 2*sqrt(F)*(xi*L)^2*exp(-(pi*(M/2)/(xi*L))^2)/(sqrt(3)*pi^2*L*(M/2));
est_small = @(M,xi,L,F) 4*sqrt(F)/(sqrt(3)*pi*L)*exp(-(pi*M/(2*xi*L))^2);
e_vec = [];
ee_vec = [];
eee_vec = [];
F = sum(norm(f.^2));
A = @(a,b,c) sqrt(a)*b*(b*c)^0;
disp(num2str(A(F,opt.xi,opt.box(1))/rms_ref))
for M = [112 MM]
	
	opt.M = M*[1,1,1];
	PP = [2:1:10,10:2:32];
	rms_err = [];
    
	for P = PP
    		opt.P = P;
    		u = SE3P_Stokes(1:N, x, f, opt);
    		rms_err = [rms_err rmse(u-ref)];	
         end
    e = est(M,opt.xi,L,F);
    e_vec = [e_vec, e];
	semilogy(PP,rms_err)
    ee = est_old(M,opt.xi,L,F);
    ee_vec = [ee_vec ee];
    eee = est_small(M,opt.xi,L,F);
    eee_vec = [eee_vec eee];
	hold on
	str = [str strcat('M = ',num2str(M))];
end

if(strcmp(opt.window,'kaiser_exact'))
semilogy(PP,10*exp(-2.5.*PP),'--')
str = [str ,'10exp(-2.5\beta)'];
else
semilogy(PP,exp(-(pi/2)*PP*opt.c),'--')
str = [str,'exp(-(\pi/2)Pc)'];
end
for i = 1:length(e_vec)
    semilogy(PP,ones(1,length(PP))*e_vec(i),('bl--'))
%    semilogy(PP,ones(1,length(PP))*ee_vec(i),'r--')
    semilogy(PP,ones(1,length(PP))*eee_vec(i),'g--')
end
str = [str,'estimate'];
opt
legend(str)
xlabel('P')
xlim([1,32])
ylim([10^-14,1])
grid on
exportgraphics(gcf,'error_P.png')
xlim([1,15])
exportgraphics(gcf,'error_P_zoom.png')

