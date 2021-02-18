%param selection
% Emanuel Ravemyr 19/11 -2020

%%setup
clear
rng(1);
N = 10; % number of source particles


%% Parameter selection
L = 10; % box side length
box = [L L L]; % periodic box
opt.box = box;

%Ewald params


M0 = 128; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*opt.box;
opt.xi = 2;
opt.betaP = 2.5;
opt.c = sqrt(0.91);
%opt.window = 'kaiser_exact';
opt.window = 'gaussian';
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
	
	opt.P = 20;
	opt.M = 196*[1,1,1];
	ref = SE3P_Stokes(1:N,x,f,opt);
%	ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
end
	%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
rms_ref = rmse(ref);
%% Compare solutions with changing P
MM = 60:4:84;
str = {};
est = @(M,xi,L,f) sqrt(f)*(xi^3*L^2/(pi^4*(M/2)^(3/2)))*exp(-(pi*(M/2)/(xi*L))^2);
e_vec = [];
F = sum(norm(f.^2));
A = @(a,b,c) sqrt(a)*b*(b*c);
disp(num2str(A(F,opt.xi,opt.box(1))/rms_ref))
for M = [112 MM]
	
	opt.M = M*[1,1,1];
	PP = [2:1:10,10:2:32];
	rms_err = [];
    
	for P = PP
    		opt.P = P;
    		u = SE3P_Stokes(1:N, x, f, opt);
    		rms_err = [rms_err rmse(u-ref)/rms_ref];	
    end
    e = est(M,opt.xi,L,F);
    e_vec = [e_vec, e./rms_ref];
	semilogy(PP,rms_err)
   
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
    semilogy(PP,ones(1,length(PP))*e_vec(i),('bl-'))
end
str = [str,'estimate'];
legend(str)
xlabel('P')
xlim([1,32])
ylim([10^-14,1])
grid on
exportgraphics(gcf,'error_P.png')
xlim([1,15])
exportgraphics(gcf,'error_P_zoom.png')

