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
opt.P = M0;
opt.betaP = 2.5;

% charge-neutral system
[x, f] = SE_charged_system(N,box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.box = box;
refv = [];
rms_ev = [];
xx = (8:4:12)*pi;

    %% Direct solution
    if(exist('refval.mat'))
            reffile = load('refval.mat');
	    if(xx==reffile.xx)
		disp('Using existing reference solution')
		refv = reffile.refv;
		rms_ev = reffile.rms_ev;
	    else
                disp('Generating reference solution')
                for xi = xx
                    ED_opt.xi = xi;
		    u = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
		    rms_e = rmse(u);
		    rms_ev =[rms_ev, rms_e];
                    refv = [refv,u];
                    save('refval.mat','refv','xx','rms_ev');
                end
	    end
	    
    else
        disp('Generating reference solution')
        for xi = xx
            ED_opt.xi = xi;
	    u = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
	    rms_e = rmse(u);
	    rms_ev = [rms_ev,rms_e];
            refv = [refv,u];
            save('refval.mat','refv','xx','rms_ev');
        end
    end
    %opt.window = 'gaussian';
    %ref = SE3P_Stokes(1:N,x,f,opt);
    %% Estimate

    F = sum(norm(f).^2);
    est = @(M,xi,L,f) 2*L^2*sqrt(f)*(2*sqrt(pi)*M/2 + 3*xi*L)*exp(-((M/2)*pi/(xi*L))^2)/sqrt(pi);
    est2 = @(M,xi,L,f) sqrt(f)*(xi^3*L^2/(pi^4*(M/2)^(3/2)))*exp(-(pi*(M/2)/(xi*L))^2);
    MM = [50,64,80,96,112,128,146];
str = {};

for xi = xx
    err = [];
    ee = [];
    eee = [];
    idx = (find(xx==xi));
    ref = refv(:,3*idx-2:3*idx);
    rms_e = rms_ev(idx);
    for i = MM
        M0 = i;
        opt.M = M0*opt.box;
        opt.P = M0/2;
        ED_opt.layers = (opt.M(1)-1)/2;
        ED_opt.xi = xi;
	%u = SE3P_Stokes(1:N, x, f, opt);
        u = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
        err = [err rmse(u-ref)/rmse(ref)];
        esti = est(opt.M(1),xi,L,F)/rms_e;
        ee = [ee, esti];
	esti2 = est2(opt.M(1),xi,L,F)/rms_e;
	eee = [eee esti2];
    end
semilogy(MM./2,err,'.-')
hold on
semilogy(MM./2,ee,'.--')
semilogy(MM./2,eee,':*')
str = [str,strcat('computed error, \xi =', num2str(xi))];
str = [str, strcat('error estimate, \xi =', num2str(xi))];
str = [str, strcat('second estimate, \xi =', num2str(xi))];
end
axis([1,MM(end)/2,10^-15,100])
legend(str{:},'Location','North East')


% disp('error rate as M increases')
% disp(err)
% disp('estimate error comparison')
% disp(ee)
% 
% legend('relative error','error estimate','Location','North East')
exportgraphics(gcf,'error_kplot.png')
%semilogy(MM./2,ee)
%exportgraphics(gcf,'error_est.png','Resolution',1500)



