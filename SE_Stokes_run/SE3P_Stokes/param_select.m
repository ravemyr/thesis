%param selection
% Emanuel Ravemyr 19/11 -2020

%%setup
clear
rng(1);
N = 20; % number of source particles


%% Parameter selection
L = 4; % box side length
box = [L L L]; % periodic box
opt.box = box;

%Ewald params


M0=128; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*[1,1,1];
opt.P = M0;
opt.betaP = 2.5;

% charge-neutral system
[x, f] = SE_charged_system(N,box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.box = box;
refv = [];
rms_ev = [];
xx = (4:4:12)*pi;

    %% Direct solution
if(exist('refval.mat'))
        reffile = load('refval.mat');
    if(size(xx)==size(reffile.xx))
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
    %% Estimate

F = sum(norm(f).^2);
%est = @(M,xi,L,f) 2*L^2*sqrt(f)*(2*sqrt(pi)*M/2 + 3*xi*L)*exp(-((M/2)*pi/(xi*L))^2)/sqrt(pi);
est = @(M,xi,L,f) sqrt(f)*(xi*(xi*L)^2)/(pi^4*(M/2)^(2))*exp(-(pi*(M/2)/(xi*L))^2);
est2 = @(M,xi,L,f) sqrt(f)*(4*exp(-(pi*M/(2*xi*L))^2)/(L*pi*3^(1/2)))*( (M>=2)*(xi*L)^2/(2*pi*M/2));
MM = [8,24,32,48,64,72,80,88,96,104];
str = {};
opt.window = 'kaiser_exact';

for xi = xx
    err = [];
    ee = [];
    eee = [];
    idx = (find(xx==xi));
    disp(idx)
    ref = refv(:,3*idx-2:3*idx);
    rms_e = rms_ev(idx);
    for i = MM
        M0 = i;
        opt.M = M0*[1,1,1];
        opt.P = M0/2;
        ED_opt.layers = (opt.M(1)-1)/2;
        ED_opt.xi = xi;
        %u = SE3P_Stokes(1:N, x, f, opt);
        u = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
        err = [err rmse(u-ref)/rms_e];
        esti = est(opt.M(1),xi,L,F)/rms_e;
        ee = [ee, esti];
        esti2 = est2(opt.M(1),xi,L,F)/rms_e;
        eee = [eee esti2];
    end
    semilogy(MM./2,err,'.-')
    hold on
    %semilogy(MM./2,ee,'--o')
    semilogy(MM./2,eee,'--*')
    str = [str,strcat('computed error, \xi =', num2str(xi))];
    %str = [str, strcat('2D-estimate, \xi =', num2str(xi))];
    str = [str, strcat('3P-estimate, \xi =', num2str(xi))];
end
ylim([10^-8,1])
xlim([0,73])
legend(str{:},'Location','Best','FontSize',14)
xlabel('k_{\infty}')

% disp('error rate as M increases')
% disp(err)
% disp('estimate error comparison')
% disp(ee)
% 
% legend('relative error','error estimate','Location','North East')
exportgraphics(gcf,'error_kplot.png')
%semilogy(MM./2,ee)
%exportgraphics(gcf,'error_est.png','Resolution',1500)



