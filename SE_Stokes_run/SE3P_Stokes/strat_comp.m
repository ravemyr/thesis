% Time comparison for strategies
N = 1000;
xi = 15;
L = 1;


opt.box = [L,L,L];
tolerance = 10.^[-4,-6,-8,-10];
opt.betaP = 2.5;
opt.xi = xi;
opt.L = L;
[x, f] = SE_charged_system(N,opt.box,'vector');
%opt.M = [1,1,1]*200;
%opt.P = 32;
%opt.window = 'kaiser_exact';
%ref = SE3P_Stokes(1:N,x,f,opt);
load('refmat.mat')
ref = refmat.ref;
opt.x =x;
opt.N = N;
opt.window = 'kaiser_poly';
times = zeros(2,length(tolerance));
errors = zeros(2,length(tolerance));
for i = 1:2
       disp(i)	
    k=0;
    for tol = tolerance
        opt.f = f;
	k = k+1;
        opt = param_select_stokes(tol, opt);
        if(i==1)
            opt.M = opt.M.*1.1;
            opt.P = opt.P+4;
        else
            opt.M = opt.M.*1.18;
            opt.P = opt.P+2;
        end
        t = tic;
        for it = 1:50
            u = SE3P_Stokes(1:N,x,f,opt);
        end
        tt = toc(t);
        errors(i,k) = rmse(u-ref)
        times(i,k) = tt/it;
    end
end
disp(times)
semilogy(times(1,:),errors(1,:))
hold on
semilogy(times(2,:),errors(2,:))
legend('Strategy 3','Strategy 4','Location','Best')
xlabel('time (s)')
exportgraphics(gcf, 'strat_comparison.png')
