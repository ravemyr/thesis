% Time comparison for strategies
N = 1000;
xi = 15;
L = 1;

opt.window = 'kaiser_poly';
opt.box = [L,L,L];
tolerance = 10.^[-4.-6,-8,-10];
timings = [];
opt.N = N;
opt.xi = xi;
opt.L = L;
[x, f] = SE_charged_system(N,opt.box,'vector');    
opt.x =x;
opt.f = f;
times = zeroes(2,length(tolerance));
opt.betaP = 2.5;
for i = 1:2 
    k=0;
    for tol = tolerance
        k = k+1;
        opt = param_select_stokes(tol, opt);
        if(i==1)
            M = M*1.1;
            P = P+4;
        else
            M = M*1.18;
            P = P+2;
        end
        t = tic;
        for it = 1:20
            u = SE3P_Stokes(1:N,x,f,opt);
        end
        tt = toc(t);
        times(i,k);
    end
end
disp(times)
semilogy(
