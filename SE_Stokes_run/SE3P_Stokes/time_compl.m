%Timing 

%%setup
clear
rng(1);
NN = [10,50,250,1250,6250]; % number of source particles


%% Parameter selection
L = 1;
opt.window = 'kaiser_poly';
opt.box = [L,L,L];
opt.xi = 5;
tolerance = 10^(-10);
timings = [];

for N = NN
    opt.N = N;
    opt.box = 5^(1/3)*opt.box;
    [x, f] = SE_charged_system(N,opt.box,'vector');    
    opt.x =x;
    opt.f = f;
    opt = param_select_stokes(tolerance, opt);
    tic; 
    SE3P_Stokes(1:N,x,f,opt);
    toc;
    t = tic;
    for it = 1:20
    	u = SE3P_Stokes(1:N, x, f, opt);
    end
    ti = toc(t);
    timings = [timings ti/it];
end
loglog(NN,timings)
hold on
loglog(NN,NN.*log(NN),'--')

exportgraphics(gcf,'complexity.png')
disp(timings)

%%
L = 1;
opt.window = 'kaiser_poly';
opt.box = [L,L,L];
opt.xi = 5;
tolerance = 10^(-10);
timings = [];
N = 100;
opt.N = N;
[x, f] = SE_charged_system(N,opt.box,'vector');    
opt.x =x;
opt.f = f;
times = [];
opt.betaP = 2.5;
MM = 10:10:200;
for M = MM
opt.M = M*[1,1,1];
opt.window = 'kaiser_exact';
opt.polynomial_degree = 9;
opt.P = 32;
t = tic;
for it = 1:20
    u = SE3P_Stokes(1:N,x,f,opt);
end
tt = toc(t);
times = [times tt];
end
disp(times)
figure
loglog(MM,times,'*')
hold on
loglog(MM,MM.^3.*log(MM.^3))
exportgraphics(gcf,'MP_complexity.png')
disp(timings)
