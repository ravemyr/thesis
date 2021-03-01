%Timing 

%%setup
clear
rng(1);
NN = [10,100,400,1000,5000,10000,50000]; % number of source particles


%% Parameter selection
L = 1;
opt.window = 'kaiser_poly';
opt.box = [L,L,L];
opt.xi = 10;
tolerance = 10^(-10);





timings = [];

for N = NN
    opt.N = N;
    [x, f] = SE_charged_system(N,opt.box,'vector');    
    opt.x =x;
    opt.f = f;
    opt = param_select_stokes(tolerance, opt);
    t = tic;
    for it = 1:10
    	u = SE3P_Stokes(1:N, x, f, opt);
    end
    ti = toc(t);
    timings = [timings ti/it];
end
loglog(NN,timings)
hold on
loglog(NN,NN.*log(NN),'--')
exportgraphics(gcf,'complexity.png')
