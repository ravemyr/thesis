%Timing 

%%setup
clear
rng(1);
NN = [10,20,50,100,200,400]; % number of source particles


%% Parameter selection
L = 1; % box side length
box = [L L L]; % periodic box
opt.box = box;

%Ewald params

xi = 10;
M0=144; % Set M0=M/L, the restu * 1+ is automatic
opt.M = M0*[1,1,1];
opt.P = 32;
opt.betaP = 2.5;
opt.xi = xi;
opt.P = 32;
opt.window = 'kaiser_poly';
timings = [];
for N = NN
    [x, f] = SE_charged_system(N,box,'vector');
    t = tic;
    u = SE3P_Stokes(1:N, x, f, opt);
    timings = [timings toc(t)];
end
plot(NN,timings)
