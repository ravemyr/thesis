%demo final

N = 10;
L = 1;
opt.box = [L,L,L];
[x, f] = SE_charged_system(N,opt.box,'vector');

opt.xi = 20;
tolerance = 10^(-10);
opt.N = N;
opt.x =x;
opt.f = f;


t = tic();
opt = param_select_stokes(tolerance, opt);
t_par = toc(t);
t = tic();
u = SE3P_Stokes(1:N, x, f, opt);
tSE = toc(t);
fprintf('parameter selection time: %f \n',t_par)
fprintf('computation time: ', tSE)
ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
opt
fprintf(tolerance>rmse(u-ref)/rmse(ref))