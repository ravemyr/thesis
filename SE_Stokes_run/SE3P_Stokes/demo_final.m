%demo final

N = 100;
L = 1;
opt.window = 'kaiser_poly';
opt.box = [L,L,L];
[x, f] = SE_charged_system(N,opt.box,'vector');

opt.xi = 10;
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
fprintf('parameter selection time: %f (s) \n',t_par)
fprintf('computation time: %f (s) \n', tSE)
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = opt.box;
%%
% Direct computation for reference
t = tic();
ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
tDirReal = toc(t);
opt
format long;
fprintf('tolerance: %.16g \n',tolerance)
fprintf('%.16g \n',(rmse(u-ref)))
fprintf('%.16g \n',rmse(u-ref)/rmse(ref))
plot(opt.xi, tolerance,'*')
hold on
plot(opt.xi, rmse(u-ref),'o')
exportgraphics(gcf,'results.png')
