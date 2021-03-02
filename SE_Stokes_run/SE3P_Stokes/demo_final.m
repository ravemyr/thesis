%demo final

N = 100000;
L = 2;
opt.window = 'kaiser_poly';
opt.box = [L,L,L];
[x, f] = SE_charged_system(N,opt.box,'vector');
F = norm(f.^2)
opt.xi = 15;
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
opt
fprintf('parameter selection time: %f (s) \n',t_par)
fprintf('computation time: %f (s) \n', tSE)
ED_opt.layers = (opt.M(1)-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = opt.box;
%%
% Direct computation for reference
t = tic();
opt.M = 196*[1,1,1];
opt.P = 32;
ref = SE3P_Stokes(1:N,x,f,opt);
tDirReal = toc(t);

format long;
fprintf('tolerance: %.16g \n',tolerance)
fprintf('Absolute error: %.16g \n',(rmse(u-ref)))
fprintf('Relative error: %.16g \n',rmse(u-ref)/rmse(ref))
fprintf('Approximate relative error: %.16g \n',rmse(u-ref)/(sqrt(F)*opt.xi))
semilogy(opt.xi, tolerance,'*')
hold on
semilogy(opt.xi, rmse(u-ref),'or')
semilogy(opt.xi, rmse(u-ref)/rmse(ref),'ob')
exportgraphics(gcf,'results.png')
