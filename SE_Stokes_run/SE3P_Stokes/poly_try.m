N = 1000;
opt.box = [1,1,1];    
opt.window = 'kaiser_poly';
[x, f] = SE_charged_system(N,opt.box,'vector');


tolerance = 10^-10;
opt.f = f;
opt.xi = 20;

t = tic();
%varg = {tolerance, opt, 'Relative'};
%varg{1} = tolerance;
%varg{2} = opt;
%varg{3} = 'Relative';
opt = param_select_stokes(tolerance, opt);
t_par = toc(t);
opt.polynomial_degree = 9;
opt.P = 16;
t = tic();
u = SE3P_Stokes(1:N, x, f, opt);
tSE = toc(t);
disp(opt)
%%
% Direct computation for reference
opt.M = 350*[1,1,1];
opt.P = 32;
t = tic();
opt.window = 'kaiser_exact';
ref = SE3P_Stokes(1:N,x,f,opt);
tDirReal = toc(t);
res = rmse(u-ref);
format long;
fprintf('tolerance: %.16g \n',tolerance)
fprintf('Absolute error: %.16g \n',(rmse(u-ref)))
fprintf('Relative error: %.16g \n',rmse(u-ref)/rmse(ref))
fprintf('Approximate relative error: %.16g \n',rmse(u-ref)/(sqrt(F)*opt.xi))
