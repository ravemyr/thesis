function [res, t_par, tSE, tDirReal] = test_function(N,xi,L,tol)
    opt.box = [L,L,L];    
    opt.window = 'kaiser_poly';
    [x, f] = SE_charged_system(N,opt.box,'vector');
    F = norm(f.^2);
    
    tolerance = tol;
    opt.N = N;
    opt.x =x;
    opt.f = f;
    opt.xi = xi;
    
    
    t = tic();
    opt = param_select_stokes(tolerance, opt);
    t_par = toc(t);
    
    t = tic();
    u = SE3P_Stokes(1:N, x, f, opt);
    tSE = toc(t);
    
    %%
    % Direct computation for reference
    opt.M = 280*[1,1,1];
    opt.P = 32;
    t = tic();

    ref = SE3P_Stokes(1:N,x,f,opt);
    tDirReal = toc(t);
    res = rmse(u-ref);
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
end
