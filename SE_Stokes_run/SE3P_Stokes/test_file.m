%SE Tests

N_vals = [100,1000,10000];
xi_vals = [15,20,30];
L_vals = [1,2,4,5];
tol_vals = [10^-6, 10^-8, 10^-10, 10^-12];
r = [];
tt = [];
for n = N_vals
    for x = xi_vals
        for L = L_vals
            for tol = tol_vals
                t = test_function(n,x,L,tol);	
                assert(t(1)<tol,'N =%g, xi = %g, L = %f, tol = %g' ,n,x,L,tol)
                r = [r, t(1)];
                tt = [tt, t(3)];
            end
        end
    end
end
semilogy(tt,r,'*')
exportgraphics(gcf,'testimg.png')
