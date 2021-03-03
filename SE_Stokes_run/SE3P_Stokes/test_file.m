%SE Tests

N_vals = [100];
xi_vals = [15,20];
L_vals = [2,4];
tol_vals = [10^-6, 10^-8, 10^-10, 10^-12];
r = [];
tt = [];
for n = N_vals
    for x = xi_vals
        for L = L_vals
            for tol = tol_vals
                [res, par, SE,~ ] = test_function(n,x,L,tol);	
                assert(res<tol,'N =%g, xi = %g, L = %f, tol = %g' ,n,x,L,tol)
                r = [r, res];
                tt = [tt, SE];
            end
        end
    end
end
figure
semilogy(tt,r,'*')
exportgraphics(gcf,'testimg.png')
