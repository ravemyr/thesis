%SE Tests

N_vals = [100,1000];
xi_vals = [5,10];
L_vals = [1,2];
tol_vals = [10^-6, 10^-8];
r = [];
tt = [];
tols = [];
for n = N_vals
    for x = xi_vals
        for L = L_vals
            for tol = tol_vals
                [res, par, SE,~ ] = test_function(n,x,L,tol);	
                assert(res<tol,'N =%g, xi = %g, L = %f, tol = %g' ,n,x,L,tol)
                r = [r, res];
                tt = [tt, SE];
                tols = [tols, tol];
            end
        end
    end
end
figure
ll = length(N_vals)*length(xi_vals)*length(L_vals)*length(tol_vals);
semilogy(1:ll,r,'*')
hold on
semilogy(1:ll,tols,'o')
exportgraphics(gcf,'testimg.png')
