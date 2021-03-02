%SE Tests

N_vals = [10,100,1000,10000];
xi_vals = [5,10,15,20,30];
L_vals = [1,2,4,5,10];
tol_vals = [10^-4, 10^-6, 10^-8, 10^-10, 10^-12];

for n = N_vals
    for x = xi_vals
        for L = L_vals
            for tol = tol_vals
                t = test_function(N,x,L,tol);
                assert(t(1)<tol)
            end
        end
    end
end