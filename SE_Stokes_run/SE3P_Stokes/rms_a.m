%rms comparison with A
clear
rng(1);

opt.betaP = 2.5;
opt.window = 'kaiser_exact';

% parameters for (reference) direct Ewald sum


%% Direct solution
	%opt.window = 'gaussian';
%ref = SE3P_Stokes(1:N,x,f,opt);
%% Compare solutions with changing P
str = {};
M = 96; 

PP = 4:2:32;

A = @(F,xi,L) F^(1/2)*(xi)^(1)*(L^(0));


for N = [10,50,100]
    for L = [1,2,3]
        opt.box = [L,L,L];
        opt.M = M*opt.box;
        [x,f] = SE_charged_system(N,opt.box,'vector');
        F = sum(norm(f).^2);
        for xi = 4:2:16
	    rms_err = [];
	    a = A(F,xi,L);
            opt.xi = xi;
            ED_opt.layers = (opt.M(1)-1)/2;
            ED_opt.xi = opt.xi;
            ED_opt.box = opt.box;
            ref = SE3P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);
            for P = PP
                opt.P = P;
                u = SE3P_Stokes(1:N, x, f, opt);
                rms_err = [rms_err rmse(u-ref)/a];
	    end
         semilogy(PP,rms_err,'.b')
         hold on
            
        end
	exportgraphics(gcf,'rms_a.png')
    end
end

exportgraphics(gcf,'rms_a.png')
