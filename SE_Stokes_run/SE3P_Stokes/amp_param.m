%% Amplitute
rng(1);
A = @(F,xi,L) (F^2*xi^2/(L)).^(1/2); %Hypothesis

M0 = 128;
opt.M = [M0,M0,M0];
opt.P = 96;
opt.window = 'gaussian';
N = 100;
LL = [0.5,1,1.5,2];
xx = [2,4,6,8,10,12];
MA = [];
fileid = fopen('outdata.txt','w');
fprintf(fileid,'%2s %3s %4s %9s %6s %5s \n','N','L','xi', 'rms','A','rms/a');
for L = LL
    opt.box = [L,L,L];
    [x,f] = SE_charged_system(N,opt.box,'vector');
    for xi = xx
        opt.xi = xi;
        for j = 0:3
	    ff = f*2^j;
            eu = rmse(SE3P_Stokes(1:N, x, ff, opt));
            F = sqrt(sum(norm(ff).^2));
            a = A(F,xi,L);
            fprintf(fileid,'%3i %2.2f %1i %6.4f %6.4f %6.4f \n',N,L,xi,eu,norm(a),eu/a);
            MA = [MA eu];
        end
    end
end
plot(1:length(MA),MA,'*')
exportgraphics(gcf,'error_amp.png');
fclose(fileid);
%disp(MA)
