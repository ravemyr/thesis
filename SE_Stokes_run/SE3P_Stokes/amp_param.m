%% Amplitute
rng(1);
%A = @(F,xi,L) F*(xi^(3/4)/(L^(9/16))); %Hypothesis for new generated f
A = @(F,xi,L) F^(1/2)*(xi)^(1)*(L^(0));

M0 = 128;
opt.M = [M0,M0,M0];
opt.P = 32;
opt.window = 'gaussian';
N = 100;
LL = [1,2,2.5,3,4,5,6,7];
xx = [8,10,12,14,18,20,22,25,30];
MA = [];
%opt.box = [1,1,1];
%[x,f] =SE_charged_system(N,opt.box,'vector');
fileid = fopen('outdata.txt','w');
fprintf(fileid,'%2s %3s %4s %9s %6s %5s \n','N','L','xi', 'rms','A','rms/a');
for L = LL
    opt.box = [L,L,L];
    [x,f] = SE_charged_system(N,opt.box,'vector');
    for xi = xx
        opt.xi = xi;
        for j = 0:2
	    ff = f*2^j;
            eu = rmse(SE3P_Stokes(1:N, x, ff, opt));
            F = sum(norm(ff).^2);
            a = A(F,xi,L);
            fprintf(fileid,'%3i %2.2f %1i %6.4f %6.4f %6.4f \n',N,L,xi,eu,norm(a),eu/a);
            MA = [MA eu/a];
        end
    end
end
plot(1:length(MA),MA,'*')
axis([0,200,0.2,0.5]);
exportgraphics(gcf,'error_amp.png');
fclose(fileid);
%disp(MA)
