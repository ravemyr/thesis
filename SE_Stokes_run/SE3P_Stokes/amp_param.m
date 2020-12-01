%% Amplitute

A = @(F,xi,L) (F*xi^2/L).^(1/2); %Hypothesis

M0 = 128;
opt.M = [M0,M0,M0];
opt.P = 96;
opt.window = 'gaussian';
NN = [25,50,100];
LL = [1,2,3];
xx = [4,6,8];
MA = [];
fileid = fopen('outdata.txt','w');
fprintf(fileid,'%2s %3s %4s %9s %11s %5s \n','N','L','xi', 'rms','A','rms/a');
for L = LL
    opt.box = [L,L,L];
    for N = NN
        [x,f] = SE_charged_system(N,opt.box,'vector');
        F = sqrt(sum(abs(f).^2));
        for xi = xx
            opt.xi = xi;
            eu = rms(SE3P_Stokes(1:N, x, f, opt));
            a = A(F,xi,L);
            fprintf(fileid,'%3i %2.2f %1i %6.4f %6.4f %6.4f %6.4f %6.4f \n',N,L,xi,eu,norm(a),eu(2)/norm(a));
            MA = [MA eu/norm(a)];
        end
    end
end
plot(1:length(MA),MA,'*')
exportgraphics(gcf,'error_amp.png');
fclose(fileid);
disp(MA)
