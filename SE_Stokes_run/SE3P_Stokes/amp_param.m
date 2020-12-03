%% Amplitute
rng(1);
A = @(F,xi,L) (F*xi/(L)).^(1/2); %Hypothesis

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



FF = [f 2*f 4*f 8*f];


for L = LL
    opt.box = [L,L,L];
    [x,f] = SE_charged_system(N,opt.box,'vector');
    FF = [f, 2*f, 4*f, 8*f];
    for xi = xx
        opt.xi = xi;
        for f = FF
            eu = rmse(SE3P_Stokes(1:N, x, f, opt));
            F = sqrt(sum(norm(f).^2));
            a = A(F,xi,L);
            fprintf(fileid,'%3i %2.2f %1i %6.4f %6.4f %6.4f \n',N,L,xi,eu,norm(a),eu/a);
            MA = [MA eu];
        end
    end
end
plot(1:length(MA),MA,'*')
exportgraphics(gcf,'error_amp.png','Resolution',1500);
fclose(fileid);
%disp(MA)
