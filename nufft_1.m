close all
%Non-uniform FFT

%Generate f and x
N = 256;
data = rand(1,N)*2*pi;
data = sort(data);
fun = @(x) sin(x);
fj = fun(data);

res_norm = [];
res_fft_slow = [];
%Dimension and tau
M = N/2;
Mr = 2*M;
t_vec = (2:2:24)/(M^2);
% Gaussian kernel - infinite sum not feasible numerically
q = -10000:10000;
l = length(q);
for t = t_vec
    gf = @(x) sum(exp(-((x*ones(1,l)-2*pi*q).^2)/(4*t)));

    %Compute the convolution f*g -> ft
    ft = zeros(1,Mr);
    for j = 0:Mr-1
        for jj = 1:N
            ft(j+1) = ft(j+1) + fj(jj)*gf(2*pi*j/Mr - data(jj));
        end 
    end

    %Perform fft on ft - Matlab built in function
    ff = fft(ft)/length(ft);
    fft_sum = zeros(1,Mr);
    for k = -Mr/2:Mr/2-1
        fft_sum(k+Mr/2+1) = sum(ft.*exp(-1i*k*2*pi*(0:Mr-1)/Mr))/Mr;
    end
    fft_sol = norm(ff-fft_sum);


    %Adjustment multiplication (negative exponential renders better results?)
    fff = sqrt(pi/t)*exp(t*(0:Mr-1).^2).*ff;
    fff_slow = sqrt(pi/t)*exp(t*(-Mr/2:Mr/2-1).^2).*fft_sum;

    %Comparative computation/Validation
    valid = zeros(1,Mr);
    for k = 1:Mr
        valid(k) = (1/N)*sum(fj.*exp(-1i*(k-1-Mr/2)*data));
    end

    %Compute 2-norm of error

    a = norm(fff-valid);
    res_norm = [res_norm a];
    a = norm(fff_slow-valid);
    res_fft_slow = [res_fft_slow a];
end
loglog(t_vec,res_fft_slow)
figure
loglog(t_vec,res_norm)