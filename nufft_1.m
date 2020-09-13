close all
%Non-uniform FFT

%Preamble, generate f and x
N = 1024;
rng(0)
data = rand(1,N)*2*pi;
data = sort(data);
fun = @(x) x+x.^2;
fj = fun(data);
res_norm = [];

%Dimension and tau
M = 4*N;
Mr = 2*M;
t_vec = [1,6,12]/(M^2);

% Gaussian kernel - infinite sum not feasible numerically
% q = -40000:40000;
% l = length(q);
%gf = @(x,p) sum(exp((-1/(4*p))*((x'*ones(1,l)-2*pi*ones(length(x),1)*q).^2)),2)';

%Comparative computation/Validation
valid = zeros(1,M);
for k = 1:M
    valid(k) = sum(fj.*exp(-1i*(k-1-M/2)*data))/N;
end

%Comppute NUFFT for different t, given fj, x and Mr
for t = t_vec
    
    %Compute the convolution f*g -> ft
    ft = zeros(1,Mr);
    for j = 0:Mr-1
        ft(j+1) = sum(fj.*gtau(2*pi*j/Mr - data,t));
    end

    %Perform fft on ft - Matlab built in function with k-shift
    ff = fft(ft.*exp(-1i*pi*2*(0:Mr-1)*(-Mr/2)/Mr))/Mr;
    
    %Matlab gives F(k) where 0<=k<N-1, adjust size due to oversampling
    p = Mr/M; %Assumes Mr>M
    ff = ff(1+(p-1)*M/2:(p+1)*M/2);
    
    %Adjustment multiplication (negative exponential renders better results?)
    fff = sqrt(pi/t)*exp(t*(-M/2:M/2-1).^2).*ff/N;

    %Compute 2-norm of error
    a = norm(fff-valid);
    res_norm = [res_norm a];
end
loglog(t_vec,res_norm)
figure
plot(0:2*pi/(Mr-1):2*pi,ft)