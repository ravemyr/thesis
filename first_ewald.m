%Fast Ewald prelim work
%Emanuel Ravemyr 21-9-2020

%Preamble
N = 64;
M = 4*N;
d = 1;
L = 1;
P = 9;
Q = zeros(N,1);
x = rand(N,d)*L;
q = (-N/2:N/2-1).*rand(1,N);
%Grid size h, length L
L = 1;
h = L/M;
% w = hP/2
%Validation computation
for j =1:N
    for i =1:N
        if(i~=j)
            Q(j) = Q(j) + q(i)/norm(x(j,:)-x(i,:));
        end
    end
end
%%
close all;
%Fast Ewald
% Decomposed into real sum and k-space sum
% real sum:
%First find neighbours within a radius r_c 
f = @(x,y) norm(x-y);
A = zeros(N,N);
rc =0.9;
for i = 1:N
    for j = 1:N
        if f(x(i,:),x(j,:)) < rc
            A(i,j) = 1;
        end
    end
end

%Compute the real sum
R = zeros(N,1);
xi = 0.01;

for i = 1:N
    for ii = 1:N
        if(i~=ii)
            R(i) = R(i) + q(ii)*erfc(xi*norm(x(i,:)-x(ii,:)))*A(i,ii)/(norm(x(i,:)-x(ii,:)));
        end
    end
end

%k-space sum
w = h*P/2;
m = sqrt(P);
eta = (2*xi*w/m)^2; %Optimal choice for Gaussian: (2wxi/m)^2
k = (1:M)'*2*pi/L;

%Uniform grid on [0,L)^d
p = zeros(M,d);
if(d==2)
    for i = 1:sqrt(M)
        for ii = 1:sqrt(M)
            p((i-1)*8 + ii,:) = [i,ii]/(sqrt(M)+1);
        end
    end
elseif(d ==1)
        for i = 1:M
           p(i) = i/(M+1);
        end
end
Ha = H(p,x,q,xi,eta,P);
if(d==2)
    FHa = fft2(Ha);
elseif(d==1)
    FHa = fft(Ha);
end
F = FHa.*exp(-(1-eta)*(k.^2)/(4*xi^2))./(k.^2);
if(d==2)
    FHb = ifft2(F);
elseif(d==1)
    FHb = ifft(F);
end
Hb = zeros(M,1);
if(d==2)
    for i = 1:sqrt(M)
        for ii = 1:sqrt(M)
            Hb(8*(i-1)+ii) = FHb(i,ii);
        end
    end
elseif(d==1)
    Hb = FHb;
end
T = zeros(N,1);

for i = 1:N
    phi = zeros(M,1);
    for ii = 1:M
        phi(ii) = phi(ii) + Hb(ii)*((2*xi^2/(pi*eta))^(3/2))*exp((-2*(xi^2)*(norm(p(ii,:)-x(i,:)))^2)/eta);
    end
    T(i) = 4*pi*h*(sum(phi(2:M-1))+(phi(M)+phi(1))/2);
end

Fin = R+T;
norm(Fin-Q)
plot(Fin)
figure
plot(Q)
figure
plot(abs(Fin-Q))
