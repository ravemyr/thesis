% ------------------------------------------------------------------------------
function [G Gr G0] = se1p_k_scaling(H, Hr, H0, opt)

% Non-zero mode (without upsampling)
% k-vectors
[k1 zidx1] = k_vectors(opt.M(1), opt.L, 1);
kappa2 = k_vectors(opt.My, opt.Ly, 1);
kappa3 = k_vectors(opt.Mz, opt.Lz, 1);
[K1 KAPPA2 KAPPA3] = ndgrid(k1,kappa2,kappa3);
Ksq = K1.^2 + KAPPA2.^2 + KAPPA3.^2;

% scale
GR = Green(Ksq, 0);
B = (1 + Ksq/(4*opt.xi^2)) .* GR .* exp(-(1-opt.eta)*Ksq/(4*opt.xi^2));
KdotH = K1.*H{1} + KAPPA2.*H{2} + KAPPA3.*H{3};
G{1} = B .* (Ksq.*H{1} - KdotH.*K1);
G{2} = B .* (Ksq.*H{2} - KdotH.*KAPPA2);
G{3} = B .* (Ksq.*H{3} - KdotH.*KAPPA3);

% Non-zero mode (with upsampling)
% k-vectors
if (numel(opt.local_pad) > 0)
    kappa2 = k_vectors(opt.My, opt.Ly, opt.s);
    kappa3 = k_vectors(opt.Mz, opt.Lz, opt.s);
    [K1 KAPPA2 KAPPA3] = ndgrid(k1(opt.local_pad),kappa2,kappa3);
    Ksq = K1.^2 + KAPPA2.^2 + KAPPA3.^2;

    % scale
    GR = Green(Ksq, 0);
    B = (1 + Ksq/(4*opt.xi^2)) .* GR .* exp(-(1-opt.eta)*Ksq/(4*opt.xi^2));
    KdotH = K1.*Hr{1} + KAPPA2.*Hr{2} + KAPPA3.*Hr{3};
    Gr{1} = B .* (Ksq.*Hr{1} - KdotH.*K1);
    Gr{2} = B .* (Ksq.*Hr{2} - KdotH.*KAPPA2);
    Gr{3} = B .* (Ksq.*Hr{3} - KdotH.*KAPPA3);
end

% Zero mode
% k-vectors
kappa2 = k_vectors(opt.My, opt.Ly, opt.s0);
kappa3 = k_vectors(opt.Mz, opt.Lz, opt.s0);
[KAPPA2 KAPPA3] = ndgrid(kappa2,kappa3);
K1 = k1(zidx1);
Ksq = K1.^2 + KAPPA2.^2 + KAPPA3.^2;

% scale
GR = -Green(Ksq, opt.R); % minus, because that's how it's done everywhere else
B = (1 + Ksq/(4*opt.xi^2)) .* GR .* exp(-(1-opt.eta)*Ksq/(4*opt.xi^2));
KdotH = K1.*H0{1} + KAPPA2.*H0{2} + KAPPA3.*H0{3};
G0{1} = B .* (Ksq.*H0{1} - KdotH.*K1);
G0{2} = B .* (Ksq.*H0{2} - KdotH.*KAPPA2);
G0{3} = B .* (Ksq.*H0{3} - KdotH.*KAPPA3);


% ------------------------------------------------------------------------------
function [ks idxz] = k_vectors(M,L,q)

% q: oversampling factor (1 = no oversampling)
Mq = q*M;

if mod(Mq,2)==0
    k = (-Mq/2):(Mq/2-1);
else
    k = -(Mq-1)/2:(Mq-1)/2;
end

k = fftshift(k); % standard reordering

idxz = 1;
ks = 2*pi*k/(L*q);
if (abs(ks(idxz))>eps)
    ks = circshift(ks,[1 1]);
end

assert(abs(ks(idxz))<eps)

% ------------------------------------------------------------------------------
function G = Green(Ksq, R)
    if R == 0 % non-zero mode
        G = 8*pi./(Ksq.*Ksq);
        G(Ksq==0) = 0;
    else % zero mode
        K = sqrt(Ksq);
        G = 8*pi*( ...
          (besselj(0,R*K) - 1)./(Ksq.*Ksq) ...
          - (R^3*(log(R)-1).*besselj(1,R*K) ./ (4*K)) ...
          + (R*log(R)*besselj(1,R*K))./(K.*Ksq) ...
          - (R^2*(2*log(R)-1)*besselj(0,R*K))./(4*Ksq) ...
        );
        G(K==0) = 8*pi*(1/64)*R^4*(5-4*log(R));
    end

