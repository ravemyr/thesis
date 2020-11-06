% ------------------------------------------------------------------------------
function [G Gr G0] = se2p_k_scaling_kaiser(H, Hr, H0, opt, pre_window)

% Non-zero mode (without upsampling)
% k-vectors
[k1 zidx1] = k_vectors(opt.M(1), opt.L, 1);
[k2 zidx2] = k_vectors(opt.M(2), opt.L, 1);
kappa = k_vectors(opt.Mz, opt.Lz, 1);
[K1 K2 KAPPA] = ndgrid(k1,k2,kappa);
Ksq = K1.^2 + K2.^2 + KAPPA.^2;

% scale
GR = Green(Ksq, 0);
B = (1 + Ksq/(4*opt.xi^2)) .* GR .* exp(-Ksq/(4*opt.xi^2)) .* pre_window.F;
KdotH = K1.*H{1} + K2.*H{2} + KAPPA.*H{3};
G{1} = B .* (Ksq.*H{1} - KdotH.*K1);
G{2} = B .* (Ksq.*H{2} - KdotH.*K2);
G{3} = B .* (Ksq.*H{3} - KdotH.*KAPPA);

% Non-zero mode (with upsampling)
% k-vectors
if (numel(opt.local_pad) > 0)
    kappa = k_vectors(opt.Mz, opt.Lz, opt.s);
    [K1 K2 KAPPA] = ndgrid(k1(opt.local_pad),k2(opt.local_pad),kappa);
    Ksq = K1.^2 + K2.^2 + KAPPA.^2;

    % scale
    GR = Green(Ksq, 0);
    B = (1 + Ksq/(4*opt.xi^2)) .* GR .* exp(-Ksq/(4*opt.xi^2)) .* pre_window.Fr;
    KdotH = K1.*Hr{1} + K2.*Hr{2} + KAPPA.*Hr{3};
    Gr{1} = B .* (Ksq.*Hr{1} - KdotH.*K1);
    Gr{2} = B .* (Ksq.*Hr{2} - KdotH.*K2);
    Gr{3} = B .* (Ksq.*Hr{3} - KdotH.*KAPPA);
end

% Zero mode
H0{1} = H0{1}(:);
H0{2} = H0{2}(:);
H0{3} = H0{3}(:);

% k-vectors
kappa = k_vectors(opt.Mz, opt.Lz, opt.s0);
[K1 K2 KAPPA] = ndgrid(k1(zidx1),k2(zidx2),kappa);
K1 = K1(:);
K2 = K2(:);
KAPPA = KAPPA(:);
Ksq = K1.^2 + K2.^2 + KAPPA.^2;

% scale
GR = -Green(Ksq, opt.R); % minus, because that's how it's done everywhere else
B = (1 + Ksq/(4*opt.xi^2)) .* GR .* exp(-Ksq/(4*opt.xi^2)) .* pre_window.F0;
G0{1} = B.*H0{1};
G0{2} = B.*H0{2};
G0{3} = 0*H0{3};

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
        G = 8*pi*(R*K.*sin(R*K) + cos(R*K) - 1) ./ Ksq;
        G(K==0) = 8*pi*R^2/2;
    end

