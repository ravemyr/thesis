function pre = SE0P_Laplace_kernel_fft_precomp(opt)
% Precomputation for resolving the free-space kernel

% Modified Green's function, corresponding to delta in corner
[k1, k2, k3] = NEW_k_vectors(opt.oversampled_M, opt.oversampled_box, 'shifted');
[K1, K2, K3] = ndgrid(k1, k2, k3);
K = sqrt(K1.^2 + K2.^2 + K3.^2);
GR = 8*pi*(sin(K*opt.R/2)./K).^2;
GR(K==0) = 2*pi*opt.R^2;
GR = ifftn(GR);
% GR now has rubbish in the center,
% truncate in real space by shifting, picking out center and shifting back.
% This is faster than centering by multiplying with e^(i*k*xc) in k-space.
n = opt.padded_M;
N = opt.oversampled_M;
GR = fftshift(GR);
start = floor((N-n)/2);
GR = GR(start(1) + (1:n(1)), ...
        start(2) + (1:n(2)), ...
        start(3) + (1:n(3)));
GR = ifftshift(GR);
pre.GR = fftn(GR);
