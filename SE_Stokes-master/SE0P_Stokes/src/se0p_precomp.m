function pre = se0p_precomp(opt, kernel)

se0p_opt = se0p_parse_params(opt);

% Modified Greens function, corresponding to delta in corner
[k1 k2 k3] = k_vectors(se0p_opt.oversampled_M, se0p_opt.oversampled_box);
[K1 K2 K3] = ndgrid(k1,k2,k3);
K = sqrt(K1.^2 + K2.^2 + K3.^2);
GR = kernel(K, se0p_opt.R);
%GR = ifftn(ifftshift(GR));
GR = ifftn(GR);
% GR now has rubbish in the center,
% truncate in real space by shifting, picking out center and shifting back.
% This is faster than centering by multiplying with e^(i*k*xc) in k-space.
n = se0p_opt.padded_M;
N = se0p_opt.oversampled_M;
GR = fftshift(GR);
start = floor( (N-n)/2 );
GR = GR( start(1) + (1:n(1)), ...
         start(2) + (1:n(2)), ...
         start(3) + (1:n(3)) );
GR = ifftshift(GR);
%pre.GR = fftshift(fftn(GR));
pre.GR = fftn(GR);
