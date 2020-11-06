function G = se0p_k_scaling_kaiser(H, opt, pre_window, pre_kernel)

G = cell(3, 1);
GR = -pre_kernel.GR; % why minus? because the paper says so!

% k-vectors
[k1 k2 k3] = k_vectors(opt.padded_M, opt.padded_box);
[K1 K2 K3] = ndgrid(k1,k2,k3);
Ksq = K1.^2 + K2.^2 + K3.^2;

c = Ksq/(4*opt.xi^2);
%d = 8*pi*(1+c)./(Ksq.*Ksq);
%d(K1==0 & K2==0 & K3==0) = 0;

% scale
B = (1+c) .* exp(-c) .* GR .* pre_window.F;
KdotH = K1.*H{1} + K2.*H{2} + K3.*H{3};
G{1} = B .* (Ksq.*H{1} - KdotH.*K1);
G{2} = B .* (Ksq.*H{2} - KdotH.*K2);
G{3} = B .* (Ksq.*H{3} - KdotH.*K3);
