function [G1,G2,G3] = SE3P_Stokes_scaling(H1,H2,H3,opt)
% SE3P_Stokes_Scaling  Fourier space scaling for the stokeslet with
% different windowfunctions
%    H = SE3P_Stokes_SCALING(H,OPT,PRE)
%
% Parameters:
% :param H: Input data in Fourier space
% :param OPT: Structure with Ewald options
% :param PRE: Optional precomputation structure:
% :param PRE.F: W^(-2) where W is the Fourier transform of the window function
%
% :returns: **H** -- Output data in Fourier space


% Compute k-vectors
[k1, k2, k3] = NEW_k_vectors(opt.M, opt.box);
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% Scale
if strcmp(opt.window, 'gaussian')
  Z = 8*pi .* (1 + KK / (4*opt.xi^2)) .* exp(-(1-opt.eta) .* KK/4*opt.xi^2)./ KK.^2;
% else
%   if ~isempty(pre)
%     F = pre.F;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    b2 = opt.beta^2;
    f1 = kaiser_exact_ft(k1.^2, b2, opt.w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(k2.^2, b2, opt.w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(k3.^2, b2, opt.w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^2, f2.^2, f3.^2);
    F = F1.*F2.*F3; % tensor product of spatial directions
    F = 1./F;
%     pre.F = F;
    Q = (KK / 4 * opt.xi^2);
    Z = F .* 8*pi.*(1 + Q) .* exp(-Q) ./ (KK .* KK); 
%     G = se3p_k_scaling_kaiser({H1,H2,H3},opt,pre);
%     G1 = G{1}; G2 = G{2}; G3 = G{3};
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
%   end
end
Z(1,1,1) = 0;
KH = K1 .* H1 + K2 .* H2 + K3 .* H3;
G1 = Z .* (KK .* H1 - KH .* K1);
G2 = Z .* (KK .* H2 - KH .* K2);
G3 = Z .* (KK .* H3 - KH .* K3);
end
% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;

end