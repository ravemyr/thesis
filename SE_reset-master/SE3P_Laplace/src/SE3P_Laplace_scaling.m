function [H, K1, K2, K3] = SE3P_Laplace_scaling(H, opt, pre)
% SE3P_LAPLACE_SCALING  Fourier-space scaling for 3-periodic electrostatics
%    H = SE3P_LAPLACE_SCALING(H,OPT,PRE)
%
% Parameters:
% :param H: Input data in Fourier space
% :param OPT: Structure with Ewald options
% :param PRE: Optional precomputation structure:
% :param PRE.F: W^(-2) where W is the Fourier transform of the window function
%
% :returns: **H** -- Output data in Fourier space

if nargin < 3
  pre = [];
end

% Compute k-vectors
[k1, k2, k3] = NEW_k_vectors(opt.M, opt.box);
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% Scale
if strcmp(opt.window, 'gaussian')
  Z = exp(-(1-opt.eta)/(4*opt.xi^2)*KK) ./ KK;
else
  if ~isempty(pre)
    F = pre.F;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    b2 = opt.beta^2;
    f1 = kaiser_exact_ft(k1.^2, b2, opt.w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(k2.^2, b2, opt.w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(k3.^2, b2, opt.w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^2, f2.^2, f3.^2);
    F = F1.*F2.*F3; % tensor product of spatial directions
    F = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Z = F .* exp(-KK/(4*opt.xi^2)) ./ KK;
end
Z(1,1,1) = 0;
H = H.*Z;

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
