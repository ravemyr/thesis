function G = SE0P_Laplace_scaling(H, opt, pre_kernel, pre_window)
% SE0P_LAPLACE_SCALING  Fourier-space scaling for 0-periodic electrostatics
%    G = SE0P_LAPLACE_SCALING(H,OPT,PRE_KERNEL,PRE_WINDOW)
%
% Parameters:
% :param H: Input data in Fourier space
% :param OPT: Structure with Ewald options
% :param PRE_KERNEL: Mandatory precomputation structure:
% :param PRE_KERNEL.GR: Green's function
% :param PRE_WINDOW: Optional precomputation structure:
% :param PRE_WINDOW.F: W^(-2) where W is the Fourier transform of the window function
%
% :returns: **G** -- Output data in Fourier space

if nargin < 4
  pre_window = [];
end

% Compute k-vectors
[k1, k2, k3] = NEW_k_vectors(opt.padded_M, opt.padded_box, 'shifted');
[K1, K2, K3] = ndgrid(k1, k2, k3);
KK = K1.^2 + K2.^2 + K3.^2;

% Scale
if strcmp(opt.window, 'gaussian')
  Z = pre_kernel.GR .* exp(-(1-opt.eta)/(4*opt.xi^2)*KK);
else
  if ~isempty(pre_window)
    F = pre_window.F;
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
  Z = pre_kernel.GR .* F .* exp(-KK/(4*opt.xi^2));
end
G = H.*Z;

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;
