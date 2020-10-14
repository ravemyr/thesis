function [G, Gr, G0] = SE1P_Laplace_scaling(G, Gr, G0, opt, pre)
% SE1P_LAPLACE_SCALING  Fourier-space scaling for 1-periodic electrostatics
%    [G,GR,G0] = SE1P_LAPLACE_SCALING(G,GR,G0,OPT,PRE)
%
% Parameters:
% :param [G,GR,G0]: Input data in Fourier space
% :param OPT: Structure with Ewald options
% :param PRE: Optional precomputation structure:
% :param PRE.F: W^(-2) where W is the Fourier transform of the window function
%
% :returns: **[G,GR,G0]** -- Output data in Fourier space

if nargin < 5
  pre = [];
end

% Compute k-vectors
[k1, zidx1] = SE1P_k_vectors(opt.M(1), opt.box(1), 1);
kappa_2     = SE1P_k_vectors(opt.My, opt.Ly, opt.sg);
kappa_3     = SE1P_k_vectors(opt.Mz, opt.Lz, opt.sg);
[K1, KAPPA2, KAPPA3] = ndgrid(k1, kappa_2, kappa_3);

% Scale the whole Fourier domain
KK = K1.^2 + KAPPA2.^2 + KAPPA3.^2;
if strcmp(opt.window, 'gaussian')
  Z = exp(-(1-opt.eta)/(4*opt.xi^2)*KK) ./ KK;
else
  if ~isempty(pre)
    F = pre.F;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    b2 = opt.beta^2;
    f1 = kaiser_exact_ft(k1.^2, b2, opt.w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(kappa_2.^2, b2, opt.w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(kappa_3.^2, b2, opt.w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^2, f2.^2, f3.^2);
    F = F1.*F2.*F3; % tensor product of spatial directions
    F = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Z = F .* exp(-KK/(4*opt.xi^2)) ./ KK;
end
Z(zidx1,1,1) = 0; % Apparently, zidx1=1 always?

% Scale the zero mode
if opt.sg == opt.s0
  KK = squeeze(KK(1,:,:));
else
  kappa_2 = SE1P_k_vectors(opt.My, opt.Ly, opt.s0);
  kappa_3 = SE1P_k_vectors(opt.Mz, opt.Lz, opt.s0);
  [KAPPA2, KAPPA3] = ndgrid(kappa_2, kappa_3);
  KK = k1(zidx1)^2 + KAPPA2.^2 + KAPPA3.^2;
end
K = sqrt(KK);
R = opt.R;
Green = (1-besselj(0,R*K))./KK - R*log(R)*besselj(1,R*K)./K;
if strcmp(opt.window, 'gaussian')
  Znum = exp(-(1-opt.eta)/(4*opt.xi^2)*KK);
  Z0 = Znum .* Green;
else
  if ~isempty(pre)
    F0 = pre.F0;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    f2 = kaiser_exact_ft(kappa_2.^2, b2, opt.w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(kappa_3.^2, b2, opt.w, opt.kaiser_scaling);
    F1 = f1(1).^2; [F2, F3] = ndgrid(f2.^2, f3.^2);
    F0 = F1*F2.*F3; % tensor product of spatial directions
    F0 = 1./F0;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Znum = exp(-KK/(4*opt.xi^2));
  Z0 = F0 .* Znum .* Green;
end
% Finite limit at k2=k3=0
Z0(1,1) = R^2/2 * (1/2 - log(R));
if opt.sg == opt.s0 && opt.sg == opt.sl
  Z(1,:,:) = Z0;
  G = G.*Z;
  return
end
G = G.*Z;
G0 = G0.*Z0;

% Scale the local pad
if numel(opt.local_pad) > 0
  n = opt.local_pad;
  kappa_2 = SE1P_k_vectors(opt.My, opt.Ly, opt.s);
  kappa_3 = SE1P_k_vectors(opt.Mz, opt.Lz, opt.s);
  [K1, KAPPA2, KAPPA3] = ndgrid(k1(n), kappa_2, kappa_3);
  KK = K1.^2 + KAPPA2.^2 + KAPPA3.^2;
  if strcmp(opt.window, 'gaussian')
    Zr = exp(-(1-opt.eta)/(4*opt.xi^2)*KK) ./ KK;
  else
    if ~isempty(pre)
      Fr = pre.Fr;
    elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
      f2 = kaiser_exact_ft(kappa_2.^2, b2, opt.w, opt.kaiser_scaling);
      f3 = kaiser_exact_ft(kappa_3.^2, b2, opt.w, opt.kaiser_scaling);
      [F1, F2, F3] = ndgrid(f1(n).^2, f2.^2, f3.^2);
      Fr = F1.*F2.*F3; % tensor product of spatial directions
      Fr = 1./Fr;
    else
      error('Precomputed Fourier transform necessary for the %s window', opt.window);
    end
    Zr = Fr .* exp(-KK/(4*opt.xi^2)) ./ KK;
  end
  Gr = Gr.*Zr;
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;

% ------------------------------------------------------------------------------
function [ks, idxz] = SE1P_k_vectors(M, L, q)
% q: oversampling factor (1 = no oversampling)
% FIXME: this function is identical for SE1P and SE2P, move it out!

Mq = round(q*M);

if mod(Mq,2) == 0
    k = (-Mq/2):(Mq/2-1);
else
    k = -(Mq-1)/2:(Mq-1)/2;
end

k = fftshift(k); % standard reordering (in 1P at least)

idxz = 1;
ks = 2*pi*k/(L*q);
if (abs(ks(idxz)) > eps)
    ks = circshift(ks, [1 1]);
end

assert(abs(ks(idxz)) < eps);
