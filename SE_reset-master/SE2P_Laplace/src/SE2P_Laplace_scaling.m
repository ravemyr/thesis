function [G, Gr, G0] = SE2P_Laplace_scaling(G, Gr, G0, opt, pre)
% SE2P_LAPLACE_SCALING  Fourier-space scaling for 2-periodic electrostatics
%    [G,GR,G0] = SE2P_LAPLACE_SCALING(G,GR,G0,OPT,PRE)
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
[k1, zidx1] = SE2P_k_vectors(opt.M(1), opt.box(1), 1);
[k2, zidx2] = SE2P_k_vectors(opt.M(2), opt.box(2), 1);
kappa       = SE2P_k_vectors(opt.Mz, opt.Lz, 1);
[K1, K2, KAPPA] = ndgrid(k1, k2, kappa);

% Scale the whole Fourier domain
KK = K1.^2 + K2.^2 + KAPPA.^2;
if strcmp(opt.window, 'gaussian')
  Z = exp(-(1-opt.eta)/(4*opt.xi^2)*KK) ./ KK;
else
  if ~isempty(pre)
    F = pre.F;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    b2 = opt.beta^2;
    f1 = kaiser_exact_ft(k1.^2, b2, opt.w, opt.kaiser_scaling);
    f2 = kaiser_exact_ft(k2.^2, b2, opt.w, opt.kaiser_scaling);
    f3 = kaiser_exact_ft(kappa.^2, b2, opt.w, opt.kaiser_scaling);
    [F1, F2, F3] = ndgrid(f1.^2, f2.^2, f3.^2);
    F = F1.*F2.*F3; % tensor product of spatial directions
    F = 1./F;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Z = F .* exp(-KK/(4*opt.xi^2)) ./ KK;
end
Z(zidx1,zidx2,1) = 0; % Apparently, zidx1=zidx2=1 always?
G = G.*Z;

% Scale the zero mode
kappa = SE2P_k_vectors(opt.Mz, opt.Lz, opt.s0)';
KK = k1(zidx1).^2 + k2(zidx2).^2 + kappa.^2;
K = sqrt(KK);
R = opt.R;
Green = -(R*K.*sin(R*K)+cos(R*K)-1)./KK;
if strcmp(opt.window, 'gaussian')
  Znum = exp(-(1-opt.eta)/(4*opt.xi^2)*KK);
  Z0 = Znum .* Green;
else
  if ~isempty(pre)
    F0 = pre.F0;
  elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
    f3 = kaiser_exact_ft(kappa.^2, b2, opt.w, opt.kaiser_scaling);
    F1 = f1(1).^2; F2 = f2(1).^2; F3 = f3.^2;
    F0 = F1*F2*F3; % tensor product of spatial directions
    F0 = 1./F0;
  else
    error('Precomputed Fourier transform necessary for the %s window', opt.window);
  end
  Znum = exp(-KK/(4*opt.xi^2));
  Z0 = F0 .* Znum .* Green;
end
% Finite limit at k3=0
Z0(1) = R^2/4 * (1-2*log(R));
G0 = G0.*Z0;

% Scale the local pad
if numel(opt.local_pad) > 0
  n = opt.local_pad;
  kappa = SE2P_k_vectors(opt.Mz, opt.Lz, opt.s);
  [K1, K2, KAPPA] = ndgrid(k1(n), k2(n), kappa);
  KK = K1.^2 + K2.^2 + KAPPA.^2;
  if strcmp(opt.window, 'gaussian')
    Zr = exp(-(1-opt.eta)/(4*opt.xi^2)*KK) ./ KK;
  else
    if ~isempty(pre)
      Fr = pre.Fr;
    elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
      f3 = kaiser_exact_ft(kappa.^2, b2, opt.w, opt.kaiser_scaling);
      [F1, F2, F3] = ndgrid(f1(n).^2, f2(n).^2, f3.^2);
      Fr = F1.*F2.*F3; % tensor product of spatial directions
      Fr = 1./Fr;
    else
      error('Precomputed Fourier transform necessary for the %s window', opt.window);
    end
    Zr = Fr .* exp(-KK/(4*opt.xi^2)) ./ KK;
  end
  % Now the elements corresponding to the zero mode are Inf.
  % This is fine since we overwrite them later on.
  Gr = Gr.*Zr;
end

% ------------------------------------------------------------------------------
function F = kaiser_exact_ft(k2, b2, w, scaling)

t = sqrt(b2 - k2*w^2);
F = 2*w*sinh(t)./t * scaling;

% ------------------------------------------------------------------------------
function [ks, idxz] = SE2P_k_vectors(M, L, q)
% q: oversampling factor (1 = no oversampling)
% FIXME: this function is identical for SE1P and SE2P, move it out!

Mq = round(q*M);

if mod(Mq,2) == 0
    k = (-Mq/2):(Mq/2-1);
else
    k = -(Mq-1)/2:(Mq-1)/2;
end

k = fftshift(k); % standard reordering (in 2P at least)

idxz = 1;
ks = 2*pi*k/(L*q);
if (abs(ks(idxz)) > eps)
    ks = circshift(ks, [1 1]);
end

assert(abs(ks(idxz)) < eps);
