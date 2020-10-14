function pre = SE3P_window_fft_precomp(opt)

if strcmp(opt.window, 'gaussian')
  error('FFT precomputing not supported for Gaussian window');
elseif strcmp(opt.window, 'expsemicirc')
  window = @(x) expsemicirc(x, opt.beta, opt.w);
elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
  window = @(x) kaiser_exact(x, opt.beta, opt.w, opt.kaiser_scaling);
else
  error('Unsupported window function');
end

L = opt.L;
h = opt.h;
x = 0:h:(L-h);
f = h*window(x-L/2);

F = fft(f);
[F1,F2,F3] = ndgrid(F.^2);
F = F1.*F2.*F3; % tensor product of spatial directions
pre.F = 1./F;

% TODO: In old se3p_window_precomp.m, pre.F(:,:,M/2)=0. Why?
% In old precomp.m pre.F(isnan(F))=0 instead, which doesn't happen.
% It seems to be fine to use what we get here directly.

% ------------------------------------------------------------------------------
function z = expsemicirc(x, beta, w)

t = sqrt(1-(x/w).^2);
z = exp(beta*(t-1)) .* (abs(x) <= w);

% ------------------------------------------------------------------------------
function z = kaiser_exact(x, beta, w, scaling)

t = sqrt(1-(x/w).^2);
I = besseli(0,beta*t) * scaling;
z = I .* (abs(x) <= w);
