function pre = SE0P_window_fft_precomp(opt)

L = opt.oversampled_box(1);
M = opt.oversampled_M(1);
h = L/M;
w = opt.w * h/opt.h; % w for oversampled grid

if strcmp(opt.window, 'gaussian')
  error('FFT precomputing not supported for Gaussian window');
elseif strcmp(opt.window, 'expsemicirc')
  window = @(x) expsemicirc(x, opt.beta, w);
elseif strcmp(opt.window, 'kaiser_exact') || strcmp(opt.window, 'kaiser_poly')
  window = @(x) kaiser_exact(x, opt.beta, w, opt.kaiser_scaling);
else
  error('Unsupported window function');
end

x = 0:h:(L-h);
f = h*window(x-L/2);
n = opt.padded_M(1);
N = opt.oversampled_M(1);
start = round((N-n)/2);
f = f(start + (1:n));
F = fft(f);
FF = 1./F.^2;
[F1,F2,F3] = ndgrid(FF, FF, FF);
F = F1.*F2.*F3; % tensor product of spatial directions
pre.F = fftshift(F); % TODO: why fftshift?

% ------------------------------------------------------------------------------
function z = expsemicirc(x, beta, w)

t = sqrt(1-(x/w).^2);
z = exp(beta*(t-1)) .* (abs(x) <= w);

% ------------------------------------------------------------------------------
function z = kaiser_exact(x, beta, w, scaling)

t = sqrt(1-(x/w).^2);
I = besseli(0,beta*t) * scaling;
z = I .* (abs(x) <= w);
