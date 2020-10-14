function pre = SE2P_window_fft_precomp(opt)

% We assume periodicity in x and y and free in z.

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
Lz = opt.Lz;
Mz = opt.Mz;
h = opt.h;

% Common in all parts
assert(L == opt.box(2), 'Code assumes Lx = Ly');
x = 0:h:(L-h);
f = h*window(x-L/2);
F = fft(f);
FF = F.^2;

% z direction, whole domain, no upsampling
x = 0:h:(Lz-h);
f = h*window(x-Lz/2);
Fz = fft(f, Mz);
FzFz = Fz.^2;
[F1,F2,F3] = ndgrid(FF, FF, FzFz);
F = F1.*F2.*F3; % tensor product of spatial directions
pre.F = 1./F;
%pre.F(:,:,Mz/2+1) = 0;
% TODO: In old se2p_window_precomp.m, pre.F(:,:,Mz/2+1)=0. Why?
% It seems to be fine to use what we get here directly.

% z direction, local pad, upsampling s
n = opt.local_pad;
s = opt.s;
x = 0:h:(s*Lz-h);
f = h*window(x-s*Lz/2);
Fz = fft(f, round(s*Mz));
FzFz = Fz.^2;
[F1,F2,F3] = ndgrid(FF(n), FF(n), FzFz);
Fr = F1.*F2.*F3; % tensor product of spatial directions
pre.Fr = 1./Fr;

% z direction, zero mode, upsampling s0
s0 = opt.s0;
x = 0:h:(s0*Lz-h);
f = h*window(x-s0*Lz/2);
Fz = fft(f', round(s0*Mz));
FzFz = Fz.^2;
F1 = FF(1); F2 = FF(1); F3 = FzFz;
F0 = F1*F2*F3; % tensor product of spatial directions
pre.F0 = 1./F0;

% ------------------------------------------------------------------------------
function z = expsemicirc(x, beta, w)

t = sqrt(1-(x/w).^2);
z = exp(beta*(t-1)) .* (abs(x) <= w);

% ------------------------------------------------------------------------------
function z = kaiser_exact(x, beta, w, scaling)

t = sqrt(1-(x/w).^2);
I = besseli(0,beta*t) * scaling;
z = I .* (abs(x) <= w);
