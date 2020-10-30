function pre = SE1P_window_fft_precomp(opt)

% We assume periodicity in x and free in y and z.

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
Ly = opt.Ly;
Lz = opt.Lz;
My = opt.My;
Mz = opt.Mz;
h = opt.h;

% Common in all parts
x = 0:h:(L-h);
f = h*window(x-L/2);
F = fft(f);
FF = F.^2;

% yz plane, whole domain, upsampling sg
sg = opt.sg;
y = 0:h:(sg*Ly-h);
fy = h*window(y-sg*Ly/2);
Fy = fft(fy, round(sg*My));
FyFy = Fy.^2;
z = 0:h:(sg*Lz-h);
fz = h*window(z-sg*Lz/2);
Fz = fft(fz, round(sg*Mz));
FzFz = Fz.^2;
[F1,F2,F3] = ndgrid(FF, FyFy, FzFz);
F = F1.*F2.*F3; % tensor product of spatial directions
pre.F = 1./F;

% yz plane, local pad, upsampling sl
n = opt.local_pad;
sl = opt.sl;
y = 0:h:(sl*Ly-h);
fy = h*window(y-sl*Ly/2);
Fy = fft(fy, round(sl*My));
FyFy = Fy.^2;
z = 0:h:(sl*Lz-h);
fz = h*window(z-sl*Lz/2);
Fz = fft(fz, round(sl*Mz));
FzFz = Fz.^2;
[F1,F2,F3] = ndgrid(FF(n), FyFy, FzFz);
Fr = F1.*F2.*F3; % tensor product of spatial directions
pre.Fr = 1./Fr;

% yz plane, zero mode, upsampling s0
s0 = opt.s0;
y = 0:h:(s0*Ly-h);
fy = h*window(y-s0*Ly/2);
Fy = fft(fy', round(s0*My));
FyFy = Fy.^2;
z = 0:h:(s0*Lz-h);
fz = h*window(z-s0*Lz/2);
Fz = fft(fz', round(s0*Mz));
FzFz = Fz.^2;
F1 = FF(1); [F2,F3] = ndgrid(FyFy, FzFz);
F0 = F1*F2.*F3; % tensor product of spatial directions
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
