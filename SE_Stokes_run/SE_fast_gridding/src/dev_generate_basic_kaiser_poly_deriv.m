% Script to generate `gen_basic_kaiser_poly_deriv.c'
% which is included by `SE_fg_windows.c' for computing forces.

% The generated code in `gen_basic_kaiser_poly_deriv.c' expects input:
% - z: the local grid offset scaled to [-1,1]
% - P: support of window function (number of subintervals)
% - degree: degree of polynomial approximation
% and writes the values of the window function to "out".

function dev_generate_basic_kaiser_poly_deriv(window)
  opt.degs = 1:13;
  opt.Ps = {
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
    [2:2:16],
  };
  opt.Ppad = true; % pad evaluation points to multiple of 4

  if nargin < 1
    window = 'kaiser';
    %window = 'expsemicirc';
  end
  if strcmp(window, 'kaiser')
    I0 = @(z) besseli(0,z);
    I1 = @(z) besseli(1,z);
    K = @(x,beta) I1(beta*sqrt(1-x.^2))./I0(beta) .* (-beta*x./sqrt(1-x.^2));
  elseif strcmp(window, 'expsemicirc')
    K = @(x,beta) exp(beta*(sqrt(1-x.^2)-1)) .* (-beta*x./sqrt(1-x.^2));
  else
    error('Invalid window function name');
  end
  opt.window = window;

  betafun = @(P) 2.5*P;
  dev_generate_polynomial_code('gen_basic_kaiser_poly_deriv.c', 'dev_generate_basic_kaiser_poly_deriv.m', ...
                               K, betafun, opt);
end
