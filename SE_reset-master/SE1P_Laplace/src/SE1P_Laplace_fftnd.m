function [x,xr,x0] = SE1P_Laplace_fftnd(x,opt)
%SE1P_LAPLACE_FFTND 3-dimensional discrete Fourier Transform.
%   [x,xr,x0] = SE1P_Laplace_fftnd(x,opt)
%      Returns the 3-dimensional DFT of the 3D array x with one periodic
%      direction and two free directions, in the arrays x, xr and
%      x0. The output xr contains the local pad result and x0
%      contains the zero mode.
%
% Note: For simplicity, the local pad contains the zero mode as
% well but it will be updated by the zero mode result later on.
%
%   INPUT:
%       x:              Input vector of size (Mx,My,Mz)
%       opt.M:          Vector containing the sizes [Mx,My,Mz]
%       opt.My:         Size of second dimension of x
%       opt.Mz:         Size of third dimension of x
%       opt.local_pad:  List of modes to apply a large oversampling to
%       opt.s0:         Oversampling factor for the zero mode
%       opt.sl:         Oversampling factor on opt.local_pad
%       opt.sg:         Oversampling factor on the rest of the domain
%
%   OUTPUT:
%       x:              Output vector of size (Mx,My*sg,Mz*sg)
%       xr:             Output vector of size (numel(opt.local_pad),My*sl,Mz*sl)
%       x0:             Output vector of size (1,My*s0,Mz*s0)
%
% The periodic direction is x, while the free directions are y and z.

if opt.sg == opt.sl && opt.sg == opt.s0
  x = fftn(x, [opt.M(1) round(opt.My*opt.sg) round(opt.Mz*opt.sg)]);
  xr = [];
  x0 = [];
else
  % Global domain
  Fx = fft(x, opt.M(1));                          % 1D fft in x
  x = fft(fft(Fx, round(opt.My*opt.sg), 2), ...
          round(opt.Mz*opt.sg), 3);               % 2D fft with global padding in y and z

  % Local pad
  xr = fft(fft(Fx(opt.local_pad,:,:), round(opt.My*opt.sl), 2), ...
           round(opt.Mz*opt.sl), 3);              % 2D fft with local padding in y and z

  % Zero mode
  F0 = squeeze(Fx(1,:,:));
  x0 = fft2(F0, round(opt.My*opt.s0), round(opt.Mz*opt.s0));
end

end
