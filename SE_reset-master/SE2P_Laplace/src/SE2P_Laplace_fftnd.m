function [x,xr,x0] = SE2P_Laplace_fftnd(x,opt)
%SE2P_LAPLACE_FFTND 3-dimensional discrete Fourier Transform.
%   [x,xr,x0] = SE2P_Laplace_fftnd(x,opt)
%      Returns the 3-dimensional DFT of the 3D array x with two periodic
%      directions and one free direction, in the arrays x, xr and
%      x0. The output xr contains the local pad result and x0
%      contains the zero mode.
%
% Note: For simplicity, the local pad contains the zero mode as
% well but it will be updated by the zero mode result later on.
%
%   INPUT:
%       x:              Input vector of size (Mx,My,Mz)
%       opt.M:          Vector containing the sizes [Mx,My,Mz]
%       opt.Mz:         Size of third dimension of x
%       opt.local_pad:  List of modes to apply a large oversampling to
%       opt.s0:         Oversampling factor for the zero mode
%       opt.s:          Oversampling factor on opt.local_pad
%
%   OUTPUT:
%       x:              Output vector of size (Mx,My,Mz)
%       xr:             Output vector of size (numel(opt.local_pad),numel(opt.local_pad),Mz*s)
%       x0:             Output vector of size (1,1,Mz*s0)
%
% The periodic directions are x and y, while the free direction is z.

% Global domain
Fxy = fft2(x, opt.M(1), opt.M(2));              % 2D fft in x and y
x   = fft(Fxy, opt.Mz, 3);                      % 1D fft with no padding in z

% Local pad
xr  = fft(Fxy(opt.local_pad,opt.local_pad,:),round(opt.Mz*opt.s),3); % 1D fft with padding in z

% Zero mode
F0 = Fxy(1,1,:); F0 = F0(:);
x0 = fft(F0,round(opt.Mz*opt.s0));

end
