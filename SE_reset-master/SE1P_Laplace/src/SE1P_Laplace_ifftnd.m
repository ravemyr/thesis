function x = SE1P_Laplace_ifftnd(x,xr,x0,opt)
%SE1P_LAPLACE_IFFTND 3-dimensional inverse discrete Fourier Transform.
%   x = SE1P_Laplace_ifftnd(x,xr,x0,opt)
%      Returns the 3-dimensional inverse DFT of the 3D array x
%      with one periodic direction and two free directions.
%      The input xr contains the local pad result and x0
%      contains the zero mode.
%
% Note: For simplicity, the local pad contains the zero mode as
% well but it will be updated by the zero mode result later on.
%
%   INPUT:
%       x:              Input vector of size (Mx,My*sg,Mz*sg)
%       xr:             Local pad vector of size (numel(opt.local_pad),My*sl,Mz*sl)
%       x0:             Zero mode vector of size (1,My*s0,Mz*s0)
%       opt.M:          Vector containing the sizes [Mx,My,Mz]
%       opt.My:         Size of second dimension of output
%       opt.Mz:         Size of third dimension of output
%       opt.local_pad:  List of modes to apply a large oversampling to
%       opt.s0:         Oversampling factor for the zero mode
%       opt.sl:         Oversampling factor on opt.local_pad
%       opt.sg:         Oversampling factor on the rest of the domain
%
%   OUTPUT:
%       x:              Output vector of size (Mx,My,Mz)
%
% The periodic direction is x, while the free directions are y and z.

if opt.sg == opt.sl && opt.sg == opt.s0
  x = ifftn(x);
  x = x(1:opt.M(1), 1:opt.My, 1:opt.Mz); % restrict
else
  global_pad = setdiff(1:opt.M(1), [1 opt.local_pad]); % global pad vector

  % Zero mode
  F0 = ifft2(x0, round(opt.My*opt.s0), round(opt.Mz*opt.s0)); % 2D ifft in y and z

  % Local pad
  Fr = ifft(ifft(xr, round(opt.My*opt.sl), 2), round(opt.Mz*opt.sl), 3);

  % Global domain
  F = ifft(ifft(x, round(opt.My*opt.sg), 2), round(opt.Mz*opt.sg), 3);

  % Put together
  Fxy = zeros(opt.M(1), opt.My, opt.Mz);
  Fxy(opt.local_pad,:,:) = Fr(:, 1:opt.My, 1:opt.Mz); % restrict to (numel(local_pad),My,Mz)
  Fxy(1,:,:) = F0(1:opt.My, 1:opt.Mz); % restrict to (1,My,Mz)
  Fxy(global_pad,:,:) = F(global_pad, 1:opt.My, 1:opt.Mz); % restrict

  x = ifft(Fxy, opt.M(1)); % 1D ifft in x
end

end
