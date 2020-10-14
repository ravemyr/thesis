function x = SE2P_Laplace_ifftnd(x,xr,x0,opt)
%SE2P_LAPLACE_IFFTND 3-dimensional inverse discrete Fourier Transform.
%   x = SE2P_Laplace_ifftnd(x,xr,x0,opt)
%      Returns the 3-dimensional inverse DFT of the 3D array x
%      with two periodic directions and one free direction.
%      The input xr contains the local pad result and x0
%      contains the zero mode.
%
% Note: For simplicity, the local pad contains the zero mode as
% well but it will be updated by the zero mode result later on.
%
%   INPUT:
%       x:              Input vector of size (Mx,My,Mz)
%       xr:             Local pad vector of size (numel(opt.local_pad),numel(opt.local_pad),Mz*s)
%       x0:             Zero mode vector of size (1,1,Mz*s0)
%       opt.M:          Vector containing the sizes [Mx,My,Mz]
%       opt.Mz:         Size of third dimension of output
%       opt.local_pad:  List of modes to apply a large oversampling to
%       opt.s0:         Oversampling factor for the zero mode
%       opt.s:          Oversampling factor on opt.local_pad
%
%   OUTPUT:
%       x:              Output vector of size (Mx,My,Mz)
%
% The periodic directions are x and y, while the free direction is z.

% Zero mode
F0 = ifft(x0, round(opt.Mz*opt.s0)); % since this is a vector we skip 3

% Local pad
Fr = ifft(xr, round(opt.Mz*opt.s), 3); % 1D ifft in z

% Global domain
F = ifft(x, opt.Mz, 3); % 1D ifft in z

% Put together
Fz = zeros(opt.M(1), opt.M(2), opt.Mz);
Fz(:,:,:) = F(:, :, 1:opt.Mz); % restrict
Fz(opt.local_pad,opt.local_pad,:) = Fr(:, :, 1:opt.Mz); % restrict
Fz(1,1,:) = F0(1:opt.Mz); % restrict to (1,1,Mz)

x = ifft2(Fz, opt.M(1), opt.M(2)); % 2D ifft in x and y

end
