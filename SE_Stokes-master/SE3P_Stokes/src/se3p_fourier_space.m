function [u varargout]  = se3p_fourier_space(eval_idx,x,f,xi,opt)
% Compute Fourier space part of Ewald sum for periodic stokeslet potential.
%
% u = se3p_fourier_space(eval_idx,x,f,xi,opt)
%   Return potential
%
% [U1, U2, U3] = se3p_fourier_space(eval_idx,x,f,xi,opt)
%   Return Fourier coefficients
%
% :param eval_idx: index of source locations where potential should be evaluated
% :param x: source locations (Nx3)
% :param f: source strengths (Nx3)
% :param xi: Ewald paramter
% :param opt: Ewald options
% :param opt.M: grid size (M1, M2, M3)
% :param opt.P: Gaussian width
% :param opt.box: Box size (L1, L2, L3)
% :param opt.window: Window function ('kaiser' or 'gaussian')
% :returns: **phi** -- Fourier space potential
% :returns: **U1,U2,U3** -- Fourier space potential

verb = false;
if isfield(opt, 'window') && strcmp(opt.window, 'kaiser')
  use_kaiser = true;
else
  use_kaiser = false;
end

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

% parameters and constants
opt = se3p_parse_params(opt);
[w m M P] = unpack_params(opt);
opt.xi = xi;
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;
% to grid function
if use_kaiser
  pre_t = tic;
  S = precomp_kaiser(x,opt);
  walltime.pre = toc(pre_t);
  grid_fcn = @(f) SE_fg_grid_split_kaiser_mex(x(S.perm,:),f(S.perm),opt,S.zx,S.zy,S.zz, ...
                                            S.idx);
  SI = S;
  iperm = @(u) u(SI.iperm,:);
  int_fcn = @(F) iperm(SE_fg_int_split_kaiser_mex(0,F,opt,SI.zx,SI.zy,SI.zz,SI.idx));
else
  pre_t = tic;
  S = SE_FGG_precomp(x,xi,opt);
  walltime.pre = toc(pre_t);
  grid_fcn = @(f) SE_fg_grid_split_thrd_mex(x(S.perm,:),f(S.perm),opt,S.zs,S.zx,S.zy,S.zz, ...
                                            S.idx);
  SI = S;
  iperm = @(u) u(SI.iperm,:);
  int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));
end
grid_t = tic;
H1 = grid_fcn(f(:,1));
H2 = grid_fcn(f(:,2));
H3 = grid_fcn(f(:,3));
walltime.grid = walltime.grid + toc(grid_t);

% transform and shift
fft_t = tic;
G1 = fftshift( fftn(H1) );
G2 = fftshift( fftn(H2) );
G3 = fftshift( fftn(H3) );
walltime.fft = walltime.fft + toc(fft_t);

% multiply with modified greens function
if isreal(G1) || isreal(G2) || isreal(G3)
    G1 = complex(G1);
    G2 = complex(G2);
    G3 = complex(G3);
end
if use_kaiser
  pre = se3p_window_precomp(opt);
  scale_t = tic;
  G = se3p_k_scaling_kaiser({G1,G2,G3}, opt, pre);
  walltime.scale = toc(scale_t);
  G1 = G{1}; G2 = G{2}; G3 = G{3};
else
  scale_t = tic;
  [G1 G2 G3] = se3p_fast_k_scaling(G1,G2,G3,xi,opt.box,eta);
  walltime.scale = toc(scale_t);
end

% inverse shift and inverse transform
fft_t = tic;
F1 = real( ifftn( ifftshift( G1 )));
F2 = real( ifftn( ifftshift( G2 )));
F3 = real( ifftn( ifftshift( G3 )));
walltime.fft = walltime.fft + toc(fft_t);

% Integrate
u = zeros(length(eval_idx),3);
int_t = tic;
u(:,1) = int_fcn(F1);
u(:,2) = int_fcn(F2);
u(:,3) = int_fcn(F3);
walltime.int = walltime.int + toc(int_t);

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end


% ------------------------------------------------------------------------------
function [w m M P] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
P = opt.P;
