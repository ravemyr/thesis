function varargout  = SE_Stokes(eval_idx,x,f,xi,opt)
% Compute Fourier space part of Ewald sum for periodic stokeslet potential.
%
% u = SE_Stokes(eval_idx,x,f,xi,opt)
%   Return potential
%
% [U1, U2, U3] = SE_Stokes(eval_idx,x,f,xi,opt)
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
% :returns: **phi** -- Fourier space potential
% :returns: **U1,U2,U3** -- Fourier space potential

verb = false;

% parameters and constants
opt = parse_params(opt);
[w m M P] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;
% to grid function
S = SE_FGG_precomp(x,xi,opt);
grid_fcn = @(f) SE_fg_grid_split_thrd_mex(x(S.perm,:),f(S.perm),opt,S.zs,S.zx,S.zy,S.zz, ...
                                          S.idx);
H1 = grid_fcn(f(:,1));
H2 = grid_fcn(f(:,2));
H3 = grid_fcn(f(:,3));

cprintf(verb, 'M = [%d %d %d] P = %d m=%d w=%f\n',M,P,m,w);
cprintf(verb, 'eta = %f\t a=%f\n', eta, pi^2/opt.c);

% transform and shift
G1 = fftshift( fftn(H1) );
G2 = fftshift( fftn(H2) );
G3 = fftshift( fftn(H3) );

% multiply with modified greens function
if isreal(G1) || isreal(G2) || isreal(G3)
    G1 = complex(G1);
    G2 = complex(G2);
    G3 = complex(G3);
end
[G1 G2 G3] = stokeslet_fast_k_scaling(G1,G2,G3,xi,opt.box,eta);

if nargout > 1
    varargout{1} = G1;
    varargout{2} = G2;
    varargout{3} = G3;
    return
end

% inverse shift and inverse transform
F1 = real( ifftn( ifftshift( G1 )));
F2 = real( ifftn( ifftshift( G2 )));
F3 = real( ifftn( ifftshift( G3 )));

% Integrate
SI = SE_FGG_precomp(x(eval_idx,:),xi,opt);
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));
u = zeros(length(eval_idx),3);
u(:,1) = int_fcn(F1);
u(:,2) = int_fcn(F2);
u(:,3) = int_fcn(F3);

varargout{1} = u;

% ------------------------------------------------------------------------------
function p = parse_params(opt)

% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% verify all assumptions on parameters

% step size
L = opt.box(1);
h = L/opt.M(1);
assert(abs(opt.box(2)/opt.M(2) - opt.box(1)/opt.M(1)) < eps)
assert(abs(opt.box(3)/opt.M(3) - opt.box(2)/opt.M(2)) < eps)

% Gaussian
P = opt.P;
if( isfield(opt,'m')), m = opt.m; else m = 0.9*sqrt(pi*P); end;
w = h*(P-1)/2;

% collect
p.M=opt.M;
p.P = P;
p.w = w;
p.m = m;
p.box = opt.box;
p.L = L;
p.h = h;


% ------------------------------------------------------------------------------
function [w m M P] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
P = opt.P;
