%Stokes_kaiser
function varargout  = SE3P_Stokes(eval_idx,x,f,opt)

%Get input and verify that it is on correct form
opt = parse_params(opt); 

%Precompute the grid and compute the convolution onto the precomputed
%points
W_precomp_gridding = @(x,opt) SE3P_Laplace_gridding_precomp(x, opt); %Needs to change to SE3P_Stokes_grid_precomp
if strcmp(opt.window, 'gaussian')
   W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_thrd_mex(x,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  W_fast_int = @(F,opt,S) SE_fg_int_split_mex(0,F,opt,...
                                S.zs,S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_mex is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.

elseif strcmp(opt.window, 'expsemicirc') || strcmp(opt.window, 'kaiser_exact') ...
    || strcmp(opt.window, 'kaiser_poly')
  W_fast_grid = @(x,F,opt,S) SE_fg_grid_split_kaiser_mex(x,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  W_fast_int = @(F,opt,S) SE_fg_int_split_kaiser_mex(0,F,opt,...
                                S.zx,S.zy,S.zz,S.idx);
  % NB: The first argument to SE_fg_int_split_kaiser_mex is unused.
  % It is not possible to specify the evaluation points, so we
  % have to pick out the correct evaluation points at the end.
else
  error('Unsupported window function');
end
S = W_precomp_gridding(x, opt);
grid_fcn = @(F) W_fast_grid(x(S.perm,:), F(S.perm), opt, S); 

%Integration functions, check if needs to be remade to suit Stokes
iperm = @(u) u(S.iperm,:);
int_fcn = @(F) iperm(W_fast_int(F, opt, S));

%To grid
  
H1 = grid_fcn(f(:,1));
H2 = grid_fcn(f(:,2));
H3 = grid_fcn(f(:,3));

%FFT step
G1 = fftshift( fftn(H1) );
G2 = fftshift( fftn(H2) );
G3 = fftshift( fftn(H3) );

if isreal(G1) || isreal(G2) || isreal(G3)
    G1 = complex(G1);
    G2 = complex(G2);
    G3 = complex(G3);
end
[G1, G2, G3] = SE3P_Stokes_scaling(G1,G2,G3,opt); 

% inverse shift and inverse transform
F1 = real( ifftn( ifftshift( G1 )));
F2 = real( ifftn( ifftshift( G2 )));
F3 = real( ifftn( ifftshift( G3 )));

% Integrate

u = zeros(length(eval_idx),3);
u(:,1) = int_fcn(F1);
u(:,2) = int_fcn(F2);
u(:,3) = int_fcn(F3);

varargout{1} = u;
end

function opt = parse_params(opt)
% Check that mandatory options are present
assert(isfield(opt, 'box'), 'cell size box must be given in opt struct');
assert(isfield(opt, 'M'), 'grid size M must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald parameter xi must be given in opt struct');
assert(isfield(opt, 'P'), 'window support P must be given in opt struct');

if ~isfield(opt,'potential'), opt.potential = true; end
if ~isfield(opt,'force'), opt.force = false; end
if ~isfield(opt,'fourier_differentiation'), opt.fourier_differentiation = false; end

% Verify assumptions on parameters
opt.L = opt.box(1);
opt.h = opt.L/opt.M(1); % step size (M contains number of subintervals)
% Check that h is the same in all directions
assert(abs(opt.h - opt.box(2)/opt.M(2)) < eps);
assert(abs(opt.h - opt.box(3)/opt.M(3)) < eps);

if ~isfield(opt,'window'), opt.window = 'gaussian'; end
if ~isfield(opt,'w'), opt.w = opt.h*opt.P/2; end
if ~isfield(opt,'m'), opt.m = 1.71*sqrt(opt.P); end
opt.eta = (2*opt.w*opt.xi/opt.m)^2;
opt.c = 2*opt.xi^2/opt.eta;
if ~isfield(opt,'fast_gridding'), opt.fast_gridding = true; end

% Options for ExpSemiCirc and Kaiser-Bessel windows
if strcmp(opt.window,'expsemicirc') || strcmp(opt.window,'kaiser_exact') ...
    || strcmp(opt.window,'kaiser_poly')
  if ~isfield(opt,'betaP'), opt.betaP = 2.5; end
  opt.beta = opt.betaP*opt.P;
end
if strcmp(opt.window,'kaiser_exact') || strcmp(opt.window,'kaiser_poly')
  opt.kaiser_scaling = 1/besseli(0,opt.beta);
end
if strcmp(opt.window,'kaiser_poly')
  if ~isfield(opt,'polynomial_degree'), opt.polynomial_degree = 1; end
end

% Half-support of window
if mod(opt.P,2)==0
  opt.p_half = opt.P/2;
else
  opt.p_half = (opt.P-1)/2;
end

% External evaluation points
if isfield(opt,'eval_x') && numel(opt.eval_x) > 0
  opt.eval_external = true;
  fprintf('WARNING: External evaluation points are not implemented\n');
else
  opt.eval_external = false;
  opt.eval_x = [];
end
% TODO: external evaluation points are not implemented yet!
end