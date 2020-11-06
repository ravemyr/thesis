function [u varargout]= se0p_fourier_space(x, f, opt, pre_kernel, k_scaling)

if isfield(opt, 'window') && strcmp(opt.window, 'kaiser')
  use_kaiser = true;
else
  use_kaiser = false;
end

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

% Check
assert(isfield(opt, 'xi'), 'xi must be given in opt struct')
assert(isfield(opt, 'oversampling'), 'oversampling rate must be given in opt struct')

% Setup vars, modify opt
se0p_opt = se0p_parse_params(opt);
x = bsxfun(@plus, x, se0p_opt.delta);
opt = se0p_opt;
opt.box = se0p_opt.extended_box; % Work with box extended for Gaussian support
opt.M = se0p_opt.extended_M;

fsize = size(f);
N = fsize(1);
dim_in = fsize(2:end);

if use_kaiser
  % Precompute Kaiser window
  pre_window = se0p_window_precomp(se0p_opt);
end

% === Use vectorized code
%opt
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
  S = SE_FGG_precomp(x,opt.xi,opt);
  walltime.pre = toc(pre_t);
  grid_fcn = @(f) SE_fg_grid_split_mex(x(S.perm,:),f(S.perm),opt,S.zs,S.zx,S.zy,S.zz, ...
                                            S.idx);
  SI = S;
  iperm = @(u) u(SI.iperm,:);
  int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));
end

% === Uncomment for direct code
if use_kaiser
    %grid_fcn = @(f) SE_fg_grid_kaiser_mex(x,f,opt);
    %int_fcn = @(f) SE_fg_int_kaiser_mex(x,f,opt);
else
    %grid_fcn = @(f) SE_fg_grid_mex(x,f,opt);
    %int_fcn = @(f) SE_fg_int_mex(x,f,opt);
end


% grid + FFT
H = cell([dim_in, 1]);
for i=1:prod(dim_in)
    grid_t = tic;
    % Grid
    tmp = grid_fcn(f(:,i));
    walltime.grid = walltime.grid + toc(grid_t);
    % transform and shift, let FFT do padding
    fft_t = tic;
    H{i} = fftshift( fftn(tmp, se0p_opt.padded_M) );
    walltime.fft = walltime.fft + toc(fft_t);
end

% scale
scale_t = tic;
if use_kaiser
  G = se0p_k_scaling_kaiser(H, se0p_opt, pre_window, pre_kernel);
else
  G = se0p_k_scaling(H, se0p_opt, opt, pre_kernel);
end
walltime.scale = toc(scale_t);
dim_out = numel(G);

% inverse shift and inverse transform
for i=1:dim_out
    fft_t = tic;
    G{i} = ifftn( ifftshift(G{i}));
    walltime.fft = walltime.fft + toc(fft_t);
end

% Option 1: Truncate grid before integration
M = se0p_opt.extended_M;

% Option 2: Integrate directly on padded grid
%opt.M = se0p_opt.padded_M;
%opt.box = se0p_opt.padded_box;

u = zeros(N, dim_out);
for i=1:dim_out
    int_t = tic;
    u(:,i) = int_fcn(G{i}(1:M(1), 1:M(2), 1:M(3)));
    walltime.int = walltime.int + toc(int_t);
end

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end

