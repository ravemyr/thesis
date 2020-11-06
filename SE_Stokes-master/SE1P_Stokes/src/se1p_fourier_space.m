function [u walltime]  = se1p_fourier_space(eval_idx,x,f,opt)

if isfield(opt, 'window') && strcmp(opt.window, 'kaiser')
  use_kaiser = true;
else
  use_kaiser = false;
end

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);
    
% parameters and constants
opt = se1p_parse_params(opt);
fsize = size(f);
N = fsize(1);
dim_in = fsize(2:end);

if use_kaiser
  % Precompute Kaiser window
  pre_window = se1p_window_precomp(opt);
end

% === Use vectorized code
% Gridder
if use_kaiser
	pre_t = tic;
	S = se1p_precomp_kaiser(x,opt);
	walltime.pre = toc(pre_t);
	grid_fcn = @(f) SE_fg_grid_split_kaiser_mex_1p(x(S.perm,:),f(S.perm), ...
						  opt,S.zx,S.zy,S.zz,S.idx);
	% Integrator
	SI = S;
	iperm = @(u) u(SI.iperm,:);
	int_fcn = @(F) iperm(SE_fg_int_split_kaiser_mex_1p(0,F, ...
                                         opt,SI.zx,SI.zy,SI.zz,SI.idx));
else
	pre_t = tic;
	S = se1p_precomp(x,opt);
	walltime.pre = toc(pre_t);
	grid_fcn = @(f) SE_fg_grid_split_mex_1p(x(S.perm,:),f(S.perm), ...
						  opt,S.zs,S.zx,S.zy,S.zz,S.idx);
	% Integrator
	SI = S;
	iperm = @(u) u(SI.iperm,:);
	int_fcn = @(F) iperm(SE_fg_int_split_mex_1p(0,F, ...
						 opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));
end

% === Uncomment for direct code
if use_kaiser
    %grid_fcn = @(f) SE_fg_grid_kaiser_mex_1p(x,f,opt);
    %int_fcn = @(f) SE_fg_int_kaiser_mex_1p(x,f,opt);
else
    %grid_fcn = @(f) SE_fg_grid_mex_1p(x,f,opt);
    %int_fcn = @(f) SE_fg_int_mex_1p(x,f,opt);
end

% grid + FFT
H = cell([dim_in, 1]); 
Hr = cell([dim_in, 1]); 
H0 = cell([dim_in, 1]); 
for i=1:prod(dim_in)
    grid_t = tic;
    % Grid
    tmp = grid_fcn(f(:,i));
    walltime.grid = walltime.grid + toc(grid_t);
    % transform and shift, let FFT do padding
    fft_t = tic;
    [H{i} Hr{i} H0{i}] = fftnd1p(tmp, opt);
    walltime.fft = walltime.fft + toc(fft_t);
end

% scale
scale_t = tic;
if use_kaiser
	[G, Gr, G0] = se1p_k_scaling_kaiser(H, Hr, H0, opt, pre_window);
else
	[G, Gr, G0] = se1p_k_scaling(H, Hr, H0, opt);
end
walltime.scale = toc(scale_t);
dim_out = numel(G);

% inverse shift and inverse transform
F = cell([dim_in, 1]); 
for i=1:dim_out
    fft_t = tic;
    F{i} = ifftnd1p(G{i}, Gr{i}, G0{i}, opt);
    walltime.fft = walltime.fft + toc(fft_t);
end

% integrate
u = zeros(N, dim_out);
for i=1:dim_out
    int_t = tic;
    u(:,i) = int_fcn(F{i});
    walltime.int = walltime.int + toc(int_t);
end
if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end
