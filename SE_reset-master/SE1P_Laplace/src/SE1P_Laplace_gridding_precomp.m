function SE_static = SE1P_Laplace_gridding_precomp(x, opt)
% SE1P_LAPLACE_GRIDDING_PRECOMP  Precomputation of fast gridding vectors
%    SE_STATIC = SE1P_LAPLACE_GRIDDING_PRECOMP(X,OPT)
%
% Parameters:
% :param X: source locations (N×3)
% :param OPT: Structure with Ewald options
%
% :returns: **SE_STATIC** -- Structure containing precomputed data

x = recenter_points(x, opt.box);

if strcmp(opt.window, 'gaussian')
  if opt.force
    [zx, zy, zz, zfx, zfy, zfz, idx] = SE_fgg_expand_all_force_mex_1p(x,opt);
  else
    [zx, zy, zz, idx] = SE_fgg_expand_all_mex_1p(x,opt);
  end
  SE_static.zs = SE_fgg_base_gaussian_mex_1p(opt);
else
  [zx, zy, zz, idx] = SE_fkg_expand_all_mex_1p(x,opt);
end

[idx, s] = sort(idx);
x = x(s,:);
zx = zx(1:opt.P,s);
zy = zy(1:opt.P,s);
zz = zz(1:opt.P,s);

SE_static.zx = zx;
SE_static.zy = zy;
SE_static.zz = zz;
SE_static.idx = idx;
SE_static.perm = s';
SE_static.iperm(s) = 1:length(s);

if opt.force
  zfx = zfx(1:opt.P,s);
  zfy = zfy(1:opt.P,s);
  zfz = zfz(1:opt.P,s);
  SE_static.zfx = zfx;
  SE_static.zfy = zfy;
  SE_static.zfz = zfz;
end
