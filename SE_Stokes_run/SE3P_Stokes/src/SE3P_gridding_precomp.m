function SE_static = SE3P_gridding_precomp(x, opt)
%
% Parameters:
% :param X: source locations (N×3)
% :param OPT: Structure with Ewald options
%
% :returns: **SE_STATIC** -- Structure containing precomputed data

x = recenter_points(x, opt.box);

if strcmp(opt.window, 'gaussian')
  [zx, zy, zz, idx] = SE_fgg_expand_all_mex(x,opt);
  SE_static.zs = SE_fgg_base_gaussian_mex(opt);
else
  [zx, zy, zz, idx] = SE_fkg_expand_all_mex(x,opt);
end

[idx, s] = sort(idx);

%Changes
x = x(s,:);
zx = zx(:,s);
zy = zy(:,s);
zz = zz(:,s);

SE_static.zx = zx;
SE_static.zy = zy;
SE_static.zz = zz;
SE_static.idx = idx;
SE_static.perm = s';
SE_static.iperm(s) = 1:length(s);
