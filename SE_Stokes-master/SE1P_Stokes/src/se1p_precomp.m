function SE_static  = se1p_precomp(x,opt)

verb = false;

% get parameters
[zx zy zz idx] = SE_fgg_expand_all_mex_1p(x,opt);

[~, s] = sort(idx);
x = x(s,:);

[zx zy zz idx] = SE_fgg_expand_all_mex_1p(x,opt);
zs = SE_fgg_base_gaussian_mex_1p(opt);

SE_static.zs=zs;
SE_static.zx=zx;
SE_static.zy=zy;
SE_static.zz=zz;
SE_static.idx=idx;
SE_static.perm=s';
SE_static.iperm(s)=1:length(s);
