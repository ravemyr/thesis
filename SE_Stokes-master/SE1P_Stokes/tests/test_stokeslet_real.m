clear
rng(1);

L = 1;
box = [L L L];
N = 10;
[x, f] = vector_system(N, box, 3);

M0 = 2000; % Set M0, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0 / 12000000;
opt.rc = 6 / opt.xi;
opt.box = box;

% parameters for (reference) direct Ewald sum
opt.layers = (opt.M-1)/2;

% compute FD Ewald sum
idx = 1:10;
%ref_fd = SE2P_Stokes_direct_fd_mex(idx,x,f,opt);
ref_self = SE1P_Stokes_direct_self_mex(idx,x,f,opt);
ref_r = SE1P_Stokes_direct_real_mex(idx,x,f,opt);
ref_rc = SE1P_Stokes_direct_rsrc_mex(idx,x,f,opt);
ref_d = SE1P_Stokes_direct_mex(idx,x,f,opt);
%ref_k0 = SE2P_Stokes_direct_k0_mex(idx,x,f,opt);
%ref = ref_fd + ref_k0;

%u = se2p_fourier_space(1:N,x,f,opt);

%err = rms(u(idx,:)-ref) / rms(ref)
ref_diff = ref_r + ref_self - ref_d;

err = rms(ref_r - ref_rc) / rms(ref_r)
err2 = rms(ref_diff) / rms(ref_d)
