clear
rng(1);

L = 1;
box = [L L L];
N = 10;
[x, f] = vector_system(N, box, 3);

M0 = 24; % Set M0, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.P = 32;
opt.s = 4;
opt.s0 = 2.5;
opt.n = ceil(opt.M(1)/10);
opt.rc = 6 / opt.xi;
opt.box = box;
opt.beta = 2.4;
opt.window = 'kaiser';

% parameters for (reference) direct Ewald sum
opt.layers = 1000000;%(opt.M-1)/2;

% compute FD Ewald sum
idx = 1:10;
%ref_fd = SE1P_Stokes_direct_fd_mex(idx,x,f,opt);
ref_self = SE1P_Stokes_direct_self_mex(idx,x,f,opt);
%ref_r = SE1P_Stokes_direct_real_mex(idx,x,f,opt);
ref_rc = SE1P_Stokes_direct_rsrc_mex(idx,x,f,opt);
ref_d = SE1P_Stokes_direct_mex(idx,x,f,opt);
ref_fd = ref_d - ref_rc - ref_self;
%ref_k0 = SE2P_Stokes_direct_k0_mex(idx,x,f,opt);
ref = ref_fd;

u = se1p_fourier_space(1:N,x,f,opt);
err = rms(u(idx,:)-ref) / rms(ref)
