clear
rng(1);

L = 1;
box = [L L L];
N = 1000;
[x, f] = vector_system(N, box, 3);

M0 = 24; % Set M0, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0 / 12;
%opt.xi = 4;
opt.P = 18;
opt.s = 4;
opt.s0 = 2;
opt.n = ceil(opt.M(1)/10);
opt.rc = 6 / opt.xi;
opt.box = box;
opt.beta = 2.4;
opt.window = 'kaiser';

% parameters for (reference) direct Ewald sum
ED_opt.layers = (opt.M-1)/2;
ED_opt.xi = opt.xi;
ED_opt.box = opt.box;
ED_opt.rc = opt.rc;

% compute FD Ewald sum
idx = 1:10;
ref_fd = SE2P_Stokes_direct_fd_mex(idx,x,f,ED_opt);
ref_self = SE2P_Stokes_direct_self_mex(idx,x,f,ED_opt);
ref_rc = SE2P_Stokes_direct_rsrc_mex(idx,x,f,ED_opt);
ref_k0 = SE2P_Stokes_direct_k0_mex(idx,x,f,ED_opt);
ref = ref_fd + ref_k0;

u = se2p_fourier_space(1:N,x,f,opt);

err = rms(u(idx,:)-ref) / rms(ref)
