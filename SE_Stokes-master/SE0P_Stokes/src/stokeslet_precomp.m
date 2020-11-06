function pre = stokeslet_precomp(opt)

pre = se0p_precomp(opt, @kernels.biharmonic);
