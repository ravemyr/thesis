function [u varargout] = stokeslet_fourier_space(x, f, opt, pre)

if nargout==2
	[u walltime] = se0p_fourier_space(x, f, opt, pre, @stokeslet_k_scaling);
    varargout{1} = walltime;
else
	u = se0p_fourier_space(x, f, opt, pre, @stokeslet_k_scaling);
end
