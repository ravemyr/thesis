function u = SE0P_Laplace_real_space(eval_idx, x, f, opt)
% SE0P_LAPLACE_REAL_SPACE
% Compute real-space part of Ewald sum for 0-periodic (free-space)
% electrostatic potential.
%
% u = SE0P_Laplace_real_space(eval_idx, x, f, opt)
%   Return potential
%
% Parameters:
% :param eval_idx: index of source locations where potential should be evaluated
% :param x: source locations (N×3)
% :param f: source charges (N×1)
% :param opt: Ewald options:
% :param opt.xi:              Ewald parameter (required)
% :param opt.rc:              Cutoff radius (required)
% :param opt.eval_x:          External points to evaluate potential in (Nex×3)
%
% :returns: **u** -- real-space potential

assert(isfield(opt, 'xi'), 'Ewald parameter xi must be given in opt struct');
assert(isfield(opt, 'rc'), 'Cutoff radius rc must be given in opt struct');

% External evaluation points
if isfield(opt,'eval_x') && numel(opt.eval_x) > 0
  opt.eval_external = true;
  fprintf('WARNING: External evaluation points are not implemented\n');
else
  opt.eval_external = false;
  opt.eval_x = [];
end
% TODO: external evaluation points are not implemented yet!

xi = opt.xi;
rc = opt.rc;

N = size(x, 1);
[idx, d] = rangesearch(x, x, rc);

u = zeros(N, 1);
for target=eval_idx
    nnb_idx = idx{target}(2:end);
    nnb_r = d{target}(2:end);
    u(target) = sum(f(nnb_idx)' .* erfc(xi*nnb_r)./nnb_r);
end
