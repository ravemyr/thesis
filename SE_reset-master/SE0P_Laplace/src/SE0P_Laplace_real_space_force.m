function du = SE0P_Laplace_real_space_force(eval_idx, x, f, opt)
% SE0P_LAPLACE_REAL_SPACE_FORCE
% Compute real-space part of Ewald sum for 0-periodic (free-space)
% electrostatic potential.
%
% du = SE0P_Laplace_real_space_force(eval_idx, x, f, opt)
%   Return force
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
% :returns: **du** -- real-space force

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
c = 2/sqrt(pi)*xi;

N = size(x, 1);
[idx, d] = rangesearch(x, x, rc);

du = zeros(N, 3);
for target=eval_idx
    nnb_idx = idx{target}(2:end);
    nnb_r = d{target}(2:end);
    r = bsxfun(@minus, x(target,:), x(idx{target}(2:end),:));
    u1 = c*exp(-(xi^2*nnb_r.^2));
    u2 = erfc(xi*nnb_r) ./ nnb_r;
    u = f(nnb_idx)' .* (u1+u2) ./ nnb_r.^2;
    du(target,:) = sum(bsxfun(@times, u', r));
end
