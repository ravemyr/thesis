function du = SE0P_Laplace_direct_full_force(eval_idx, x, f)
% SE0P_Laplace_direct_full_force  Direct evaluation for free-space Laplace
%   du = SE0P_Laplace_direct_full_force(eval_idx, x, f)
%
%   Returns the full free-space Laplace force.
%
%   Parameters:
%   :param eval_idx: index of source locations where force should be evaluated
%   :param x: source locations (N×3)
%   :param f: source charges (N×1)
%
%   :returns: **du** -- full force

N = size(x, 1);
du = zeros(N, 3);
for target=eval_idx
  source = [1:target-1 target+1:N];
  rvec = bsxfun(@minus, x(target,:), x(source,:));
  dist = sqrt(sum(rvec.^2, 2));
  u0 = f(source)./dist.^3;
  du(target,:) = sum(bsxfun(@times, u0, rvec));
end
