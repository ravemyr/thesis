function u = SE0P_Laplace_direct_full(eval_idx, x, f)
% SE0P_Laplace_direct_full  Direct evaluation for free-space Laplace
%   u = SE0P_Laplace_direct_full(eval_idx, x, f)
%
%   Returns the full free-space Laplace potential.
%
%   Parameters:
%   :param eval_idx: index of source locations where potential should be evaluated
%   :param x: source locations (N×3)
%   :param f: source charges (N×1)
%
%   :returns: **u** -- full potential

N = size(x, 1);
u = zeros(N, 1);
for target=eval_idx
  source = [1:target-1 target+1:N];
  rvec = bsxfun(@minus, x(target,:), x(source,:));
  dist = sqrt(sum(rvec.^2, 2));
  u(target) = sum(f(source)./dist);
end
