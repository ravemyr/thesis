function [x,f] = vector_system(N, box, varargin)
% VECTOR_SYSTEM  Generate particle systems of different dimensions.
%    [x,f] = VECTOR_SYSTEM(N,[L1 L2 L3]) or
%    [x,f] = VECTOR_SYSTEM(N,[L1 L2 L3],1) generates a
%    one-dimensional random system on a box of size [L1 L2 L3].
%
%    [x,f] = VECTOR_SYSTEM(N,[L1 L2 L3],DIM) generates a
%    DIM-dimensional random system.
%
%    The generated system will be charge neutral, unless the
%    function is called as VECTOR_SYSTEM(...,'not_charge_neutral').

% Default values
dim = 1;
charge_neutral = true;

for k=1:numel(varargin)
  if isa(varargin{k}, 'double')
    dim = varargin{k};
  elseif strcmp(varargin{k}, 'not_charge_neutral')
    charge_neutral = false;
  end
end

x = bsxfun(@times, rand(N, 3), box);
f = 1 - 2*rand(N, dim);

if charge_neutral
  f = f - repmat(mean(f,1), N, 1); % neutrality
end
