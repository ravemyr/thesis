function init()
root = getenv('SE_reset_ROOT');
if isempty(root)
  % Run main initialization
  path = fileparts(mfilename('fullpath')); % path to the directory containing this file
  run([path, '/../init.m']);
  root = getenv('SE_reset_ROOT');
end

addpath([root, '/SE3P_Laplace/src']);
addpath([root, '/util']);
addpath([root, '/bin']);

% TEMPORARY PATHS (TODO FIXME)
%addpath([root, '/../SE_unified.git/SE_kaiser/src']);
%addpath([root, '/../SE_unified.git/util']);
