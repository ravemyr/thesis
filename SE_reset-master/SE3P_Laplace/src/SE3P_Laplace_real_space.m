function [pot, force] = SE3P_Laplace_real_space(eval_idx, x, f, opt)
% SE3P_LAPLACE_REAL_SPACE
% Compute real-space part of Ewald sum for 3-periodic electrostatic potential.
%
% pot = SE3P_Laplace_real_space(eval_idx, x, f, opt)
%   Return potential
%
% [pot, force] = SE3P_Laplace_real_space(...)
%   Return potential and force
%
% Unless opt.force is true, the force output will be empty.
%
% Parameters:
% :param eval_idx: index of source locations where potential should be evaluated
% :param x: source locations (N×3)
% :param f: source charges (N×1)
% :param opt: Ewald options:
% :param opt.potential:     Compute potential (default: true)
% :param opt.force:         Compute force (default: false)
% :param opt.box:           Size of periodic cell [L1,L2,L3] (required)
% :param opt.xi:            Ewald parameter (required)
% :param opt.rc:            Cutoff radius (required)
% :param opt.eval_x:        External points to evaluate potential in (Nex×3)
%
% :returns: **pot** -- real-space potential
% :returns: **force** -- real-space force

assert(isfield(opt, 'box'), 'cell size box must be given in opt struct');
assert(isfield(opt, 'xi'), 'Ewald parameter xi must be given in opt struct');
assert(isfield(opt, 'rc'), 'Cutoff radius rc must be given in opt struct');

if ~isfield(opt,'potential'), opt.potential = true; end
if ~isfield(opt,'force'), opt.force = false; end

% External evaluation points
if isfield(opt,'eval_x') && numel(opt.eval_x) > 0
  opt.eval_external = true;
  fprintf('WARNING: External evaluation points are not implemented\n');
else
  opt.eval_external = false;
  opt.eval_x = [];
end
% TODO: external evaluation points are not implemented yet!

pot = []; force = [];
if opt.potential
  pot = SE3P_Laplace_real_rc_cell_mex(x, f, opt.rc, opt.xi, opt.box);
  pot = pot(eval_idx,:);
end
if opt.force
  force = SE3P_Laplace_real_rc_cell_force_mex(x, f, opt.rc, opt.xi, opt.box);
  force = force(eval_idx,:);
end
