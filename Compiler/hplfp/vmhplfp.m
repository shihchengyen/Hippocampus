function vmhplfp(varargin)
% rh = rplhighpass('auto',varargin{:});
% rl = rpllfp('auto',varargin{:});
% rp = rplparallel('auto',varargin{:});
vh = vmhighpass('auto','SaveLevels',2,varargin{:});
vl = vmlfp('auto','SaveLevels',2,varargin{:});
