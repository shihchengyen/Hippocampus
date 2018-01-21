function [obj, varargout] = vmlfp(varargin)
%@vmlfp Constructor function for vmlfp class
%   OBJ = vmlfp(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = vmlfp('auto') attempts to create a vmlfp object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmlfp %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmlfp('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmlfp';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vp';

% To decide the method to create or load the object
command = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
    l = load(Args.matname);
    obj = eval(['l.' Args.matvarname]);
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever is needed here
    rp = rpllfp('auto',modvarargin{:});
    rl = rplparallel('auto',modvarargin{:});
    if( isempty(rp) | isempty(rl) )
    	obj = createEmptyObject(Args);
    else
	    obj = createObject(rp,rl,Args,modvarargin{:});
	end
end

function obj = createObject(rp,rl,Args,varargin)

data = rp.data;
data = arrangeMarkers(data,rl);
data.numSets = size(data.trialIndices,1);
	
% create nptdata so we can inherit from it    
data.Args = Args;
n = nptdata(data.numSets,0,pwd);
d.data = data;
obj = class(d,Args.classname,n);
saveObject(obj,'ArgsC',Args);	

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
