function [obj, varargout] = vmhighpass(varargin)
%@vmhighpass Constructor function for vmhighpass class
%   OBJ = vmhighpass(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = vmhighpass('auto') attempts to create a vmhighpass object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmhighpass %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmhighpass('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, 'Raw',0);
Args.flags = {'Auto','ArgsOnly','Raw'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmhighpass';
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
    % whatever needed here
    if(Args.Raw)
    	rh = rplraw('auto',modvarargin{:});
    else
    	rh = rplhighpass('auto',modvarargin{:});
    end
    rl = rplparallel('auto',modvarargin{:});
    if( isempty(rh) | isempty(rl) )
    	obj = createEmptyObject(Args);
    else
	    obj = createObject(rh,rl,Args,modvarargin{:});
	end
end

function obj = createObject(rh,rl,Args,varargin)

data = rh.data;
% compute analogTime to make it easier to plot the data
% replace analogTime in rplhighpass, which is single precision, with double 
% precision variable to avoid plotting errors after 40 or so trials
data.analogTime = (0:(data.analogInfo.NumberSamples-1))' ./ data.analogInfo.SampleRate;
% data = arrangeMarkers(data,rl);
data.markers = rl.data.markers;
data.timeStamps = rl.data.timeStamps;
data.trialIndices = rl.data.trialIndices;
data.numSets = size(data.trialIndices,1);

% create nptdata so we can inherit from it    
data.Args = Args;
n = nptdata(data.numSets,0,pwd);
d.data = data;
obj = class(d,Args.classname,n);
saveObject(obj,'ArgsC',Args);	

function obj = createEmptyObject(Args)

% useful fields for most objects
data.analogTime = [];
data.markers = [];
data.timeStamps = [];
data.trialIndices = [];
data.numSets = 0;

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
