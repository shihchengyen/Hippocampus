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

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
    'OldMarkerFormat',0,'OldMarkerFormat2',0);
Args.flags = {'Auto','ArgsOnly','OldMarkerFormat','OldMarkerFormat2'};
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
    rh = rplhighpass('auto',varargin{:});
    rl = rplparallel('auto',varargin{:});
    if( isempty(rh) | isempty(rl) )
    	obj = createEmptyObject(Args);
    else
	    obj = createObject(rh,rl,Args,modvarargin{:});
	end
end

function obj = createObject(rh,rl,Args,varargin)

data = rh.data;
if(Args.OldMarkerFormat)
    data.markers = reshape(rl.data.markers,2,[])';
    % get start time for each trial
    data.trialIndices = floor(reshape(rl.data.timeStamps*data.analogInfo.SampleRate,2,[])');
elseif(Args.OldMarkerFormat2)
    rawMarkers = rl.data.markers;
    % format is marker followed by 0, so we will reshape into 2 rows
    % to make it easy to remove the 0's
    rm1 = reshape(rawMarkers,2,[]);
    % remove the 2nd row, which is all 0's
    % format is 1 marker for start of trial, and another for either
    % correct or incorrect, so reshape into 2 columns
    data.markers = reshape(rm1(1,:),2,[])';
    % get start time for each trial
    rtime = rl.data.timeStamps;
    % reshape to remove timestamps for 0's
    rt1 = reshape(rtime,2,[]);
    % select 1st row to remove 0's, and the reshape timestamps into 
    % 2 columns
    data.timeStamps = reshape(rt1(1,:),2,[])';
    % compute the data indices corresponding to the marker timestamps
    data.trialIndices = floor(data.timeStamps*data.analogInfo.SampleRate);
else
    rawMarkers = rl.data.markers;
    % format is marker followed by 0, so we will reshape into 2 rows
    % to make it easy to remove the 0's
    % remove the 2nd row, which is all 0's
    % format is: 1 for cue onset/start of trial; 2 for cue offset; 
    % and 3 or 4 for reward or error/timeout
    % so we will reshape markers into 3 columns
    rm1 = reshape(rawMarkers,6,[]);
    data.markers = rm1([1 3 5],:)';
    % get start time for each trial
    rtime = rl.data.timeStamps;
    rt1 = reshape(rtime,6,[]);
    data.timeStamps = rt1([1 3 5],:)';    
    data.trialIndices = floor(data.timeStamps*data.analogInfo.SampleRate);
end
data.numSets = length(data.trialIndices);
	
% create nptdata so we can inherit from it    
data.Args = Args;
n = nptdata(data.numSets,0,pwd);
d.data = data;
obj = class(d,Args.classname,n);
saveObject(obj,'ArgsC',Args);	

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.setNames = '';

% these are object specific fields
data.dlist = [];
data.setIndex = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
