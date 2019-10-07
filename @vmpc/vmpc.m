function [obj, varargout] = vmpc(varargin)
%@vmpc Constructor function for vmpc class
%   OBJ = vmpc(varargin)
%
%   OBJ = vmpc('auto') attempts to create a vmpc object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmpc %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmpc('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell', 'RequiredFile','spiketrain.mat', ...
				'MinTrials',5, 'GridSteps',5, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4,'AdaptiveSmooth',1, 'FiltLowOcc',0);
            
Args.flags = {'Auto','ArgsOnly','FRSIC','UseAllTrials','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'MaxTimeDiff','MinTrials','GridSteps','ShuffleLimits', ...
	'NumShuffles'};                           

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmpc';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vmpc';

% To decide the method to create or load the object
[command,robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
    % l = load(Args.matname);
    % obj = eval(['l.' Args.matvarname]);
	obj = robj;
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

% example object
dlist = nptDir;
% get entries in directory
dnum = size(dlist,1);

% check if the right conditions were met to create object
if(~isempty(dir(Args.RequiredFile)))
	
	um = umaze('auto',varargin{:});
	rp = rplparallel('auto',varargin{:});
	spiketrain = load(Args.RequiredFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spiketimes = spiketrain.timestamps/1000; % now in seconds
    maxTime = rp.data.timeStamps(end,3);
    tShifts = [0 ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
    full_arr = repmat(spiketimes, Args.NumShuffles+1, 1);
    full_arr = full_arr + tShifts';
    keepers = length(spiketimes) - sum(full_arr>maxTime, 2);
    for row = 1:size(full_arr,1)
        if keepers(row) < length(spiketimes) 
            full_arr(row,:) = [full_arr(row,1+keepers:end)-maxTime-1 full_arr(row,1:keepers)];
        end
    end
    
    % calculating proportion of occupied time in each grid position across
    % entire session.
    
    gpsessiondur = sum(um.data.gpDurations, 2);
    Pi = gpsessiondur/sum(gpsessiondur);

    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% create nptdata so we can inherit from it    
	data.numSets = 1;
	data.Args = Args;
	n = nptdata(1,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);
else
	% create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% these are object specific fields
data.dlist = [];
data.setIndex = [];

% create nptdata so we can inherit from it
% useful fields for most objects
data.numSets = 0;
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
