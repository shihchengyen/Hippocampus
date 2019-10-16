function [obj, varargout] = spatialview(varargin)
%@vmplacecell Constructor function for spatialview class
%   OBJ = spatialview(varargin)
%
%   OBJ = spatialview('auto') attempts to create a spatialview object
%   using an unitymaze object, a rplparallel object, and the spiketrain.mat
%   in a Cell directory.
%   
%Example, to create a spatialview object in a cell directory:
%   sv = spatialview('auto');
%
%To look at the fields in the object:
%   sv.data
%
%To create a spatialview object, and save it in the current directory:
%   sv = spatialview('auto','save');
%
%To create a spatialview object with different GridSteps, which will be
%passed to the unitymaze object, and to save both the new unitymaze object
%and the spatialview objects:
%   sv = spatialview('auto','GridSteps',10,'SaveLevels',2);
%
%To create spatialview from a channel directory:
%   sv = ProcessLevel(vmplacecell,'Levels','Channel','save');
%
%To create spatialview from an array directory:
%   sv = ProcessLevel(spatialview,'Levels','Array','save');
%
%dependencies: unitymaze, rplparallel

Args = struct('RedoLevels',0, 'SaveLevels',1, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell', 'RequiredFile','spiketrain.mat', ...
				'MaxTimeDiff',0.002,'MinTrials',5, 'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'UseAllTrials',0, 'AdaptiveSmooth',1, 'FiltLowOcc',0);
Args.flags = {'Auto','ArgsOnly','HPC','FRSIC','UseAllTrials','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'MinTrials','GridSteps','ShuffleLimits', ...
	'NumShuffles','UseAllTrials'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'spatialview';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'sv';

% To decide the method to create or load the object
[command, robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
	obj = robj;
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

if(~isempty(dir(Args.RequiredFile)))
    cwd = pwd;
	% load gaze object
	gz = gaze('auto',varargin{:});
	% load rplparallel object
	rp = rplparallel('auto',varargin{:});
    % load unity maze object
    um = umaze('auto',varargin{:});
    cd(cwd);
	% load spike train file
	spiketrain = load(Args.RequiredFile); 

    % run in Matlab
    data = spatialview_shuffle(gz,rp,um,spiketrain,Args);

	% create nptdata so we can inherit from it    
	data.numSets = 1;
	data.Args = Args;
	n = nptdata(1,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);

else
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.gridSteps = [];
data.meanFRs = [];
data.semFRs = [];
data.SIC = [];
data.SICsh = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);