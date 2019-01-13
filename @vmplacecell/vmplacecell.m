function [obj, varargout] = vmplacecell(varargin)
%@vmplacecell Constructor function for vmplacecell class
%   OBJ = vmplacecell(varargin)
%
%   OBJ = vmplacecell('auto') attempts to create a vmplacecell object
%   using an unitymaze object, a rplparallel object, and the spiketrain.mat
%   in a Cell directory.
%   
%Example, to create a vmplacecell object in a cell directory:
%   vpc = vmplacecell('auto');
%
%To look at the fields in the object:
%   vpc.data
%
%To create a vmplacecell object, and save it in the current directory:
%   vpc = vmplacecell('auto','save');
%
%To create a vmplacecell object with different GridSteps, which will be
%passed to the unitymaze object, and to save both the new unitymaze object
%and the vmplacecell objects:
%   vpc = vmplacecell('auto','GridSteps',10,'SaveLevels',2);
%
%To create vmplacecells from a channel directory:
%   vpc = ProcessLevel(vmplacecell,'Levels','Channel','save');
%
%To create vmplacecells from an array directory:
%   vpc = ProcessLevel(vmplacecell,'Levels','Array','save');
%
%dependencies: unitymaze, rplparallel

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell', 'RequiredFile','spiketrain.mat', ...
				'HPC',0, 'HPCInputFilename','vpData.mat', ...
				'HPCCmd','condor_submit placeselect_submit.txt', ...
				'MaxTimeDiff',0.002, 'MinTrials',5);
Args.flags = {'Auto','ArgsOnly','HPC'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'MinTrials'};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmplacecell';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vp';

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
	% load unitymaze object
	um = unitymaze('auto',varargin{:});
	% load rplparallel object
	rp = rplparallel('auto',varargin{:});
	% save(Args.HPCInputFilename,'Args','um','rp','varargin');

	if(Args.HPC)
		% use HTCondor to run rplsplitcreateObject
		[s,w] = system(Args.HPCCmd);
		% make sure shell script ran without problems
		if(s~=0)
			error([Args.classname ': Error splitting files!'])
		else
			% display output
			fprintf('%s\n',w);
		end	
	else
		% run in Matlab
		% data = placeselect(pwd);
		data = placeselect(um,rp,Args);
	end

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

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
