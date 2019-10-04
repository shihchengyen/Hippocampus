function [obj, varargout] = vmspiketrain(varargin)
%@vmspiketrain Constructor function for VMSPIKETRAIN class
%   OBJ = vmspiketrain(varargin)
%
%   OBJ = vmspiketrain('auto') attempts to create a vmspiketrain object by checking for
%   the existence of 'spiketrain.mat' (the name can be changed using the 'DataFile'
%   argument. If the file exists, the directory is saved in the object.
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmspiketrain %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Example [as, Args] = vmspiketrain('save','redo')
%           vs = ProcessLevel(vmspiketrain,'Levels','Days', ...
%              'Include',{'20180810','20180823','20180824','20181101','20181102','20181105'}, ...
%              'Exclude',{'sessioneye'})
%   Dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, 'ObjectLevel', 'Cell', ...
				'DataFile','spiketrain.mat');
Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmspiketrain';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vs';

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
dlist = nptDir(Args.DataFile);
% get entries in directory
dnum = size(dlist,1);

% try to create unitymaze object
um = unitymaze('auto',varargin{:});

% check if the right conditions were met to create object
if(dnum>0 && ~isempty(um))	
	data.horGridBound = um.data.horGridBound;
	data.vertGridBound = um.data.vertGridBound;

	% load spike train
	st = load(Args.DataFile);
	ts = (st.timestamps/1000)';

	% perform histogram on spike times using unity time bins in order to get the closest xy
	% position
	[~,~,ubin] = histcounts(ts,um.data.unityTime);
	% since ts may contain spike times after the unity program has ended, we need to look only
	% at the non-zero values in ubin
	spike_xy = um.data.unityData(ubin(ubin>0),3:4);

	% perform histogram on spike times to see which grid position the above information should
	% go to
	% [n,tsbin] = histc(ts,sTime(:,1));
	[~,~,~,binH,binV] = histcounts2(spike_xy(:,1),spike_xy(:,2),data.horGridBound,data.vertGridBound);
	% compute grid position number
	spike_gridposition = binH + ((binV - 1) * um.data.gridSteps);

	% get the grid position for each spike
	% spike_gridposition = sTime(tsbin(tsbin>0),2);

	% sort according to grid position
	[sorted_spike_gridposition,sorted_sgpi] = sort(spike_gridposition);

	% find the points where grid position changes
	% first set of values in sorted_spike_gridposition is for 0, which is where there is no 
	% grid position, e.g. during cue presentation, etc., so the first change is where we want
	% to start. 
	change_ssgpindex = find(diff(sorted_spike_gridposition));
	change_ssgpindex(end+1) = size(spike_xy,1);

	data.spike_xy = {spike_xy};
	data.sorted_sgpi = {sorted_sgpi};
	data.change_ssgpindex = {change_ssgpindex};
	
	% create nptdata so we can inherit from it
	data.numSets = 1;    
    data.Args = Args;
	n = nptdata(data.numSets,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);
else
	% create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% create nptdata so we can inherit from it
% useful fields for most objects
data.horGridBound = [];
data.vertGridBound = [];
data.spike_xy = {};
data.sorted_sgpi = {};
data.change_sgpindex = {};
data.numSets = 0;
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
