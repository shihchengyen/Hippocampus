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
%   Example: [as, Args] = vmspiketrain('save','redo');
%
%   Example: vs = ProcessLevel(vmspiketrain,'Levels','Days', ...
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

% try to create umaze and unityfile objects
ufile = unityfile('auto',varargin{:});
um = umaze('auto',varargin{:});

% check if the right conditions were met to create object
if( dnum>0 && ~isempty(um) && ~isempty(ufile) )	
	data.horGridBound = um.data.horGridBound;
	data.vertGridBound = um.data.vertGridBound;
	sTime = um.data.sessionTime;
	sortedGPindinfo = um.data.sortedGPindinfo;
	unityTime = ufile.data.unityTime;

	% load spike train
	st = load(Args.DataFile);
	ts1 = (st.timestamps/1000)';
	% remove any spike times before Unity started and after Unity ended
	ts = ts1(ts1<sTime(end,1) & ts1>unityTime(1));

	% find spikes that occurred during navigation
	% e.g.
%	row		sessionTime			ts		bins	zero_indices
%	1	0		0		7.8217	6.7649	1	1
%	2	7.8217	3.0000	3.6656	6.7917	1	6
%	3	11.4873	4.0000	1.5644	8.1015	2	14
%	4	13.0517	5.0000	0.3790	9.3976	2	20
%	5	13.4307	10.0000	1.1994	16.4411	6	26
%	6	14.6301	0		3.1062	16.4855	6	34
%	7	17.7363	10.0000	2.4343	16.5127	6	40
%	8	20.1706	15.0000	0.4862	16.5400	6	46
%	9	20.6568	14.0000	1.7571	16.8922	6	54
%	10	22.4139	13.0000	0.8260	16.9383	6	61
	% value in bins match up with the values in zero_indices, so we can just find the values
	% in bins that are not members of zero_indices
	[~,~,bins] = histcounts(ts,sTime(:,1));
	% compute the repetition number within each grid position each spike occurs in
	% first find the row in sortedGPindices that corresponds to the occurence of the grid
	% position that the spike occurred in
	[~, spike_index] = ismember(bins,um.data.sortedGPindices);
	% get indices that correspond to non-zero grid positions
	n0bins = spike_index>sortedGPindinfo(1,2);
	% get the grid position for each spike in non-zero grid positions
	spike_gridposition = sTime(bins(n0bins),2);
	% compute the repetition number for each spike by finding the row that corresponds to
	% the start of the appropriate grid position
	rep_num = mod(spike_index(n0bins),sortedGPindinfo(um.data.gp2ind(spike_gridposition),2)) + 1;
	
	% perform histogram on spike times that don't correspond to zero bins using unity time 
	% bins in order to get the closest xy position
	[~,~,ubin] = histcounts(ts(n0bins),ufile.data.unityTime);
	spike_xy = ufile.data.unityData(ubin,3:4);

	% sort according to grid position and compute number of spikes at each position
    [sorted_sgpi,sorted_sgpi_info] = groupdata(spike_gridposition);
    % in order to extract the number of observations for the grid positions
    % that have spikes, we find the indices for the grid positions that
    % have spikes
    lia = ismember(sortedGPindinfo(:,1),sorted_sgpi_info);

	data.spike_xy = {spike_xy};
	data.sorted_sgpi = {sorted_sgpi};
    % we add a 5th column that stores the number of observations for that
    % grid position
	data.sorted_sgpi_info = {[sorted_sgpi_info sortedGPindinfo(lia,4)]};
	data.rep_num = rep_num;
	
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
