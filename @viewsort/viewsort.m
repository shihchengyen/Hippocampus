function [obj, varargout] = viewsort(varargin)
%@viewsort Constructor function for viewsort class
%   OBJ = viewsort(varargin)
%
%   OBJ = viewsort('auto') attempts to create a viewsort object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on viewsort %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = viewsort('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'FileName','hmmsort.mat','HMMDir','hmmsort', ...
                'HMMFile','spike_templates.hdf5','HMMNoise','cinv');
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
Args.classname = 'viewsort';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vs';

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
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

% example object
dlist = dir(Args.FileName);
% get entries in directory
dnum = size(dlist,1);

% check if the right conditions were met to create object
if(dnum>0)
	% these are fields that are useful for most objects
	data.numChannels = 1;
	% this is a valid object

    l = load(Args.FileName);	
	[sf1,sf2,sf3] = size(l.spikeForms);
	data.spikeForms = squeeze(l.spikeForms);
	if(sf1==1)
		data.spikeForms = data.spikeForms';
    end
    
    % get noise data from spike templates
    cwd = pwd;
    cd(Args.HMMDir);
    data.Noise = 1/(hdf5read(Args.HMMFile,Args.HMMNoise));
    cd(cwd);
    
    if (isfield(l,'mlseq') > 0) % skip incomplete files
        % get and plot ISI coefficient of variation for waveform
        spikeSort = l.mlseq;
        for ind = 1:sf1
            [M,I] = min(data.spikeForms(ind,:)); % get index of peak
            spikeTimes = find(spikeSort(ind,:) == I); % get index of peaks
            spikeTimes = spikeTimes/30; spike_ISI = diff(spikeTimes);
            data.coeffV_ISI(ind) = std(spike_ISI)/mean(spike_ISI);
        end
    else
        data.coeffV_ISI(1:sf1) = NaN;
    end
    data.coeffV_ISI = data.coeffV_ISI';
    
    if (sf1>1)
        % get spike similarities (dot product)
        perms = nchoosek(1:size(data.spikeForms,1),2); % number of possible pair-wise comparisons
        for a = 1:size(perms,1)
            perms(a,3) = (dot(data.spikeForms(perms(a,1),:),data.spikeForms(perms(a,2),:)))/(norm(data.spikeForms(perms(a,1),:))*norm(data.spikeForms(perms(a,2),:))); % dot product divided by the magnitude of each vector
            perms(a,3) = round(perms(a,3),2); % round to 2 decimal places
            % amp1 = abs(min(data.spikeForms(perms(a,1),:)))+abs(max(data.spikeForms(perms(a,1),:))); % peak-to-peak amplitude of first waveform in comparison
            % amp2 = abs(min(data.spikeForms(perms(a,2),:)))+abs(max(data.spikeForms(perms(a,2),:))); % peak-to-peak amplitude of second waveform in comparison
            % perms(a,4) = round(abs(amp2-amp1),1); % abs difference in peak-to-peak amplitude
        end
        data.spikesim = perms(:,3);
        numperms = nchoosek(sf1,2);
    else
        data.spikesim = NaN;
        numperms = 1;
    end
    
	% set index to keep track of which data goes with which directory
    data.ChannelIndex = [0; sf1];
    data.spikesimIndex = [0; numperms];
	data.ArrayIndex = [0; 1];
	data.SessionIndex = [0; 1];
	data.DayIndex = [0; 1];
	% get channel string
	[data.arrstr(1,:), chnstr] = nptFileParts(pwd);
	% get array string
	[data.sesstr(1,:), arrstr] = nptFileParts(data.arrstr(1,:));
	% get session string
	[data.daystr(1,:), sesstr] = nptFileParts(data.sesstr(1,:));
	% get day string
	% [p4, data.daystr(1,:)] = nptFileParts(p3);
    data.Args = Args;
		
	% create nptdata so we can inherit from it
    data.Args = Args;
	n = nptdata(data.numChannels,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);
else
	% create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numChannels = 0;
data.numArrays = 0;
data.numSessions = 0;
data.numDays = 0;
data.spikeForms = [];


% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
