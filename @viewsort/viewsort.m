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
                'HMMFile','spike_templates.hdf5','HMMNoise','cinv', ...
                'SavedFileName','spiketrain.mat', 'sd', 4, 'Saved', 0);
Args.flags = {'Auto','ArgsOnly','Saved'};
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
    
    if (Args.Saved==1)
    	if(exist('cell01', 'dir')==7)
			% load saved spiketrains
			cwd = pwd;
			data.spikeForms = [];
			celldir = 'cell0';
			cellname = 'cell01';
			i = 1;
       	    while exist(cellname, 'dir')
				cd(cellname)
				s = load(Args.SavedFileName);
            	data.spikeForms(end+1,:) = s.spikeForm;
            	i = i + 1;
            	cellname = strcat(celldir, int2str(i));
            	cd(cwd)
			end  % while exist(cellname, 'dir')
        
	    	[sf1,sf2] = size(data.spikeForms);
		else  % if(exist('cell01', 'dir')==7)
			data.spikeForms = [];
			sf1 = 0;
		end  % if(exist('cell01', 'dir')==7)        
    else %if (Args.Saved==1) && (exist('cell01', 'dir')==7)
        % load waveforms
        l = load(Args.FileName);
        [sf1,sf2,sf3] = size(l.spikeForms);
        data.spikeForms = squeeze(l.spikeForms);
        if(sf1==1)
            data.spikeForms = data.spikeForms';
        end
        
        % get and plot ISI coefficient of variation for waveform
        data.coeffV_ISI = zeros(sf1,1);
        if (isfield(l,'mlseq') > 0) % skip incomplete files
            spikeSort = l.mlseq;
            for ind = 1:sf1
                [M,I] = min(data.spikeForms(ind,:)); % get index of peak
                spikeTimes = find(spikeSort(ind,:) == I); % get index of peaks
                spikeTimes = spikeTimes/30; 
                spike_ISI = diff(spikeTimes);
                data.coeffV_ISI(ind) = std(spike_ISI)/mean(spike_ISI);
            end
        else
            data.coeffV_ISI = repmat(NaN,sf1,1);
        end
        
        if (sf1>1)
            % get spike similarities (correlation coefficient)
            % number of possible pair-wise comparisons
            perms = nchoosek(1:size(data.spikeForms,1),2); 
            perms_size = size(perms,1);
            % create memory
            corrcoefs = zeros(perms_size,1);
            p2pdiffs = corrcoefs;
            for a = 1:perms_size
            	% compute correlation coefficient
            	cc = corrcoef(data.spikeForms(perms(a,1),:),data.spikeForms(perms(a,2),:));
                corrcoefs(a) = cc(2,1);
                % peak-to-peak amplitude of first waveform in comparison
                amp1 = abs(min(data.spikeForms(perms(a,1),:)))+abs(max(data.spikeForms(perms(a,1),:))); 
                % peak-to-peak amplitude of second waveform in comparison
                amp2 = abs(min(data.spikeForms(perms(a,2),:)))+abs(max(data.spikeForms(perms(a,2),:))); 
                % abs difference in peak-to-peak amplitude
                p2pdiffs(a) = round(abs(amp2-amp1),1); 
            end
            data.spikesim = corrcoefs;
            data.spikesp2pdiffs = p2pdiffs;
            numperms = nchoosek(sf1,2);
        else
            data.spikesim = NaN;
            data.spikesp2pdiffs = NaN;
            numperms = 1;
        end
        data.spikesimIndex = [0; numperms];
    end %if (Args.Saved==1) && (exist('cell01', 'dir')==7)
    
    % get noise data from spike templates
    cwd = pwd;
    cd(Args.HMMDir);
    data.Noise = Args.sd*sqrt(1/(hdf5read(Args.HMMFile,Args.HMMNoise)));
    cd(cwd);
    
	% set index to keep track of which data goes with which directory
    data.ChannelIndex = [0; sf1];
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
