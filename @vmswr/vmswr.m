function [obj, varargout] = vmswr(varargin)
%@vmswr Constructor function for vmswr class
%   OBJ = vmswr(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = vmswr('auto') attempts to create a vmswr object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmswr %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmswr('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'Data',[], 'LowpassFreqs',[150 250], 'LPFOrder',8);
Args.flags = {'Auto','ArgsOnly'};
% The arguments that are critical when loading saved data
Args.DataCheckArgs = {'LowpassFreqs'};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmswr';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vr';

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

if(~isempty(Args.Data))
	data = Args.Data;
	data.analogTime = (0:(data.analogInfo.NumberSamples-1))' ./ data.analogInfo.SampleRate;
	data.numSets = 1;
	% clear Data in Args so it is not saved
	Args.Data = [];
	Args.LowpassFreqs = [data.analogInfo.HighFreqCorner/1000 ...
		data.analogInfo.LowFreqCorner/1000];
		
	% create nptdata so we can inherit from it   
	data.Args = Args; 
	n = nptdata(data.numSets,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);	
else
	rw = rplraw('auto',varargin{:});
	if(~isempty(rw))
%         rlndata = nptRemoveLineNoise(rw.data.analogData, 50, rw.data.analogInfo.SampleRate); % remove 50Hz line noise
		[bpdata,resampleRate] = nptLowPassFilter(rw.data.analogData,rw.data.analogInfo.SampleRate, ...
					Args.LowpassFreqs(1),Args.LowpassFreqs(2));
		data.analogData = bpdata;
        data.analogInfo = rw.data.analogInfo;
		data.analogInfo.SampleRate = resampleRate;
		data.analogInfo.MinVal = min(bpdata);
		data.analogInfo.MaxVal = max(bpdata);
		data.analogInfo.HighFreqCorner = Args.LowpassFreqs(1)*1000;
		data.analogInfo.LowFreqCorner = Args.LowpassFreqs(2)*1000;
		data.analogInfo.NumberSamples = length(bpdata);
		data.analogInfo.HighFreqOrder = Args.LPFOrder;
		data.analogInfo.LowFreqOrder = Args.LPFOrder;
		data.analogInfo.ProbeInfo = strrep(data.analogInfo.ProbeInfo,'raw','lfp');
		data.analogTime = (0:(data.analogInfo.NumberSamples-1))' ./ data.analogInfo.SampleRate;
        
        data.analogRmsData = nptRms(bpdata,10,5,1); % Calculate RMS every 5ms using a moving 10ms window
        data.analogRmsInfo.Swr = nptSwr(data.analogRmsData);
        data.analogRmsInfo.Mean = mean(data.analogRmsData);
        data.analogRmsInfo.Std = std(data.analogRmsData);
        
        rl = rplparallel('auto',varargin{:});
		data = arrangeMarkers(data,rl);
        data.numSets = size(data.trialIndices,1);
        
		% create nptdata so we can inherit from it    
		data.Args = Args; 
		n = nptdata(data.numSets,0,pwd);
		d.data = data;
		obj = class(d,Args.classname,n);
		saveObject(obj,'ArgsC',Args);			
	else
		% create empty object
		obj = createEmptyObject(Args);
	end
end

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
