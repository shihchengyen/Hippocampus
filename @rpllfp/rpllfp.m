function [obj, varargout] = rpllfp(varargin)
%@rpllfp Constructor function for rpllfp class
%   OBJ = rpllfp(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = rpllfp('auto') attempts to create a rpllfp object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rpllfp %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = rpllfp('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'Data',[], 'LowpassFreqs',[1 150], 'LPFOrder',8, ...
				'ResampleRate', 1000);
Args.flags = {'Auto','ArgsOnly'};
% The arguments that are critical when loading saved data
Args.DataCheckArgs = {'LowpassFreqs'};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'rpllfp';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'df';

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
    oname = fieldnames(l);
    tmpobj = eval(['l.' oname{:}]);
    if(isstruct(tmpobj))
        obj = createEmptyObject(Args);
        for fn = fieldnames(tmpobj)'    %enumerat fields
            if(~strcmp(fn,'nptdata'))
                try
                    obj.(fn{1}) = tmpobj.(fn{1});   %and copy
                catch
                    warning('Could not copy field %s', fn{1});
                end
            else
               obj = set(obj,'Number',obj.data.numSets);
               obj = set(obj,'SessionDirs',tmpobj.nptdata.sessiondirs);
            end
        end
    elseif(isa(tmpobj,Args.classname))
        obj = tmpobj;
    end
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

if(~isempty(Args.Data))
	data = Args.Data;
	resampleRate = Args.ResampleRate;
	data.numSets = 1;
	% clear Data in Args so it is not saved
	Args.Data = [];
	Args.LowpassFreqs = [data.analogInfo.HighFreqCorner/resampleRate ...
		data.analogInfo.LowFreqCorner/resampleRate];
		
	% create nptdata so we can inherit from it   
	data.Args = Args; 
	n = nptdata(data.numSets,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);	
else
	rw = rplraw('auto',varargin{:});
	if(~isempty(rw))
		% analogData in rplraw should be in single precision format, so we have to
		% convert to double to avoid errors in nptLowPassFilter
        [lpdata,resampleRate] = nptLowPassFilter(double(rw.data.analogData),rw.data.analogInfo.SampleRate, ...
					Args.LowpassFreqs(1),Args.LowpassFreqs(2));
		% convert to single precision float to save disk space, and to make loading the files faster
		data.analogData = single(lpdata);
		data.analogInfo = rw.data.analogInfo;
		data.analogInfo.SampleRate = resampleRate;
		data.analogInfo.MinVal = min(lpdata);
		data.analogInfo.MaxVal = max(lpdata);
		data.analogInfo.HighFreqCorner = Args.LowpassFreqs(1)*resampleRate;
		data.analogInfo.LowFreqCorner = Args.LowpassFreqs(2)*resampleRate;
		data.analogInfo.NumberSamples = length(lpdata);
		data.analogInfo.HighFreqOrder = Args.LPFOrder;
		data.analogInfo.LowFreqOrder = Args.LPFOrder;
		data.analogInfo.ProbeInfo = strrep(data.analogInfo.ProbeInfo,'raw','lfp');
		data.numSets = 1;
		
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
