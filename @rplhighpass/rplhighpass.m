function [obj, varargout] = rplhighpass(varargin)
%@rplhighpass Constructor function for rplhighpass class
%   OBJ = rplhighpass(varargin) extracts highpass signals from a RIPPLE recording
%
%   OBJ = rplhighpass('auto') attempts to create a rplhighpass object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rplhighpass %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = rplhighpass('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'Data',[], 'HighpassFreqs',[500 7500], 'HPFOrder',8);
Args.flags = {'Auto','ArgsOnly'};
% The arguments that are critical when loading saved data
Args.DataCheckArgs = {'HighpassFreqs'};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'rplhighpass';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'rh';

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
	data.numSets = 1;
	% clear Data in Args so it is not saved
	Args.Data = [];
	Args.HighpassFreqs = [data.analogInfo.HighFreqCorner/1000 ...
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
		% Matlab Compiler appears to be calling the wrong subsref when accessing 
		% fields inside the data structure inside objects, so we are going to
		% have to get the data structure first and then access the fields inside it
		rwdata = rw.data;
		% analogData in rplraw should be in single precision format, so we have to
		% convert to double to avoid errors in nptHighPassFilter
		hpdata = nptHighPassFilter(double(rwdata.analogData),rwdata.analogInfo.SampleRate, ...
					Args.HighpassFreqs(1),Args.HighpassFreqs(2));
		% convert to single precision float to save disk space, and to make loading the files faster
		data.analogData = single(hpdata);
		data.analogInfo = rwdata.analogInfo;
		data.analogInfo.MinVal = min(hpdata);
		data.analogInfo.MaxVal = max(hpdata);
		data.analogInfo.HighFreqCorner = Args.HighpassFreqs(1)*1000;
		data.analogInfo.LowFreqCorner = Args.HighpassFreqs(2)*1000;
		data.analogInfo.NumberSamples = length(hpdata);
		data.analogInfo.HighFreqOrder = Args.HPFOrder;
		data.analogInfo.LowFreqOrder = Args.HPFOrder;
		data.analogInfo.ProbeInfo = strrep(data.analogInfo.ProbeInfo,'raw','hp');
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
