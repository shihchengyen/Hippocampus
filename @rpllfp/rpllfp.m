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

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, 'Channel',1,'StartMarker',1);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.UnimportantArgs = {'RedoLevels','SaveLevels'};                            

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
    obj = eval(['l.' Args.matvarname]);
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,varargin{:});
end

function obj = createObject(Args,varargin)

% look for ripple file
dfile = dir('*.ns2');
dnum = size(dfile,1);

% check if the right conditions were met to create object
if(dnum>0)
	% this is a valid object
	if(dnum>1)
		% there is more than 1 ns2 file in the current directory
		% print a warning that only the 1st file will be used to create the object
		warning('More than 1 ns2 file found')
		warning(['Creating object only from ' dfile(1).name])
	end
	
	% these are object specific fields
	data.lfpfname = dfile(1).name;
	% open the file, and read the markers
	[ns_status, hFile] = ns_OpenFile(dfile(1).name); 
	entityID = find(cellfun(@strcmpi, {hFile.Entity.Reason},...
		repmat({'Parallel Input'}, size({hFile.Entity.Reason})))); % Note: original 'Digital Input'
	% Extract channel info
	[ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID(end));
	% Get events and time stamps
	numCount = entityInfo.ItemCount;
	ddata = NaN(1, numCount); timeStamps = NaN(1, numCount); sz = NaN(1, numCount);
	for i = 1:numCount
		[~, timeStamps(i), ddata(i), dataSize(i)] = ns_GetEventData(hFile, entityID, i);
	end 
	% get LFP data
	EntityIndices = find([hFile.Entity(:).ElectrodeID] == Args.Channel);         
	for i = 1:length(EntityIndices)       
		fileTypeNum = hFile.Entity(EntityIndices(i)).FileType;
		fileType = hFile.FileInfo(fileTypeNum).Type;            
		if strcmp('ns2', fileType); entityID = EntityIndices(i); break; end 
	end
	% Extract channel info
	% analog info contains things like range and sampling rate	
	[ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID(end));     
	TimeStamps = hFile.FileInfo(hFile.Entity(entityID).FileType).TimeStamps;
	numSamples = sum(TimeStamps(:,end));
	analogInputData = zeros(1,numSamples);
	startIndex = 1; 
	indexCount = TimeStamps(2,1);
	for i = 1:size(TimeStamps,2)                
	    [~, ~, tempData] = ns_GetAnalogData(hFile, entityID, startIndex, indexCount);
	    dataRange = TimeStamps(1,i) + (1:TimeStamps(2,i));
	    analogInputData(dataRange) = tempData';
	    clear tempData
	    if i ~= size(TimeStamps,2) 
	        startIndex = startIndex + TimeStamps(2,i);
	        indexCount = TimeStamps(2,i+1);
	    end
	end
	data.sampleRate = analogInfo.SampleRate;
	data.ltime = (0:numSamples-1)' ./ data.sampleRate;
	data.lfp = analogInputData;
	data.markers = ddata;
	% get start time for each trial
	data.trialIndices = floor(reshape(timeStamps*data.sampleRate,2,[])');
	data.numSets = length(data.trialIndices);
	ns_status = ns_CloseFile(hFile);
		
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
