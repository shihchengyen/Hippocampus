function [obj, varargout] = rplsplit(varargin)
%@rplsplit Constructor function for rplsplit class
%   OBJ = rplsplit(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = rplsplit('auto') attempts to create a rplsplit object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rplsplit %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = rplsplit('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, 'SkipRaw', 0, 'SkipLFP', 0, ...
				'SkipParallel',0, 'Channels',[], 'ChannelsPerArray',32);
Args.flags = {'Auto','ArgsOnly','SkipRaw','SkipLFP','SkipParallel'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'rplsplit';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'rs';

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
dfile = dir('*.ns5');
dnum = size(dfile,1);

% check if the right conditions were met to create object
if(dnum>0)
	% this is a valid object
	if(dnum>1)
		% there is more than 1 ns5 file in the current directory
		% print a warning that only the 1st file will be used to create the object
		warning('More than 1 ns5 file found')
		warning(['Creating object only from ' dfile(1).name])
	end
	
	% these are object specific fields
	data.rawfname = dfile(1).name;
	% open the file, and read the information
	[ns_status, hFile] = ns_OpenFile(dfile(1).name); 
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
	% get number of EntityCount
	nec = nsFileInfo.EntityCount;
	% go through and create the appropriate subdirectory for entity
	for ni = 1:nec
		% get info on the type of entity
		[ns_status, nsEI] = ns_GetEntityInfo(hFile, ni);
		numSamples = nsEI.ItemCount;
		if(nsEI.EntityType==2)
			% get label
			eLabel = nsEI.EntityLabel;
			% check if it is raw data
			if(~isempty(strfind(eLabel,'raw')))
				chan_num = sscanf(eLabel,'raw %d');
				b_raw = 1;
				b_lfp = 0;
			elseif(~isempty(strfind(eLabel,'lfp')))
				chan_num = sscanf(eLabel,'lfp %d');
				b_lfp = 1;
				b_raw = 0;
			end
			% check if we should process this channel
			chanArgs = isempty(Args.Channels);
			if( chanArgs | (~chanArgs && ~isempty(find(Args.Channels==chan_num)) ) )
				if( (b_raw * ~Args.SkipRaw) | (b_lfp * ~Args.SkipLFP) )
					% entity is raw data, so create a channel directory
					% read data
					[ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, ni);
					[ns_RESULT, ~, analogData] = ns_GetAnalogData(hFile, ni, 1, numSamples);
					% check how channels are arranged with respect to arrays
					array_num = floor((chan_num-1)/Args.ChannelsPerArray)+1;
					array_dir = sprintf('array%02d',array_num);
					cwd = pwd;
					nptMkDir(array_dir);
					cd(array_dir);
					chan_dir = sprintf('channel%03d',rem(chan_num-1,Args.ChannelsPerArray)+1);
					nptMkDir(chan_dir);
					cd(chan_dir);
					tData.analogInfo = analogInfo;
					tData.analogInfo.NumberSamples = numSamples;
					tData.analogData = analogData;
					% create and save obj
					if(b_raw)
						rplraw('auto','Data',tData,varargin{:});				
					else
						rpllfp('auto','Data',tData,varargin{:});
					end
					cd(cwd);
				end % if( (b_raw * ~Args.SkipRaw) | (b_lfp * ~Args.SkipLFP) )
			end % if( chanArgs && (~chanArgs && ~isempty(find(Args.Channels==chan_num)) )
		elseif(nsEI.EntityType==1 && ~Args.SkipParallel) % if(nsEI.EntityType==2)
			% parallel input
			% Get events and time stamps
			ddata = NaN(1, numSamples); timeStamps = NaN(1, numSamples);
			for i = 1:numSamples
				[~, timeStamps(i), ddata(i)] = ns_GetEventData(hFile, ni, i);
			end 
			tData.markers = ddata;
			tData.timeStamps = timeStamps;
			% create and save obj
			rplparallel('auto','Data',tData,varargin{:});
		end % if(nsEI.EntityType==2)
	end
	ns_status = ns_CloseFile(hFile);
		
	% create nptdata so we can inherit from it    
    data.Args = Args;
    data.numSets = 1;
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
