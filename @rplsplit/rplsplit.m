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

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
	'SkipRaw',0, 'SkipLFP',0, 'SkipParallel',0, 'SkipAnalog',0, ...
	'Channels',[], 'ChannelsPerArray',32, 'HPCCmd','', ...
    'HPCInputFilename','rsData.mat');
Args.flags = {'Auto','ArgsOnly','SkipRaw','SkipLFP','SkipParallel','SkipAnalog'};
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
	
		% save Args and modvarargin so that compiled program can retrieve them
		% look for ripple file
		save(Args.HPCInputFilename,'Args','modvarargin','dfile');
		if(~isempty(Args.HPCCmd))
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
			rsCreateObject();
		end
	
		% create nptdata so we can inherit from it    
		data.Args = Args;
		data.rawfname = dfile(1).name;
		data.numSets = 1;
		n = nptdata(data.numSets,0,pwd);
		d.data = data;
		obj = class(d,Args.classname,n);
		saveObject(obj,'ArgsC',Args);	
		
	else  % if(dnum>0)
		% create empty object
		obj = createEmptyObject(Args);
	end
end

function obj = createEmptyObject(Args)

% create nptdata so we can inherit from it
data.Args = Args;
data.rawfname = '';
data.numSets = 0;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
