function [obj, varargout] = rplview(varargin)
%@rplview Constructor function for rplview class
%   OBJ = rplview(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = rplview('auto') attempts to create a rplview object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rplview %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = rplview('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'rplview';
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
    obj = createObject(Args,modvarargin{:});
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
	ns_status = ns_CloseFile(hFile);
		
	% create nptdata so we can inherit from it    
    data.Args = Args;
    data.numSets = nec;
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
