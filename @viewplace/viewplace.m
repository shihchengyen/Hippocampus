function [obj, varargout] = viewplace(varargin)
%@viewplace Constructor function for viewplace class
%   OBJ = viewplace(varargin)
%
%   OBJ = viewplace('auto') attempts to create a viewplace object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on viewplace %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = viewplace('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'FileName','unit*_placefield.mat', 'ObjectLevel','Channel');
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
Args.classname = 'viewplace';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vp';

% To decide the method to create or load the object
[command, robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
	obj = robj;
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

	data.locSI = zeros(dnum,1);
	for didx = 1:dnum
		l = load(dlist(didx).name);
		data.locSI(didx) = l.unit_locSI;
	end
	% set index to keep track of which data goes with which directory
	data.ChannelIndex = [0; dnum];
	data.ArrayIndex = [0; dnum];
	data.SessionIndex = [0; dnum];
	data.DayIndex = [0; dnum];
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
data.locSI = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
