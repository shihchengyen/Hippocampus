function [obj, varargout] = rplparallel(varargin)
%@rplparallel Constructor function for rplparallel class
%   OBJ = rplparallel(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = rplparallel('auto') attempts to create a rplparallel object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rplparallel %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = rplparallel('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'Data',[], 'ObjectLevel','Session');
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'rplparallel';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'df';

% change to the right directory
% [p,cwd] = getDataOrder('session','CDNow');

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
    % l = load(Args.matname);
    % obj = eval(['l.' Args.matvarname]);
	obj = robj;
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,modvarargin{:});
end

% change back to previous directory
% cd(cwd)

function obj = createObject(Args,varargin)

if(~isempty(Args.Data))
	data = Args.Data;
	% call function to figure out the data indices for different trials,
	% which will create the markers, timeStamps, and trialIndices fields
    try 
        data = arrangeMarkers(data);
    catch
        fprintf('Problem with ArrangeMarkers\n');
    end
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

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;

% these are object specific fields
data.markers = [];
data.timeStamps = [];
data.SampleRate = [];
data.trialIndices = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
