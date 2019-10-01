function [obj, varargout] = umaze(varargin)
%@umaze Constructor function for umaze class
%   OBJ = umaze(varargin)
%
%   OBJ = umaze('auto') attempts to create a umaze object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on umaze %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = umaze('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Session', ...
				'GridSteps',5, ...
                'MinObs',5);

Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps', 'MinObs'};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'umaze';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'uma';

% To decide the method to create or load the object
[command,robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

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

function obj = createObject(Args,varargin)

% move to correct directory
[pdir,cwd] = getDataOrder('Session','relative','CDNow');

% example object
dlist = dir('unityfile.mat');
% get entries in directory
dnum = size(dlist,1);

% check if the right conditions were met to create object
if(dnum>0)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% umaze specific calculations start from here  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
	disp('found unity file!');
    
    ufdata = unityfile('auto');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% umaze specific calculations ends here  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% create nptdata so we can inherit from it
	data.numSets = 1;    
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

% these are object specific fields
data.dlist = [];
data.setIndex = [];

% create nptdata so we can inherit from it
% useful fields for most objects
data.numSets = 0;
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
