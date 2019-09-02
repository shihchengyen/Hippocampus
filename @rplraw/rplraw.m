function [obj, varargout] = rplraw(varargin)
%@rplraw Constructor function for rplraw class
%   OBJ = rplraw(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = rplraw('auto') attempts to create a rplraw object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rplraw %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = rplraw('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, 'Channel',1, ...
				'StartMarker',1, 'Data',[]);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'rplraw';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'rw';

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
