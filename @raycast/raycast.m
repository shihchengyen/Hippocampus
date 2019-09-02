function [obj, varargout] = raycast(varargin)
%@raycast Constructor function for raycast class
%   OBJ = raycast(varargin)
%
%   OBJ = raycast('auto') attempts to create a raycast object by 
%   extracting data from a csv file.
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on raycast %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = raycast('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0);
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
Args.classname = 'raycast';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'el';

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

% example object
dlist = nptDir;
% get entries in directory
dnum = size(dlist,1);

% check if the right conditions were met to create object
if(dnum>0)
	% these are object specific fields
	%data.dlist = dlist;
	% set index to keep track of which data goes with which directory
	%data.setIndex = [0; dnum];
    
    rd = pwd;
	cd ~/Documents/CSV 
    file = readtable('session_1_5112018105323_out.csv');
    cd (rd);
    el = load ('eyelink.mat');
    %assuming the eel file is completed 
    trialTimestamps = el.el.data.trial_timestamps + double(el.el.data.expTime);
    allTimes = file.Var2; %stores all the times in the csv file 
    
    %look for all the timestamps for the start, cue and end/timeout 
    %store the index they appear at in the file in the structure "index"
    index(:,1) = arrayfun(@(x) find(allTimes==x,1), trialTimestamps(:,1));
    index(:,2) = arrayfun(@(x) find(allTimes==x,1), trialTimestamps(:,2));
    index(:,3) = arrayfun(@(x) find(allTimes==x,1), trialTimestamps(:,3));
    index(2:end, 1) = index(2:end, 1)+1;
    index(:, 2:3) = index(:, 2:3) + 1;

    %Stores the indices for fixations as well to be abe to access relevant
    %information at that event
    fixTimes = sortrows(vertcat(el.el.data.fix_times(:,1), el.el.data.fix_times(:,2)));
    fixIndex = find(ismember(allTimes,fixTimes));
    
    %Stores the names to check for any missing trials or markers 
    trialNames(:,1) = file.Var3(index(:,1));
    trialNames (:,3)= file.Var3(index(:,3));
    trialNames (:,2)= file.Var3(index(:,2));
    
    % create nptdata so we can inherit from it
	data.numSets = size(trialNames,1); %no of trials in the session    
    data.Args = Args;
    
    rawGazeData(:,1) = file.Var4;
    rawGazeData(:,2) = file.Var5;
    
    data.timestamps = file.Var2;
    data.trialTimestamps = trialTimestamps;
    data.fixatedObj = file.Var3;
    data.index = index;
    data.fixIndex = fixIndex;
    data.rawGazeData =  rawGazeData;
    data.playerLocation = horzcat (file.Var6, file.Var7, file.Var8, file.Var9);  
    data.playerGazeLocation = horzcat(file.Var10, file.Var11, file.Var12);
    data.fixatedObjLoc = horzcat(file.Var13, file.Var14, file.Var15);
    data.RelativeToFixdObjGaze = horzcat(file.Var16, file.Var17);
    
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
