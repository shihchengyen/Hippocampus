function [obj, varargout] = rosdata(varargin)
%@rosdata Constructor function for rosdata class
%   OBJ = rosdata(varargin)
%
%   OBJ = rosdata('auto') attempts to create a rosdata object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rosdata %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Example [as, Args] = rosdata('save','redo')
%
%   Dependencies: 
%     The rosbag name saved in the data folder has to changed to 'posebag.bag'
    Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
                    'BagName','posebag*bag','MapName','maptest0.pgm', ...                
                    'ObjectLevel', 'Session');
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
    Args.classname = 'rosdata';
    Args.matname = [Args.classname '.mat'];
    Args.matvarname = 'rd';

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
end

function obj = createObject(Args,varargin)
    %% Import Data
    % Import RosBag files, and extract X,Y Coordinates
    poseBag=rosbag(Args.BagName);
    data.mapimage = Args.MapName;

    bSelpose = select(poseBag,'Topic','/amcl_pose');
    data.poseStructs = readMessages(bSelpose,'DataFormat','struct');
    bSeltrigg = select(poseBag,'Topic','/trigger_msgs');
    data.triggStructs = readMessages(bSeltrigg,'DataFormat','struct');
    
    %% Splitting pose into trials using trigger time

    % take the data of trigger_msgs out as a vector of double
    numMsg = cellfun(@(m) double(m.Data),data.triggStructs);

    poseTime = bSelpose.MessageList.Time;
    poseTime_diff = poseTime-poseTime(1);
    
    % extract the index of the starting trigger_msgs' time
    indexTg1 = find(numMsg>10 & numMsg<20);
    % extract the index of the ending trigger_msgs' time
    index = find(numMsg>30);    
    triggTime = bSeltrigg.MessageList.Time(index);
    triggTime_diff = triggTime-poseTime(1);
    triggStartTime = bSeltrigg.MessageList.Time(indexTg1);
    triggStartTime = triggStartTime-poseTime(1);
    % extract the two digits stopped at the end of each trial
    poster_stopped = numMsg(index);
    % extract the animal cues
    poster_indx = find(numMsg>20 & numMsg<30);
    animal_cues = rem(numMsg(poster_indx),20);
    data.animal_cues = animal_cues;
    data.poster_stopped = poster_stopped;

    eots = [];  % eots[] saves the timeStamp of the trials' end points
    heads = [];
    for i=1:length(triggTime)
        
        tg1= find(poseTime_diff <= triggStartTime(i)); 
        heads = [heads tg1(end)];

        % triggTime(i) is used to mark the end of the ith trial
        eot = find(poseTime_diff <= triggTime_diff(i));  %eot: end of trial
        eots = [eots eot(end)];
    end
     
    totalTS = size(poseTime_diff);
    if totalTS(1)>eots(end)
    	eots = [eots totalTS(1)];
    end
    data.eots = eots;
    data.heads = heads;
    
    % inherit from DPV
    data.numSets = length(eots);
	data.Args = Args;
	n = nptdata(data.numSets,0,pwd);
    r.data = data;
	obj = class(r,Args.classname,n);
	saveObject(obj,'ArgsC',Args);
end

function obj = createEmptyObject(Args)

    % these are object specific fields
    data.eots = [];
    data.poseStructs = {};
    data.triggStructs = {};
    data.mapimage = strings;

    % create nptdata so we can inherit from it
    % useful fields for most objects
    data.numSets = 0;
    data.Args = Args;
    n = nptdata(0,0);
    r.data = data;
    obj = class(r,Args.classname,n);
end
% ceck object reco. as a DPV object