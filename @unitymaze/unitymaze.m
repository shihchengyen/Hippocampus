function [obj, varargout] = unitymaze(varargin)
%@unitymaze Constructor function for UNITYMAZE class
%   OBJ = unitymaze(varargin)
%
%   OBJ = unitymaze('auto') attempts to create a DIRFILES object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on unitymaze %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = unitymaze('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'FileLineOfffset',15, 'DirName','RawData_*', ...
				'FileName','session*txt', 'TriggerVal1',10, 'TriggerVal2',20, ...
				'TriggerVal3',30);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'unitymaze';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'um';

% To decide the method to create or load the object
% change to proper directory to check for saved object
% [pdir,cdir] = getDataOrder('session','relative','CDNow');
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

% return to previous directory
% if(~isempty(cdir))
%     cd(cdir);
% end

function obj = createObject(Args,varargin)

% save current directory
cwd = pwd;
	
% look for session_1_*.txt in RawData_T*
rd = dir(Args.DirName);
if(~isempty(rd))
	cd(rd(1).name)
	dlist = nptDir(Args.FileName);
	% get entries in directory
	dnum = size(dlist,1);

	% check if the right conditions were met to create object
	if(dnum>0)
		% this is a valid object
		% these are fields that are useful for most objects
		data.numSets = dnum;
	
		% these are object specific fields
	    % Concatenate unity data into one matrix (accommodates data structure of trial-by-trial training sessions)
	    length = 0;
	    for i = 1:dnum % go through all 'session_' files and extract data
	        temp = dlmread(dlist(i).name,'',Args.FileLineOfffset,0); % start reading at row 15 (skip logged parameters)
	        unityData(length+1:length+size(temp,1),1:5) = temp;
	        length = length + size(temp,1);
	    end

	    unityTriggers(:,1) = find(unityData(:,1) > Args.TriggerVal1 & unityData(:,1) < Args.TriggerVal2); % look for triggers '1' (start trial, cue onset) in first column
	    unityTriggers(:,2) = find(unityData(:,1) > Args.TriggerVal2 & unityData(:,1) < Args.TriggerVal3); % look for triggers '2' (start trial, cue offset) in first column
	    unityTriggers(:,3) = find(unityData(:,1) > Args.TriggerVal3); % look for triggers '3' (reward) or '4' (timeout) in last column

		%% DIJKSTRA 
		% Initialize adjacency matrix for graph (21 vertices, distance weights)
		A = [0 5 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
		     5 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
		     0 5 0 5 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
		     0 0 5 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
		     0 0 0 5 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0;
		     5 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0;
		     0 0 5 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0;
		     0 0 0 0 5 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0;
		     0 0 0 0 0 5 0 0 0 5 0 0 0 5 0 0 0 0 0 0 0;
		     0 0 0 0 0 0 0 0 5 0 5 0 0 0 0 0 0 0 0 0 0;
		     0 0 0 0 0 0 5 0 0 5 0 5 0 0 5 0 0 0 0 0 0;
		     0 0 0 0 0 0 0 0 0 0 5 0 5 0 0 0 0 0 0 0 0;
		     0 0 0 0 0 0 0 5 0 0 0 5 0 0 0 5 0 0 0 0 0;
		     0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 5 0 0 0 0;
		     0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 5 0 0;
		     0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 5;
		     0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 5 0 0 0;
		     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 5 0 0;
		     0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 5 0 5 0;
		     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 5;
		     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 5 0;
		     ];

		% Vertices coordinates: 
		vertices = [-10 10; -5 10; 0 10; 5 10; 10 10; -10 5; 0 5; 10 5; -10 0; -5 0; 0 0; 5 0; 10 0; -10 -5; 0 -5; 10 -5; -10 -10; -5 -10; 0 -10; 5 -10; 10 -10;];

		% Poster coordinates
		posterpos = [-5 -7.55; -7.55 5; 7.55 -5; 5 7.55; 5 2.45; -5 -2.45]; % Posters 1 to 6 respectively

		% Plot boundaries
		xBound = [-12.5,12.5,12.5,-12.5,-12.5]; % (clockwise from top-left corner) outer maze wall
		zBound = [12.5,12.5,-12.5,-12.5,12.5];
		x1Bound =[-7.5,-2.5,-2.5,-7.5,-7.5]; % yellow pillar 
		z1Bound =[7.5,7.5,2.5,2.5,7.5];
		x2Bound =[2.5,7.5,7.5,2.5,2.5]; % red pillar
		z2Bound =[7.5,7.5,2.5,2.5,7.5];
		x3Bound =[-7.5,-2.5,-2.5,-7.5,-7.5]; % blue pillar
		z3Bound =[-2.5,-2.5,-7.5,-7.5,-2.5];
		x4Bound =[2.5,7.5,7.5,2.5,2.5]; % green pillar
		z4Bound =[-2.5,-2.5,-7.5,-7.5,-2.5];

		totTrials = size(unityTriggers,1); % total trials inluding repeated trials due to error/timeout
		trialCounter = 0; % set up trial counter

		for a = 1:totTrials
			trialCounter = trialCounter + 1;

			% get target identity
			target = rem(unityData(unityTriggers(a,3),1),10); % target identity coded by digit in ones position  
    
			% (starting position) get nearest neighbour vertex
			S = pdist2(vertices,unityData(unityTriggers(a,2),3:4),'euclidean'); % get distance from all vertices
		    [M1,I1] = min(S); % M is the minimum value, I is the index of the minimum value
			startPos = I1;
        
			% (destination, target) get nearest neighbour vertex
			D = pdist2(vertices,posterpos(target,:),'euclidean'); % get distance from all vertices
			[M2,I2] = min(D); % M is the minimum value, I is the index of the minimum value
			destPos = I2;
        
			% get least distance (ideal route)
			[idealCost,idealroute] = dijkstra(A,destPos,startPos);
         
		    % get actual route taken (match to vertices)
			for b = 0:unityTriggers(a,3)-unityTriggers(a,2), % go through unity file line by line
		        currPos = unityData(unityTriggers(a,2)+b,3:4);

		        % (current position) get nearest neighbour vertex
		        CP = pdist2(vertices,currPos,'euclidean'); % get distance from all vertices
		        [M3,I3] = min(CP); % M is the minimum value, I is the index of the minimum value
		        path(b+1,1) = I3; % store nearest vertice for each location (per frame) into matrix 'path'
			end
        
			pathdiff = diff(path); change = [1; pathdiff(:)]; index = find(abs(change)>0); % get index of vertex change
			actualRoute = path(index); % get actual route
			actualCost = (size(actualRoute,1)-1)*5; % get cost of actual route
			clear path; clear pathdiff; clear change; clear index;         
        
			% Store summary
			sumCost(a,1) = idealCost; 
		    sumCost(a,2) = actualCost; 
		    sumCost(a,3) = actualCost - idealCost; % get difference between least cost and actual cost incurred        
			sumCost(a,4) = target; % mark out current target identity
			sumCost(a,5) = unityData(unityTriggers(a,3),1)-target; % mark out correct/incorrect trials        
			sumRoute(a,1:size(idealroute,2)) = idealroute; % store shortest route
			sumActualRoute(a,1:size(actualRoute,1)) = actualRoute; % store actual route
    
		    if sumCost(a,3) <= 0, % least distance taken
		        sumCost(a,6) = 1; % mark out trials completed via shortest route
			elseif sumCost(a,3) > 0 && sumCost(a,5) == 30, 
		        pathdiff = diff(actualRoute); % check if there's a one grid change of mind. If so, enable user to inspect trajectory 
		        for c = 1:size(pathdiff,1)-1,
		            if pathdiff(c) == pathdiff(c+1)*(-1);

		                plot(xBound,zBound,'k','LineWidth',1.5);
						hold on
		                plot(x1Bound,z1Bound,'k','LineWidth',1);
		                plot(x2Bound,z2Bound,'k','LineWidth',1);
		                plot(x3Bound,z3Bound,'k','LineWidth',1);
		                plot(x4Bound,z4Bound,'k','LineWidth',1);
		                plot(unityData(unityTriggers(a,2):unityTriggers(a,3)-1,3),unityData(unityTriggers(a,2):unityTriggers(a,3)-1,4),'b','LineWidth',1); % plot current trial trajectory
		                plot(unityData(unityTriggers(a,3),3),unityData(unityTriggers(a,3),4),'k.','MarkerSize',20); % plot end point identifier
						hold off

		                userInput = input('Enter (1) to ACCEPT and (2) to REJECT as shortest route: ');

		                    if userInput == 1,
		                        sumCost(a,6) = 1;
		                    end

		                    % close all % close figure     
		                break % exit for loop
		            end % if pathdiff(c) == pathdiff(c+1)*(-1);
		        end % for c = 1:size(pathdiff,1)-1,
		        clear pathdiff 
		    end % if sumCost(a,3) <= 0,
		end % for a = 1:totTrials

		cd(cwd)

		% Calculate performance
		errorInd = find(sumCost(:,5) == 40); 
		sumCost(errorInd,6) = 0; sumCost(errorInd+1,6) = 0; % exclude rewarded trials that were preceded by a timeout
		perf = sum(sumCost(:,6))/50; disp(strcat('% trials completed via shortest path = ', num2str(perf))); % get percentage of correct trials completed via the shortest route (calculate as a percentage of correct trials preceded by a correct trial)
		processTrials = find(sumCost(:,6) == 1); % Analyse only one-hit trials (comment out to plot all trials indiscriminately)

		% clearvars -except day session directory unityTriggers unityData sumCost sumRoute sumActualRoute perf processTrials
		% save('processTrials','processTrials'); save('unityData','unityData'); save('unityTriggers','unityTriggers');

		data.unityData = unityData;
		data.unityTriggers = unityTriggers;
		data.sumCost = sumCost;
		data.sumRoute = sumRoute;
		data.sumActualRoute = sumActualRoute;
		data.perf = perf;
		data.processTrials = processTrials;
		data.setIndex = [1; totTrials];

		% create nptdata so we can inherit from it
	    data.Args = Args;
		n = nptdata(data.numSets,0,pwd);
		d.data = data;
		obj = class(d,Args.classname,n);
		saveObject(obj,'ArgsC',Args);
	else % if(dnum>0)
		% create empty object
		obj = createEmptyObject(Args);
	end % if(dnum>0)

else % if(~isempty(rd))
	% create empty object
	obj = createEmptyObject(Args);
end % if(~isempty(rd))


function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;

% these are object specific fields
data.unityData = [];
data.unityTriggers = [];
data.sumCost = [];
data.sumRoute = [];
data.sumActualRoute = [];
data.perf = [];
data.processTrials = [];
data.setIndex = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
