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
				'ObjectLevel','Session', 'FileLineOfffset',15, 'DirName','RawData*', ...
				'FileName','session*txt', 'TriggerVal1',10, 'TriggerVal2',20, ...
				'TriggerVal3',30, 'GridSteps',5);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {'GridSteps'};

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
% [pdir,cdir] = getDataOrder(Args.ObjectLevel,'relative','CDNow');
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

% return to previous directory
% if(~isempty(cdir))
%     cd(cdir);
% end

function obj = createObject(Args,varargin)

% move to correct directory
[pdir,cwd] = getDataOrder(Args.ObjectLevel,'relative','CDNow');

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

        % in case of an aborted session we might have to toss out the markers in the last
        % trial. So check the number of rows in each column before creating an array
		% with the smallest number of rows, which means we throw out the incomplete
		% set of markers in the last trial
		utRows = zeros(1,3);
	    uT1 = find(unityData(:,1) > Args.TriggerVal1 & unityData(:,1) < Args.TriggerVal2); % look for triggers '1' (start trial, cue onset) in first column
        utRows(1,1) = size(uT1,1);
	    uT2 = find(unityData(:,1) > Args.TriggerVal2 & unityData(:,1) < Args.TriggerVal3); % look for triggers '2' (start trial, cue offset) in first column
        utRows(1,2) = size(uT2,1);
        uT3 = find(unityData(:,1) > Args.TriggerVal3);  % look for triggers '3' (reward) or '4' (timeout) in last column
        utRows(1,3) = size(uT3,1);
		% find smallest number of rows
		mutrows = min(utRows);
		% find largest number of rows
		incompletetrials = max(utRows) - mutrows;
		% check if there were any incomplete markers
		if(incompletetrials>0)
			% print warning to say that there were incomplete markers
			disp(['Incomplete session! Last ' num2str(incompletetrials) ' trial discarded\n'])
		end
		unityTriggers = zeros(mutrows,3);
		utIndices = 1:mutrows;
		unityTriggers(:,1) = uT1(utIndices);
		unityTriggers(:,2) = uT2(utIndices);
		unityTriggers(:,3) = uT3;
		totTrials = size(unityTriggers,1); % total trials inluding repeated trials due to error/timeout

		% default grid is 5 x 5 but can be changed 
		gridSteps = Args.GridSteps;
		% these should be read from the unity file
        gridBins = gridSteps * gridSteps;
		overallGridSize = 25;
		oGS2 = overallGridSize/2;
		gridSize = overallGridSize/gridSteps;
		horGridBound = -oGS2:gridSize:oGS2;
		vertGridBound = horGridBound;
		% need to add one more bin to the end as histcounts counts by doing: edges(k) ≤ X(i) < edges(k+1)
		gpEdges = 1:(gridBins+1);

		% get gridpositions
		[h2counts,horGridBound,vertGridBound,binH,binV] = histcounts2(unityData(:,3),unityData(:,4),horGridBound,vertGridBound);

		% compute grid position number
		gridPosition = binH + ((binV - 1) * gridSteps);

		% create memory for arrays
		gpDurations = zeros(gridBins,totTrials);
		% add 1 to have the correct number of rows because of the subtraction
		% e.g. if uT3 was 50 and unityTriggers(:,2) was 1, there should be 50-1+1 rows
		% add another 1 to have 0 as the first bin for histcounts
		unityTrialTime = nan(max(uT3-unityTriggers(:,2))+2,totTrials);

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

		trialCounter = 0; % set up trial counter
        % initialize index for setting non-navigating gridPositions to 0
        gpreseti = 1;

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
		        mpath(b+1,1) = I3; % store nearest vertice for each location (per frame) into matrix 'mpath'
			end
        
			pathdiff = diff(mpath); 
			change = [1; pathdiff(:)]; 
			index = find(abs(change)>0); % get index of vertex change
			actualRoute = mpath(index); % get actual route
			actualCost = (size(actualRoute,1)-1)*5; % get cost of actual route
            actualTime = index; 
			actualTime(end+1) = size(mpath,1); 
			actualTime = diff(actualTime)+1; % get time in each space bin
            % clear mpath; 
			clear pathdiff; 
			clear change; 
			clear index;         
        
			% Store summary
			sumCost(a,1) = idealCost; 
		    sumCost(a,2) = actualCost; 
		    sumCost(a,3) = actualCost - idealCost; % get difference between least cost and actual cost incurred        
			sumCost(a,4) = target; % mark out current target identity
			sumCost(a,5) = unityData(unityTriggers(a,3),1)-target; % mark out correct/incorrect trials        
			sumRoute(a,1:size(idealroute,2)) = idealroute; % store shortest route
			sumActualRoute(a,1:size(actualRoute,1)) = actualRoute; % store actual route
            sumActualTime(a,1:size(actualTime,1)) = actualTime; % store actual time spent in each space bin
            
		    if sumCost(a,3) <= 0, % least distance taken
		        sumCost(a,6) = 1; % mark out trials completed via shortest route
			elseif sumCost(a,3) > 0 && sumCost(a,5) == 30, 
		        pathdiff = diff(actualRoute); % check if there's a one grid change of mind. If so, enable user to inspect trajectory 
		        
                for c = 1:size(pathdiff,1)-1,
                    if pathdiff(c) == pathdiff(c+1)*(-1);
                        timeingrid = size(find(mpath == actualRoute(c+1)),1);

                        if timeingrid > 165, % greater than 5 seconds spent in grid (unlikely to have entered by accident due to poor steering)
                            break
                        else
                            % disp('im here');
                            sumCost(a,6) = 1;
                        end
                    else                
                    end
                end
               
		        clear mpath pathdiff 
		    end % if sumCost(a,3) <= 0,

	        % get indices for this trial
	        % unityTriggers(a,2) will be the index where the 2nd marker was found
	        % and the time recorded on that line should be time since the last update
	        % so we really want the time for the next update instead
			uDidx = (unityTriggers(a,2)+1):unityTriggers(a,3);
			
			% get cumulative time from cue offset to end of trial for this specific trial
			% this will make it easier to correlate to spike times
			% add zero so that the edges for histcount will include 0 at the beginning
			unityTrialTime(1:(size(uDidx,2)+1),a) = [0; cumsum(unityData(uDidx,2))]; 
	
			% get grid positions for this trial
			tgp = gridPosition(uDidx);
	
			% get unique positions
			utgp = unique(tgp);
	
			for pidx = 1:size(utgp,1)
				tempgp = utgp(pidx);
				% find indices that have this grid position
				utgpidx = find(tgp==tempgp);
				gpDurations(tempgp,a) = sum(unityData(utgpidx,2));
			end  % for pidx = 1:size(utgp,1)
			
			% set gridPositions when not navigating to 0
			gridPosition(gpreseti:uDidx(1)) = 0;
            gpreseti = unityTriggers(a,3)+1;
		end % for a = 1:totTrials

		% Calculate performance
		errorInd = find(sumCost(:,5) == 40); 
		sumCost(errorInd,6) = 0; 
		sumCost(errorInd+1,6) = 0; % exclude rewarded trials that were preceded by a timeout
		perf = sum(sumCost(:,6))/50; 
		disp(strcat('% trials completed via shortest path = ', num2str(perf))); % get percentage of correct trials completed via the shortest route (calculate as a percentage of correct trials preceded by a correct trial)
		processTrials = find(sumCost(:,6) == 1); % Analyse only one-hit trials (comment out to plot all trials indiscriminately)

		data.gridSteps = gridSteps; 
		data.overallGridSize = overallGridSize;
		data.oGS2 = oGS2;
		data.gridSize = gridSize;
		data.horGridBound = horGridBound;
		data.vertGridBound = vertGridBound;
		data.gpEdges = gpEdges;
		data.unityData = unityData;
		data.unityTriggers = unityTriggers;
		data.sumCost = sumCost;
		data.sumRoute = sumRoute;
		data.sumActualRoute = sumActualRoute;
        data.sumActualTime = sumActualTime;
		data.perf = perf;
		data.processTrials = processTrials;
		data.gridPosition = gridPosition;
		data.gpDurations = gpDurations;
		data.unityTrialTime = unityTrialTime;
		data.setIndex = [1; totTrials];
        % compute cumulative sum of unity time to make it easy for
        % placeselect.m to compute histograms for shuffled data
        % add a zero at the beginning to avoid spike from being missed
        data.unityTime = [0; cumsum(unityData(:,2))];

		% move back to session directory from RawData directory to make
		% the object is created and saved in the correct directory
		cd ..
		% create nptdata so we can inherit from it
	    data.Args = Args;
		n = nptdata(data.numSets,0,pwd);
		d.data = data;
		obj = class(d,Args.classname,n);
		saveObject(obj,'ArgsC',Args);
	else % if(dnum>0)
		% move back to session directory from RawData directory
		cd ..
		% create empty object
		obj = createEmptyObject(Args);
	end % if(dnum>0)

else % if(~isempty(rd))
	% create empty object
	obj = createEmptyObject(Args);
end % if(~isempty(rd))

% move back to previous directory
cd(cwd)


function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;

% these are object specific fields
data.gridSteps = []; 
data.overallGridSize = [];
data.oGS2 = [];
data.gridSize = [];
data.horGridBound = [];
data.vertGridBound = [];
data.gpEdges = [];
data.unityData = [];
data.unityTriggers = [];
data.sumCost = [];
data.sumRoute = [];
data.sumActualRoute = [];
data.sumActualTime = [];
data.perf = [];
data.processTrials = [];
data.gridPosition = [];
data.gpDurations = [];
data.unityTrialTime = [];
data.setIndex = [];
data.unityTime = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
