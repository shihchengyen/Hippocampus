
function [obj, varargout] = unityfile(varargin)
%@unityfile Constructor function for unityfile class
%   OBJ = unityfile(varargin)
%
%   OBJ = unityfile('auto') attempts to create a DIRFILES object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on unityfile %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = unityfile('save','redo')
%
%dependencies: 
Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Session', 'FileLineOfffset',15, 'DirName','RawData*', ...
				'FileName','session*txt', 'TriggerVal1',10, 'TriggerVal2',20, ...
				'TriggerVal3',30, 'GridSteps',5, 'MaxTimeDiff',0.002);
Args.flags = {'Auto','ArgsOnly'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {'GridSteps'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'unityfile';
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
[pdir,cwd] = getDataOrder('Session','relative','CDNow');

% need to correct timestamps by comparing to the marker timestamps
% in rplparallel
rp = rplparallel('auto',varargin{:});
% compute trial durations from the timestamps in rplparallel
rplTrialDur = diff(rp.data.timeStamps(:,2:3),1,2);

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
	    
		% move back to session directory from RawData directory to make
		% the object is created and saved in the correct directory
		cd ..

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
		% need to add one more bin to the end as histcounts counts by doing: edges(k) ² X(i) < edges(k+1)
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

		% create memory for array for histogram time bins, grid position,
		% and interval
		sessionTime = zeros(size(gridPosition,1),3);	
		% sessionTime will have the start and end timestamps of each trial
		% along with when the grid position changed in the 1st column and 
		% the grid positions in the 2nd column. The 3rd column will be
		% filled in later with the interval between successive timestamps
		% in the 1st column. 
		% e.g.
		% 0			0
		% 7.9724	181
		% 8.1395	221
		% ...
		% 10.6116	826
		% 10.7876	0
		% 13.8438	826
		% ...
		% 22.7543	693
		% 23.5624	0
		% 27.2817	0
		% 36.2458	0
		% 39.4584	701
		% The 1st trial started at 7.9724 s, with the animal at grid position
		% 181. The animal then moved to position 221 at 8.1395 s. Near
		% the end of the trial, the animal moved to position 826 at 10.6116 s,
		% and stayed there until the end of the 1st trial at 10.7876 s, which 
		% was marked by a 0 in the 2nd column. The 2nd trial started at 
		% 13.8438 s at position 826, and ended at 23.5624 s at position 693.
		% The next 2 lines illustrates what happens if the mismatch in the 
		% duration between Unity and Ripple is too large. The start of the
		% trial at 27.2817 s is marked with position 0, and is immediately
		% followed by the end of the trial at 36.2458 s, marked as before
		% with a 0 in the 2nd column.
		% 
		% start the array with 0 to make sure any spike times before the 
		% first trigger	are captured 
		% increment the index for sessionTime
		sTi = 2;

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

			% check how much of a discrepancy there is between Unity and Ripple
	        % get indices for this trial
	        % unityTriggers(a,2) will be the index where the 2nd marker was found
	        % and the time recorded on that line should be time since the last update
	        % so we really want the time for the next update instead
			uDidx = (unityTriggers(a,2)+1):unityTriggers(a,3);
			% get number of frames in this trial
			numUnityFrames = size(uDidx,2);
			
			% get cumulative time from cue offset to end of trial for this specific trial
			% this will make it easier to correlate to spike times
			% add zero so that the edges for histcount will include 0 at the beginning
			% get indices for this trial
			tindices = 1:(numUnityFrames+1);
			tempTrialTime = [0; cumsum(unityData(uDidx,2))]; 
			% unityTrialTime(tindices,a) = [0; cumsum(unityData(uDidx,2))]; 
			
			% get Unity end time for this trial
			uet = tempTrialTime(end);
			% get Ripple end time for this trial
			ret = rplTrialDur(a);
			% compute the difference
			tdiff = uet - ret;
	
			% get grid positions for this trial
			tgp = gridPosition(uDidx);
	
			% get the starting Ripple timestamp
			tstart = rp.data.timeStamps(a,2);
			tend = rp.data.timeStamps(a,3);
			
			% compare trial durations
			% if the difference is acceptable, we will add the timestamps
			% to the histogram bin limits. Otherwise, we will skip to the
			% end of the trial
			if(abs(tdiff) < Args.MaxTimeDiff)
				sessionTime(sTi,:) = [tstart tgp(1) 0];
				sTi = sTi + 1;
				% shift the Unity timestamps to match the Ripple timestamps
				% by distributing the difference over all the frames
				unityTrialTime(tindices,a) = tempTrialTime - [0; cumsum(repmat(tdiff/numUnityFrames,numUnityFrames,1))];

				% find the timepoints where grid positions changed
				gpc = find(diff(tgp)~=0);
				ngpc = size(gpc,1);
			
				% add the Unity frame intervals to the starting timestamp to
				% create corrected version of unityTime, which will also be the
				% bin limits for the histogram function call
				sessionTime(sTi:(sTi+ngpc-1),1:2) = [unityTrialTime(gpc+2,a)+tstart tgp(gpc+1)];	
				sTi = sTi + ngpc;		
			else
				unityTrialTime(tindices,a) = tempTrialTime;
				% leave the 2nd column as 0 to indicate this was a skipped trial
				sessionTime(sTi,1) = tstart;
				sTi = sTi + 1;			
			end			
				
			% add an entry for the end of the trial
			sessionTime(sTi,1:2) = [tend 0];
			sTi = sTi + 1;

			% get unique positions
			utgp = unique(tgp);
	
			for pidx = 1:size(utgp,1)
				tempgp = utgp(pidx);
				% find indices that have this grid position
				utgpidx = find(tgp==tempgp);
				gpDurations(tempgp,a) = sum(unityData(utgpidx,2));
			end  %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             for pidx = 1:size(utgp,1)
			
			% set gridPositions when not navigating to 0
			% subtract 1 from uDidx(1) as we set the start of uDidx to 1
			% row after unityTrigger(a,2)
			gridPosition(gpreseti:(uDidx(1)-1)) = 0;
            gpreseti = unityTriggers(a,3)+1;
		end % for a = 1:totTrials

		% Calculate performance
		errorInd = find(sumCost(:,5) == 40); 
		sumCost(errorInd,6) = 0; 
		sumCost(errorInd+1,6) = 0; % exclude rewarded trials that were preceded by a timeout
		perf = sum(sumCost(:,6))/50; 
		disp(strcat('% trials completed via shortest path = ', num2str(perf))); % get percentage of correct trials completed via the shortest route (calculate as a percentage of correct trials preceded by a correct trial)
		processTrials = find(sumCost(:,6) == 1); % Analyse only one-hit trials (comment out to plot all trials indiscriminately)

		% get number of rows in sessionTime
		snum = sTi - 1;
		% reduce memory for sessionTime
		sTime = sessionTime(1:snum,:);
		% fill in 3rd column with time interval so it will be easier to compute
		% firing rate
		sTime(1:(snum-1),3) = diff(sTime(:,1));
		% sort the 2nd column so we can extract the firing rates by position
		[sTP,sTPi] = sort(sTime(:,2));
		% find the number of observations per position by looking for 
		% the indices when position changes. The first change should be from
		% position 0 to 1st non-zero position. Add 1 to adjust for the change
		% in index when using diff. These will be the starting indices for 
		% the unique positions.
		sTPsi = find(diff(sTP)~=0) + 1;
		% find the ending indices by subtracting 1 from sTPsi, and adding
		% snum at the end
		sTPind = [sTPsi [[sTPsi(2:end)-1]; snum]];
		sTPin = diff(sTPind,1,2);
		% then we find the differences between consecutive sTPi2
		% this means that we are counting the number of points from the 1st
		% non-zero position, so we don't need to drop the 1st value.
		% But we need to add the last value, which will go to the end of snum.
		% sTPn = [diff(sTPsi); (snum-sTPsi(end))];
		% find the largest number of observations, which can be used to 
		% pre-allocate memory for matrix to store firing rates
		% sTPnmax = max(sTPn);
		% find unique positions, which will include 0
		tempsTPu = unique(sTP);
		% remove 0
		sTPu = tempsTPu(2:end);
		% find number of unique positions
		nsTPu = size(sTPu,1);
		% when we need to pre-allocate memory, we can do something like: 
		% fri = nan(sTPnmax,nsTPu);
		% compute occupancy proportion
		ou_i = zeros(nsTPu,1);
		for pi = 1:nsTPu
			ou_i(pi) = sum(sTime(sTPi(sTPind(pi,1):sTPind(pi,2)),3));
        end
        
        %After creating unityData, you need to correct the times generated
        %using the eyelink object.
        unityData = synchronise(unityData, unityTriggers);
        
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
		data.setIndex = [0; totTrials];
        % compute cumulative sum of unity time to make it easy for
        % placeselect.m to compute histograms for shuffled data
        % add a zero at the beginning to avoid spike from being missed
        data.unityTime = [0; cumsum(unityData(:,2))];
        data.sTime = sTime;
        data.sTPi = sTPi;
        data.sTPind = sTPind;
        data.sTPin = sTPin;
        data.sTPu = sTPu;
        data.nsTPu = nsTPu;
        data.ou_i = ou_i;
        data.P_i = ou_i / sum(ou_i);

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
