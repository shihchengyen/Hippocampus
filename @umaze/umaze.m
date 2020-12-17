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
				'GridSteps',40, 'overallGridSize',25, ...
                'MinObs',5);
            
% GridSteps: number of bins to cut the maze into
% overallGridSize: length/width of entire maze

Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps', 'MinObs', 'overallGridSize'};                            

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
    unityTime = ufdata.data.unityTime;
    unityData = ufdata.data.unityData;
	unityTriggers = ufdata.data.unityTriggers;
    unityTrialTime = ufdata.data.unityTrialTime;
    totTrials = size(unityTriggers,1); % total trials inluding repeated trials due to error/timeout

	% default grid is 5 x 5 but can be changed 
	gridSteps = Args.GridSteps;
    overallGridSize = Args.overallGridSize;  
    
    gridBins = gridSteps * gridSteps;
	oGS2 = overallGridSize/2;
    gridSize = overallGridSize/gridSteps; % size of each partition
    gpEdges = 1:(gridBins+1);
    
	horGridBound = -oGS2:gridSize:oGS2;
	vertGridBound = horGridBound;

	% get gridpositions (resolution of unityData lowered and allocated to bins)
	[h2counts,horGridBound,vertGridBound,binH,binV] = histcounts2(unityData(:,3),unityData(:,4),horGridBound,vertGridBound);

    % binH is the col number, binV is the row
	% compute grid position number for each timestamp (something like flattened index)
	gridPosition = binH + ((binV - 1) * gridSteps); 
		
	% will be filled with duration spent in each bin, for each trial
	gpDurations = zeros(gridBins,totTrials);
        

    trialCounter = 0; % set up trial counter
    % initialize index for setting non-navigating gridPositions to 0
    % (marks the cue onset index), to be used in conjunction with cue
    % offset later to form range to set to 0.
    gpreseti = 1;

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
    sessionTime = zeros(size(gridPosition,1),3); % was 5, with the coarse direction

    % start the array with 0 to make sure any spike times before the 
    % first trigger	are captured 
    % tracks row number for sessionTime
    sTi = 2;

    for a = 1:totTrials
        trialCounter = trialCounter + 1;

        % unityTriggers(a,2) will be the index where the 2nd marker was found
        % and the time recorded on that line should be time since the last update
        % so we really want the time for the next update instead
        uDidx = (unityTriggers(a,2)+1):unityTriggers(a,3);
        % get number of frames in this trial (for the cue offset till
        % end trial portion only
        numUnityFrames = size(uDidx,2);

        % get cumulative time from cue offset to end of trial for this specific trial
        % this will make it easier to correlate to spike times
        % add zero so that the edges for histcount will include 0 at the beginning
        % get indices for this trial
        tindices = 1:(numUnityFrames+1);
        tempTrialTime = [0; cumsum(unityData(uDidx,2))]; 
        tstart = unityTime(uDidx(1));
        tend = unityTime(uDidx(end)+1); % new mod 12062020

%             disp('trial:');
%             disp(a);
%             disp(tempTrialTime(end) - tempTrialTime(1));

        % get grid positions for this trial
        tgp = gridPosition(uDidx);

        % if tempTrialTime's last-first is 0, it means that the discrepency between
        % unity trial duration and ripple trial duration was too large,
        % and has already been flagged out within the unityfile process
        % by setting timestamps in the trial to the initial timestamp.

        if tempTrialTime(end)-tempTrialTime(1) ~= 0

            sessionTime(sTi,1:3) = [tstart tgp(1) 0];
            sTi = sTi + 1;

            % find the timepoints where grid positions changed
            gpc = find(diff(tgp)~=0);
            ngpc = size(gpc,1);
            
            % add the Unity frame intervals to the starting timestamp to
            % create corrected version of unityTime, which will also be the
            % bin limits for the histogram function call
            sessionTime(sTi:(sTi+ngpc-1),1:2) = [unityTrialTime(gpc+2,a)+tstart tgp(gpc+1)];            
%             sessionTime(sTi:(sTi+ngpc-1),4:5) = [binVt(gpc+1)-binVt(gpc) binHt(gpc)-binHt(gpc+1)];
            sTi = sTi + ngpc;

            % occasionally we will get a change in grid position in the frame interval
            % when we get the end of trial message. In that case, we will get an entry
            % in sessionTime that is the time of the end of the trial. Since the end
            % of the trial is added later, and because this will be a very brief visit
            % to the new position, we are going to remove it.
            if( (~isempty(gpc)) && (gpc(end) == (numUnityFrames-1)) )
                sTi = sTi - 1;
            end
        else
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
            utgpidx = uDidx(utgpidx);
            gpDurations(tempgp,a) = sum(unityData(utgpidx+1,2)); % kw mod to utgpidx+1 was utgpidx alone
        end
        
        % set gridPositions when not navigating to 0
        % subtract 1 from uDidx(1) as we set the start of uDidx to 1
        % row after unityTrigger(a,2)
        gridPosition(gpreseti:(uDidx(1)-1)) = 0;
        gpreseti = unityTriggers(a,3)+1;

    end % for a = 1:totTrials    
 
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
    % the unique positions not including 0.
    sTPsi = find(diff(sTP)~=0) + 1;
    if sTP(1) == -1
        sTPsi = sTPsi(2:end);
    end
    % find the ending indices by subtracting 1 from sTPsi, and adding
    % snum at the end
    sTPind = [sTPsi [ [sTPsi(2:end)-1]; size(sTP,1)]];
    % compute the number of observations per position
    sTPin = diff(sTPind,1,2) + 1;
    % arrange the information into a matrix for easier access
    sortedGPindinfo = [sTP(sTPsi) sTPind sTPin];
    % set up the conversion from grid position to index to make accessing the information
    % easier
    [~,gp2ind] = ismember(1:gridBins,sortedGPindinfo(:,1));
    % find positions (excluding 0) with more than MinReps observations
    sTPinm = find(sTPin>(Args.MinObs-1));
    
    
    
    % create temporary variable that will be used for calculations below
    % this is not a variable that we will save
    sTPsi2 = sTPsi(sTPinm);
    % save the number of observations for the positions that exceed the
    % minimum number of observations
    sTPin2 = sTPin(sTPinm);
    % save the position numbers
    sTPu = sTP(sTPsi2);
    % find number of positions
    nsTPu = size(sTPu,1);	
    % save the subset of starting and ending indices for sTime
    sTPind2 = sTPind(sTPinm,:);
    % compute occupancy proportion
    ou_i = zeros(nsTPu,1);
    for pi = 1:nsTPu
        ou_i(pi) = sum(sTime(sTPi(sTPind2(pi,1):sTPind2(pi,2)),3));
    end    
    
		data.gridSteps = gridSteps; 
		data.overallGridSize = overallGridSize;
		data.oGS2 = oGS2;
		data.gridSize = gridSize;
		data.horGridBound = horGridBound;
		data.vertGridBound = vertGridBound;
		data.gpEdges = gpEdges;

		data.gridPosition = gridPosition;
		data.gpDurations = gpDurations;
		data.setIndex = [0; totTrials];

        %data.processTrials = find(ufdata.data.sumCost(:,6)==1);
        data.processTrials = 1:totTrials; % should be top, but now unityfile.m still doesn't do shortest path calculation, placeholder variable
        
        data.sessionTime = sTime;
        data.sortedGPindices = sTPi;
        data.sortedGPindinfo = sortedGPindinfo;
        data.sGPi_minobs = sTPinm;
        data.sTPu = sTPu;
        data.nsTPu = nsTPu;
        data.ou_i = ou_i;
        data.P_i = ou_i / sum(ou_i); 
        data.gp2ind = gp2ind;
        
        %%% temporary, for sfn generation %%%
        data.unityTriggers = ufdata.data.unityTriggers;
        data.unityTrialTime = ufdata.data.unityTrialTime;

    data.unityData = unityData;
    data.unityTime = unityTime;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% umaze specific calculations ends here  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% create nptdata so we can inherit from it
    data.numSets = ufdata.data.numSets; 
    data.Args = Args;
	n = nptdata(data.numSets,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);
else
	% create empty object
	obj = createEmptyObject(Args);
end

% move back to previous directory
cd(cwd)

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
