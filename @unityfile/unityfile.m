
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
				'TriggerVal3',30, 'MaxTimeDiff',0.002);
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
        
		unityTrialTime = nan(max(uT3-unityTriggers(:,2))+2,totTrials);


		trialCounter = 0; % set up trial counter

		sTi = 2;

		for a = 1:totTrials
			trialCounter = trialCounter + 1;

			% get target identity
			target = rem(unityData(unityTriggers(a,3),1),10); % target identity coded by digit in ones position  

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

			% get the starting Ripple timestamp
			tstart = rp.data.timeStamps(a,2);
			tend = rp.data.timeStamps(a,3);
			
			% compare trial durations
			% if the difference is acceptable, we will add the timestamps
			% to the histogram bin limits. Otherwise, we will skip to the
			% end of the trial
			if(abs(tdiff) < Args.MaxTimeDiff)

				% shift the Unity timestamps to match the Ripple timestamps
				% by distributing the difference over all the frames
				unityTrialTime(tindices,a) = tempTrialTime - [0; cumsum(repmat(tdiff/numUnityFrames,numUnityFrames,1))];
                
			else
				unityTrialTime(tindices,a) = tempTrialTime;
	
			end			

		end % for a = 1:totTrials

        
        %After creating unityData, you need to correct the times generated
        %using the eyelink object.
        unityData = synchronise(unityData, unityTriggers);
        
		data.unityData = unityData;
		data.unityTriggers = unityTriggers;
		data.unityTrialTime = unityTrialTime;

        % compute cumulative sum of unity time to make it easy for
        % placeselect.m to compute histograms for shuffled data
        % add a zero at the beginning to avoid spike from being missed
        data.unityTime = [0; cumsum(unityData(:,2))];

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
data.unityData = [];
data.unityTriggers = [];

data.unityTrialTime = [];
data.unityTime = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
