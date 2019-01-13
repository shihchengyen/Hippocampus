function data = placeselect(um,rp,Args)
% function data = placeselect(dirstr)

% load arguments
% we are doing this so this program can run in compiled form
% load('vpData');

% um = unitymaze('auto');
% rp = rplparallel('auto');

gridBins = um.data.gridSteps * um.data.gridSteps;
% triggers that indicate cue onset, offset, and end of trial
unityTriggers = um.data.unityTriggers;
% trial numbers that should be processed, i.e. the animal took the shortest 
% route to the destination
processTrials = um.data.processTrials;
% raw Unity data that indicates position of the animal
unityData = um.data.unityData;
% get grid position number
gridPosition = um.data.gridPosition;
% get duration spent at each grid position
gpDurations = um.data.gpDurations;
unityTrialTime = um.data.unityTrialTime;
ntrials = size(unityTriggers,1);
nptrials = size(processTrials,1);

sampleRate = rp.data.SampleRate;
% compute trial durations from the timestamps in rplparallel
rpTrialDur = diff(rp.data.timeStamps(:,2:3),1,2);

% Get maze coordinates for each spike during trial (cue offset to end trial)
trialcount = 0; % count total number of trials
spikecountSess = 0; % count number of spikes (per trial)

spiketrain = load('spiketrain.mat');
sTimes = spiketrain.timestamps;

% create memory
spikeLocTrial = zeros(gridBins,nptrials);
anovaMatrix = spikeLocTrial;
tDiff = nan(nptrials,1);

% loop over number of trials            
for npi = 1:nptrials
    g = processTrials(npi);
    
	tstart = rp.data.timeStamps(g,2);
	tend = rp.data.timeStamps(g,3);

	% get indices for this trial
	uDidx = unityTriggers(g,2):unityTriggers(g,3);
	
	% get grid positions for this trial
	tgp = gridPosition(uDidx);

	temp = find(sTimes > tstart & sTimes < tend);
	% get spike times aligned to cue offset time
	trialSpkTimes = sTimes(temp) - tstart; 
	valIndices = ~isnan(unityTrialTime(:,g));
	uTrialTime = unityTrialTime(valIndices,g);
	tDiff(npi) = uTrialTime(end) - rpTrialDur(g);

	% only process trials with trigger discrepancies < 2 ms
	if (abs(tDiff(npi)) < Args.MaxTimeDiff)
		% compute number of spikes at each grid position
		% e.g. trialSpikeTimes is 0.3885    0.4156    6.1729    6.1768    6.1810    6.3559    6.3575    6.4199    8.6585    8.6661    8.9505
		% bins will be 9     9   150   150   150   156   156   158   240   240   248
		% get non-nan values
		[hcounts,uTT,bins] = histcounts(trialSpkTimes,uTrialTime);
		% now we use bins to get grid position
		% e.g. tgp(bins) will be  3     3    18    18    18    18    18    23    24    24    24
		% we then call histcounts again to count number of times each grid position appears
		% to get spikeLocTrial to be 0     0     2     0     0     0     0     0     0     0     0     0     0     0     0     0     0     5     0     0     0     0     1     3
		% need to remove 0’s from bins in case there were any spikes that were just outside the edges of the histcounts
		spikeLocTrial(:,npi) = histcounts(tgp(bins(bins~=0)),um.data.gpEdges);
	end
end 

% if there are any grid positions with less than 5 observations
% set them to 0 to avoid having outliers caused by a single observation
spikeLocTrial( (sum(spikeLocTrial>0,2)<Args.MinTrials) ,:) = 0;
% we divide by gpDurations to get number of spikes per duration
anovaMatrix = spikeLocTrial./gpDurations(:,processTrials);
meanFRs = nanmean(anovaMatrix,2);
semFRs = nanstd(anovaMatrix,0,2)./sqrt(sum(~isnan(anovaMatrix),2));

data.gridSteps = um.data.gridSteps;
data.meanFRs = meanFRs;
data.semFRs = semFRs;

if(isdeployed)
	% save data into a mat file
	save vmplacecelldata.mat data
end
