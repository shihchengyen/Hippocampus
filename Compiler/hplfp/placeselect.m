function data = placeselect(um,rp,spiketrain,Args)
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
if(Args.UseAllTrials)
	nptrials = ntrials;
else
	nptrials = size(processTrials,1);
end

sampleRate = rp.data.SampleRate;
% compute trial durations from the timestamps in rplparallel
rpTrialDur = diff(rp.data.timeStamps(:,2:3),1,2);

% Get maze coordinates for each spike during trial (cue offset to end trial)
trialcount = 0; % count total number of trials
spikecountSess = 0; % count number of spikes (per trial)

% spiketrain = load('spiketrain.mat');
% convert spike times to seconds to make it consistent with Ripple
% markers and Unity timestamps
sTimes = spiketrain.timestamps/1000;

% create memory
spikeLocTrial = zeros(gridBins,nptrials);
anovaMatrix = spikeLocTrial;
tDiff = nan(nptrials,1);

% loop over number of trials            
for npi = 1:nptrials
	if(Args.UseAllTrials)
		g = npi;
	else
	    g = processTrials(npi);
	end
    fprintf(1,'Processing Trial %i\n',g);
    
	tstart = rp.data.timeStamps(g,2);
	tend = rp.data.timeStamps(g,3);

	% get indices for this trial
	uDidx = unityTriggers(g,2):unityTriggers(g,3);
	
	% get grid positions for this trial
	tgp = gridPosition(uDidx);

	temp = find(sTimes > tstart & sTimes < tend);
    fprintf(1,'Found %i spikes between %f and %f\n',size(temp,2),tstart,tend);
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
		% compute histogram in parallel, and then extract the 1st column for data and the rest for shuffle
		spikeLocTrial(:,npi) = histcounts(tgp(bins(bins~=0)),um.data.gpEdges);
	end
end 

% if there are any grid positions with less than 5 observations
% set them to 0 to avoid having outliers caused by a single observation
spikeLocTrial( (sum(spikeLocTrial>0,2)<Args.MinTrials) ,:) = 0;
% we divide by gpDurations to get number of spikes per duration
if(Args.UseAllTrials)
	anovaMatrix = spikeLocTrial./gpDurations;
else
	anovaMatrix = spikeLocTrial./gpDurations(:,processTrials);
end
meanFRs = nanmean(anovaMatrix,2);
semFRs = nanstd(anovaMatrix,0,2)./sqrt(sum(~isnan(anovaMatrix),2));

% compute SIC
% compute total duration for each position
if(Args.UseAllTrials)
	o_i = sum(gpDurations,2);
else
	o_i = sum(gpDurations(:,processTrials),2);
end
% compute total duration across all positions
So_i = sum(o_i);
% compute proportion of occupied time
P_i = o_i/So_i;
% compute mean firing rate over all positions weighed by proportion of
% occupied time
% replace NaN with 0 in meanFRs;
lambda_i = meanFRs;
lambda_i(isnan(lambda_i)) = 0;
lambda_bar = P_i' * lambda_i;
% divide firing for each position by the overall mean
FRratio = lambda_i/lambda_bar;
% compute first term in SIC
SIC1 = P_i .* FRratio;
SIC2 = log2(FRratio);
SIC2(isinf(SIC2)) = 0;
SIC = SIC1' * SIC2;

% compare SIC to SIC from shuffled spike trains
% remove spikes that occurred after the unity program ended to avoid
% problems with the histogram function
% get max time, i.e. last marker in the session
maxTime = um.data.unityTime(end);
si = find(sTimes>maxTime);
% get number of spikes
if(isempty(si))
    % no spikes greater than unityTime so we will just use all the spikes
    nSpikes = size(sTimes,2);
else
    % ignore spikes after unityTime
    nSpikes = si(1) - 1;
end
% create shuffled spikes array
shuffleSize = Args.NumShuffles + 1;
spTimes = repmat(sTimes(1:nSpikes)',1,shuffleSize);

% generate circularly shifted spike trains
% first seed the random number generator to make sure we get a different
% sequence of numbers
rng('shuffle');
% generate 1000 random time shifts between 0.1 and 0.9 of maxTime
tShifts = ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime;
shuffleInd = 2:shuffleSize;
spTimes(:,shuffleInd) = sort(mod(spTimes(:,shuffleInd) + tShifts,maxTime));
% compute histogram of data and shuffled data
[n,bin] = histcie(spTimes,um.data.unityTime);
% compute histogram for location
[n2,bin2] = histcie(gridPosition(bin),um.data.gpEdges);
% compute SIC
% compute total duration for each position
o_i = sum(gpDurations,2);
% compute total duration across all positions
So_i = sum(o_i);
% compute proportion of occupied time
P_i = o_i/So_i;
% compute mean firing rate over all positions weighed by proportion of
% occupied time
% replace NaN with 0 in meanFRs;
lambda_i = n2(1:gridBins,:) ./ repmat(o_i,1,shuffleSize);
lambda_i(isnan(lambda_i)) = 0;
lambda_bar = P_i' * lambda_i;
% divide firing for each position by the overall mean
FRratio = lambda_i ./ repmat(lambda_bar,gridBins,1);
% compute first term in SIC
SIC1 = repmat(P_i,1,shuffleSize) .* FRratio;
SIC2 = log2(FRratio);
SIC2(isinf(SIC2)) = 0;
SICsh = sum(SIC1 .* SIC2);

data.gridSteps = um.data.gridSteps;
if(Args.FRSIC)
	data.meanFRs = lambda_i;
	data.semFRs = [];
else
	data.meanFRs = meanFRs;
	data.semFRs = semFRs;
end
data.SIC = SIC;
data.SICsh = SICsh';

if(isdeployed)
	% save data into a mat file
	save vmplacecelldata.mat data
end
