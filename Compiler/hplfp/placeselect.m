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

% This is the 1st of 2 ways that we are going to compute place selectivity.
% For this method, we loop through trials to compute spike counts, and
% the mean and stderr of the firing rate across trials. 
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
		% need to remove 0?s from bins in case there were any spikes that were just outside the edges of the histcounts
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

% We take the mean firing rates computed above to compute the SIC.
% compute SIC
% compute total duration for each position
if(Args.UseAllTrials)
	o_i = sum(gpDurations,2);
else
	o_i = sum(gpDurations(:,processTrials),2);
end
% Compute spike count per grid position
spikeLoc = sum(spikeLocTrial,2);
if Args.AdaptiveSmooth
    % Scaling factor
    alpha = 1e5; 
    % Restructure bins from linear to square grid
    o_iGrid = nan(Args.GridSteps);
    spikeLocGrid = nan(Args.GridSteps);
    for ii = 1:Args.GridSteps
        o_iGrid(ii,:) = o_i( (ii-1)*Args.GridSteps+1:ii*Args.GridSteps );
        spikeLocGrid(ii,:) = spikeLoc( (ii-1)*Args.GridSteps+1:ii*Args.GridSteps );
    end
    % Smooth firing rate maps
    [meanFRsGrid,spikeLocGrid_smoothed,o_iGrid_smoothed] = adaptivesmooth(o_iGrid,spikeLocGrid,alpha);
    semFRs = nanstd(anovaMatrix,0,2)./sqrt(sum(~isnan(anovaMatrix),2)); %%% CHECK???
    % Restructure bins back into linear array
    meanFRs = nan(Args.GridSteps*Args.GridSteps,1);
    for ii = 1:Args.GridSteps
        meanFRs((ii-1)*Args.GridSteps+1:ii*Args.GridSteps,1) = meanFRsGrid(ii,:);
        o_i((ii-1)*Args.GridSteps+1:ii*Args.GridSteps,1) = o_iGrid_smoothed(ii,:);
        o_i(isnan(o_i)) = 0;
    end
elseif (Args.UseMedian)
    meanFRs = nanmedian(anovaMatrix,2);
	semFRs = [prctile(anovaMatrix,25,2) prctile(anovaMatrix,75,2)];
else
    meanFRs = nanmean(anovaMatrix,2);
    semFRs = nanstd(anovaMatrix,0,2)./sqrt(sum(~isnan(anovaMatrix),2));
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

% This is the 2nd method for computing place selectivity, which compares
% the results we get from the data to shuffled surrogates. This method
% attempts to count spikes at 1 go for the entire session, instead of 
% looping over trials, to make computing place selectivity for the data and 
% the shuffled spikes more efficient.
% compare SIC to SIC from shuffled spike trains
% remove spikes that occurred after the unity program ended to avoid
% problems with the histogram function
% get max time, i.e. last marker in the session
% maxTime = um.data.unityTime(end);
% use rplparallel timestamps instead of unity timestamps
maxTime = rp.data.timeStamps(end,3);
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
% rng('shuffle');
% generate 1000 random time shifts between 0.1 and 0.9 of maxTime
tShifts = ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime;
shuffleInd = 2:shuffleSize;
spTimes(:,shuffleInd) = sort(mod(spTimes(:,shuffleInd) + tShifts,maxTime));

% get variables stored in unitymaze
sTime = um.data.sessionTime;
% get indices to sTime sorted by position
sTPi = um.data.sTPi;
% get positions with a minimum number of observations
sTPu = um.data.sTPu;
% get number of positions
nsTPu = um.data.nsTPu;
% get starting and ending indices to sTPi for each of those positions
sTPind = um.data.sTPind;
% get number of observations for each of those positions
sTPin = um.data.sTPin;
% get occupancy for each position
ou_s = um.data.ou_i;
% divide ou_s by sum of ou_s
P_i = um.data.P_i;
nFRBins = Args.NumFRBins;

% compute SIC
% compute histogram of data and shuffled data
% get corrected timeline for entire session
n = histcie(spTimes,sTime(:,1));
% divide by interval to get firing rates
fr = n ./ repmat(sTime(:,3),1,size(n,2));
% get maximum firing rate across data and surrogates
frmax = max(max(fr));
% divide range into nFRbins
frbins = 0:frmax/nFRBins:frmax;
% compute histogram of fr using frbins
[frn,frni] = histcie(fr,frbins,'DropLast');
% compute P_j from frn
P_j = frn / sum(frn);
% create array for storing mean and stderr firing rate
lambda_i = zeros(nsTPu,shuffleSize);
lambda_ise = lambda_i;

% create array for P_ij
P_ij = zeros(nFRBins,nsTPu,shuffleSize);
% bin indices should be 1 to nFRBins, so shift the limits by 0.5
frnbins = 0.5:1:(nFRBins+0.5);
for fri = 1:nsTPu
	% get relevant indices
	frindices = sTPi(sTPind(fri,1):sTPind(fri,2));
	% get firing rates for this position for data and surrogates
	frncounts = frni(frindices,:);
	% count number of rbins in each category
	P_ij(:,fri,:) = histcie(frncounts,frnbins,'DropLast','DataCols');
	frvals = fr(frindices,:);
	lambda_i(fri,:) = mean(frvals);
	lambda_ise(fri,:) = std(frvals)/sqrt(sTPin(fri));
end

% compute total duration across all positions
% So_i = sum(o_i);
% compute proportion of occupied time
% P_i = o_i/So_i;

% divide P_ij by its sum to get probability
P_ijs = sum(sum(P_ij));
P_ij = P_ij ./ P_ijs;

% create matrices for computing denominator in MI calculation
P_i_mat = repmat(P_i',[nFRBins,1,shuffleSize]);
P_j_mat = repmat(P_j,[1,nsTPu,shuffleSize]);
MI = sum(nansum(P_ij .* log2 ( P_ij ./ (P_i_mat .* P_j_mat) )));

% compute mean firing rate over all positions weighed by proportion of
% occupied time
% replace NaN with 0 in meanFRs;
% lambda_i = n2(1:gridBins,:) ./ repmat(o_i,1,shuffleSize);
% lambda_i(isnan(lambda_i)) = 0;

lambda_bar = P_i' * lambda_i;
% divide firing for each position by the overall mean
FRratio = lambda_i ./ repmat(lambda_bar,nsTPu,1);
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
data.MI = MI;

if(isdeployed)
	% save data into a mat file
	save vmplacecelldata.mat data
end


function [smoothedRate,smoothedSpk,smoothedPos] = adaptivesmooth(pos,spk,alpha)
% Adapted from rates_adaptivesmooth.m (Wills et al)
% pos = occupancy map/dwell time in each position bin (in seconds)
% spk = spike map/spike count in each position bin
% alpha = scaling parameter (1e6 for Skaggs et al 1996, 1e5 for Wills et al 2010)

% Check for empty spk maps %
if sum(sum(spk))==0
    smoothedPos=pos;    smoothedPos(pos==0)=nan;
    smoothedSpk=spk;    smoothedSpk(pos==0)=nan;
    smoothedRate=spk;   smoothedRate(pos==0)=nan;
    return
end
% Pre-assign output %
smoothedPos=zeros(size(pos));
smoothedSpk=zeros(size(pos));
% Visited env template: use this to get numbers of visited bins in filter at edge of environemnt %
vis=zeros(size(pos));
vis(pos>0)=1;
% Pre-assign map which records which bins have passed %
smoothedCheck=false(size(pos));
smoothedCheck(pos==0)=true; % Disregard unvisited - mark as already done.
% Pre-assign list of radii used (this is for reporting purposes, not used for making maps) %
radiiUsedList=nan(1,sum(sum(pos>0)));
radiiUsedCount=1;
% These parameters depend on place or dir mode %
if size(pos,2)>1
    boundary=0;             % IMFILTER boundary condition
    rBump=0.5;              % Increase radius in 0.5 bin steps.
elseif size(pos,2)==1
    boundary='circular';
    rBump=1;                % Increase radius in 1 bin steps.
end

%%% Run increasing radius iterations %%%
r=1; % Circle radius
while any(any(~smoothedCheck))
    % Check radius isn't getting too big (if >map/2, stop running) %
    if r>max(size(pos))/2
        smoothedSpk(~smoothedCheck)=nan;
        smoothedPos(~smoothedCheck)=nan;
        break
    end
    % Construct filter kernel ...
    if size(pos,2)>1
        % Place: Flat disk, where r>=distance to bin centre %
        f=fspecial('disk',r); 
        f(f>=(max(max(f))/3))=1;
        f(f~=1)=0;
    elseif size(pos,2)==1 
        % Direction: boxcar window, r bins from centre symmetrically %
        f=ones(1+(r*2),1);
    end     
    % Filter maps (get N spikes and pos sum within kernel) %
    fSpk=imfilter(spk,f,boundary);
    fPos=imfilter(pos,f,boundary);
    fVis=imfilter(vis,f,boundary);
    % Which bins pass criteria at this radius? %
    warning('off', 'MATLAB:divideByZero');
    binsPassed=alpha./(sqrt(fSpk).*fPos) <= r;
    warning('on', 'MATLAB:divideByZero');
    binsPassed=binsPassed & ~smoothedCheck; % Only get the bins that have passed in this iteration.
    % Add these to list of radii used %
    nBins=sum(binsPassed(:));
    radiiUsedList(radiiUsedCount:radiiUsedCount+nBins-1)=r;
    radiiUsedCount=radiiUsedCount+nBins;
    % Assign values to smoothed maps %
    smoothedSpk(binsPassed)=fSpk(binsPassed)./fVis(binsPassed);
    smoothedPos(binsPassed)=fPos(binsPassed)./fVis(binsPassed);
    % Record which bins were smoothed this iteration %
    smoothedCheck(binsPassed)=true;
    % Increase circle radius (half-bin steps) %
    r=r+rBump;
end

% Assign Output %
warning('off', 'MATLAB:divideByZero');
smoothedRate=smoothedSpk./smoothedPos;
warning('on', 'MATLAB:divideByZero');
smoothedRate(pos==0)=nan;
smoothedPos(pos==0)=nan;
smoothedSpk(pos==0)=nan;

% report radii sizes?