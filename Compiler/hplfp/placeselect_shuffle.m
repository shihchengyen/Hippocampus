function data = placeselect_shuffle(um,rp,spiketrain,Args)
% function data = placeselect(dirstr)

% load arguments
% we are doing this so this program can run in compiled form
% load('vpData');

% um = unitymaze('auto');
% rp = rplparallel('auto');

gridBins = um.data.gridSteps * um.data.gridSteps;
% triggers that indicate cue onset, offset, and end of trial
unityTriggers = um.data.unityTriggers;
% get grid position number
gridPosition = um.data.gridPosition;
% get duration spent at each grid position
gpDurations = um.data.gpDurations;
unityTrialTime = um.data.unityTrialTime;
ntrials = size(unityTriggers,1);
if(Args.UseAllTrials)
    processTrials = (1:ntrials)';
    nptrials = ntrials;
else
    % trial numbers that should be processed, i.e. the animal took the shortest 
    % route to the destination
    processTrials = um.data.processTrials;
    nptrials = size(processTrials,1);
end

sampleRate = rp.data.SampleRate;
% compute trial durations from the timestamps in rplparallel
rpTrialDur = diff(rp.data.timeStamps(:,2:3),1,2);

% Get maze coordinates for each spike during trial (cue offset to end trial)
trialcount = 0; % count total number of trials
spikecountSess = 0; % count number of spikes (per trial)

% This is the 1st of 2 ways that we are going to compute place selectivity.
% For this method, we loop through trials to compute spike counts, and
% the mean and stderr of the firing rate across trials. 
% loop over number of trials. Loop through this for every shuffle

% create shuffled spikes array
shuffleSize = Args.NumShuffles + 1;
% generate circularly shifted spike trains
% first seed the random number generator to make sure we get a different
% sequence of numbers
% rng('shuffle');
% generate 1000 random time shifts between 0.1 and 0.9 of maxTime
maxTime = rp.data.timeStamps(end,3);
tShifts = [0 ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
% shuffleInd = 2:shuffleSize;
% spikeLocTrialsh(shuffleInd,:,:) = sort(mod(spTimes(:,shuffleInd) + tShifts,maxTime));

% convert spike times to seconds to make it consistent with Ripple
% markers and Unity timestamps
sTimesOrig = spiketrain.timestamps/1000;
% If any spike times exceed maxTime, trim
indExceed = sTimesOrig > maxTime;
sTimesOrig(indExceed) = [];
sTimes = zeros(shuffleSize,size(sTimesOrig,2));
sTimes(1,:) = sTimesOrig;

for kk = 1:3 % For either full session (1), 1st half (2) or 2nd half (3)
    % Set up trial numbers to include depending on extent of session
    switch kk
        case 1 % Full session
            subset = processTrials;
            reps = shuffleSize;
        case 2 % First half of sesson
            x = floor((ntrials-1)/2);
            subset = processTrials(ismember(processTrials,1:x));
            reps = 1;
        case 3
            x = floor((ntrials-1)/2);
            subset = processTrials(ismember(processTrials,ntrials-x+1:ntrials));
            reps = 1;
    end
    SICsh = zeros(reps,1);
    maps_raw = zeros(gridBins,reps);
    maps_rawmedian = maps_raw;
    maps_adsmooth = maps_raw;
    maps_boxsmooth = maps_raw;
    mapsSEM_raw = maps_raw;
    mapsSEM_rawmedian = maps_raw;
    lambda_ish = maps_raw;
    % Get spatial data
    for shi = 1:reps
        if shi > 1 % Shift spike train by random amount of time
                sTimetemp = sTimes(1,:) + tShifts(1,shi);
                ind = find(sTimetemp > maxTime,1);
                if isempty(ind)
                    sTimes(shi,:) = sTimetemp;
                else
                    sTimes(shi,:) = [sTimetemp(1,ind:end)-maxTime-1 sTimetemp(1,1:ind-1)];
                end
        end
        % Show progress of shuffle
        if shi == 1 || mod(shi,100) == 0 
            fprintf(1,'Processing shuffle %i\n',shi);
        end

        % create memory
        spikeLocTrial = zeros(gridBins,size(subset,1));
        anovaMatrix = spikeLocTrial;
        tDiff = nan(size(subset,1),1);
        for npi = 1:size(subset,1)
            g = subset(npi);
            tstart = rp.data.timeStamps(g,2);
            tend = rp.data.timeStamps(g,3);

            % get indices for this trial
            uDidx = unityTriggers(g,2):unityTriggers(g,3);

            % get grid positions for this trial
            tgp = gridPosition(uDidx);

            temp = find(sTimes(shi,:) > tstart & sTimes(shi,:) < tend);
%             fprintf(1,'Found %i spikes between %f and %f\n',size(temp,2),tstart,tend);
            % get spike times aligned to cue offset time
            trialSpkTimes = sTimes(shi,temp) - tstart; % 
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
        
        % compute total duration for each position
        o_i = sum(gpDurations(:,subset),2);
        
        % Filter for low occupancy bins 
        zero_occ = o_i == 0; % unvisited grid positions
        low_occ = false(size(o_i));
        if Args.FiltLowOcc
            low_occ = (sum(gpDurations>0,2) < Args.MinTrials); % grid positions with less than 5 observations
            gpDurations(low_occ,:) = 0;
            spikeLocTrial(low_occ,:) = 0;
            o_i(low_occ,:) = 0;
        end
        
        % we divide by gpDurations to get number of spikes per duration
        anovaMatrix = spikeLocTrial./gpDurations(:,subset);
        mapTrial = nanmean(anovaMatrix,2);

        % We take the mean firing rates computed above to compute the SIC.
        % compute SIC

        % Compute spike count per grid position
        spikeLoc = sum(spikeLocTrial,2);
        
        % Restructure bins from linear to square grid
        o_iGrid = nan(Args.GridSteps);
        spikeLocGrid = nan(Args.GridSteps);
        mapGrid = nan(Args.GridSteps);
        for ii = 1:Args.GridSteps
            o_iGrid(ii,:) = o_i( (ii-1)*Args.GridSteps+1:ii*Args.GridSteps );
            spikeLocGrid(ii,:) = spikeLoc( (ii-1)*Args.GridSteps+1:ii*Args.GridSteps );
            mapGrid(ii,:) = mapTrial( (ii-1)*Args.GridSteps+1:ii*Args.GridSteps );
        end
        % Adaptive smooth scaling factor
        alpha = 1e4; % Scaling factor
        % Boxcar filter
        boxfilt = [0.0025 0.0125 0.0200 0.0125 0.0025;...
               0.0125 0.0625 0.1000 0.0625 0.0125;...
               0.0200 0.1000 0.1600 0.1000 0.0200;...
               0.0125 0.0625 0.1000 0.0625 0.0125;...
               0.0025 0.0125 0.0200 0.0125 0.0025;];
        
        % COMPUTE RATE MAPS
        % Unsmoothed
        map_raw = mapTrial;
        mapSEM_raw = nanstd(anovaMatrix,0,2)./sqrt(sum(~isnan(anovaMatrix),2)); 
%         % Unsmoothed median
%         map_rawmedian = nanmedian(anovaMatrix,2);
%         mapSEM_rawmedian = [prctile(anovaMatrix,25,2) prctile(anovaMatrix,75,2)];
        % Adaptive smoothing
        [map_adsmoothGrid,spikeLoc_adsmoothGrid,o_i_adsmoothGrid] = adaptivesmooth(o_iGrid,spikeLocGrid,alpha);
        % Restructure grid bins back into linear array
        map_adsmooth = nan(Args.GridSteps*Args.GridSteps,1);
        for ii = 1:Args.GridSteps
            map_adsmooth((ii-1)*Args.GridSteps+1:ii*Args.GridSteps,1) = map_adsmoothGrid(ii,:);
            o_i_adsmooth((ii-1)*Args.GridSteps+1:ii*Args.GridSteps,1) = o_i_adsmoothGrid(ii,:);
            o_i_adsmooth(isnan(o_i_adsmooth)) = 0;
        end
        % Boxcar smoothing
        map_boxsmooth = boxcarsmooth(mapTrial, boxfilt, low_occ);
        
        % COMPUTE SIC
        % compute total duration across all positions
        So_i = sum(o_i);
        % compute proportion of occupied time
        P_i = o_i/So_i;
        % compute mean firing rate over all positions weighed by proportion of
        % occupied time
        % replace NaN with 0 in meanFRs;
        lambda_i = map_adsmooth;
        lambda_i(isnan(lambda_i)) = 0;
        lambda_bar = P_i' * lambda_i;
        % divide firing for each position by the overall mean
        FRratio = lambda_i/lambda_bar;
        % compute first term in SIC
        SIC1 = P_i .* lambda_i; 
        ind = find(SIC1 > 0);
        SIC1 = SIC1(ind);
        SIC2 = log2(FRratio(ind));
        bits_per_sec = SIC1' * SIC2;
        if lambda_bar > 0
            bits_per_spike = bits_per_sec/lambda_bar;
        else
            bits_per_spike = NaN;
        end
        SICsh(shi,1) = bits_per_spike;
        maps_adsmooth(:,shi) = map_adsmooth;
        maps_raw(:,shi) = map_raw;
        maps_boxsmooth(:,shi) = map_boxsmooth;
        mapsSEM_raw(:,shi) = mapSEM_raw;
%         maps_rawmedian(:,shi) = map_rawmedian;
%         mapsSEM_rawmedian(:,shi) = mapSEM_rawmedian;
        lambda_ish(:,shi) = lambda_i;
    end
    % Store data
    data.gridSteps = um.data.gridSteps;
    switch kk
        case 1
            if(Args.FRSIC)
                data.maps_raw = lambda_ish(:,1);
                data.mapsSEM_raw = [];
            else
                data.maps_raw = maps_raw(:,1);
                data.maps_boxsmooth = maps_boxsmooth(:,1);
                data.maps_adsmooth = maps_adsmooth(:,1);
                data.mapsSEM_raw = mapsSEM_raw(:,1);
            end
            data.SIC = SICsh(1,1);
            data.SICsh = SICsh;
        case 2
            if(Args.FRSIC)
                data.maps_raw = lambda_ish(:,1);
                data.mapsSEM_raw = [];
            else
                data.maps_raw1sthalf = maps_raw(:,1);
                data.maps_boxsmooth1sthalf = maps_boxsmooth(:,1);
                data.maps_adsmooth1sthalf = maps_adsmooth(:,1);
                data.mapsSEM_raw1sthalf = mapsSEM_raw(:,1);
            end
            data.SIC1sthalf = SICsh(1,1);
        case 3
            if(Args.FRSIC)
                data.maps_raw = lambda_ish(:,1);
                data.mapsSEM_raw = [];
            else
                data.maps_raw2ndhalf = maps_raw(:,1);
                data.maps_boxsmooth2ndhalf = maps_boxsmooth(:,1);
                data.maps_adsmooth2ndhalf = maps_adsmooth(:,1);
                data.mapsSEM_raw2ndhalf = mapsSEM_raw(:,1);
            end
            data.SIC2ndhalf = SICsh(1,1);
    end
end

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


function [map_smooth]=boxcarsmooth(map,filt,low_occ)
% Smooth a place map. map=map, k=boxcar kernel, unvis=index of unvisited bins.
if max(size(filt))==1;  filt=ones(filt);  end   % Expand single parameter to flat k-by-k square
map(low_occ)=0;
visTemplate=ones(size(map));
visTemplate(low_occ)=0;
filtMap=imfilter(map,filt);
filtVis=imfilter(visTemplate,filt);
warning('off', 'MATLAB:divideByZero');
map_smooth=filtMap./filtVis;
warning('on', 'MATLAB:divideByZero');
map_smooth(low_occ)=nan;