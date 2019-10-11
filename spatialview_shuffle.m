function data = spatialview_shuffle(gz,rp,um,spiketrain,Args)

binDepths = [gz.data.binDepths(:,2) gz.data.binDepths(:,1)]; % x and y are flipped in Unity coordinates
ngridBins = sum(binDepths(:,1).*binDepths(:,2));

% get duration spent at each grid position per trial
trial_dur = gz.data.binnedRelGazeDur/1000;
trialTime = gz.data.trialTimestamps;
ntrials = gz.data.numTrials;
if(Args.UseAllTrials)
    processTrials = (1:ntrials)';
    nptrials = ntrials;
else
    % trial numbers that should be processed, i.e. the animal took the shortest 
    % route to the destination
    processTrials = um.data.processTrials;
    nptrials = size(processTrials,1);
end

% compute trial durations from the timestamps in rplparallel
rpTrialDur = diff(rp.data.timeStamps(:,2:3),1,2);

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
% Trim out any spike times that exceed maxTime from ripple timestamps
indExceed = sTimesOrig > maxTime;
sTimesOrig(indExceed) = [];
sTimes = zeros(shuffleSize,size(sTimesOrig,2));
sTimes(1,:) = sTimesOrig;

% Loop through shuffling process for full session, 1st half or 2nd half of
% session
for kk = 1:3 % full session (1), 1st half (2) or 2nd half (3)
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
    % Initialise output variables
    linsh_SIC = zeros(reps,1);
    linsh_map_raw = zeros(ngridBins,reps);
    linsh_map_adsmooth = linsh_map_raw;
    linsh_map_boxsmooth = linsh_map_raw;
    linsh_lambda_i = linsh_map_raw;
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

        % Bin spikes according to where gaze bins
        trial_spikeLoc = zeros(ngridBins,size(subset,1));
        tDiff = nan(size(subset,1),1);
        for npi = 1:size(subset,1)
            g = subset(npi);
            tstart = rp.data.timeStamps(g,2);
            tend = rp.data.timeStamps(g,3);

            % get indices for this trial
            uDidx = 1:gz.data.binnedRelGazeTrialInds(g);
            % get grid positions for this trial
            tgp = gz.data.binnedRelGazeTrial(uDidx,g);
            
            % get spike times aligned to cue offset time
            temp = find(sTimes(shi,:) > tstart & sTimes(shi,:) < tend);
            trialSpkTimes = sTimes(shi,temp) - tstart; 
            
            % get difference in timing between unity and ripple
            uTrialTime = (1:(trialTime(g,3)-trialTime(g,2)))'/1000; % Convert ms to s to match ripple time
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
                trial_spikeLoc(:,npi) = histcounts(tgp(bins(bins~=0)),1:sum(binDepths(:,1).*binDepths(:,2))+1);
            end
        end 
        
        % compute total duration for each position
        lin_o_i = nansum(trial_dur(:,subset),2);
        
        % Filter for low occupancy bins 
        lin_zero_occ = lin_o_i == 0; % unvisited grid positions
        lin_low_occ = false(size(lin_o_i));
        if Args.FiltLowOcc
            lin_low_occ = (sum(trial_dur>0,2) < Args.MinTrials); % grid positions with less than 5 observations
            trial_dur(lin_low_occ,:) = 0;
            trial_spikeLoc(lin_low_occ,:) = 0;
            lin_o_i(lin_low_occ,:) = 0;
        end
        
        % we divide by durTrial to get number of spikes per duration
        trial_anovaMatrix = trial_spikeLoc./trial_dur(:,subset);
        lin_map = nanmean(trial_anovaMatrix,2);

        % Compute spike count per grid position
        lin_spikeLoc = sum(trial_spikeLoc,2);
        
        % Restructure bins from linear to separate grids
        grid_o_i = cell(size(binDepths,1),1);
        grid_spikeLoc = grid_o_i;
        grid_map = grid_o_i;
        grid_low_occ = grid_o_i;
        for jj = 1:size(binDepths,1) % for each grid
            % Initialise empty matrices
            o_i = nan(binDepths(jj,1),binDepths(jj,2));
            spikeLoc = o_i;
            map = o_i;
            low_occ = o_i;
            % Assign linear bin to grid bin
            for mm = 1:binDepths(jj,1)*binDepths(jj,2) % For every point in linear map
                if mod(mm,binDepths(jj,2)) == 0
                    y = binDepths(jj,2);
                else
                    y = mod(mm,binDepths(jj,2));
                end
                x = ceil(mm/binDepths(jj,2));
                indbins_lin = mm + sum(binDepths(1:jj-1,1).*binDepths(1:jj-1,2));
                % Assign
                o_i(x,y) = lin_o_i(indbins_lin);
                spikeLoc(x,y) = lin_spikeLoc(indbins_lin);
                map(x,y) = lin_map(indbins_lin);
                low_occ(x,y) = lin_low_occ(indbins_lin);
            end
            % Collect output 
            grid_o_i{jj} = o_i;
            grid_spikeLoc{jj} = spikeLoc;
            grid_map{jj} = map;
            grid_low_occ{jj} = logical(low_occ);
        end
        
        % SMOOTH MAPS
        % Adaptive smooth scaling factor
        alpha = 1e3;
        % Boxcar filter
        boxfilt = [0.0025 0.0125 0.0200 0.0125 0.0025;...
               0.0125 0.0625 0.1000 0.0625 0.0125;...
               0.0200 0.1000 0.1600 0.1000 0.0200;...
               0.0125 0.0625 0.1000 0.0625 0.0125;...
               0.0025 0.0125 0.0200 0.0125 0.0025;];
       % Initialise output matrices
        grid_map_adsm = cell(size(grid_map));
        grid_map_boxsm = grid_map_adsm;
        if shi == 1
            grid_sm_radii = grid_map_adsm;
        end
        % Smooth
        for jj = 1:size(grid_map,1) % For each separate grid
            
            if binDepths(jj,1)*binDepths(jj,2) > 2 % For non-cue/non-hint grids
                o_i = grid_o_i{jj};
                spikeLoc = grid_spikeLoc{jj};
                map = grid_map{jj};
                % Pad each grid map with adjoining bins from other grids
                % Pad with <<5>> extra bin rows
                n = 5;
                [retrievemap,o_i,spikeLoc,map] = padgrids(n,o_i,spikeLoc,map,grid_o_i,grid_spikeLoc,grid_map,gz.data.gazeSections,jj);
                % Adaptive smooth
                [map_adsm,~,~,smooth_r] = adaptivesmooth(o_i,spikeLoc,alpha);
                % Boxcar smooth  
                map_boxsm = boxcarsmooth(map, boxfilt, grid_low_occ{jj});
                % Store smoothing radii for posthoc analysis
                if shi == 1
                    if isempty(smooth_r)
                        grid_sm_radii{jj} = NaN;
                    else
                        grid_sm_radii{jj} = smooth_r;
                    end
                end
                % Remove padding from map
                map_adsm = map_adsm(retrievemap(1,1):retrievemap(1,2),retrievemap(2,1):retrievemap(2,2));
                map_boxsm = map_boxsm(retrievemap(1,1):retrievemap(1,2),retrievemap(2,1):retrievemap(2,2));

            else
                map_adsm = grid_map{jj};
                map_boxsm = grid_map{jj};
                if shi == 1
                    grid_sm_radii{jj} = NaN;
                end
            end
            
            % Collect output
            grid_map_adsm{jj} = map_adsm;
            grid_map_boxsm{jj} = map_boxsm;
        end
        
        % Restructure grid bins back into linear array
        lin_map_adsm = nan(size(lin_map));
        lin_map_boxsm = lin_map_adsm;
        for ii = 1:size(binDepths,1)
            if ii > 1
                % index of corresponding linear bin
                ind_lin = sum(binDepths(1:ii-1,1).*binDepths(1:ii-1,2))+1:binDepths(ii,1)*binDepths(ii,2)+sum(binDepths(1:ii-1,1).*binDepths(1:ii-1,2));
                % Reshape grid so as to concatenate each row from left to right
                a = permute(grid_map_adsm{ii},[2,1]);
                % Fill in adaptive smooth map
                lin_map_adsm(ind_lin,1) = a(:);
                % Reshape grid so as to concatenate each row from left to right
                b = permute(grid_map_boxsm{ii},[2,1]);
                % Fill in boxsmooth map
                lin_map_boxsm(ind_lin,1) = b(:);
            else
                lin_map_adsm(ii,1) = grid_map_adsm{ii};
            end
        end
        
        % COMPUTE SIC
        % compute total duration across all positions
        So_i = sum(lin_o_i);
        % compute proportion of occupied time
        P_i = lin_o_i/So_i;
        % compute mean firing rate over all positions weighed by proportion of
        % occupied time
        % replace NaN with 0 in meanFRs;
        lambda_i = lin_map_adsm;
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
        linsh_SIC(shi,1) = bits_per_spike;
        linsh_map_adsmooth(:,shi) = lin_map_adsm;
        linsh_map_raw(:,shi) = lin_map;
        linsh_map_boxsmooth(:,shi) = lin_map_boxsm;
        linsh_lambda_i(:,shi) = lambda_i;
    end
    % Store data
    data.binDepths = binDepths;
    switch kk
        case 1
            
            data.maps_raw = linsh_map_raw(:,1);
            data.maps_boxsmooth = linsh_map_boxsmooth(:,1);
            data.maps_adsmooth = linsh_map_adsmooth(:,1);

            data.SIC = linsh_SIC(1,1);
            data.SICsh = linsh_SIC;
            data.smoothingradii = grid_sm_radii';
        case 2
            
            data.maps_raw1sthalf = linsh_map_raw(:,1);
            data.maps_boxsmooth1sthalf = linsh_map_boxsmooth(:,1);
            data.maps_adsmooth1sthalf = linsh_map_adsmooth(:,1);
                
            data.SIC1sthalf = linsh_SIC(1,1);
            data.smoothingradii1sthalf = grid_sm_radii';
        case 3
            
            data.maps_raw2ndhalf = linsh_map_raw(:,1);
            data.maps_boxsmooth2ndhalf = linsh_map_boxsmooth(:,1);
            data.maps_adsmooth2ndhalf = linsh_map_adsmooth(:,1);
                
            data.SIC2ndhalf = linsh_SIC(1,1);
            data.smoothingradii2ndhalf = grid_sm_radii';
    end
end

if(isdeployed)
	% save data into a mat file
	save spatialviewdata.mat data
end


function [smoothedRate,smoothedSpk,smoothedPos,radiiUsedList] = adaptivesmooth(pos,spk,alpha)
% Adapted from rates_adaptivesmooth.m (Wills et al)
% pos = occupancy map/dwell time in each position bin (in seconds)
% spk = spike map/spike count in each position bin
% alpha = scaling parameter (1e6 for Skaggs et al 1996, 1e5 for Wills et al 2010)

% Check for empty spk maps %
if sum(sum(spk))==0
    smoothedPos=pos;    smoothedPos(pos==0)=nan;
    smoothedSpk=spk;    smoothedSpk(pos==0)=nan;
    smoothedRate=spk;   smoothedRate(pos==0)=nan;
    radiiUsedList=nan(1,sum(sum(pos>0)));
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

function [retrievemap,o_i,spikeLoc,map] = padgrids(n,o_i,spikeLoc,map,grid_o_i,grid_spikeLoc,grid_map,gazeSections,jj)

% Pad maps with adjoining bins from adjacent maps

switch gazeSections{jj}
    case 'Ground'
        wallsection_ind = strcmp(gazeSections,'Walls');
        wall_o_i = grid_o_i{wallsection_ind};
        wall_spikeLoc = grid_spikeLoc{wallsection_ind};
        wall_map = grid_map{wallsection_ind};
        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = spikeLoc;
        map_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        map_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = map;
        % Pad with wall data
        o_i_temp(1:n,n+1:n+size(o_i,1)) = wall_o_i(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1)); % top
        o_i_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1)),-1); % right
        o_i_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1)),-2); % bottom
        o_i_temp(n+1:size(o_i,1)+n,1:n) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1)),1); % left
        spikeLoc_temp(1:n,n+1:n+size(o_i,1)) = wall_spikeLoc(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1)); % top
        spikeLoc_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1)),-1); % right
        spikeLoc_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1)),-2); % bottom
        spikeLoc_temp(n+1:size(o_i,1)+n,1:n) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1)),1); % left
        map_temp(1:n,n+1:n+size(o_i,1)) = wall_map(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1)); % top
        map_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end) = rot90(wall_map(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1)),-1); % right
        map_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n) = rot90(wall_map(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1)),-2); % bottom
        map_temp(n+1:size(o_i,1)+n,1:n) = rot90(wall_map(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1)),1); % left
        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [n+1 n+size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;
        map = map_temp;
    case 'Ceiling'
        wallsection_ind = strcmp(gazeSections,'Walls');
        wall_o_i = grid_o_i{wallsection_ind};
        wall_spikeLoc = grid_spikeLoc{wallsection_ind};
        wall_map = grid_map{wallsection_ind};
        % Flip walldata upside down
        wall_o_i = flipud(wall_o_i);
        wall_spikeLoc = flipud(wall_spikeLoc);
        wall_map = flipud(wall_map);
        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = spikeLoc;
        map_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        map_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = map;
        % Pad with wall data
        o_i_temp(1:n,n+1:n+size(o_i,1)) = fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1))); % top
        o_i_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1))),-1); % right
        o_i_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1))),-2); % bottom
        o_i_temp(n+1:size(o_i,1)+n,1:n) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1))),1); % left
        spikeLoc_temp(1:n,n+1:n+size(o_i,1)) = fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1))); % top
        spikeLoc_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1))),-1); % right
        spikeLoc_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1))),-2); % bottom
        spikeLoc_temp(n+1:size(o_i,1)+n,1:n) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1))),1); % left
        map_temp(1:n,n+1:n+size(o_i,1)) = fliplr(wall_map(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1))); % top
        map_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end) = rot90(fliplr(wall_map(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1))),-1); % right
        map_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n) = rot90(fliplr(wall_map(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1))),-2); % bottom
        map_temp(n+1:size(o_i,1)+n,1:n) = rot90(fliplr(wall_map(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1))),1); % left
        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [n+1 n+size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;
        map = map_temp;
    case 'Walls'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};
        ground_map = grid_map{groundsection_ind};
        ceilingsection_ind = strcmp(gazeSections,'Ceiling');
        ceiling_o_i = grid_o_i{ceilingsection_ind};
        ceiling_spikeLoc = grid_spikeLoc{ceilingsection_ind};
        ceiling_map = grid_map{ceilingsection_ind};
        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = spikeLoc;
        map_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n);
        map_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2)) = map;
        % Pad with ground data
        o_i_temp(n+size(o_i,1)+1:end,n+1:size(ground_o_i,2)+n) = rot90(ground_o_i(:,1:n),-1);
        o_i_temp(n+size(o_i,1)+1:end,n+size(ground_o_i,2)+1:n+2*size(ground_o_i,2)) = ground_o_i(1:n,:);
        o_i_temp(n+size(o_i,1)+1:end,n+2*size(ground_o_i,2)+1:n+3*size(ground_o_i,2)) = rot90(ground_o_i(:,size(ground_o_i,1)-n+1:end),1);
        o_i_temp(n+size(o_i,1)+1:end,n+3*size(ground_o_i,1)+1:n+4*size(ground_o_i,1)) = rot90(ground_o_i(size(ground_o_i,1)-n+1:end,:),2);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+1:size(ground_o_i,2)+n) = rot90(ground_spikeLoc(:,1:n),-1);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+size(ground_o_i,2)+1:n+2*size(ground_o_i,2)) = ground_spikeLoc(1:n,:);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+2*size(ground_o_i,2)+1:n+3*size(ground_o_i,2)) = rot90(ground_spikeLoc(:,size(ground_spikeLoc,1)-n+1:end),1);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+3*size(ground_o_i,1)+1:n+4*size(ground_o_i,1)) = rot90(ground_spikeLoc(size(ground_spikeLoc,1)-n+1:end,:),2);
        map_temp(n+size(o_i,1)+1:end,n+1:size(ground_o_i,2)+n) = rot90(ground_map(:,1:n),-1);
        map_temp(n+size(o_i,1)+1:end,n+size(ground_o_i,2)+1:n+2*size(ground_o_i,2)) = ground_map(1:n,:);
        map_temp(n+size(o_i,1)+1:end,n+2*size(ground_o_i,2)+1:n+3*size(ground_o_i,2)) = rot90(ground_map(:,size(ground_map,1)-n+1:end),1);
        map_temp(n+size(o_i,1)+1:end,n+3*size(ground_o_i,1)+1:n+4*size(ground_o_i,1)) = rot90(ground_map(size(ground_map,1)-n+1:end,:),2);
        % Pad with ceiling data
        o_i_temp(1:n,n+1:size(ceiling_o_i,1)+n) = fliplr(rot90(ceiling_o_i(:,size(ceiling_o_i,1)-n+1:end),1));
        o_i_temp(1:n,n+size(ceiling_o_i,1)+1:n+2*size(ceiling_o_i,1)) = fliplr(ceiling_o_i(1:n,:));
        o_i_temp(1:n,n+2*size(ceiling_o_i,1)+1:n+3*size(ceiling_o_i,1)) = fliplr(rot90(ceiling_o_i(:,1:n),-1));
        o_i_temp(1:n,n+3*size(ceiling_o_i,1)+1:n+4*size(ceiling_o_i,1)) = fliplr(rot90(ceiling_o_i(size(ceiling_o_i,1)-n+1:end,:),2));
        spikeLoc_temp(1:n,n+1:size(ceiling_o_i,1)+n) = fliplr(rot90(ceiling_spikeLoc(:,size(ceiling_spikeLoc,1)-n+1:end),1));
        spikeLoc_temp(1:n,n+size(ceiling_o_i,1)+1:n+2*size(ceiling_o_i,1)) = fliplr(ceiling_spikeLoc(1:n,:));
        spikeLoc_temp(1:n,n+2*size(ceiling_o_i,1)+1:n+3*size(ceiling_o_i,1)) = fliplr(rot90(ceiling_spikeLoc(:,1:n),-1));
        spikeLoc_temp(1:n,n+3*size(ceiling_o_i,1)+1:n+4*size(ceiling_o_i,1)) = fliplr(rot90(ceiling_spikeLoc(size(ceiling_spikeLoc,1)-n+1:end,:),2));
        map_temp(1:n,n+1:size(ceiling_o_i,1)+n) = fliplr(rot90(ceiling_map(:,size(ceiling_map,1)-n+1:end),1));
        map_temp(1:n,n+size(ceiling_o_i,1)+1:n+2*size(ceiling_o_i,1)) = fliplr(ceiling_map(1:n,:));
        map_temp(1:n,n+2*size(ceiling_o_i,1)+1:n+3*size(ceiling_o_i,1)) = fliplr(rot90(ceiling_map(:,1:n),-1));
        map_temp(1:n,n+3*size(ceiling_o_i,1)+1:n+4*size(ceiling_o_i,1)) = fliplr(rot90(ceiling_map(size(ceiling_map,1)-n+1:end,:),2));
        % Pad with wall data on either end
        o_i_temp(n+1:n+size(o_i,1),1:n) = o_i(:,size(o_i,2)-n+1:end);
        o_i_temp(n+1:n+size(o_i,1),size(o_i_temp,2)-n+1:end) = o_i(:,1:n);
        spikeLoc_temp(n+1:n+size(o_i,1),1:n) = spikeLoc(:,size(o_i,2)-n+1:end);
        spikeLoc_temp(n+1:n+size(o_i,1),size(o_i_temp,2)-n+1:end) = spikeLoc(:,1:n);
        map_temp(n+1:n+size(o_i,1),1:n) = map(:,size(map,2)-n+1:end);
        map_temp(n+1:n+size(o_i,1),size(o_i_temp,2)-n+1:end) = map(:,1:n);
        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [n+1 n+size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;
        map = map_temp;
    case 'Pillar1'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};
        ground_map = grid_map{groundsection_ind};
        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = spikeLoc;
        map_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        map_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = map;
        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_o_i(25:32,25-n:24),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_o_i(25-n:24,25:32);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_o_i(25:32,33:32+n),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_o_i(33:32+n,25:32),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_spikeLoc(25:32,25-n:24),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_spikeLoc(25-n:24,25:32);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_spikeLoc(25:32,33:32+n),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_spikeLoc(33:32+n,25:32),2);
        map_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_map(25:32,25-n:24),-1);
        map_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_map(25-n:24,25:32);
        map_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_map(25:32,33:32+n),1);
        map_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_map(33:32+n,25:32),2);
        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n) = o_i(:,size(o_i,2)-n+1:end);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = o_i(:,1:n);
        spikeLoc_temp(1:size(o_i,1),1:n) = spikeLoc(:,size(o_i,2)-n+1:end);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = spikeLoc(:,1:n);
        map_temp(1:size(o_i,1),1:n) = map(:,size(map,2)-n+1:end);
        map_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = map(:,1:n);
        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;
        map = map_temp;
    case 'Pillar2'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};
        ground_map = grid_map{groundsection_ind};
        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = spikeLoc;
        map_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        map_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = map;
        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_o_i(25:32,9-n:8),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_o_i(25-n:24,9:16);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_o_i(25:32,17:16+n),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_o_i(33:32+n,9:16),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_spikeLoc(25:32,9-n:8),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_spikeLoc(25-n:24,9:16);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_spikeLoc(25:32,17:16+n),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_spikeLoc(33:32+n,9:16),2);
        map_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_map(25:32,9-n:8),-1);
        map_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_map(25-n:24,9:16);
        map_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_map(25:32,17:16+n),1);
        map_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_map(33:32+n,9:16),2);
        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n) = o_i(:,size(o_i,2)-n+1:end);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = o_i(:,1:n);
        spikeLoc_temp(1:size(o_i,1),1:n) = spikeLoc(:,size(o_i,2)-n+1:end);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = spikeLoc(:,1:n);
        map_temp(1:size(o_i,1),1:n) = map(:,size(map,2)-n+1:end);
        map_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = map(:,1:n);
        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;
        map = map_temp;
    case 'Pillar3'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};
        ground_map = grid_map{groundsection_ind};
        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = spikeLoc;
        map_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        map_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = map;
        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_o_i(9:16,25-n:24),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_o_i(9-n:8,25:32);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_o_i(9:16,33:32+n),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_o_i(17:16+n,25:32),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_spikeLoc(9:16,25-n:24),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_spikeLoc(9-n:8,25:32);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_spikeLoc(9:16,33:32+n),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_spikeLoc(17:16+n,25:32),2);
        map_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_map(9:16,25-n:24),-1);
        map_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_map(9-n:8,25:32);
        map_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_map(9:16,33:32+n),1);
        map_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_map(17:16+n,25:32),2);
        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n) = o_i(:,size(o_i,2)-n+1:end);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = o_i(:,1:n);
        spikeLoc_temp(1:size(o_i,1),1:n) = spikeLoc(:,size(o_i,2)-n+1:end);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = spikeLoc(:,1:n);
        map_temp(1:size(o_i,1),1:n) = map(:,size(map,2)-n+1:end);
        map_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = map(:,1:n);
        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;
        map = map_temp;
    case 'Pillar4'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};
        ground_map = grid_map{groundsection_ind};
        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = spikeLoc;
        map_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n);
        map_temp(1:size(o_i,1), n+1:n+size(o_i,2)) = map;
        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_o_i(9:16,9-n:8),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_o_i(9-n:8,9:16);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_o_i(9:16,17:16+n),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_o_i(17:16+n,9:16),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_spikeLoc(9:16,9-n:8),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_spikeLoc(9-n:8,9:16);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_spikeLoc(9:16,17:16+n),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_spikeLoc(17:16+n,9:16),2);
        map_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n) = rot90(ground_map(9:16,9-n:8),-1);
        map_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4)) = ground_map(9-n:8,9:16);
        map_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4)) = rot90(ground_map(9:16,17:16+n),1);
        map_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4)) = rot90(ground_map(17:16+n,9:16),2);
        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n) = o_i(:,size(o_i,2)-n+1:end);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = o_i(:,1:n);
        spikeLoc_temp(1:size(o_i,1),1:n) = spikeLoc(:,size(o_i,2)-n+1:end);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = spikeLoc(:,1:n);
        map_temp(1:size(o_i,1),1:n) = map(:,size(map,2)-n+1:end);
        map_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end) = map(:,1:n);
        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;
        map = map_temp;
end