function [obj, varargout] = vmpc(varargin)
%@vmpc Constructor function for vmpc class
%   OBJ = vmpc(varargin)
%
%   OBJ = vmpc('auto') attempts to create a vmpc object by ...
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmpc %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmpc('save','redo')
%
%dependencies:

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
                'ObjectLevel','Cell', 'RequiredFile','spiketrain.mat', ...
                'GridSteps',40, 'pix',1,...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4,'SmoothType','Adaptive', 'UseMinObs',0, 'ThresVel',1, 'UseAllTrials',1,...
                'SelectiveCriteria','SIC','Alpha', 10000);
            
Args.flags = {'Auto','ArgsOnly','FRSIC','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps','NumShuffles','UseMinObs','SmoothType','ThresVel','UseAllTrials', 'Alpha'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
    'subtract',{'RedoLevels','SaveLevels'}, ...
    'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
    'remove',{'Auto'});


Args.classname = 'vmpc';
filename = hashFileName(Args);
Args.matname = sprintf('%s_%s.mat', Args.classname, filename);
Args.matvarname = 'vmp';
% testing hash functions
% test_hashFileName;
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
    % NEW: if object did not exist, create the object instead
    if isfile(Args.matname)
        obj = robj; 
    else
        disp(['Object not existed for' Args.matname, ', create object instead']);
        obj = createObject(Args,modvarargin{:});
    end
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!!
    % If there is additional requirements for creating the object, add
    % whatever needed here
    % NEW: if the object already existed, avoid creating a new one and load
    % instead
    if isfile(Args.matname)
        disp(['Object already existed for' Args.matname]);
        obj = robj; 
    else
        obj = createObject(Args,modvarargin{:});
    end
end

function obj = createObject(Args,varargin)

% example object
dlist = nptDir;
% get entries in directory
dnum = size(dlist,1);

% added code here to create a hash map for the required arguments
% and check if the new object has the arguments
argsMap = containers.Map();
argFields = Args.DataCheckArgs;
for i = 1:length(argFields)
     argsMap(argFields{i}) = Args.(argFields{i});
end

requiredKeys = Args.DataCheckArgs;
missingKeys = {};

for i = 1:length(requiredKeys)
    if ~isKey(argsMap, requiredKeys{i})
        missingKeys{end+1} = requiredKeys{i};
    end
end

if ~isempty(missingKeys)
    error(['Missing required arguments in Args: ' strjoin(missingKeys, ',')]);
end
% check if the right conditions were met to create object
if(~isempty(dir(Args.RequiredFile)))
    
    ori = pwd;

    data.origin = {pwd};
%     pv = vmpv('auto', varargin{:});

    %%%% PATCH
    cd ..; cd ..; cd ..;
    pv = load([num2str(Args.pix) 'vmpv.mat']);
    % pv = load('mac_vmpv.mat');
    pv = pv.pv;
    %%%%%%%

    cd(ori);
    [d,fn,ext] = fileparts(Args.RequiredFile);
    if strcmp(ext,'.mat')
        spiketrain = load(Args.RequiredFile);
    elseif strcmp(ext,'.csv')
        spiketrain.timestamps = load(Args.RequiredFile)';
    else
        error(['Unkown file type ' + Args.RequiredFile]);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    NumShuffles_saved = Args.NumShuffles;

    for repeat = 1:3 % 1 = full trial, 2 = 1st half, 3 = 2nd half
        
        if repeat == 1
            disp('Full session:');
        elseif repeat == 2
            disp('First half:');
        elseif repeat == 3
            disp('Second half:');
        end
        
        if repeat > 1
            Args.NumShuffles = 0;
        end

        if repeat == 1
            stc = pv.data.sessionTimeC;
        end
        
        % spike shuffling

        spiketimes = spiketrain.timestamps/1000; % now in seconds
        maxTime = pv.data.rplmaxtime;
        % spiketimes(spiketimes>maxTime) = [];
        tShifts = [0 ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
        full_arr = repmat(spiketimes, Args.NumShuffles+1, 1);
        full_arr = full_arr + tShifts';
        keepers = length(spiketimes) - sum(full_arr>maxTime, 2);
        for row = 2:size(full_arr,1)
            full_arr(row,:) = [full_arr(row,1+keepers(row):end)-maxTime full_arr(row,1:keepers(row))];
        end
        flat_spiketimes = NaN(2,size(full_arr,1)*size(full_arr,2));
        temp = full_arr';
        flat_spiketimes(1,:) = temp(:);
        flat_spiketimes(2,:) = repelem(1:size(full_arr,1), size(full_arr,2));
        flat_spiketimes = flat_spiketimes';
        flat_spiketimes = sortrows(flat_spiketimes);

        flat_spiketimes(flat_spiketimes(:,1) < stc(1,1),:) = [];
        
        % selecting rows from sessionTimeC
        if repeat == 1
            disp('      Filtering...');
            stc(:,5) = [diff(stc(:,1)); 0];
            stc(:,6) = zeros(size(stc,1),1); % For spike binning
        end
        
        conditions = ones(size(stc,1),1);

        if Args.UseAllTrials == 0
            conditions = conditions & pv.data.good_trial_markers;
        end
        
        if repeat == 2
            conditions = conditions & (pv.data.halving_markers==1);
        elseif repeat == 3
            conditions = conditions & (pv.data.halving_markers==2);
        end

        if Args.ThresVel > 0
            conditions = conditions & get(pv,'SpeedLimit',Args.ThresVel);
        end
        
        if Args.UseMinObs
            bins_sieved = pv.data.place_good_bins;
            conditions = conditions & (pv.data.pv_good_rows); % Make sure minobs take into account both place and view
        else
            bins_sieved = 1:(Args.GridSteps * Args.GridSteps);
        end

        if repeat == 1
            % Group into intervals those consecutive rows where same place bin is occupied
            dstc = diff(stc(:,1));
            stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
            stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
            stc_changing_ind = stc_changing_ind(1:end-1,:);
        end
        
        consol_arr = zeros(Args.GridSteps * Args.GridSteps,Args.NumShuffles + 1);
        
        if repeat == 1
            disp(['Assigning '  num2str(size(flat_spiketimes,1)) ' spikes to bins...']);
        end
        interval = 1;
        for sp = 1:size(flat_spiketimes,1)

            while interval < size(stc_changing_ind,1)
                if flat_spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && ... % > start timestamp of this interval but < start timestamp of next
                        flat_spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
                    break;
                end
                interval = interval + 1; % didn't fall in this interval, search in the next interval
            end

            % Bin all spikes into stc, unfiltered. If > 1 row for same time sample (i.e. large view cone), add spike to last row, backfill later
            if flat_spiketimes(sp,2) == 1
                stc(stc_changing_ind(interval,2),6) = stc(stc_changing_ind(interval,2),6) + 1;
            end
            % Keep only bins that meet filter criteria and have all of place, view, and hd data
            bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),[2 3 4]); % find the relevant place and view bin
            bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))),:); % take out bins that don't satisfy filters
            bins_hit(~(bins_hit(:,1)>0),:) = []; % take out bins where place bin = 0
            bins_hit(~(bins_hit(:,3)>0),:) = []; % take out bins where view bin = nan
            bins_hit(~(bins_hit(:,2)>0),:) = []; % take out bins where HD bin = 0
            consol_arr(bins_hit(:,1),flat_spiketimes(sp,2)) = consol_arr(bins_hit(:,1),flat_spiketimes(sp,2)) + 1;

        end
        
        spike_count_full = consol_arr';

        %% This portion for place-related calculations

            % Remove non-place and non-view rows for duration
            stc_filt = stc(find(conditions==1),:);
            stc_filt(~(stc_filt(:,2) > 0),:) = []; % remove place bin = 0
            stc_filt(isnan(stc_filt(:,4)),:) = []; % remove NaN view bins
            stc_filt(~(stc_filt(:,3) > 0),:) = []; % remove hd bin = 0
            stc_ss = stc_filt(:,[2 5]); % [place dur];
            stc_ss = [stc_ss; [1600 0]];
    
            gpdurfull = accumarray(stc_ss(:,1),stc_ss(:,2))';

        % %% This portion for combined sessionTimeC with view to output for later mixed sel calculations
        %
        %     fillindex = false(size(stc,1),1);
        %     % back-filling spikes for view bins that occupy the same time bin
        %     stcfill = stc;
        %     stcfill(stcfill(:,6)==0,6) = nan;
        %     stcfill(:,7) = stcfill(:,5)~=0;
        %     stcfill(isnan(stcfill(:,6)) & stcfill(:,7), 6) = 0;
        %     stcfill(:,7) = [];
        %     fillindex(isnan(stcfill(1:end-1,6)),1) = true;
        %     stcfill(:,6) = fillmissing(stcfill(:,6), 'next');
        %     stcfill(isnan(stcfill(:,6)),6) = 0;
        %     % back-filling duration for view bins that occupy the same time bin
        %     stc_lasttime = stcfill(end,5); % Make sure if last duration sample is zero, it remains zero, not nan
        %     stcfill(stcfill(:,5)==0,5) = nan;
        %     stcfill(end,5) = stc_lasttime;
        %     % fillindex(isnan(stcfill(1:end-1,5)),2) = true;
        %     stcfill(:,5) = fillmissing(stcfill(:,5), 'next'); % [timestamp place hd view dur spk]
        %
        %     % Remove non-place and non-view rows for duration
        %     stcfill_filt = stcfill(find(conditions==1),:);
        %     stcfill_filt(~(stcfill_filt(:,2) > 0),:) = []; % remove place bin = 0
        %     stcfill_filt(isnan(stcfill_filt(:,4)),:) = []; % remove NaN view bins
        %     stcfill_filt(~(stcfill_filt(:,3) > 0),:) = []; % remove hd bin = 0
        
        % Remove low observation bins
        spikes_count = zeros(Args.NumShuffles+1,Args.GridSteps*Args.GridSteps);
        dur_raw = zeros(1,Args.GridSteps*Args.GridSteps);
        spikes_count(:,bins_sieved) = spike_count_full(:,bins_sieved);
        dur_raw(1,bins_sieved) = gpdurfull(1,bins_sieved);
        
        maps_raw = spikes_count./repmat(dur_raw,size(spikes_count,1),1);
        
        % Save raw maps
        map_raw = maps_raw(1,:);
        spk_raw = spikes_count(1,:);
        if repeat == 1
            data.sessionTimeC = stc; % Full unfiltered
            % data.sessionTimeCfill = stcfill; % Full unfiltered, backfilled for view
            data.stcfilt = stc_filt; % Filtered for conditions, must have place, view, and hd data
            % data.stcfillfilt = stcfill_filt; % Filtered for conditions, must have place, view, and hd data, backfilled for view
            data.maps_raw = map_raw;
            data.dur_raw = dur_raw;
            data.spk_raw = spk_raw;
            data.filtspknum = sum(spk_raw);
            % data.fillindex = fillindex; % rows that were backfilled
        elseif repeat == 2
            data.maps_raw1 = map_raw;
            data.dur_raw1 = dur_raw;
            data.spk_raw1 = spk_raw;
        elseif repeat == 3
            data.maps_raw2 = map_raw;
            data.dur_raw2 = dur_raw;
            data.spk_raw2 = spk_raw;
        end
        
            if 1 % Smoothing
                
                if repeat == 1
                    disp('Adaptive smoothing...');
                end
                nan_track = isnan(map_raw);

                alpha = Args.Alpha;

                % smoothing part here, need to reshape to 3d matrix
                % 1. add in nan values for pillar positions (variables with ones suffix)
                % 2. reshape each row to 5x5
                % after permute step, now structured 5x5x10001, with each grid in a
                % slice as following:
                %
                % 1 6 11 16 21
                % 2 - 12 -  22
                % 3 8 13 18 23
                % 4 - 14 -  24
                % 5 10 15 20 25
                %
                % but will be reverted back to usual linear representation by the
                % end of the smoothing chunk
                
               

                [maps_adsm, durs_adsm, rad_adsm, maps_bcsm, maps_dksm, durs_bcsm, durs_dksm,rad_adsm_grid] = smoothMaps(maps_raw, dur_raw, spk_raw, spikes_count, Args);

               

                % smoothing part ends
                switch Args.SmoothType
                    case 'Adaptive'
                        maps_sm = maps_adsm;
                    case 'Boxcar'
                        maps_sm = maps_bcsm;
                    case 'Disk'
                        maps_sm = maps_dksm;
                end
                
                if repeat == 1
                    if data.filtspknum < 100
                        data.discard = true;
                    else
                        data.discard = false;
                    end
                    if max(maps_sm(1,:),[],'omitnan') < 0.7
                        data.rateok = false;
                    else
                        data.rateok = true;
                    end
                    data.maps_adsm = maps_adsm(1,:);
                    data.maps_adsmsh = maps_adsm(2:end,:);
                    data.dur_adsm = durs_adsm(1,:);
                    data.dur_adsmsh = durs_adsm(2:end,:);
                    data.radii = rad_adsm_grid(1,:);
                    data.radiish = rad_adsm_grid(2:end,:);
                    data.maps_bcsm = maps_bcsm(1,:);
                    data.maps_bcsmsh = maps_bcsm(2:end,:);
                    data.maps_dksm = maps_dksm(1,:);
                    data.maps_dksmsh = maps_dksm(2:end,:);
                    data.maps_sm = maps_sm(1,:);
                    data.maps_smsh = maps_sm(2:end,:);
                elseif repeat == 2
                    data.maps_adsm1 = maps_adsm(1,:);
                    data.dur_adsm1 = durs_adsm(1,:);
                    data.radii1 = rad_adsm_grid(1,:);
                    data.maps_bcsm1 = maps_bcsm(1,:);
                    data.maps_dksm1 = maps_dksm(1,:);
                    data.maps_sm1 = maps_sm(1,:);
                    data.maps_smsh1 = maps_sm(2:end,:);
                elseif repeat == 3
                    data.maps_adsm2 = maps_adsm(1,:);
                    data.dur_adsm2 = durs_adsm(1,:);
                    data.radii2 = rad_adsm_grid(1,:);
                    data.maps_bcsm2 = maps_bcsm(1,:);
                    data.maps_dksm2 = maps_dksm(1,:);
                    data.maps_sm2 = maps_sm(1,:);
                    data.maps_smsh2 = maps_sm(2:end,:);
                end

            else
                maps_adsm = maps_raw;
                durs_adsm = repmat(dur_raw,Args.NumShuffles+1,1); % HM added

            end

            % Calculating SIC for adaptive smoothing
            disp('Calculating SIC...');
            sic_adsm = skaggs_sic(maps_adsm',durs_adsm');
            sic_adsm = sic_adsm';
            
            % Calculating SIC for boxcar smoothing
            sic_bcsm = skaggs_sic(maps_bcsm',durs_bcsm');
            sic_bcsm = sic_bcsm';
            
            % Calculating SIC for disk smoothing
            sic_dksm = skaggs_sic(maps_dksm',durs_dksm');
            sic_dksm = sic_dksm';

%             % ISE part
%             lambda_i = firing_rates_full;
%
%             if repeat == 1
%                 ise_adsm = ise(lambda_i(1,:), lambda_i(2:end,:), Args.GridSteps, Args.GridSteps);
%                 data.ISE_sm = ise_adsm(1);
%                 data.ISEsh_sm = ise_adsm(2:end,1);
%             elseif repeat == 2
%                 ise_sm = ise(lambda_i, [], Args.GridSteps, Args.GridSteps);
%                 data.ISE_sm1 = ise_adsm;
%             elseif repeat == 3
%                 ise_sm = ise(lambda_i, [], Args.GridSteps, Args.GridSteps);
%                 data.ISE_sm2 = ise_adsm;
%             end
        
        switch Args.SmoothType
            case 'Adaptive'
                sic_sm = sic_adsm;
            case 'Boxcar'
                sic_sm = sic_bcsm;
            case 'Disk'
                sic_sm = sic_dksm;
        end
        
        switch Args.SelectiveCriteria
            case 'SIC'
                crit_sm = sic_sm;
            case 'ISE'
                
        end

        % Calculate sparsity
        sparsity = spatial_sparsity(dur_raw,map_raw);

        % Calculate selectivity (signal-to-noise)
        sig2noise = spatial_sig2noise(map_raw);

        % Calculate coherence using raw map

        coherence = spatial_coherence('place',[Args.GridSteps Args.GridSteps],map_raw,1); % raw
        coherence_sm = spatial_coherence('place',[Args.GridSteps Args.GridSteps],maps_adsm(1,:),1); % adsm, cheat
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if repeat == 1
            data.SIC_adsm = sic_adsm(1);
            data.SICsh_adsm = sic_adsm(2:end,1);
            data.SIC_bcsm = sic_bcsm(1);
            data.SICsh_bcsm = sic_bcsm(2:end,1);
            data.SIC_dksm = sic_dksm(1);
            data.SICsh_dksm = sic_dksm(2:end,1);
            data.crit_sm = crit_sm(1);
            data.critsh_sm = crit_sm(2:end,1);
            data.critthrcell = prctile(crit_sm(2:end,1),95);
            data.sparsity = sparsity;
            data.sig2noise = sig2noise;
            data.coherence = coherence;
            data.coherence_sm = coherence_sm;
        %     data.median_occ_firings = median_stats';
        %     data.variance_occ_firings = var_stats';
        %     data.perc_occ_firings = perc_stats';
        %     data.occ_data = occ_data;
        elseif repeat == 2
            data.SIC_adsm1 = sic_adsm;
            data.SIC_bcsm1 = sic_bcsm;
            data.SIC_dksm1 = sic_dksm;
            data.crit_sm1 = crit_sm;
            data.sparsity1 = sparsity;
            data.sig2noise1 = sig2noise;
            data.coherence1 = coherence;
            data.coherence_sm1 = coherence_sm;
        elseif repeat == 3
            data.SIC_adsm2 = sic_adsm;
            data.SIC_bcsm2 = sic_bcsm;
            data.SIC_dksm2 = sic_dksm;
            data.crit_sm2 = crit_sm;
            data.sparsity2 = sparsity;
            data.sig2noise = sig2noise;
            data.coherence2 = coherence;
            data.coherence_sm2 = coherence_sm;
        end

        %     data.median_occ_firings = median_stats';
        %     data.variance_occ_firings = var_stats';
        %     data.perc_occ_firings = perc_stats';

        %     data.occ_data = occ_data;

            
    end

    % Calculate intra-session stability
    map1 = data.maps_bcsm1;
    map2 = data.maps_bcsm2;
    vis1 = ~isnan(map1);
    vis2 = ~isnan(map2);
    vis = vis1 & vis2; % Correlate only visited bins;
    intracorr = corr2(map1(vis), map2(vis));
    map1z = zscore(map1(vis));
    map2z = zscore(map2(vis));
    intracorrz = corr2(map1z, map2z);
    % Store stability data
    data.intracorr = intracorr;
    data.intracorrz = intracorrz;
    
    % create nptdata so we can inherit from it
    data.gridSteps = Args.GridSteps;
    Args.NumShuffles = NumShuffles_saved;
    data.numSets = 1;
    data.Args = Args;
    n = nptdata(1,0,pwd);
    d.data = data;
    obj = class(d,Args.classname,n);
    saveObject(obj,'ArgsC',Args);
    % NEW: add the new matname to the filemap
    % vmpcFileMap(Args.matname) = true;

else
    % create empty object
    obj = createEmptyObject(Args);
end



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

function hash = javaHash(inputStr, algorithm)
    % Convert string to bytes (uint8)
    data = uint8(inputStr);
    % Create MessageDigest object for the specified algorithm (e.g., SHA-256)
    md = java.security.MessageDigest.getInstance(algorithm);
    % Update the MessageDigest with the data
    md.update(data);
    % Get the resulting hash (as a byte array)
    hashBytes = md.digest();
    % Convert hash to a hex string
    hash = '';
    for i = 1:length(hashBytes)
        hash = [hash, sprintf('%02x', hashBytes(i))];  % Format each byte as 2-digit hex
    end



function filename = hashFileName(Args)
    hash_input = '';
    for i = 1:length(Args.DataCheckArgs)
        key = Args.DataCheckArgs{i};
        val = Args.(key);
        if isnumeric(val)
            valStr = mat2str(val);
        elseif ischar(val)
            valStr = val;
        else
            error(['Unsupported argument type: ' key]);
        end
        hash_input = [hash_input '_'  valStr];
   
    end
    hashString = javaHash(hash_input, 'SHA-256');

    % Create final short and unique filename
    filename = sprintf('%s_%s.mat', Args.classname, hashString(1:16));  % Shorten to 16 chars

    


function test_hashFileName
    % --- Setup ---
    % Define a mock Args struct with required fields
    Args = struct();
    Args.GridSteps = 40;
    Args.NumShuffles = 10000;
    Args.UseMinObs = 0;
    Args.SmoothType = 'Adaptive';
    Args.ThresVel = 1;
    Args.UseAllTrials = 1;
    Args.Alpha = 10000;
    
    % Specify which fields to include in the hash
    Args.DataCheckArgs = {'GridSteps','NumShuffles','UseMinObs','SmoothType','ThresVel','UseAllTrials','Alpha'};
    
    % --- Act ---
    hash1 = hashFileName(Args);

    % --- Test for consistency: same input → same hash ---
    hash2 = hashFileName(Args);
    assert(strcmp(hash1, hash2), 'Hash should be consistent for same Args');

    % --- Test for change: modifying one field should change the hash ---
    Args.NumShuffles = 50;
    hash3 = hashFileName(Args);
    assert(~strcmp(hash1, hash3), 'Hash should change when Arg.NumShuffles changes');

    Args.NumShuffles = 100;
    hash4 = hashFileName(Args);
    assert(~strcmp(hash4, hash3), 'Hash should change when Args.NumShuffles changes');
    
    Args.UseMinObs = 40;
    hash5 = hashFileName(Args);
    assert(~strcmp(hash4, hash5), 'Hash should change when Args.UseMinObs changes');

    Args.GridSteps = 41;
    hash6 = hashFileName(Args);
    assert(~strcmp(hash5, hash6), 'Hash should change when Arg.GridSteps changes');

    Args.GridSteps = 45;
    hash7 = hashFileName(Args);
    assert(~strcmp(hash6, hash7), 'Hash should change when Arg.GridSteps changes');

    disp('All hashFileName tests passed.');