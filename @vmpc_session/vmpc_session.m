function [obj, varargout] = vmpc_session(varargin)

%@vmpc_session Constructor function for session-level vmpc class
%   OBJ = vmpc_session(varargin)
%
%   This is a SESSION-LEVEL constructor that creates vmpc objects without
%   spike shuffling. Shuffling should be performed at the day-level after
%   combining multiple session objects.
%
%   Key differences from day-level vmpc:
%   - No spike shuffling (NumShuffles forced to 0)
%   - No shuffle-based statistics (SICsh, critsh, etc.)
%   - Stores raw spike counts and durations for later aggregation
%   - Still computes session-level metrics (SIC, sparsity, coherence)
%
%example [as, Args] = vmpc_session('save','redo')

% ---------------ADDITIONAL DESCRIPTION ----------
% VMPC_DAY and VMPC_session are functions built to construct vmpc objects
% for individual cells on days with multiple sessions
% VMPC_Session create lightweight .mat objects contain session level raw
% maps
% VMPC_DAY aggregates all the individual vmpc_session for a specific cell
% selected (based on array__/channel__/cell__ path used in argument

% if there is only ONE SESSION for a day, you can run VMPC('auto') at the
% cell-level directory

% if there are multiple sessions:

%   1) Run vmpc_session on all individual session/array/cell folders to generate lightweight .mat files 
%   containing raw maps (this is where the spiketrains.mat exists)

%   2) Run vmpc_day at the day-level targeting a specific cell 
%       using the arugments in the format mentioned above 
%       It will automatically load the previously created session objects, 
%       aggregate them, and run the shuffling statistics to determine spatial info.
%


Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
                'ObjectLevel','Session', 'RequiredFile','spiketrain.mat', ...
                'GridSteps',40, 'pix',1,...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',0, ... % FORCED TO 0
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4,'SmoothType','Adaptive', 'UseMinObs',0, ...
                'ThresVel',1, 'UseAllTrials',1,...
                'SelectiveCriteria','SIC','Alpha', 10000,'UseFileHash',0);
            
Args.flags = {'Auto','ArgsOnly','FRSIC','UseMedian'};
Args.DataCheckArgs = {'GridSteps','UseMinObs','SmoothType','ThresVel','UseAllTrials', 'Alpha'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
    'subtract',{'RedoLevels','SaveLevels'}, ...
    'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
    'remove',{'Auto'});

Args.classname = 'vmpc_session';
if Args.UseFileHash
    filename = hashFileName(Args);
else
    filename = [Args.classname '.mat'];
end
Args.matname = filename;
Args.matvarname = 'vmp_sess';

[command,robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
    if exist(Args.matname, 'file')
        obj = robj; 
    elseif exist([Args.classname '.mat'],'file')
        obj = robj;
    else
        disp(['Object not existed for ' Args.matname ', create object instead']);
        obj = createObject(Args,modvarargin{:});
    end
elseif(strcmp(command,'createObj'))
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

dlist = nptDir;
dnum = size(dlist,1);

% Validate required arguments
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

% Check if required file exists
if(~isempty(dir(Args.RequiredFile)))
    
    ori = pwd;
    data.origin = {pwd};

    % Load SESSION-LEVEL vmpv object
    cd ..; cd ..; cd ..;
    % pv = load([num2str(Args.pix) 'vmpv.mat']);
    pv = load('vmpv.mat');
    pv = pv.pv;
    cd(ori);
    
    % Load spike train
    [d,fn,ext] = fileparts(Args.RequiredFile);
    if strcmp(ext,'.mat')
        spiketrain = load(Args.RequiredFile);
    elseif strcmp(ext,'.csv')
        spiketrain.timestamps = load(Args.RequiredFile)';
    else
        error(['Unknown file type ' Args.RequiredFile]);
    end
    
    for repeat = 1:3
        
        if repeat == 1
            disp('Session-level VMPC: Full session');
        elseif repeat == 2
            disp('Session-level VMPC: First half');
        elseif repeat == 3
            disp('Session-level VMPC: Second half');
        end
        
        if repeat == 1
            stc = pv.data.sessionTimeC;
        end
        
        %% NO SPIKE SHUFFLING - Just process actual spikes
        spiketimes = spiketrain.timestamps/1000; % convert to seconds
        maxTime = pv.data.rplmaxtime;
        sessionStart = stc(1,1);
        sessionEnd = sessionStart + maxTime; 
        spiketimes = spiketimes(spiketimes >= sessionStart & spiketimes <= sessionEnd);
        % Create spike time array (only 1 row, no shuffles)
        flat_spiketimes = [spiketimes(:), ones(length(spiketimes),1)];
        flat_spiketimes = sortrows(flat_spiketimes);
        flat_spiketimes(flat_spiketimes(:,1) < stc(1,1),:) = [];
        
        %% Filter sessionTimeC based on conditions
        if repeat == 1
            disp('      Filtering...');
            stc(:,5) = [diff(stc(:,1)); 0];  % Duration
            stc(:,6) = zeros(size(stc,1),1); % Spike count
        end
        
        conditions = ones(size(stc,1),1);

        % Apply filters
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
            conditions = conditions & (pv.data.pv_good_rows);
        else
            bins_sieved = 1:(Args.GridSteps * Args.GridSteps);
        end

        % Group consecutive rows with same timestamp
        if repeat == 1
            dstc = diff(stc(:,1));
            stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
            stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
            stc_changing_ind = stc_changing_ind(1:end-1,:);
        end
        
        % Initialize spike count array (only 1 column, no shuffles)
        consol_arr = zeros(Args.GridSteps * Args.GridSteps, 1);
        
        %% Assign spikes to bins
        if repeat == 1
            disp(['      Assigning ' num2str(size(flat_spiketimes,1)) ' spikes to bins...']);
        end
        
        interval = 1;
        for sp = 1:size(flat_spiketimes,1)
            % Find time interval for this spike
            while interval < size(stc_changing_ind,1)
                if flat_spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && ...
                        flat_spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
                    break;
                end
                interval = interval + 1;
            end

            % Record spike in sessionTimeC (last row of interval)
            stc(stc_changing_ind(interval,2),6) = stc(stc_changing_ind(interval,2),6) + 1;
            
            % Get bins that meet filter criteria
            bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),[2 3 4]);
            bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))),:);
            bins_hit(~(bins_hit(:,1)>0),:) = []; % Remove invalid place bins
            bins_hit(~(bins_hit(:,3)>0),:) = []; % Remove invalid view bins
            bins_hit(~(bins_hit(:,2)>0),:) = []; % Remove invalid hd bins
            
            % Accumulate spikes
            consol_arr(bins_hit(:,1)) = consol_arr(bins_hit(:,1)) + 1;
        end
        
        spike_count_full = consol_arr';

        %% Calculate duration per bin
        stc_filt = stc(find(conditions==1),:);
        stc_filt(~(stc_filt(:,2) > 0),:) = [];
        stc_filt(isnan(stc_filt(:,4)),:) = [];
        stc_filt(~(stc_filt(:,3) > 0),:) = [];
        stc_ss = stc_filt(:,[2 5]);
        stc_ss = [stc_ss; [1600 0]];
        gpdurfull = accumarray(stc_ss(:,1),stc_ss(:,2))';

        % Apply bin sieving
        spikes_count = zeros(1, Args.GridSteps*Args.GridSteps);
        dur_raw = zeros(1, Args.GridSteps*Args.GridSteps);
        spikes_count(bins_sieved) = spike_count_full(bins_sieved);
        dur_raw(bins_sieved) = gpdurfull(bins_sieved);
        
        % Calculate raw firing rate map
        map_raw = spikes_count ./ dur_raw;
        spk_raw = spikes_count;
        
        % Store raw data
        if repeat == 1
            data.sessionTimeC = stc;
            data.stcfilt = stc_filt;
            data.maps_raw = map_raw;
            data.dur_raw = dur_raw;
            data.spk_raw = spk_raw;
            data.filtspknum = sum(spk_raw);
            data.filtered_spiketimes = spiketimes; 
            data.session_duration = maxTime;
            data.filter_conditions = conditions;
        elseif repeat == 2
            data.maps_raw1 = map_raw;
            data.dur_raw1 = dur_raw;
            data.spk_raw1 = spk_raw;
        elseif repeat == 3
            data.maps_raw2 = map_raw;
            data.dur_raw2 = dur_raw;
            data.spk_raw2 = spk_raw;
        end
        
        %% Smoothing
        if repeat == 1
            disp('      Adaptive smoothing...');
        end

        [maps_adsm, durs_adsm, rad_adsm, maps_bcsm, maps_dksm, ~, ~, rad_adsm_grid] = ...
            smoothMaps(map_raw, dur_raw, spk_raw, spikes_count, Args);

        % Select smoothing type
        switch Args.SmoothType
            case 'Adaptive'
                maps_sm = maps_adsm;
            case 'Boxcar'
                maps_sm = maps_bcsm;
            case 'Disk'
                maps_sm = maps_dksm;
        end
        
        % Quality checks
        if repeat == 1
            if data.filtspknum < 100
                data.discard = true;
            else
                data.discard = false;
            end
            if max(maps_sm,[],'omitnan') < 0.7
                data.rateok = false;
            else
                data.rateok = true;
            end
        end
        
        % Store smoothed maps
        if repeat == 1
            data.maps_adsm = maps_adsm;
            data.dur_adsm = durs_adsm;
            data.radii = rad_adsm_grid;
            data.maps_bcsm = maps_bcsm;
            data.maps_dksm = maps_dksm;
            data.maps_sm = maps_sm;
        elseif repeat == 2
            data.maps_adsm1 = maps_adsm;
            data.dur_adsm1 = durs_adsm;
            data.radii1 = rad_adsm_grid;
            data.maps_bcsm1 = maps_bcsm;
            data.maps_dksm1 = maps_dksm;
            data.maps_sm1 = maps_sm;
        elseif repeat == 3
            data.maps_adsm2 = maps_adsm;
            data.dur_adsm2 = durs_adsm;
            data.radii2 = rad_adsm_grid;
            data.maps_bcsm2 = maps_bcsm;
            data.maps_dksm2 = maps_dksm;
            data.maps_sm2 = maps_sm;
        end

        %% Calculate spatial metrics (NO SHUFFLE-BASED STATISTICS)
        disp('      Calculating SIC...');
        sic_adsm = skaggs_sic(maps_adsm', durs_adsm');
        sic_bcsm = skaggs_sic(maps_bcsm', durs_adsm');
        sic_dksm = skaggs_sic(maps_dksm', durs_adsm');

        switch Args.SmoothType
            case 'Adaptive'
                sic_sm = sic_adsm;
            case 'Boxcar'
                sic_sm = sic_bcsm;
            case 'Disk'
                sic_sm = sic_dksm;
        end

        % Calculate other spatial metrics
        sparsity = spatial_sparsity(dur_raw, map_raw);
        sig2noise = spatial_sig2noise(map_raw);
        coherence = spatial_coherence('place', [Args.GridSteps Args.GridSteps], map_raw, 1);
        coherence_sm = spatial_coherence('place', [Args.GridSteps Args.GridSteps], maps_sm, 1);

        % Store metrics
        if repeat == 1
            data.SIC_adsm = sic_adsm;
            data.SIC_bcsm = sic_bcsm;
            data.SIC_dksm = sic_dksm;
            data.crit_sm = sic_sm;  % No shuffle threshold available
            data.sparsity = sparsity;
            data.sig2noise = sig2noise;
            data.coherence = coherence;
            data.coherence_sm = coherence_sm;
        elseif repeat == 2
            data.SIC_adsm1 = sic_adsm;
            data.SIC_bcsm1 = sic_bcsm;
            data.SIC_dksm1 = sic_dksm;
            data.crit_sm1 = sic_sm;
            data.sparsity1 = sparsity;
            data.sig2noise1 = sig2noise;
            data.coherence1 = coherence;
            data.coherence_sm1 = coherence_sm;
        elseif repeat == 3
            data.SIC_adsm2 = sic_adsm;
            data.SIC_bcsm2 = sic_bcsm;
            data.SIC_dksm2 = sic_dksm;
            data.crit_sm2 = sic_sm;
            data.sparsity2 = sparsity;
            data.sig2noise2 = sig2noise;
            data.coherence2 = coherence;
            data.coherence_sm2 = coherence_sm;
        end
    end

    %% Calculate intra-session stability
    map1 = data.maps_bcsm1;
    map2 = data.maps_bcsm2;
    vis1 = ~isnan(map1);
    vis2 = ~isnan(map2);
    vis = vis1 & vis2;
    intracorr = corr2(map1(vis), map2(vis));
    map1z = zscore(map1(vis));
    map2z = zscore(map2(vis));
    intracorrz = corr2(map1z, map2z);
    
    data.intracorr = intracorr;
    data.intracorrz = intracorrz;
    
    % Store metadata
    data.gridSteps = Args.GridSteps;
    data.numSets = 1;
    data.Args = Args;
    data.sessionLevel = true;  % Flag indicating this is session-level
    
    % Create object
    n = nptdata(1,0,pwd);
    d.data = data;
    obj = class(d, Args.classname, n);
    saveObject(obj,'ArgsC',Args);

else
    obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)
data.dlist = [];
data.setIndex = [];
data.numSets = 0;
data.Args = Args;
data.sessionLevel = true;
n = nptdata(0,0);
d.data = data;
obj = class(d, Args.classname, n);

function hash = javaHash(inputStr, algorithm)
    data = uint8(inputStr);
    md = java.security.MessageDigest.getInstance(algorithm);
    md.update(data);
    hashBytes = md.digest();
    hash = '';
    for i = 1:length(hashBytes)
        hash = [hash, sprintf('%02x', hashBytes(i))];
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
        hash_input = [hash_input '_' valStr];
    end
    hashString = javaHash(hash_input, 'SHA-256');
    filename = sprintf('%s_%s.mat', Args.classname, hashString(1:16));
