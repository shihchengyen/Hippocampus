function [obj, varargout] = vmpc_day(varargin)
%@vmpc_day Aggregates session-level vmpc objects into day-level with shuffling
%   OBJ = vmpc_day(varargin)
%
%   This function:
%   1. Loads multiple session-level vmpc_session objects
%   2. Aggregates spike counts and durations across sessions
%   3. Performs spike shuffling at the DAY level
%   4. Computes day-level spatial metrics with shuffle-based statistics
%
%   Usage:
%   obj = vmpc_day('auto', 'CellPath', 'array01/channel01/cell01')
%   obj = vmpc_day('Sessions', {sess1, sess2, sess3})
%   obj = vmpc_day('SessionDirs', {'session01/array01/channel01/cell01', ...})
%
%   Arguments:
%   'CellPath' - Relative path from session dir to cell (e.g., 'array01/channel01/cell01')
%                Used with 'auto' to find the same cell across all sessions
%
%example [as, Args] = vmpc_day('save','redo','CellPath','array01/channel01/cell01')

Args = struct('RedoLevels',0, 'SaveLevels',1, 'Auto',0, 'ArgsOnly',0, ...
                'ObjectLevel','Day', ...
                'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'SmoothType','Adaptive', 'Alpha', 10000, ...
                'SelectiveCriteria','SIC', ...
                'GaussSigma', 2.0, ... %for gaussian smoothing
                'Sessions',[], ...  % Pass session objects directly
                'SessionDirs',[], ...  % Or specify full paths to cells
                'CellPath','', ...  % Relative path from session to cell
                'UseFileHash',0);
            
Args.flags = {'Auto','ArgsOnly'};
Args.DataCheckArgs = {'GridSteps','NumShuffles','SmoothType','Alpha'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
    'subtract',{'RedoLevels','SaveLevels'}, ...
    'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
    'remove',{'Auto'});

Args.classname = 'vmpc_day';
if Args.UseFileHash
    filename = hashFileName(Args);
else
    filename = [Args.classname '.mat'];
end
Args.matname = filename;
Args.matvarname = 'vmp_day';



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
    else
        disp(['Day-level object not found, creating new one']);
        obj = createObject(Args,modvarargin{:});
    end
elseif(strcmp(command,'createObj'))
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

data.origin = {pwd};

%% STEP 1: Load or receive session objects
disp('==================================================');
disp('DAY-LEVEL VMPC AGGREGATION');
disp('==================================================');

sessionObjects = [];

if ~isempty(Args.Sessions)
    % Sessions passed directly
    sessionObjects = Args.Sessions;
    disp(['Received ' num2str(length(sessionObjects)) ' session objects']);
    
elseif ~isempty(Args.SessionDirs)
    % Load from specified directories
    disp(['Loading sessions from ' num2str(length(Args.SessionDirs)) ' directories']);
    sessionObjects = cell(length(Args.SessionDirs), 1);
    for i = 1:length(Args.SessionDirs)
        cd(Args.SessionDirs{i});
        if exist('vmpc_session.mat', 'file')
            temp = load('vmpc_session.mat');
            sessionObjects{i} = temp.vmp_sess;
            disp(['  Loaded: ' Args.SessionDirs{i}]);
            
        else
            error(['vmpc_session.mat not found in ' Args.SessionDirs{i}]);
        end
        cd(data.origin{1});
    end
    
elseif ~isempty(Args.CellPath)
    % Auto-detect using CellPath: search session*/CellPath/vmpc_session.mat
    disp(['Auto-detecting sessions for cell: ' Args.CellPath]);
    dlist = dir;
    sessionObjects = {};
    sessionNames = {};
    
    for i = 1:length(dlist)
        % Look for directories starting with 'session'
        if dlist(i).isdir && strncmp(dlist(i).name, 'session', 7)
            cellFullPath = fullfile(pwd, dlist(i).name, Args.CellPath);
            vmpcFile = fullfile(cellFullPath, 'vmpc_session.mat');
            
            if exist(vmpcFile, 'file')
                cd(cellFullPath);
                temp = load('vmpc_session.mat');
                sessionObjects{end+1} = temp.vmp_sess;
                sessionNames{end+1} = dlist(i).name;
                disp(['  Found: ' dlist(i).name '/' Args.CellPath]);
                cd(data.origin{1});
           
            end
        end
    end
    
    if isempty(sessionObjects)
        error(['No vmpc_session.mat files found for cell path: ' Args.CellPath]);
    end
    tempCellPath = Args.CellPath;
    if iscell(tempCellPath)
        tempCellPath = tempCellPath{1}; % Take the first string from the cell
    end
    
    cellPathNaming = replace(tempCellPath, filesep, "-"); % Use filesep instead of "/"
    Args.matname = [Args.classname '_' cellPathNaming '.mat'];
    disp(['Day-level object will be located/saved as: ' Args.matname]);
    %if Args.RedoLevels == 0 && exist(Args.matname, 'file')
    %    disp(['File exists, loading and skipping create: ' Args.matname]);
        % We must load the object so the batch script gets the data
    %    l = load(Args.matname);
    %    obj = eval(['l.' Args.matvarname]);
    %    return; % End the function immediately
    % end
    % Store session names for reference
    data.sessionNames = sessionNames;
    
else
    % No specification provided
    error(['Must specify one of: ''CellPath'', ''Sessions'', or ''SessionDirs''. ' ...
           'Example: vmpc_day(''auto'', ''CellPath'', ''array01/channel01/cell01'')']);
end

numSessions = length(sessionObjects);
disp(['Total sessions to aggregate: ' num2str(numSessions)]);
disp(' ');

%% STEP 2: Validate session objects
disp('Validating sessions...');
for i = 1:numSessions
    sess = sessionObjects{i};
    
    % Check it's a session-level object
    if ~isfield(sess.data, 'sessionLevel') || ~sess.data.sessionLevel
        warning(['Session ' num2str(i) ' does not appear to be session-level']);
    end
    
    % Check grid compatibility
    if sess.data.gridSteps ~= Args.GridSteps
        error(['Session ' num2str(i) ' has incompatible GridSteps: ' ...
               num2str(sess.data.gridSteps) ' vs ' num2str(Args.GridSteps)]);
    end
    
    % Check for required fields
    requiredFields = {'spk_raw', 'dur_raw', 'maps_raw', ...
                      'spk_raw1', 'dur_raw1', 'maps_raw1', ...
                      'spk_raw2', 'dur_raw2', 'maps_raw2'};
    for j = 1:length(requiredFields)
        if ~isfield(sess.data, requiredFields{j})
            error(['Session ' num2str(i) ' missing field: ' requiredFields{j}]);
        end
    end
    %disp(['  Session ' num2str(i) ': ' num2str(sess.data.filtered_spiketimes) ' spiketimes, OK']);
    disp(['  Session ' num2str(i) ': ' num2str(sess.data.filtspknum) ' spikes, OK']);
end

%% STEP 3: Aggregate spike counts and durations across sessions
disp(' ');
disp('Aggregating data across sessions...');

% Initialize aggregated arrays
session_intracorrs = NaN(numSessions, 1);
spk_agg = zeros(1, Args.GridSteps * Args.GridSteps);
dur_agg = zeros(1, Args.GridSteps * Args.GridSteps);
spk_agg1 = zeros(1, Args.GridSteps * Args.GridSteps);
dur_agg1 = zeros(1, Args.GridSteps * Args.GridSteps);
spk_agg2 = zeros(1, Args.GridSteps * Args.GridSteps);
dur_agg2 = zeros(1, Args.GridSteps * Args.GridSteps);

% Sum across sessions
for i = 1:numSessions
    sess = sessionObjects{i};
    if isfield(sess.data, 'intracorr')
        session_intracorrs(i) = sess.data.intracorr;
    end
    
    % Full session aggregation
    spk_agg = spk_agg + sess.data.spk_raw;
    dur_agg = dur_agg + sess.data.dur_raw;
    
    
    
end

totalSpikes = sum(spk_agg);
disp(['Total spikes aggregated: ' num2str(totalSpikes)]);

% Calculate aggregated firing rate maps
map_agg = spk_agg ./ dur_agg;

% Store aggregated raw data
data.spk_raw = spk_agg;
data.dur_raw = dur_agg;
data.maps_raw = map_agg;
% (We no longer store the flawed aggregated halves)
data.filtspknum = totalSpikes;
data.numSessions = numSessions;

%% STEP 4: Perform DAY-LEVEL spike shuffling
disp(' ');
disp(['Performing ' num2str(Args.NumShuffles) ' spike shuffles at day level...']);

% Load spike trains and session info from all sessions
allSpikeTimes = [];
sessionBoundaries = [0];
cellPaths = {};
allSessionTimeC = [];  % NEW: Store all sessionTimeC data
sessionTimeOffsets = [0];  % NEW: Track time offsets for each session
all_conditions = [];

for i = 1:numSessions
    % Get the session object (already in memory)
    sess = sessionObjects{i};
    
    % 1. Get pre-processed data from the session object
    spikeTimes = sess.data.filtered_spiketimes;
    sessionDur = sess.data.session_duration;
    stc = sess.data.sessionTimeC; % Full sessionTimeC
    conditions = sess.data.filter_conditions; % Pre-made velocity/trial filter
    
    cellPaths{i} = sess.data.origin{1};
    sessionStart = stc(1,1);

    % 2. Offset spike times
    if i > 1
        spikeTimes = spikeTimes - sessionStart + sessionBoundaries(i);
    else
        spikeTimes = spikeTimes - sessionStart;  % Start from 0
    end
    
    % 3. Offset sessionTimeC timestamps
    stc_offset = stc;
    stc_offset(:,1) = stc(:,1) - sessionStart + sessionBoundaries(i);
    
    % 4. Append to aggregated arrays
    allSpikeTimes = [allSpikeTimes; spikeTimes(:)];
    allSessionTimeC = [allSessionTimeC; stc_offset];
    all_conditions = [all_conditions; conditions(:)]; % Append the filter
    
    sessionBoundaries(i+1) = sessionBoundaries(i) + sessionDur;
    sessionTimeOffsets(i+1) = sessionBoundaries(i+1);
end
% --- END OF THE NEW LOOP ---

totalDuration = sessionBoundaries(end);
disp(['Total day duration: ' num2str(totalDuration) ' seconds']);
disp(['Total spikes for shuffling: ' num2str(length(allSpikeTimes))]);

% Store session info
data.sessionBoundaries = sessionBoundaries;
data.cellPaths = cellPaths;

disp(['VMPC_DAY: Number of spikes before shuffle: ' num2str(length(allSpikeTimes))]);
disp(['VMPC_DAY: totalDuration: ' num2str(totalDuration)]);
disp(['VMPC_DAY: allSessionTimeC(1,1): ' num2str(allSessionTimeC(1,1))]);

% Generate shuffles using the same approach as vmpc
disp('Generating shuffle distributions...');
tShifts = [0 ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+...
           Args.ShuffleLimits(1))*totalDuration];
full_arr = repmat(allSpikeTimes', Args.NumShuffles+1, 1);
full_arr = full_arr + tShifts';

% Circular wrap
keepers = size(full_arr,2) - sum(full_arr>totalDuration, 2);
for row = 2:size(full_arr,1)
    full_arr(row,:) = [full_arr(row,1+keepers(row):end)-totalDuration ...
                       full_arr(row,1:keepers(row))];
end

% Create flat_spiketimes structure like in vmpc
flat_spiketimes = NaN(2, size(full_arr,1)*size(full_arr,2));
temp = full_arr';
flat_spiketimes(1,:) = temp(:);  % All spike times
flat_spiketimes(2,:) = repelem(1:size(full_arr,1), size(full_arr,2));  % Shuffle ID
flat_spiketimes = flat_spiketimes';
flat_spiketimes = sortrows(flat_spiketimes);  % Sort by timestamp

% Remove spikes before first sessionTimeC timestamp
flat_spiketimes(flat_spiketimes(:,1) < allSessionTimeC(1,1),:) = [];

disp('Assigning shuffled spikes to spatial bins...');

% Prepare sessionTimeC data (similar to vmpc)
stc = allSessionTimeC;
stc(:,5) = [diff(stc(:,1)); 0];  % Duration column
stc(:,6) = zeros(size(stc,1),1);  % Spike binning column

% Set up filtering conditions (adapt based on your needs)
conditions = all_conditions;

% Filter based on place, view, and hd data availability
% You may need to adjust this based on how session data is aggregated
bins_sieved = 1:(Args.GridSteps * Args.GridSteps);

% Group into intervals (consecutive rows with same timestamp)
dstc = diff(stc(:,1));
stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
stc_changing_ind = stc_changing_ind(1:end-1,:);

% Initialize spike count array
consol_arr = zeros(Args.GridSteps * Args.GridSteps, Args.NumShuffles + 1);

disp(['Assigning ' num2str(size(flat_spiketimes,1)) ' spikes to bins...']);
interval = 1;

for sp = 1:size(flat_spiketimes,1)
    
    % Find the interval this spike belongs to
    while interval < size(stc_changing_ind,1)
        if flat_spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && ...
                flat_spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
            break;
        end
        interval = interval + 1;
    end
    
    % Bin spike into stc (for original data only)
    if flat_spiketimes(sp,2) == 1
        stc(stc_changing_ind(interval,2),6) = stc(stc_changing_ind(interval,2),6) + 1;
    end
    
    % Get bins hit during this interval
    bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),[2 3 4]);
    bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))),:);
    bins_hit(~(bins_hit(:,1)>0),:) = [];  % Remove place bin = 0
    bins_hit(~(bins_hit(:,3)>0),:) = [];  % Remove view bin = nan
    bins_hit(~(bins_hit(:,2)>0),:) = [];  % Remove hd bin = 0
    
    % Assign spike to bins for this shuffle
    consol_arr(bins_hit(:,1), flat_spiketimes(sp,2)) = ...
        consol_arr(bins_hit(:,1), flat_spiketimes(sp,2)) + 1;
    
    % Progress indicator
    if mod(sp, 100000) == 0
        disp(['  Processed ' num2str(sp) '/' num2str(size(flat_spiketimes,1)) ' spikes']);
    end
end

% Extract spike counts
spike_count_shuffled = consol_arr';

disp('Computing firing rate maps for shuffles...');
maps_shuffled = spike_count_shuffled ./ repmat(dur_agg, Args.NumShuffles+1, 1);
maps_shuffled(1,:) = map_agg;
spike_count_shuffled(1,:) = spk_agg;
%% STEP 5: Smooth aggregated maps

disp('Smoothing aggregated firing rate maps...');

[maps_adsm, durs_adsm, ~, maps_bcsm, maps_dksm, durs_bcsm, durs_dksm, rad_adsm_grid] = ...
    smoothMaps(maps_shuffled, dur_agg, spk_agg, spike_count_shuffled, Args);

%%section on gaussian smoothing
%% ---------------------------- 
disp('Applying Gaussian smoothing...');
% 1. Create the 2D Gaussian kernel
if license('checkout', 'image_toolbox')
    sigma_bins = Args.GaussSigma;
    kernel_size = ceil(sigma_bins * 3) * 2 + 1; % Kernel size ~3x sigma
    gauss_kernel = fspecial('gaussian', kernel_size, sigma_bins);
else
    error('Image Processing Toolbox license not found. Required for fspecial(''gaussian'').');
end

% 2. Reshape raw maps to 2D
grid_dim = Args.GridSteps;
num_all_maps = Args.NumShuffles + 1;
spk_maps_2d = reshape(spike_count_shuffled', grid_dim, grid_dim, num_all_maps);
dur_map_2d = reshape(dur_agg, grid_dim, grid_dim); 

% 3. Initialize output arrays
maps_gasm = zeros(num_all_maps, grid_dim * grid_dim);
durs_gasm = zeros(num_all_maps, grid_dim * grid_dim); 

% 4. Smooth the 2D duration map (only needs to be done once)
smoothed_occupancy = imfilter(dur_map_2d, gauss_kernel, 'conv', 'replicate');
durs_gasm(1:end, :) = repmat(reshape(smoothed_occupancy, 1, []), num_all_maps, 1);

% 5. Loop and smooth each spike map (real + shuffles)
for i = 1:num_all_maps
    smoothed_spikes = imfilter(spk_maps_2d(:,:,i), gauss_kernel, 'conv', 'replicate');
    gasm_map_2d = smoothed_spikes ./ smoothed_occupancy;
    gasm_map_2d(smoothed_occupancy < 1e-5) = NaN; 
    maps_gasm(i, :) = reshape(gasm_map_2d, 1, []);
end
%% ---------------------------- 

switch Args.SmoothType
    case 'Adaptive'
        maps_sm = maps_adsm;
    case 'Boxcar'
        maps_sm = maps_bcsm;
    case 'Disk'
        maps_sm = maps_dksm;
    case 'Gaussian' 
        maps_sm = maps_gasm; 
end

% Store smoothed data

data.maps_adsm = maps_adsm(1,:);
data.maps_adsmsh = maps_adsm(2:end,:);
data.dur_adsm = durs_adsm(1,:);
data.radii = rad_adsm_grid(1,:);
data.radiish = rad_adsm_grid(2:end,:);
data.maps_bcsm = maps_bcsm(1,:);
data.maps_bcsmsh = maps_bcsm(2:end,:);
data.maps_dksm = maps_dksm(1,:);
data.maps_dksmsh = maps_dksm(2:end,:);
data.maps_sm = maps_sm(1,:);
data.maps_smsh = maps_sm(2:end,:);

% extra: Store Gaussian maps
data.maps_gasm = maps_gasm(1,:);
data.maps_gasmsh = maps_gasm(2:end,:);
data.durs_gasm = durs_gasm(1,:);

%% STEP 6: Calculate day-level spatial metrics with shuffle statistics
disp('Calculating spatial information content...');

sic_adsm = skaggs_sic(maps_adsm', durs_adsm');
sic_bcsm = skaggs_sic(maps_bcsm', durs_bcsm');
sic_dksm = skaggs_sic(maps_dksm', durs_dksm');
sic_gasm = skaggs_sic(maps_gasm', durs_gasm');

sic_adsm = sic_adsm';

sic_bcsm = sic_bcsm';

sic_dksm = sic_dksm';

sic_gasm = sic_gasm';

switch Args.SmoothType
    case 'Adaptive'
        sic_sm = sic_adsm;
    case 'Boxcar'
        sic_sm = sic_bcsm;
    case 'Disk'
        sic_sm = sic_dksm;

    case 'Gaussian'
        sic_sm = sic_gasm; 
end

switch Args.SelectiveCriteria
    case 'SIC'
        crit_sm = sic_sm;
end

% Calculate other metrics
sparsity = spatial_sparsity(dur_agg, map_agg);
sig2noise = spatial_sig2noise(map_agg);
coherence = spatial_coherence('place', [Args.GridSteps Args.GridSteps], map_agg, 1);
coherence_sm = spatial_coherence('place', [Args.GridSteps Args.GridSteps], maps_sm(1,:), 1);

% Store metrics with shuffle statistics
%%new gaussian
data.SIC_gasm = sic_gasm(1);
data.SICsh_gasm = sic_gasm(2:end);

data.SIC_adsm = sic_adsm(1);
data.SICsh_adsm = sic_adsm(2:end);
data.SIC_bcsm = sic_bcsm(1);
data.SICsh_bcsm = sic_bcsm(2:end);
data.SIC_dksm = sic_dksm(1);
data.SICsh_dksm = sic_dksm(2:end);
data.crit_sm = crit_sm(1);
data.critsh_sm = crit_sm(2:end);
data.critthrcell = prctile(crit_sm(2:end), 95);  % 95th percentile threshold
data.sparsity = sparsity;
data.sig2noise = sig2noise;
data.coherence = coherence;
data.coherence_sm = coherence_sm;

% Quality checks
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

% Statistical significance
if data.crit_sm > data.critthrcell
    data.significant = true;
    disp(' ');
    disp(['*** Cell is SPATIALLY SELECTIVE ***']);
    disp(['    SIC = ' num2str(data.crit_sm, '%.4f')]);
    disp(['    Threshold (95th percentile) = ' num2str(data.critthrcell, '%.4f')]);
else
    data.significant = false;
    disp(' ');
    disp(['Cell is NOT spatially selective']);
    disp(['    SIC = ' num2str(data.crit_sm, '%.4f')]);
    disp(['    Threshold (95th percentile) = ' num2str(data.critthrcell, '%.4f')]);
end

%% STEP 7: [NEW] Calculate True Intra-Day Halves
disp(' ');
disp('Processing true intra-day halves for stability...');

Args_half = Args;
Args_half.NumShuffles = 0; % No shuffling for halves

if numSessions == 0 % Should be impossible, but safe
    midpointTime = 0;
else
    midpointTime = totalDuration / 2;
end

% --- First Half of the Day ---
disp('  Calculating first half of the day...');
[map_agg1, dur_agg1, spk_agg1] = bin_day_half(midpointTime, 'first', ...
    allSessionTimeC, all_conditions, flat_spiketimes, Args);

% --- Second Half of the Day ---
disp('  Calculating second half of the day...');
[map_agg2, dur_agg2, spk_agg2] = bin_day_half(midpointTime, 'second', ...
    allSessionTimeC, all_conditions, flat_spiketimes, Args);

% --- Smooth the new block-half maps ---
% First half
[maps_adsm1, durs_adsm1, ~, maps_bcsm1, maps_dksm1, ~, ~, ~] = ...
    smoothMaps(map_agg1, dur_agg1, spk_agg1, spk_agg1, Args_half);

% Second half
[maps_adsm2, durs_adsm2, ~, maps_bcsm2, maps_dksm2, ~, ~, ~] = ...
    smoothMaps(map_agg2, dur_agg2, spk_agg2, spk_agg2, Args_half);

% Store half maps (we don't store metrics for halves)
data.maps_adsm1 = maps_adsm1;
data.maps_bcsm1 = maps_bcsm1;
data.maps_dksm1 = maps_dksm1;

data.maps_adsm2 = maps_adsm2;
data.maps_bcsm2 = maps_bcsm2;
data.maps_dksm2 = maps_dksm2;

% Store raw half data (for reference)
data.maps_raw1 = map_agg1;
data.dur_raw1 = dur_agg1;
data.spk_raw1 = spk_agg1;
data.maps_raw2 = map_agg2;
data.dur_raw2 = dur_agg2;
data.spk_raw2 = spk_agg2;

%% STEP 8: Calculate day-level stability (True Intra-Day)
% This step now correlates the two *correct* halves from STEP 7.
disp('Calculating true intra-day stability...');

map1 = data.maps_bcsm1;
map2 = data.maps_bcsm2;
vis1 = ~isnan(map1);
vis2 = ~isnan(map2);
vis = vis1 & vis2;

if sum(vis) > 0
    intracorr = corr2(map1(vis), map2(vis));
    map1z = zscore(map1(vis));
    map2z = zscore(map2(vis));
    intracorrz = corr2(map1z, map2z);
else
    intracorr = NaN;
    intracorrz = NaN;
end

data.intracorr = intracorr;
data.intracorrz = intracorrz;

disp(['Day-level stability (intra-day corr): ' num2str(intracorr, '%.4f')]);
%% Store metadata
data.gridSteps = Args.GridSteps;
data.numSets = 1;
data.Args = Args;
data.sessionLevel = false;  % This is day-level
data.sessionObjects = sessionObjects;  % Keep references to sessions
data.session_intracorrs = session_intracorrs; % Nx1 vector of session values

% Create object

n = nptdata(1,0,pwd);
d.data = data;
obj = class(d, Args.classname, n);
saveObject(obj,'ArgsC',Args);


disp(' ');
disp('==================================================');
disp('DAY-LEVEL VMPC OBJECT CREATED SUCCESSFULLY');
disp('==================================================');

function obj = createEmptyObject(Args)
data.dlist = [];
data.setIndex = [];
data.numSets = 0;
data.Args = Args;
data.sessionLevel = false;
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

function [map_agg, dur_agg, spk_agg] = bin_day_half(midpointTime, half, ...
    allSessionTimeC, all_conditions, flat_spiketimes, Args)
    
    % 1. Filter all data by the time midpoint
    if strcmp(half, 'first')
        % Get all data *before* the midpoint
        time_idx = allSessionTimeC(:,1) < midpointTime;
        spike_idx = flat_spiketimes(:,1) < midpointTime;
    else
        % Get all data *at or after* the midpoint
        time_idx = allSessionTimeC(:,1) >= midpointTime;
        spike_idx = flat_spiketimes(:,1) >= midpointTime;
    end
    
    stc_half = allSessionTimeC(time_idx, :);
    cond_half = all_conditions(time_idx);
    spikes_half = flat_spiketimes(spike_idx, :);

    % 2. Calculate Duration Map (dur_agg) for this half
    stc_filt = stc_half(find(cond_half==1),:);
    stc_filt(~(stc_filt(:,2) > 0),:) = [];
    stc_filt(isnan(stc_filt(:,4)),:) = [];
    stc_filt(~(stc_filt(:,3) > 0),:) = [];
    stc_ss = stc_filt(:,[2 5]); % [place dur]
    stc_ss = [stc_ss; [Args.GridSteps*Args.GridSteps 0]]; % Ensure full size
    dur_agg = accumarray(stc_ss(:,1),stc_ss(:,2))';
    
    % Ensure dur_agg is correct size
    if length(dur_agg) < (Args.GridSteps*Args.GridSteps)
        dur_agg(Args.GridSteps*Args.GridSteps) = 0;
    end

    % 3. Calculate Spike Map (spk_agg) for this half
    % (This is the binning logic from STEP 4, simplified for only real spikes)
    
    stc = stc_half;
    conditions = cond_half;
    
    dstc = diff(stc(:,1));
    stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
    stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
    stc_changing_ind = stc_changing_ind(1:end-1,:);

    spk_agg = zeros(1, Args.GridSteps * Args.GridSteps);
    interval = 1;

    % Filter for real spikes only
    real_spikes_half = spikes_half(spikes_half(:,2) == 1, :);

    for sp = 1:size(real_spikes_half,1)
        while interval < size(stc_changing_ind,1)
            if real_spikes_half(sp,1) >= stc(stc_changing_ind(interval,1),1) && ...
                    real_spikes_half(sp,1) < stc(stc_changing_ind(interval+1,1),1)
                break;
            end
            interval = interval + 1;
        end
        
        bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),[2 3 4]);
        bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))),:);
        bins_hit(~(bins_hit(:,1)>0),:) = [];
        bins_hit(~(bins_hit(:,3)>0),:) = [];
        bins_hit(~(bins_hit(:,2)>0),:) = [];
        
        % We add 1 spike to all valid bins in that time interval
        spk_agg(bins_hit(:,1)) = spk_agg(bins_hit(:,1)) + 1;
    end

    % 4. Calculate Rate Map
    map_agg = spk_agg ./ dur_agg;