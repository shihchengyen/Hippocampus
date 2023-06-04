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
Args.DataCheckArgs = {'GridSteps','NumShuffles','UseMinObs','AdaptiveSmooth','ThresVel','UseAllTrials', 'Alpha'};                           

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmpc';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vmp';

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

% example object
dlist = nptDir;
% get entries in directory
dnum = size(dlist,1);

% check if the right conditions were met to create object
if(~isempty(dir(Args.RequiredFile)))
    
    ori = pwd;

    data.origin = {pwd}; 
%     pv = vmpv('auto', varargin{:});


    %%%% PATCH
    cd ..; cd ..; cd ..;
    pv = load([num2str(Args.pix) 'vmpv.mat']);
    % pv = load('vmpv_YY_1px.mat');
    pv = pv.pv;
    %%%%%%%

    cd(ori);
    spiketrain = load(Args.RequiredFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     nTrials = size(pv.trial_intervals,1);
% %     midTrial = ceil(nTrials/2);
%     midTrial = pv.last_trial_first_half(1);
%     sessionTimeInds = find(um.data.sessionTime(:,2) == 0); % Look for trial end markers
%     sessionTimeInds(1) = []; % Remove first marker which marks first cue-off, instead of trial end
    
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
        flat_spiketimes = round(flat_spiketimes); % there are occasional errors where these are not integers and prevent indexing using these 

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

        % back-filling spikes for view bins that occupy the same time bin
        stc(stc(:,6)==0,6) = nan;
        stc(:,7) = stc(:,5)~=0;
        stc(isnan(stc(:,6)) & stc(:,7), 6) = 0;
        stc(:,7) = [];
        stc(:,6) = fillmissing(stc(:,6), 'next');
        stc(isnan(stc(:,6)),6) = 0;
        % back-filling duration for view bins that occupy the same time bin
        stc_lasttime = stc(end,5); % Make sure if last duration sample is zero, it remains zero, not nan
        stc(stc(:,5)==0,5) = nan;
        stc(end,5) = stc_lasttime;
        stc(:,5) = fillmissing(stc(:,5), 'next'); % [timestamp place hd view dur spk]
        
        % Remove non-place and non-view rows for duration
        stc_filt = stc(find(conditions==1),:); 
        stc_filt(~(stc_filt(:,2) > 0),:) = []; % remove place bin = 0
        stc_filt(isnan(stc_filt(:,4)),:) = []; % remove NaN view bins
        stc_filt(~(stc_filt(:,3) > 0),:) = []; % remove hd bin = 0
        stc_ss = stc_filt(:,[2 5]); % [place dur];
        stc_ss = [stc_ss; [1600 0]];

        gpdurfull = accumarray(stc_ss(:,1),stc_ss(:,2))';
        
        % Remove low observation bins
        spike_count = zeros(Args.NumShuffles+1,Args.GridSteps*Args.GridSteps);
        gpdur = zeros(1,Args.GridSteps*Args.GridSteps);
        spike_count(:,bins_sieved) = spike_count_full(:,bins_sieved);
        gpdur(1,bins_sieved) = gpdurfull(1,bins_sieved);
        
        maps_raw = spike_count./repmat(gpdur,size(spike_count,1),1);
        
        % Save raw maps
        to_save = maps_raw(1,:);
        dur_raw = gpdur;
        spk_raw = spike_count(1,:);
        if repeat == 1
            data.sessionTimeC = stc; % Full unfiltered 
            data.stcfilt = stc_filt; % Filtered for conditions and having all place, view, hd data
            data.maps_raw = to_save;
            data.dur_raw = dur_raw;
            data.spk_raw = spk_raw;
            data.filtspknum = sum(spk_raw);
        elseif repeat == 2
            data.maps_raw1 = to_save;
            data.dur_raw1 = dur_raw;
            data.spk_raw1 = spk_raw;
        elseif repeat == 3
            data.maps_raw2 = to_save;
            data.dur_raw2 = dur_raw;
            data.spk_raw2 = spk_raw;
        end
        
            if 1 % Smoothing
                
                if repeat == 1
                    disp('Adaptive smoothing...');
                end
                nan_track = isnan(to_save);

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
                
                gpdur1 = repmat(gpdur',1,Args.NumShuffles+1);
                preset_to_zeros = gpdur1 == 0; 
                
                % Switch from linear maps to grid maps
                gpdur1 = cell2mat(lineartogrid(gpdur1,'place',[Args.GridSteps Args.GridSteps]));
                firing_counts_full1 = cell2mat(lineartogrid(spike_count','place',[Args.GridSteps Args.GridSteps]));
                preset_to_zeros = logical(cell2mat(lineartogrid(preset_to_zeros,'place',[Args.GridSteps Args.GridSteps])));
                firing_rates_full_raw1 = cell2mat(lineartogrid(maps_raw','place',[Args.GridSteps Args.GridSteps]));
                
                unvis = ~(gpdur1>0);
                % Boxcar smoothing
                maps_bcsm_grid = smooth(firing_rates_full_raw1,5,unvis,'boxcar');
                pos_bcsm_grid = smooth(gpdur1,5,unvis,'boxcar');
                % Disk smoothing
                maps_dksm_grid = smooth(firing_rates_full_raw1,5,unvis,'disk');
                pos_dksm_grid = smooth(gpdur1,5,unvis,'disk');

                % Set up adaptive smoothing parameters and output vars
                to_compute = 1:0.5:Args.GridSteps/2;
                possible = NaN(length(to_compute),2,Args.GridSteps,Args.GridSteps,Args.NumShuffles + 1);
                to_fill = NaN(size(possible,3), size(possible,4), size(possible,5));
                to_fill(preset_to_zeros) = 0;
                to_fill_time = NaN(size(possible,3), size(possible,4), size(possible,5));
                to_fill_time(preset_to_zeros) = 0; 
                to_fill_radius = NaN(size(possible,3), size(possible,4), size(possible,5));
                to_fill_radius(preset_to_zeros) = 0;
                
                wip = ones(Args.NumShuffles+1,1);
                % Adaptive smoothing
                for idx = 1:length(to_compute)

                    f=fspecial('disk',to_compute(idx));
                    f(f>=(max(max(f))/3))=1;
                    f(f~=1)=0;

                    possible(idx,1,:,:,:) = repmat(imfilter(gpdur1(:,:,1), f, 'conv'), 1,1,Args.NumShuffles+1);   %./scaler;
                    possible(idx,2,:,:,find(wip)) = imfilter(firing_counts_full1(:,:,find(wip)), f, 'conv');   %./scaler;

                    logic1 = squeeze(alpha./(possible(idx,1,:,:,:).*sqrt(possible(idx,2,:,:,:))) <= to_compute(idx));
                    slice1 = squeeze(possible(idx,1,:,:,:));
                    slice2 = squeeze(possible(idx,2,:,:,:));

                    to_fill(logic1 & isnan(to_fill)) = slice2(logic1 & isnan(to_fill))./slice1(logic1 & isnan(to_fill));
                    to_fill_time(logic1 & isnan(to_fill_time)) = slice1(logic1 & isnan(to_fill_time));
                    to_fill_radius(logic1 & isnan(to_fill_radius)) = to_compute(idx);

%                     disp('smoothed with kernel size:');
%                     disp(to_compute(idx));
%                     disp('grids left');
%                     disp(sum(sum(sum(isnan(to_fill(:,:,:))))));

                    check = squeeze(sum(sum(isnan(to_fill),2),1));
                    wip(check==0) = 0;

                end
                
                % Reshape from grid to linear maps
                to_fill(preset_to_zeros) = nan; % unvisited bins should be nan
                to_fill = gridtolinear({to_fill},'place',[Args.GridSteps Args.GridSteps]);
                to_fill = to_fill';
                to_fill_time(isnan(to_fill_time) | preset_to_zeros) = 0;
                to_fill_time = gridtolinear({to_fill_time},'place',[Args.GridSteps Args.GridSteps]);
                to_fill_time = to_fill_time';
                to_fill_radius(preset_to_zeros) = nan;
                to_fill_radius = gridtolinear({to_fill_radius},'place',[Args.GridSteps Args.GridSteps]);
                to_fill_radius = to_fill_radius';
                maps_bcsm = gridtolinear({maps_bcsm_grid},'place',[Args.GridSteps Args.GridSteps]);
                maps_bcsm = maps_bcsm';
                maps_dksm = gridtolinear({maps_dksm_grid},'place',[Args.GridSteps Args.GridSteps]);
                maps_dksm = maps_dksm';
                pos_bcsm_grid(isnan(pos_bcsm_grid) | preset_to_zeros) = 0;
                pos_bcsm = gridtolinear({pos_bcsm_grid},'place',[Args.GridSteps Args.GridSteps]);
                pos_bcsm = pos_bcsm';
                pos_dksm_grid(isnan(pos_dksm_grid) | preset_to_zeros) = 0;
                pos_dksm = gridtolinear({pos_dksm_grid},'place',[Args.GridSteps Args.GridSteps]);
                pos_dksm = pos_dksm';
                
                firing_rates_full = to_fill;

                % smoothing part ends
                switch Args.SmoothType
                    case 'Adaptive'
                        maps_sm = firing_rates_full;
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
                    data.maps_adsm = firing_rates_full(1,:);
                    data.maps_adsmsh = firing_rates_full(2:end,:);
                    data.dur_adsm = to_fill_time(1,:);
                    data.dur_adsmsh = to_fill_time(2:end,:);
                    data.radii = to_fill_radius(1,:);
                    data.radiish = to_fill_radius(2:end,:);
                    data.maps_bcsm = maps_bcsm(1,:);
                    data.maps_bcsmsh = maps_bcsm(2:end,:);
                    data.maps_dksm = maps_dksm(1,:);
                    data.maps_dksmsh = maps_dksm(2:end,:);
                    data.maps_sm = maps_sm(1,:);
                    data.maps_smsh = maps_sm(2:end,:);
                elseif repeat == 2
                    data.maps_adsm1 = firing_rates_full(1,:);
                    data.dur_adsm1 = to_fill_time(1,:); 
                    data.radii1 = to_fill_radius(1,:);
                    data.maps_bcsm1 = maps_bcsm(1,:);
                    data.maps_dksm1 = maps_dksm(1,:);
                    data.maps_sm1 = maps_sm(1,:);
                    data.maps_smsh1 = maps_sm(2:end,:);
                elseif repeat == 3
                    data.maps_adsm2 = firing_rates_full(1,:);
                    data.dur_adsm2 = to_fill_time(1,:);
                    data.radii2 = to_fill_radius(1,:);
                    data.maps_bcsm2 = maps_bcsm(1,:);
                    data.maps_dksm2 = maps_dksm(1,:);
                    data.maps_sm2 = maps_sm(1,:);
                    data.maps_smsh2 = maps_sm(2:end,:);
                end

            else
                firing_rates_full = maps_raw;
                to_fill_time = repmat(gpdur,Args.NumShuffles+1,1); % HM added

            end

            % Calculating SIC for adaptive smoothing
            disp('Calculating SIC...');
            sic_adsm = skaggs_sic(firing_rates_full',to_fill_time');
            sic_adsm = sic_adsm';
            
            % Calculating SIC for boxcar smoothing
            sic_bcsm = skaggs_sic(maps_bcsm',pos_bcsm');
            sic_bcsm = sic_bcsm';
            
            % Calculating SIC for disk smoothing
            sic_dksm = skaggs_sic(maps_dksm',pos_dksm');
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
        %     data.median_occ_firings = median_stats';
        %     data.variance_occ_firings = var_stats';
        %     data.perc_occ_firings = perc_stats';
        %     data.occ_data = occ_data;
        elseif repeat == 2
            data.SIC_adsm1 = sic_adsm;
            data.SIC_bcsm1 = sic_bcsm;
            data.SIC_dksm1 = sic_dksm;
            data.crit_sm1 = crit_sm;
        elseif repeat == 3
            data.SIC_adsm2 = sic_adsm;
            data.SIC_bcsm2 = sic_bcsm;
            data.SIC_dksm2 = sic_dksm;
            data.crit_sm2 = crit_sm;
        end

        %     data.median_occ_firings = median_stats';
        %     data.variance_occ_firings = var_stats';
        %     data.perc_occ_firings = perc_stats';

        %     data.occ_data = occ_data;

            
    end
    
    % create nptdata so we can inherit from it    
    data.gridSteps = Args.GridSteps;
    Args.NumShuffles = NumShuffles_saved;
    data.numSets = 1;
    data.Args = Args;
    n = nptdata(1,0,pwd);
    d.data = data;
    obj = class(d,Args.classname,n);
    saveObject(obj,'ArgsC',Args);

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
