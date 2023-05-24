function [obj, varargout] = vmsv(varargin)
%@vmsv Constructor function for vmsv class
%   OBJ = vmsv(varargin)
%
%   OBJ = vmsv('auto') attempts to create a vmsv object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmsv %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmsv('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell', 'RequiredFile','spiketrain.mat', ...
				'GridSteps',40, 'pix',100,...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4, 'ThresVel',1, 'UseMinObs', 0, 'SmoothType', 'Adaptive', 'SelectiveCriteria','SIC',...
                'UseAllTrials',1,'Alpha',1000);
            
Args.flags = {'Auto','ArgsOnly','FRSIC','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps','NumShuffles','AdaptiveSmooth','UseMinObs','ThresVel','UseAllTrials','Alpha'};                          

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmsv';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vms';

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
% %     pv = vmpv('auto','save','MinObsView',5,'MinDurView',0.01);
    %%%%%%%
    % Patch %%%%%%
    cd ..; cd ..; cd ..;
    pv = load([num2str(Args.pix) 'vmpv.mat']);
    % pv = load('vmpv_YY_1px.mat');
    pv = pv.pv;
    %%%%%%%%%%
    %
    cd(ori);
    spiketrain = load(Args.RequiredFile);   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NumShuffles_saved = Args.NumShuffles;
    for repeat = 1:3 % 1 = full trial, 2 = 1st half, 3 = 2nd half
        
        if repeat > 1
            Args.NumShuffles = 0;
        end

        if repeat == 1
            stc = pv.data.sessionTimeC;
            disp('Full session:');
        elseif repeat == 2
            disp('1st half:');
        elseif repeat == 3
            disp('2nd half:');
        end
        
        % spike shuffling

        spiketimes = spiketrain.timestamps/1000; % now in seconds
        maxTime = pv.data.rplmaxtime;
        spiketimes(spiketimes>maxTime) = [];
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
        flat_spiketimes = flat_spiketimes'; % row 1: spiketimes concatenated in shuffle batches, serially; row 2: identity of shuffle, in ascending order
        flat_spiketimes = sortrows(flat_spiketimes); % spiketimes sorted serially across all data, shuffle identity mixed

        flat_spiketimes(flat_spiketimes(:,1) < stc(1,1),:) = [];      
        
        % selecting rows from sessionTimeC
        if repeat == 1
            disp('    Filtering...')
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
            bins_sieved = pv.data.view_good_bins;
            conditions = conditions & (pv.data.pv_good_rows); % Make sure maps take into account both place and view filters
        else
            bins_sieved = 1:5122;
        end

        if repeat == 1
            % Group into intervals those consecutive rows where same place bin is occupied
            dstc = diff(stc(:,1));
            stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
            stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
            stc_changing_ind = stc_changing_ind(1:end-1,:);
        end

        consol_arr = zeros(5122,Args.NumShuffles + 1);
        consol_arr_p = zeros(Args.GridSteps * Args.GridSteps,Args.NumShuffles + 1);

        if repeat == 1
            disp(['    Assigning ' num2str(size(flat_spiketimes,1)) ' spikes to bins...']);
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
            consol_arr(bins_hit(:,3),flat_spiketimes(sp,2)) = consol_arr(bins_hit(:,3),flat_spiketimes(sp,2)) + 1;
            consol_arr_p(bins_hit(:,1),flat_spiketimes(sp,2)) = consol_arr_p(bins_hit(:,1),flat_spiketimes(sp,2)) + 1; % place spikes

        end        

        spike_count_full = consol_arr;

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
        stc_ss = stc_filt(:,[4 5]); % [view dur];
        stc_ss = [stc_ss; [5122 0]];
        
        gpdur1_full = accumarray(stc_ss(:,1),stc_ss(:,2))';
        
        % Remove low obs bins from rest of analysis
        gpdur1 = zeros(size(gpdur1_full));
        gpdur1(bins_sieved) = gpdur1_full(bins_sieved);
        spike_count = zeros(size(spike_count_full));
        spike_count(bins_sieved,:) = spike_count_full(bins_sieved,:);
        
        % Store raw map output
        lin_o_i_Gaze = repmat(gpdur1', 1, Args.NumShuffles + 1);
        lin_spikeLoc_Gaze = spike_count;
        maps_raw = lin_spikeLoc_Gaze./lin_o_i_Gaze;
        maps_raw_out = maps_raw(:,1);

        if repeat == 1
            data.sessionTimeC = stc; % Full unfiltered 
            data.stcfilt = stc_filt; % Filtered for conditions and leaving only place and view
            data.maps_raw = maps_raw_out';
            data.dur_raw = gpdur1;
            data.spk_raw = spike_count(:,1)';
            data.filtspknum = sum(consol_arr_p(:,1));
        elseif repeat == 2
            data.maps_raw1 = maps_raw_out';
            data.dur_raw1 = gpdur1;
            data.spk_raw1 = spike_count(:,1)';
        else
            data.maps_raw2 = maps_raw_out';
            data.dur_raw2 = gpdur1;
            data.spk_raw2 = spike_count(:,1)';
        end
        
        if repeat == 1
            disp('    Adaptive smoothing...');
        end
        if 1
            
            gazeSections = {'Cue', 'Hint', 'Ground', 'Ceiling', 'Walls', 'Pillar1', 'Pillar2', 'Pillar3', 'Pillar4'};
            binDepths = [1 1;
                1 1;
                40 40;
                40 40;
                8 160;
                5 32;
                5 32;
                5 32;
                5 32];
                
            % Assign linear bin to grid
%             unvis = ~(lin_o_i_Gaze>0);
            [grid_o_i_Gaze] = lineartogrid(lin_o_i_Gaze,'view',binDepths);
            [grid_spikeBin_Gaze] = lineartogrid(lin_spikeLoc_Gaze,'view',binDepths);
            [grid_rate_Gaze] = lineartogrid(maps_raw,'view',binDepths);
%             [grid_unvis] = lineartogrid(unvis,'view',binDepths);
            emptyfloorref = grid_o_i_Gaze;
            
            % Pad grids with << 5 >> extra rows of adjoining grids
            n = 5;
            padpillar = false;
            [emptyfloorref_pad,~] = padsvmap(n,grid_o_i_Gaze,gazeSections,padpillar);
            padpillar = true;
            [grid_o_i_Gaze,retrievemap] = padsvmap(n,grid_o_i_Gaze,gazeSections,padpillar);
            [grid_spikeBin_Gaze,~] = padsvmap(n,grid_spikeBin_Gaze,gazeSections,padpillar);
            [grid_rate_Gaze,~] = padsvmap(n,grid_rate_Gaze,gazeSections,padpillar);
            
            
%             [grid_unvis,~] = padsvmap(n,grid_unvis,gazeSections);
            
            % Boxcar and disk smoothing of padded grids
            grid_bcsm_rate = cell(size(binDepths,1),1);
            grid_bcsm_rate{1} = grid_rate_Gaze{1}; % No need to smooth cue
            grid_bcsm_rate{2} = grid_rate_Gaze{2}; % No need to smooth hint
            grid_bcsm_pos = cell(size(binDepths,1),1);
            grid_bcsm_pos{1} = grid_o_i_Gaze{1};
            grid_bcsm_pos{2} = grid_o_i_Gaze{2};
            grid_dksm_rate = cell(size(binDepths,1),1);
            grid_dksm_rate{1} = grid_rate_Gaze{1}; % No need to smooth cue
            grid_dksm_rate{2} = grid_rate_Gaze{2}; % No need to smooth hint
            grid_dksm_pos = cell(size(binDepths,1),1);
            grid_dksm_pos{1} = grid_o_i_Gaze{1};
            grid_dksm_pos{2} = grid_o_i_Gaze{2};
            for jj = 3:size(binDepths,1)
                unvis = ~(emptyfloorref_pad{jj}>0) | isnan(grid_rate_Gaze{jj});
                % Boxcar smooth
                grid_bcsm_rate{jj} = smooth(grid_rate_Gaze{jj},5,unvis,'boxcar');
                grid_bcsm_pos{jj} = smooth(grid_o_i_Gaze{jj},5,unvis,'boxcar');
                % Disk smooth
                grid_dksm_rate{jj} = smooth(grid_rate_Gaze{jj},5,unvis,'disk');
                grid_dksm_pos{jj} = smooth(grid_o_i_Gaze{jj},5,unvis,'disk');
            end

            % Adaptive smoothing of padded grids
            alpha = Args.Alpha;
            grid_smoothed_Gaze = cell(size(binDepths,1),1);
            grid_smoothed_Gaze{1} = grid_spikeBin_Gaze{1}./grid_o_i_Gaze{1}; % No need to smooth cue
            grid_smoothed_Gaze{2} = grid_spikeBin_Gaze{2}./grid_o_i_Gaze{2}; % No need to smooth hint
            grid_smoothed_dur = cell(size(grid_o_i_Gaze,1),1);
            grid_smoothed_dur{1} = grid_o_i_Gaze{1};
            grid_smoothed_dur{2} = grid_o_i_Gaze{2};
            grid_ad_size = cell(size(grid_o_i_Gaze,1),1);
            grid_ad_size{1} = nan(size(grid_o_i_Gaze{1}));
            grid_ad_size{2} = nan(size(grid_o_i_Gaze{1}));
            for jj = 3:size(grid_o_i_Gaze,1) % for each grid, floor onwards
                
%                 disp(['    ...Smoothing grid ' num2str(jj)]);
                wip = ones(Args.NumShuffles+1,1);
                gpdur1 = grid_o_i_Gaze{jj};
                preset_to_zeros = gpdur1(:,:,1);
                preset_to_zeros(find(preset_to_zeros>0)) = 1;
                preset_to_zeros(find(preset_to_zeros~=1)) = 0;
                preset_to_zeros = ~preset_to_zeros;
                preset_to_zeros = repmat(preset_to_zeros, [1,1,size(gpdur1,3)]);
                
                firing_counts_full1 = grid_spikeBin_Gaze{jj};
                gpdur1(isnan(gpdur1)) = 0;
                firing_counts_full1(isnan(firing_counts_full1)) = 0;
                
%                 to_compute = 1:0.5:Args.GridSteps/2; % unit bin is actually fspecial(...0.5)
                to_compute = 1:0.5:(max(size(grid_o_i_Gaze{jj}(:,:,1))))/2;
                
                possible = NaN(2,size(firing_counts_full1,1),size(firing_counts_full1,2),Args.NumShuffles + 1);
                to_fill = NaN(size(possible,2), size(possible,3), size(possible,4));
                to_fill(preset_to_zeros) = 0;
                to_fill_smoothed_duration = NaN(size(possible,2), size(possible,3), size(possible,4));
                to_fill_smoothed_duration(preset_to_zeros) = 0;
                to_fill_size = NaN(size(possible,2), size(possible,3), size(possible,4));
                to_fill_size(preset_to_zeros) = 0;
                
                for idx = 1:length(to_compute)
                    
                    f=fspecial('disk',to_compute(idx));
                    f(f>=(max(max(f))/3))=1;
                    f(f~=1)=0;
                    
                    possible(1,:,:,:) = repmat(imfilter(gpdur1(:,:,1), f, 'conv'), 1,1,Args.NumShuffles+1);
                    possible(2,:,:,find(wip)) = imfilter(firing_counts_full1(:,:,find(wip)), f, 'conv');
                    
                    logic1 = squeeze(alpha./(possible(1,:,:,:).*sqrt(possible(2,:,:,:))) <= to_compute(idx));
                    
                    %debug
                    %                         logic1(~logic1) = 1;
                    
                    slice1 = squeeze(possible(1,:,:,:));
                    slice2 = squeeze(possible(2,:,:,:));
                    
                    to_fill(logic1 & isnan(to_fill)) = slice2(logic1 & isnan(to_fill))./slice1(logic1 & isnan(to_fill));
                    to_fill_smoothed_duration(logic1 & isnan(to_fill_smoothed_duration)) = slice1(logic1 & isnan(to_fill_smoothed_duration));
                    to_fill_size(logic1 & isnan(to_fill_size)) = to_compute(idx);
                    
                    
                    remaining = sum(sum(sum(isnan(to_fill(:,:,:)))));
%                     disp(['smoothed grid ' num2str(jj) ' with kernel size ' num2str(to_compute(idx)) ', leaving ' num2str(remaining) ' grids undone']);
                    
                    check = squeeze(sum(sum(isnan(to_fill),2),1));
                    wip(check==0) = 0;
                    
                    if remaining == 0
%                         disp('done');
                        break;
                    end
                end
                
                to_fill(preset_to_zeros) = nan;
                to_fill_size(preset_to_zeros) = nan;
                grid_smoothed_Gaze{jj} = to_fill;
                grid_smoothed_dur{jj} = to_fill_smoothed_duration;
                grid_ad_size{jj} = to_fill_size;
                
            end
            
            % Unpad smoothed maps
            grid_smoothed_Gaze = unpadsvmap(grid_smoothed_Gaze,retrievemap,emptyfloorref);
            grid_smoothed_dur = unpadsvmap(grid_smoothed_dur,retrievemap,emptyfloorref);
            grid_ad_size = unpadsvmap(grid_ad_size,retrievemap,emptyfloorref);
            grid_bcsm_rate = unpadsvmap(grid_bcsm_rate,retrievemap,emptyfloorref);
            grid_bcsm_pos = unpadsvmap(grid_bcsm_pos,retrievemap,emptyfloorref);
            grid_dksm_rate = unpadsvmap(grid_dksm_rate,retrievemap,emptyfloorref);
            grid_dksm_pos = unpadsvmap(grid_dksm_pos,retrievemap,emptyfloorref);
            
            % Convert from grid maps back to linear maps
            gpdur1 = gridtolinear(grid_smoothed_dur,'view',binDepths);
            gpdur1(isnan(gpdur1)) = 0;
            gpdur1 = gpdur1';
            lambda_i = gridtolinear(grid_smoothed_Gaze,'view',binDepths);
            lambda_i = lambda_i';
            radii = gridtolinear(grid_ad_size,'view',binDepths);
            radii = radii';
            maps_bcsm = gridtolinear(grid_bcsm_rate,'view',binDepths);
            maps_bcsm = maps_bcsm';
            pos_bcsm = gridtolinear(grid_bcsm_pos,'view',binDepths);
            pos_bcsm(isnan(pos_bcsm)) = 0;
            pos_bcsm = pos_bcsm';
            maps_dksm = gridtolinear(grid_dksm_rate,'view',binDepths);
            maps_dksm = maps_dksm';
            pos_dksm = gridtolinear(grid_dksm_pos,'view',binDepths);
            pos_dksm(isnan(pos_dksm)) = 0;
            pos_dksm = pos_dksm';
            
            switch Args.SmoothType
                case 'Adaptive'
                    maps_sm = lambda_i;
                case 'Boxcar'
                    maps_sm = maps_bcsm;
                case 'Disk'
                    maps_sm = maps_dksm;
            end
            
            if repeat == 1
                if data.filtspknum < 100 && max(maps_sm(1,:),[],'omitnan') < 0.7
                    data.discard = true;
                else 
                    data.discard = false;
                end
                data.maps_adsm = lambda_i(1,:);
                data.maps_adsmsh = lambda_i(2:end,:);
                data.radii = radii(1,:);
                data.radiish = radii(2:end,:);
                data.dur_adsm = gpdur1(1,:);
                data.dur_adsmsh = gpdur1(2:end,:);
                data.maps_bcsm = maps_bcsm(1,:);
                data.maps_bcsmsh = maps_bcsm(2:end,:);
                data.maps_dksm = maps_dksm(1,:);
                data.maps_dksmsh = maps_dksm(2:end,:);
                data.maps_sm = maps_sm(1,:);
                data.maps_smsh = maps_sm(2:end,:);
            elseif repeat == 2
                data.maps_adsm1 = lambda_i(1,:);
                data.radii1 = radii(1,:);
                data.dur_adsm1 = gpdur1(1,:);
                data.maps_bcsm1 = maps_bcsm(1,:);
                data.maps_dksm1 = maps_dksm(1,:);
                data.maps_sm1 = maps_sm(1,:);
                data.maps_smsh1 = maps_sm(2:end,:);
            elseif repeat == 3
                data.maps_adsm2 = lambda_i(1,:);
                data.radii2 = radii(1,:);
                data.dur_adsm2 = gpdur1(1,:);
                data.maps_bcsm2 = maps_bcsm(1,:);
                data.maps_dksm2 = maps_dksm(1,:);
                data.maps_sm2 = maps_sm(1,:);
                data.maps_smsh2 = maps_sm(2:end,:);
            end            
            
        else
            
            gpdur1 = repmat(gpdur1, Args.NumShuffles + 1, 1);
            lambda_i = spike_count'./ gpdur1;
            
        end

  
            
            % SIC portion
            if repeat == 1
                disp('    Calculating SIC...');
            end
            % Adaptive smoothing SIC
            sic_adsm = skaggs_sic(lambda_i',gpdur1');
            sic_adsm = sic_adsm';
            
            % Boxcar smoothing SIC
            sic_bcsm = skaggs_sic(maps_bcsm',pos_bcsm');
            sic_bcsm = sic_bcsm';
            
            % Disk smoothing SIC
            sic_dksm = skaggs_sic(maps_dksm',pos_dksm');
            sic_dksm = sic_dksm';

%             % ISE portion
%             % create overall map and insert padded portions in, to account for
%             % cross-portion pairs
%             tic;
%             firing_rates = lambda_i;
%             if repeat == 1
%                 disp('    Calculating ISE');
%             end
%             canvas = nan(51, 161, Args.NumShuffles + 1);
%             % flooring
%             floor_padded = nan(42,42,Args.NumShuffles+1);
%             floor_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);
%             floor_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,3203:3203+39),[2 1]), 40, 1, Args.NumShuffles+1),1);
%             floor_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,3243:3243+39),[2 1]), 1, 40, Args.NumShuffles+1);
%             floor_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,3283:3283+39),[2 1]), 40, 1, Args.NumShuffles+1);
%             floor_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,3323:3323+39),[2 1]), 1, 40, Args.NumShuffles+1), 2);
%             canvas(10:end,1:42,:) = floor_padded;
% 
%             % ceiling
%             ceiling_padded = nan(42,42,Args.NumShuffles+1);
%             ceiling_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,1603:3202),size(firing_rates,1),40,40), [3 2 1]), 1);
%             ceiling_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,4323:4323+39),[2 1]), 40, 1, Args.NumShuffles+1),1);
%             ceiling_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,4363:4363+39),[2 1]), 1, 40, Args.NumShuffles+1);
%             ceiling_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,4403:4403+39),[2 1]), 40, 1, Args.NumShuffles+1);
%             ceiling_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,4443:4443+39),[2 1]), 1, 40, Args.NumShuffles+1), 2);
%             canvas(10:end,44:85,:) = ceiling_padded;
% 
%             % walls
%             walls_padded = nan(8,161,Args.NumShuffles+1);
%             walls_padded(:,1:end-1,:) = flip(permute(reshape(firing_rates(:,3203:3203+1280-1), Args.NumShuffles+1, 40*4, 8),[3 2 1]), 1);
%             walls_padded(:,end,:) = walls_padded(:,1,:);
%             canvas(1:8,:,:) = walls_padded;
% 
%             % used to pad pillar base more easily
%             floor_base = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);
% 
%             % pillars
%             PTL_padded = nan(6,33,Args.NumShuffles+1);
%             PTL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4963:4963+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
%             % small diagonal issue here, diagonal floor bins at the corners are put
%             % side by side, only 16 such occurrences in total, neglected for now.
%             PTL_padded(end,1:8,:) = flip(permute(floor_base(9:16,8,:),[2 1 3]),2);
%             PTL_padded(end,9:16,:) = floor_base(8,9:16,:);
%             PTL_padded(end,17:24,:) = permute(floor_base(9:16,17,:),[2 1 3]);
%             PTL_padded(end,25:32,:) = flip(floor_base(17,9:16,:),2);
%             PTL_padded(:,end,:) = PTL_padded(:,1,:);
%             canvas(10:10+6-1,87:87+32,:) = PTL_padded;
% 
%             PTR_padded = nan(6,33,Args.NumShuffles+1);
%             PTR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4803:4803+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
%             PTR_padded(end,1:8,:) = flip(permute(floor_base(9:16,24,:),[2 1 3]),2);
%             PTR_padded(end,9:16,:) = floor_base(8,25:32,:);
%             PTR_padded(end,17:24,:) = permute(floor_base(9:16,33,:),[2 1 3]);
%             PTR_padded(end,25:32,:) = flip(floor_base(17,25:32,:),2);
%             PTR_padded(:,end,:) = PTR_padded(:,1,:);
%             canvas(10:10+6-1,121:121+32,:) = PTR_padded;
% 
%             PBL_padded = nan(6,33,Args.NumShuffles+1);
%             PBL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4643:4643+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
%             PBL_padded(end,1:8,:) = flip(permute(floor_base(25:32,8,:),[2 1 3]),2);
%             PBL_padded(end,9:16,:) = floor_base(24,9:16,:);
%             PBL_padded(end,17:24,:) = permute(floor_base(25:32,17,:),[2 1 3]);
%             PBL_padded(end,25:32,:) = flip(floor_base(33,9:16,:),2);
%             PBL_padded(:,end,:) = PBL_padded(:,1,:);
%             canvas(17:17+6-1,87:87+32,:) = PBL_padded;
% 
%             PBR_padded = nan(6,33,Args.NumShuffles+1);
%             PBR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4483:4483+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
%             PBR_padded(end,1:8,:) = flip(permute(floor_base(25:32,24,:),[2 1 3]),2);
%             PBR_padded(end,9:16,:) = floor_base(24,25:32,:);
%             PBR_padded(end,17:24,:) = permute(floor_base(25:32,33,:),[2 1 3]);
%             PBR_padded(end,25:32,:) = flip(floor_base(33,25:32,:),2);
%             PBR_padded(:,end,:) = PBR_padded(:,1,:);
%             canvas(17:17+6-1,121:121+32,:) = PBR_padded;
% 
%             actual_image = canvas(:,:,1);
%             actual_image = actual_image(:)';
%             shuffled_images = canvas(:,:,2:end);
%             shuffled_images = reshape(shuffled_images, size(shuffled_images,3),size(shuffled_images,1)*size(shuffled_images,2));
% 
%             disp(['time taken to pad map for ISE: ' num2str(toc)]);
%             tic;
% 
%             ise_out = ise(actual_image, shuffled_images, 51, 161);
%             disp(['time taken to compute ISE: ' num2str(toc)]);
            
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
%                 case 'ISE'
%                     crit_sm = ise_out;
            end
            
            if repeat == 1
%                 data.flattened = squeeze(canvas(:,:,1));
                data.SIC_adsm = sic_adsm(1);
                data.SICsh_adsm = sic_adsm(2:end,1);
%                 data.ISE_adsm = ise_out(1);
%                 data.ISEsh_adsm = ise_out(2:end,1);
                data.SIC_bcsm = sic_bcsm(1);
                data.SICsh_bcsm = sic_bcsm(2:end,1);
                data.SIC_dksm = sic_dksm(1);
                data.SICsh_dksm = sic_dksm(2:end,1);
                data.crit_sm = crit_sm(1);
                data.critsh_sm = crit_sm(2:end,1);
                data.critthrcell = prctile(crit_sm(2:end,1),95);
            elseif repeat == 2
                data.SIC_adsm1 = sic_adsm;
%                 data.ISE_adsm1 = ise_out;
                data.SIC_bcsm1 = sic_bcsm;
                data.SIC_dksm1 = sic_dksm;
                data.crit_sm1 = crit_sm;
            elseif repeat == 3
                data.SIC_adsm2 = sic_adsm;
%                 data.ISE_adsm2 = ise_out;
                data.SIC_bcsm2 = sic_bcsm;
                data.SIC_dksm2 = sic_dksm;
                data.crit_sm2 = crit_sm;
            end            
            
    end
    
    % create nptdata so we can inherit from it 
    Args.NumShuffles = NumShuffles_saved;
    data.binDepths = binDepths;
    data.gridSteps = Args.GridSteps;
    data.gazeSections = gazeSections;
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

% % relic code
% % 
% function [retrievemap,o_i,spikeLoc,map] = padgrids(n,o_i,spikeLoc,grid_o_i,grid_spikeLoc,gazeSections,jj)
% 
% % Pad maps with adjoining bins from adjacent maps
% 
% switch gazeSections{jj}
%     case 'Ground'
%         wallsection_ind = strcmp(gazeSections,'Walls');
%         wall_o_i = grid_o_i{wallsection_ind};
%         wall_spikeLoc = grid_spikeLoc{wallsection_ind};
% 
%         % Move original map to middle
%         o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
%         o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
%         spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
%         spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;
% 
%         % Pad with wall data
%         o_i_temp(1:n,n+1:n+size(o_i,1),:) = wall_o_i(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:); % top
%         o_i_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:),-1); % right
%         o_i_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:),-2); % bottom
%         o_i_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:),1); % left
%         spikeLoc_temp(1:n,n+1:n+size(o_i,1),:) = wall_spikeLoc(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:); % top
%         spikeLoc_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:),-1); % right
%         spikeLoc_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:),-2); % bottom
%         spikeLoc_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:),1); % left
% 
%         % Save indices of original grid [from_x to_x; from_y to_y]
%         retrievemap = [n+1 n+size(o_i,1); ...
%                        n+1 n+size(o_i,2)];
%         % Send vars for adaptive smoothing
%         o_i = o_i_temp;
%         spikeLoc = spikeLoc_temp;
% 
%     case 'Ceiling'
%         wallsection_ind = strcmp(gazeSections,'Walls');
%         wall_o_i = grid_o_i{wallsection_ind};
%         wall_spikeLoc = grid_spikeLoc{wallsection_ind};
% 
%         % Flip walldata upside down
%         wall_o_i = flipud(wall_o_i);
%         wall_spikeLoc = flipud(wall_spikeLoc);
% 
%         % Move original map to middle
%         o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
%         o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
%         spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
%         spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;
% 
%         % Pad with wall data
%         o_i_temp(1:n,n+1:n+size(o_i,1),:) = fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:)); % top
%         o_i_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:)),-1); % right
%         o_i_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:)),-2); % bottom
%         o_i_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:)),1); % left
%         spikeLoc_temp(1:n,n+1:n+size(o_i,1),:) = fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:)); % top
%         spikeLoc_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:)),-1); % right
%         spikeLoc_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:)),-2); % bottom
%         spikeLoc_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:)),1); % left
% 
%         % Save indices of original grid [from_x to_x; from_y to_y]
%         retrievemap = [n+1 n+size(o_i,1); ...
%                        n+1 n+size(o_i,2)];
%         % Send vars for adaptive smoothing
%         o_i = o_i_temp;
%         spikeLoc = spikeLoc_temp;
% 
%     case 'Walls'
%         groundsection_ind = strcmp(gazeSections,'Ground');
%         ground_o_i = grid_o_i{groundsection_ind};
%         ground_spikeLoc = grid_spikeLoc{groundsection_ind};
% 
%         ceilingsection_ind = strcmp(gazeSections,'Ceiling');
%         ceiling_o_i = grid_o_i{ceilingsection_ind};
%         ceiling_spikeLoc = grid_spikeLoc{ceilingsection_ind};
% 
%         % Move original map to middle
%         o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
%         o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
%         spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
%         spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;
% 
%         % Pad with ground data
%         o_i_temp(n+size(o_i,1)+1:end,n+1:size(ground_o_i,2)+n,:) = rot90(ground_o_i(:,1:n,:),-1);
%         o_i_temp(n+size(o_i,1)+1:end,n+size(ground_o_i,2)+1:n+2*size(ground_o_i,2),:) = ground_o_i(1:n,:,:);
%         o_i_temp(n+size(o_i,1)+1:end,n+2*size(ground_o_i,2)+1:n+3*size(ground_o_i,2),:) = rot90(ground_o_i(:,size(ground_o_i,1)-n+1:end,:),1);
%         o_i_temp(n+size(o_i,1)+1:end,n+3*size(ground_o_i,1)+1:n+4*size(ground_o_i,1),:) = rot90(ground_o_i(size(ground_o_i,1)-n+1:end,:,:),2);
%         spikeLoc_temp(n+size(o_i,1)+1:end,n+1:size(ground_o_i,2)+n,:) = rot90(ground_spikeLoc(:,1:n,:),-1);
%         spikeLoc_temp(n+size(o_i,1)+1:end,n+size(ground_o_i,2)+1:n+2*size(ground_o_i,2),:) = ground_spikeLoc(1:n,:,:);
%         spikeLoc_temp(n+size(o_i,1)+1:end,n+2*size(ground_o_i,2)+1:n+3*size(ground_o_i,2),:) = rot90(ground_spikeLoc(:,size(ground_spikeLoc,1)-n+1:end,:),1);
%         spikeLoc_temp(n+size(o_i,1)+1:end,n+3*size(ground_o_i,1)+1:n+4*size(ground_o_i,1),:) = rot90(ground_spikeLoc(size(ground_spikeLoc,1)-n+1:end,:,:),2);
% 
%         % Pad with ceiling data
%         o_i_temp(1:n,n+1:size(ceiling_o_i,1)+n,:) = fliplr(rot90(ceiling_o_i(:,size(ceiling_o_i,1)-n+1:end,:),1));
%         o_i_temp(1:n,n+size(ceiling_o_i,1)+1:n+2*size(ceiling_o_i,1),:) = fliplr(ceiling_o_i(1:n,:,:));
%         o_i_temp(1:n,n+2*size(ceiling_o_i,1)+1:n+3*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_o_i(:,1:n,:),-1));
%         o_i_temp(1:n,n+3*size(ceiling_o_i,1)+1:n+4*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_o_i(size(ceiling_o_i,1)-n+1:end,:,:),2));
%         spikeLoc_temp(1:n,n+1:size(ceiling_o_i,1)+n,:) = fliplr(rot90(ceiling_spikeLoc(:,size(ceiling_spikeLoc,1)-n+1:end,:),1));
%         spikeLoc_temp(1:n,n+size(ceiling_o_i,1)+1:n+2*size(ceiling_o_i,1),:) = fliplr(ceiling_spikeLoc(1:n,:,:));
%         spikeLoc_temp(1:n,n+2*size(ceiling_o_i,1)+1:n+3*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_spikeLoc(:,1:n,:),-1));
%         spikeLoc_temp(1:n,n+3*size(ceiling_o_i,1)+1:n+4*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_spikeLoc(size(ceiling_spikeLoc,1)-n+1:end,:,:),2));
% 
%         % Pad with wall data on either end
%         o_i_temp(n+1:n+size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
%         o_i_temp(n+1:n+size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
%         spikeLoc_temp(n+1:n+size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
%         spikeLoc_temp(n+1:n+size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);
% 
%         % Save indices of original grid [from_x to_x; from_y to_y]
%         retrievemap = [n+1 n+size(o_i,1); ...
%                        n+1 n+size(o_i,2)];
%         % Send vars for adaptive smoothing
%         o_i = o_i_temp;
%         spikeLoc = spikeLoc_temp;
% 
%     case 'Pillar1'
%         groundsection_ind = strcmp(gazeSections,'Ground');
%         ground_o_i = grid_o_i{groundsection_ind};
%         ground_spikeLoc = grid_spikeLoc{groundsection_ind};
% 
%         % Move original map to middle
%         o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
%         spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;
% 
%         % Pad with ground data
%         o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(25:32,25-n:24,:),-1);
%         o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(25-n:24,25:32,:);
%         o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(25:32,33:32+n,:),1);
%         o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(33:32+n,25:32,:),2);
%         spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(25:32,25-n:24,:),-1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(25-n:24,25:32,:);
%         spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(25:32,33:32+n,:),1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(33:32+n,25:32,:),2);
% 
%         % Pad with pillar data on either end
%         o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
%         o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
%         spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
%         spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);
% 
%         % Save indices of original grid [from_x to_x; from_y to_y]
%         retrievemap = [1 size(o_i,1); ...
%                        n+1 n+size(o_i,2)];
%         % Send vars for adaptive smoothing
%         o_i = o_i_temp;
%         spikeLoc = spikeLoc_temp;
% 
%     case 'Pillar2'
%         groundsection_ind = strcmp(gazeSections,'Ground');
%         ground_o_i = grid_o_i{groundsection_ind};
%         ground_spikeLoc = grid_spikeLoc{groundsection_ind};
% 
%         % Move original map to middle
%         o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
%         spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;
% 
%         % Pad with ground data
%         o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(25:32,9-n:8,:),-1);
%         o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(25-n:24,9:16,:);
%         o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(25:32,17:16+n,:),1);
%         o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(33:32+n,9:16,:),2);
%         spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(25:32,9-n:8,:),-1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(25-n:24,9:16,:);
%         spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(25:32,17:16+n,:),1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(33:32+n,9:16,:),2);
% 
%         % Pad with pillar data on either end
%         o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
%         o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
%         spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
%         spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);
% 
%         % Save indices of original grid [from_x to_x; from_y to_y]
%         retrievemap = [1 size(o_i,1); ...
%                        n+1 n+size(o_i,2)];
%         % Send vars for adaptive smoothing
%         o_i = o_i_temp;
%         spikeLoc = spikeLoc_temp;
% 
%     case 'Pillar3'
%         groundsection_ind = strcmp(gazeSections,'Ground');
%         ground_o_i = grid_o_i{groundsection_ind};
%         ground_spikeLoc = grid_spikeLoc{groundsection_ind};
% 
%         % Move original map to middle
%         o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
%         spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;
% 
%         % Pad with ground data
%         o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(9:16,25-n:24,:),-1);
%         o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(9-n:8,25:32,:);
%         o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(9:16,33:32+n,:),1);
%         o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(17:16+n,25:32,:),2);
%         spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(9:16,25-n:24,:),-1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(9-n:8,25:32,:);
%         spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(9:16,33:32+n,:),1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(17:16+n,25:32,:),2);
% 
%         % Pad with pillar data on either end
%         o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
%         o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
%         spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
%         spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);
% 
%         % Save indices of original grid [from_x to_x; from_y to_y]
%         retrievemap = [1 size(o_i,1); ...
%                        n+1 n+size(o_i,2)];
%         % Send vars for adaptive smoothing
%         o_i = o_i_temp;
%         spikeLoc = spikeLoc_temp;
% 
%     case 'Pillar4'
%         groundsection_ind = strcmp(gazeSections,'Ground');
%         ground_o_i = grid_o_i{groundsection_ind};
%         ground_spikeLoc = grid_spikeLoc{groundsection_ind};
% 
%         % Move original map to middle
%         o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
%         spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
%         spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;
% 
%         % Pad with ground data
%         o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(9:16,9-n:8,:),-1);
%         o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(9-n:8,9:16,:);
%         o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(9:16,17:16+n,:),1);
%         o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(17:16+n,9:16,:),2);
%         spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(9:16,9-n:8,:),-1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(9-n:8,9:16,:);
%         spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(9:16,17:16+n,:),1);
%         spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(17:16+n,9:16,:),2);
% 
%         % Pad with pillar data on either end
%         o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
%         o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
%         spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
%         spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);
% 
%         % Save indices of original grid [from_x to_x; from_y to_y]
%         retrievemap = [1 size(o_i,1); ...
%                        n+1 n+size(o_i,2)];
%         % Send vars for adaptive smoothing
%         o_i = o_i_temp;
%         spikeLoc = spikeLoc_temp;
% 
% end

% 
% function [smoothedRate,smoothedSpk,smoothedPos,radiiUsedList] = adaptivesmooth(pos,spk,alpha)
% % Adapted from rates_adaptivesmooth.m (Wills et al)
% % pos = occupancy map/dwell time in each position bin (in seconds)
% % spk = spike map/spike count in each position bin
% % alpha = scaling parameter (1e6 for Skaggs et al 1996, 1e5 for Wills et al 2010)
% 
% % Check for empty spk maps %
% if sum(sum(spk))==0
%     smoothedPos=pos;    smoothedPos(pos==0)=nan;
%     smoothedSpk=spk;    smoothedSpk(pos==0)=nan;
%     smoothedRate=spk;   smoothedRate(pos==0)=nan;
%     radiiUsedList=nan(1,sum(sum(pos>0)));
%     return
% end
% % Pre-assign output %
% smoothedPos=zeros(size(pos));
% smoothedSpk=zeros(size(pos));
% % Visited env template: use this to get numbers of visited bins in filter at edge of environemnt %
% vis=zeros(size(pos));
% vis(pos>0)=1;
% % Pre-assign map which records which bins have passed %
% smoothedCheck=false(size(pos));
% smoothedCheck(pos==0)=true; % Disregard unvisited - mark as already done.
% % Pre-assign list of radii used (this is for reporting purposes, not used for making maps) %
% radiiUsedList=nan(1,sum(sum(pos>0)));
% radiiUsedCount=1;
% % These parameters depend on place or dir mode %
% if size(pos,2)>1
%     boundary=0;             % IMFILTER boundary condition
%     rBump=0.5;              % Increase radius in 0.5 bin steps.
% elseif size(pos,2)==1
%     boundary='circular';
%     rBump=1;                % Increase radius in 1 bin steps.
% end
% 
% %%% Run increasing radius iterations %%%
% r=1; % Circle radius
% while any(any(~smoothedCheck))
%     % Check radius isn't getting too big (if >map/2, stop running) %
%     if r>max(size(pos))/2
%         smoothedSpk(~smoothedCheck)=nan;
%         smoothedPos(~smoothedCheck)=nan;
%         break
%     end
%     % Construct filter kernel ...
%     if size(pos,2)>1
%         % Place: Flat disk, where r>=distance to bin centre %
%         f=fspecial('disk',r); 
%         f(f>=(max(max(f))/3))=1;
%         f(f~=1)=0;
%     elseif size(pos,2)==1 
%         % Direction: boxcar window, r bins from centre symmetrically %
%         f=ones(1+(r*2),1);
%     end     
%     % Filter maps (get N spikes and pos sum within kernel) %
%     fSpk=imfilter(spk,f,boundary);
%     fPos=imfilter(pos,f,boundary);
%     fVis=imfilter(vis,f,boundary);
%     % Which bins pass criteria at this radius? %
%     warning('off', 'MATLAB:divideByZero');
%     binsPassed=alpha./(sqrt(fSpk).*fPos) <= r;
%     warning('on', 'MATLAB:divideByZero');
%     binsPassed=binsPassed & ~smoothedCheck; % Only get the bins that have passed in this iteration.
%     % Add these to list of radii used %
%     nBins=sum(binsPassed(:));
%     radiiUsedList(radiiUsedCount:radiiUsedCount+nBins-1)=r;
%     radiiUsedCount=radiiUsedCount+nBins;
%     % Assign values to smoothed maps %
%     smoothedSpk(binsPassed)=fSpk(binsPassed)./fVis(binsPassed);
%     smoothedPos(binsPassed)=fPos(binsPassed)./fVis(binsPassed);
%     % Record which bins were smoothed this iteration %
%     smoothedCheck(binsPassed)=true;
%     % Increase circle radius (half-bin steps) %
%     r=r+rBump;
% end
% 
% % Assign Output %
% warning('off', 'MATLAB:divideByZero');
% smoothedRate=smoothedSpk./smoothedPos;
% warning('on', 'MATLAB:divideByZero');
% smoothedRate(pos==0)=nan;
% smoothedPos(pos==0)=nan;
% smoothedSpk(pos==0)=nan;
% % report radii sizes?


            
            %             %%%% NEW
%             map = cell(9,1);
%             for jj = 3:size(grid_o_i_Gaze,1)
%                 [map{jj},smoothedSpk,smoothedDur]=adsmooth(grid_o_i_Gaze{jj},grid_spikeBin_Gaze{jj},alpha);
%             end
%             map_unpad = cell(9,1);
%             maplin = nan(5122,1);
%             for jj = 3:size(grid_o_i_Gaze,1)
%                 map_unpad{jj} = map{jj}(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2));
%                 % Put grid map back into linear map
%                 set = reshape(flipud(rot90(map_unpad{jj})),size(map_unpad{jj},1)*size(map_unpad{jj},2),1);
%                 lin_inds = sum(binDepths(1:jj-1,1).*binDepths(1:jj-1,2))+1:sum(binDepths(1:jj,1).*binDepths(1:jj,2));
%                 maplin(lin_inds,:) = reshape(set,1,binDepths(jj,1)*binDepths(jj,2));
%             end
%             %%% END OF NEW %%%
