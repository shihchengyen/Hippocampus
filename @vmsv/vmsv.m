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
				'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4, 'ThresVel',1, 'UseMinObs', 1, 'AdaptiveSmooth', 1, 'UseAllTrials',1,'Alpha',1000);
            
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
    pv = vmpv('auto', varargin{:});
%     pv = vmpv('auto','save','MinObsView',5,'MinDurView',0.01);
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
        end
        
        % spike shuffling

        spiketimes = spiketrain.timestamps/1000; % now in seconds
        maxTime = pv.data.rplmaxtime;
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
            stc(:,4) = [diff(stc(:,1)); 0];
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


        disp('conditioning done');
        if repeat == 1
            dstc = diff(stc(:,1));
            stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
            stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
            stc_changing_ind = stc_changing_ind(1:end-1,:);
        end

        consol_arr = zeros(5122,Args.NumShuffles + 1);

        interval = 1;

        for sp = 1:size(flat_spiketimes,1)

%             if rem(sp, 10000000) == 0
%                 disp(100*sp/size(flat_spiketimes,1))
%             end

            while interval < size(stc_changing_ind,1)
                if flat_spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && flat_spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
                    break;
                end
                interval = interval + 1;
            end   

            bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),3);
            bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))));

            bins_hit(~(bins_hit>0)) = [];

            consol_arr(bins_hit,flat_spiketimes(sp,2)) = consol_arr(bins_hit,flat_spiketimes(sp,2)) + 1;

        end        

        spike_count = consol_arr;
        
        stc(stc(:,4)==0,4) = nan;
        stc(:,4) = fillmissing(stc(:,4),'next');
        stc_ss = stc(find(conditions==1),[3 4]);
        stc_ss(~(stc_ss(:,1) > 0),:) = [];
        stc_ss(isnan(stc_ss(:,2)),:) = [];
        stc_ss = [stc_ss; [5122 0]];
        
        gpdur1 = accumarray(stc_ss(:,1),stc_ss(:,2))';
        
        if Args.AdaptiveSmooth
            
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
            
            lin_o_i_Gaze = repmat(gpdur1', 1, Args.NumShuffles + 1);
            lin_spikeLoc_Gaze = spike_count;
            maps_raw = lin_spikeLoc_Gaze./lin_o_i_Gaze;
            if Args.NumShuffles > 0
                maps_raw = maps_raw(:,1);
            end
            maps_raw_out = nan(size(maps_raw));
            maps_raw_out(bins_sieved) = maps_raw(bins_sieved);
            if repeat == 1
                data.maps_raw = maps_raw_out';
            elseif repeat == 2
                data.maps_raw1 = maps_raw_out';
            else
                data.maps_raw2 = maps_raw_out';
            end
                
            grid_o_i_Gaze = cell(size(binDepths,1),1);
            grid_spikeBin_Gaze = grid_o_i_Gaze;
            grid_smoothed_Gaze = grid_o_i_Gaze;
            
            for jj = 1:size(binDepths,1) % for each grid
                % Initialise empty matrices
                o_i = nan(binDepths(jj,1),binDepths(jj,2),Args.NumShuffles+1);
                spikeBin = o_i;
                map = o_i;
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
                    o_i(x,y,:) = lin_o_i_Gaze(indbins_lin,:);
                    spikeBin(x,y,:) = lin_spikeLoc_Gaze(indbins_lin,:);
                    
                end
                % Collect output
                grid_o_i_Gaze{jj} = o_i;
                grid_spikeBin_Gaze{jj} = spikeBin;
            end
            
            
            retrievemap = cell(size(grid_o_i_Gaze,1),1);
            grid_o_i_Gaze_temp = grid_o_i_Gaze; % keep original for reference during padding
            grid_spikeBin_Gaze_temp = grid_spikeBin_Gaze; % keep original for reference during padding
            for jj = 1:size(grid_o_i_Gaze,1) % For each separate grid
                
                if binDepths(jj,1)*binDepths(jj,2) > 2 % For non-cue/non-hint grids
                    % Pad each grid map with adjoining bins from other grids
                    % Pad with <<5>> extra bin rows
                    n = 5;
                    [retrievemap{jj},grid_o_i_Gaze{jj},grid_spikeBin_Gaze{jj}] = padgrids(n,grid_o_i_Gaze{jj},grid_spikeBin_Gaze{jj},grid_o_i_Gaze,grid_spikeBin_Gaze,gazeSections,jj);
                    
                end
            end
            
            alpha = Args.Alpha;
            %         alpha = 1;
            grid_smoothed_Gaze{1} = grid_spikeBin_Gaze{1}./grid_o_i_Gaze{1};
            grid_smoothed_Gaze{2} = grid_spikeBin_Gaze{2}./grid_o_i_Gaze{2};
            grid_smoothed_dur = cell(size(grid_o_i_Gaze,1),1);
            grid_smoothed_dur{1} = grid_o_i_Gaze{1};
            grid_smoothed_dur{2} = grid_o_i_Gaze{2};
            grid_ad_size = cell(size(grid_o_i_Gaze,1),1);
            grid_ad_size{1} = NaN;
            grid_ad_size{2} = NaN;
            
            for jj = 3:size(grid_o_i_Gaze,1) % for each grid
                
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
                to_compute = 1:0.5:(max(size(grid_o_i_Gaze{jj})))/2;
                
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
                    disp(['smoothed grid ' num2str(jj) ' with kernel size ' num2str(to_compute(idx)) ', leaving ' num2str(remaining) ' grids undone']);
                    
                    check = squeeze(sum(sum(isnan(to_fill),2),1));
                    wip(check==0) = 0;
                    
                    if remaining == 0
                        disp('done');
                        break;
                    end
                end
                
                to_fill(isnan(to_fill)) = 0;
                to_fill = to_fill(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2),:);
                grid_smoothed_Gaze{jj} = to_fill;
                to_fill_smoothed_duration(isnan(to_fill_smoothed_duration)) = 0;
                to_fill_smoothed_duration = to_fill_smoothed_duration(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2),:);
                grid_smoothed_dur{jj} = to_fill_smoothed_duration;
                to_fill_size = to_fill_size(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2),:);
                grid_ad_size{jj} = to_fill_size;
                
            end
            
            
            total_grids = 0;
            for jj = 1:size(grid_smoothed_Gaze,1)
                total_grids = total_grids + size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2);
            end
            gpdur1 = zeros(Args.NumShuffles+1,total_grids);
            lambda_i = NaN(Args.NumShuffles+1,total_grids);
            radii = zeros(Args.NumShuffles+1,total_grids);
            filling_index = 0;
            for jj = 1:size(grid_o_i_Gaze,1)
                temp4 = reshape(permute(grid_smoothed_dur{jj},[2 1 3]), [size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2) Args.NumShuffles+1]);
                temp3 = reshape(permute(grid_smoothed_Gaze{jj},[2 1 3]), [size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2) Args.NumShuffles+1]);
                gpdur1(:,filling_index+1:filling_index+size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2)) = temp4';
                lambda_i(:,filling_index+1:filling_index+size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2)) = temp3';
                
                if jj > 2
                    temp5 = reshape(permute(grid_ad_size{jj},[2 1 3]), [size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2) Args.NumShuffles+1]);
                    radii(:,filling_index+1:filling_index+size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2)) = temp5';
                end
                filling_index = filling_index + size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2);
            end
            lambda_i(gpdur1==0) = nan;
            
            if repeat == 1
                data.radii = radii(1,:);
                data.radiish = radii(2:end,:);
                data.dur_adsm = gpdur1(1,:);
                data.dur_adsmsh = gpdur1(2:end,:);
            end            
            
        else
            
            gpdur1 = repmat(gpdur1, Args.NumShuffles + 1, 1);
            lambda_i = spike_count'./ gpdur1;
            
        end

  
            
            % SIC portion
            
            Pi1 = gpdur1./sum(gpdur1,2);
            lambda_bar = nansum(Pi1 .* lambda_i,2);
            % divide firing for each position by the overall mean
            FRratio = lambda_i./repmat(lambda_bar,1,5122);
            % compute first term in SIC
            SIC1 = Pi1 .* lambda_i; 
            SIC2 = log2(FRratio);
            zeros_placing = SIC1==0;  

            bits_per_sec = SIC1 .* SIC2 ./ lambda_bar;
            bits_per_sec(zeros_placing) = NaN;
            lambda_bar_ok = lambda_bar>0;
            lambda_bar_bad = ~lambda_bar_ok;
            sic_out = nansum(bits_per_sec, 2);
            sic_out(lambda_bar_bad) = NaN;

            tic;

            % ISE portion
            % create overall map and insert padded portions in, to account for
            % cross-portion pairs
            canvas = nan(51, 161, Args.NumShuffles + 1);
            firing_rates = lambda_i;

            % flooring
            floor_padded = nan(42,42,Args.NumShuffles+1);
            floor_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);
            floor_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,3203:3203+39),[2 1]), 40, 1, Args.NumShuffles+1),1);
            floor_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,3243:3243+39),[2 1]), 1, 40, Args.NumShuffles+1);
            floor_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,3283:3283+39),[2 1]), 40, 1, Args.NumShuffles+1);
            floor_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,3323:3323+39),[2 1]), 1, 40, Args.NumShuffles+1), 2);
            canvas(10:end,1:42,:) = floor_padded;

            % ceiling
            ceiling_padded = nan(42,42,Args.NumShuffles+1);
            ceiling_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,1603:3202),size(firing_rates,1),40,40), [3 2 1]), 1);
            ceiling_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,4323:4323+39),[2 1]), 40, 1, Args.NumShuffles+1),1);
            ceiling_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,4363:4363+39),[2 1]), 1, 40, Args.NumShuffles+1);
            ceiling_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,4403:4403+39),[2 1]), 40, 1, Args.NumShuffles+1);
            ceiling_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,4443:4443+39),[2 1]), 1, 40, Args.NumShuffles+1), 2);
            canvas(10:end,44:85,:) = ceiling_padded;

            % walls
            walls_padded = nan(8,161,Args.NumShuffles+1);
            walls_padded(:,1:end-1,:) = flip(permute(reshape(firing_rates(:,3203:3203+1280-1), Args.NumShuffles+1, 40*4, 8),[3 2 1]), 1);
            walls_padded(:,end,:) = walls_padded(:,1,:);
            canvas(1:8,:,:) = walls_padded;

            % used to pad pillar base more easily
            floor_base = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);

            % pillars
            PTL_padded = nan(6,33,Args.NumShuffles+1);
            PTL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4963:4963+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
            % small diagonal issue here, diagonal floor bins at the corners are put
            % side by side, only 16 such occurrences in total, neglected for now.
            PTL_padded(end,1:8,:) = flip(permute(floor_base(9:16,8,:),[2 1 3]),2);
            PTL_padded(end,9:16,:) = floor_base(8,9:16,:);
            PTL_padded(end,17:24,:) = permute(floor_base(9:16,17,:),[2 1 3]);
            PTL_padded(end,25:32,:) = flip(floor_base(17,9:16,:),2);
            PTL_padded(:,end,:) = PTL_padded(:,1,:);
            canvas(10:10+6-1,87:87+32,:) = PTL_padded;

            PTR_padded = nan(6,33,Args.NumShuffles+1);
            PTR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4803:4803+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
            PTR_padded(end,1:8,:) = flip(permute(floor_base(9:16,24,:),[2 1 3]),2);
            PTR_padded(end,9:16,:) = floor_base(8,25:32,:);
            PTR_padded(end,17:24,:) = permute(floor_base(9:16,33,:),[2 1 3]);
            PTR_padded(end,25:32,:) = flip(floor_base(17,25:32,:),2);
            PTR_padded(:,end,:) = PTR_padded(:,1,:);
            canvas(10:10+6-1,121:121+32,:) = PTR_padded;

            PBL_padded = nan(6,33,Args.NumShuffles+1);
            PBL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4643:4643+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
            PBL_padded(end,1:8,:) = flip(permute(floor_base(25:32,8,:),[2 1 3]),2);
            PBL_padded(end,9:16,:) = floor_base(24,9:16,:);
            PBL_padded(end,17:24,:) = permute(floor_base(25:32,17,:),[2 1 3]);
            PBL_padded(end,25:32,:) = flip(floor_base(33,9:16,:),2);
            PBL_padded(:,end,:) = PBL_padded(:,1,:);
            canvas(17:17+6-1,87:87+32,:) = PBL_padded;

            PBR_padded = nan(6,33,Args.NumShuffles+1);
            PBR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4483:4483+160-1), Args.NumShuffles+1, 8*4, 5),[3 2 1]), 1);
            PBR_padded(end,1:8,:) = flip(permute(floor_base(25:32,24,:),[2 1 3]),2);
            PBR_padded(end,9:16,:) = floor_base(24,25:32,:);
            PBR_padded(end,17:24,:) = permute(floor_base(25:32,33,:),[2 1 3]);
            PBR_padded(end,25:32,:) = flip(floor_base(33,25:32,:),2);
            PBR_padded(:,end,:) = PBR_padded(:,1,:);
            canvas(17:17+6-1,121:121+32,:) = PBR_padded;

            actual_image = canvas(:,:,1);
            actual_image = actual_image(:)';
            shuffled_images = canvas(:,:,2:end);
            shuffled_images = reshape(shuffled_images, size(shuffled_images,3),size(shuffled_images,1)*size(shuffled_images,2));
%for testing
            data.actual_image=actual_image;
            data.shuffled_images=shuffled_images(:);
%             save('/Users/DataDesktop')
            
            disp(['time taken to pad map for ISE: ' num2str(toc)]);
            tic;
% %     for testing
%             bin_resolution=[0.5; 0.05; 0.005; 0.0005; 0.1; 0.01; 0.001; 0.0001; 1; 2; 4; 5; 10; 15; 20];
% 
% for i=1:size(bin_resolution,1)
% ise_out = ise(actual_image, shuffled_images, 51, 161,bin_resolution(i))
% name=['ise' num2str(bin_resolution(i)) '.mat'];
% ISE = ise_out(1);
% ISEsh = ise_out(2:end,1);
% save(name,'ISE','ISEsh');
% end
%             ise_out = ise(actual_image, shuffled_images, 51, 161);
%             ise_out = ise_ted(actual_image, shuffled_images, 51, 161);
%             ise_out = ise_matlab(actual_image, shuffled_images, 51, 161);
            ise_out = ise_ted_2(actual_image, shuffled_images, 51, 161);
            disp(['time taken to compute ISE: ' num2str(toc)]);
            
            if repeat == 1
                if ~Args.AdaptiveSmooth
                    if Args.NumShuffles > 0
                        maps_raw = firing_rates(1,:);
                    else
                        maps_raw = firing_rates;
                    end
                    maps_raw_out = nan(size(maps_raw));
                    maps_raw_out(bins_sieved) = maps_raw(bins_sieved);
                    data.maps_raw = maps_raw_out;
                else
                    if Args.NumShuffles > 0
                        maps_raw = firing_rates(1,:);
                    else
                        maps_raw = firing_rates;
                    end
                    maps_raw_out = nan(size(maps_raw));
                    maps_raw_out(bins_sieved) = maps_raw(bins_sieved);
                    data.maps_adsm = maps_raw_out;                    
                end
                maps_all = firing_rates(2:end,:);
                maps_all_out = nan(size(maps_all));
                maps_all_out(:,bins_sieved) = maps_all(:,bins_sieved);
                data.maps_adsmsh = maps_all_out;
                
                data.flattened = squeeze(canvas(:,:,1));
                data.SIC = sic_out(1);
                data.SICsh = sic_out(2:end,1);
                data.ISE = ise_out(1);
                data.ISEsh = ise_out(2:end,1);
            elseif repeat == 2
                if ~Args.AdaptiveSmooth
                    if Args.NumShuffles > 0
                        maps_raw = firing_rates(1,:);
                    else
                        maps_raw = firing_rates;
                    end
                    maps_raw_out = nan(size(maps_raw));
                    maps_raw_out(bins_sieved) = maps_raw(bins_sieved);
                    data.maps_raw1 = maps_raw_out;
                else
                    if Args.NumShuffles > 0
                        maps_raw = firing_rates(1,:);
                    else
                        maps_raw = firing_rates;
                    end
                    maps_raw_out = nan(size(maps_raw));
                    maps_raw_out(bins_sieved) = maps_raw(bins_sieved);
                    data.maps_adsm1 = maps_raw_out;                    
                end
                data.SIC1 = sic_out;
                data.ISE1 = ise_out;
            elseif repeat == 3
                if ~Args.AdaptiveSmooth
                    if Args.NumShuffles > 0
                        maps_raw = firing_rates(1,:);
                    else
                        maps_raw = firing_rates;
                    end
                    maps_raw_out = nan(size(maps_raw));
                    maps_raw_out(bins_sieved) = maps_raw(bins_sieved);
                    data.maps_raw2 = maps_raw_out;
                else
                    if Args.NumShuffles > 0
                        maps_raw = firing_rates(1,:);
                    else
                        maps_raw = firing_rates;
                    end
                    maps_raw_out = nan(size(maps_raw));
                    maps_raw_out(bins_sieved) = maps_raw(bins_sieved);
                    data.maps_adsm2 = maps_raw_out;                    
                end
                data.SIC2 = sic_out;
                data.ISE2 = ise_out;
            end            
            
    end
    
    % create nptdata so we can inherit from it 
    Args.NumShuffles = NumShuffles_saved;
    data.gridSteps = Args.GridSteps;
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

% relic code
% 
function [retrievemap,o_i,spikeLoc,map] = padgrids(n,o_i,spikeLoc,grid_o_i,grid_spikeLoc,gazeSections,jj)

% Pad maps with adjoining bins from adjacent maps

switch gazeSections{jj}
    case 'Ground'
        wallsection_ind = strcmp(gazeSections,'Walls');
        wall_o_i = grid_o_i{wallsection_ind};
        wall_spikeLoc = grid_spikeLoc{wallsection_ind};

        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
        o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
        spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;

        % Pad with wall data
        o_i_temp(1:n,n+1:n+size(o_i,1),:) = wall_o_i(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:); % top
        o_i_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:),-1); % right
        o_i_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:),-2); % bottom
        o_i_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(wall_o_i(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:),1); % left
        spikeLoc_temp(1:n,n+1:n+size(o_i,1),:) = wall_spikeLoc(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:); % top
        spikeLoc_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:),-1); % right
        spikeLoc_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:),-2); % bottom
        spikeLoc_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(wall_spikeLoc(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:),1); % left

        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [n+1 n+size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;

    case 'Ceiling'
        wallsection_ind = strcmp(gazeSections,'Walls');
        wall_o_i = grid_o_i{wallsection_ind};
        wall_spikeLoc = grid_spikeLoc{wallsection_ind};

        % Flip walldata upside down
        wall_o_i = flipud(wall_o_i);
        wall_spikeLoc = flipud(wall_spikeLoc);

        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
        o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
        spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;

        % Pad with wall data
        o_i_temp(1:n,n+1:n+size(o_i,1),:) = fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:)); % top
        o_i_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:)),-1); % right
        o_i_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:)),-2); % bottom
        o_i_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(fliplr(wall_o_i(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:)),1); % left
        spikeLoc_temp(1:n,n+1:n+size(o_i,1),:) = fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,1*size(o_i,1)+1:2*size(o_i,1),:)); % top
        spikeLoc_temp(n+1:n+size(o_i,1),size(o_i,1)+n+1:end,:) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,2*size(o_i,1)+1:3*size(o_i,1),:)),-1); % right
        spikeLoc_temp(size(o_i,1)+n+1:end,n+1:size(o_i,1)+n,:) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,3*size(o_i,1)+1:4*size(o_i,1),:)),-2); % bottom
        spikeLoc_temp(n+1:size(o_i,1)+n,1:n,:) = rot90(fliplr(wall_spikeLoc(size(wall_o_i,1)-n+1:end,0*size(o_i,1)+1:1*size(o_i,1),:)),1); % left

        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [n+1 n+size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;

    case 'Walls'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};

        ceilingsection_ind = strcmp(gazeSections,'Ceiling');
        ceiling_o_i = grid_o_i{ceilingsection_ind};
        ceiling_spikeLoc = grid_spikeLoc{ceilingsection_ind};

        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
        o_i_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+2*n,size(o_i,2)+2*n,size(o_i,3));
        spikeLoc_temp(n+1:n+size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;

        % Pad with ground data
        o_i_temp(n+size(o_i,1)+1:end,n+1:size(ground_o_i,2)+n,:) = rot90(ground_o_i(:,1:n,:),-1);
        o_i_temp(n+size(o_i,1)+1:end,n+size(ground_o_i,2)+1:n+2*size(ground_o_i,2),:) = ground_o_i(1:n,:,:);
        o_i_temp(n+size(o_i,1)+1:end,n+2*size(ground_o_i,2)+1:n+3*size(ground_o_i,2),:) = rot90(ground_o_i(:,size(ground_o_i,1)-n+1:end,:),1);
        o_i_temp(n+size(o_i,1)+1:end,n+3*size(ground_o_i,1)+1:n+4*size(ground_o_i,1),:) = rot90(ground_o_i(size(ground_o_i,1)-n+1:end,:,:),2);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+1:size(ground_o_i,2)+n,:) = rot90(ground_spikeLoc(:,1:n,:),-1);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+size(ground_o_i,2)+1:n+2*size(ground_o_i,2),:) = ground_spikeLoc(1:n,:,:);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+2*size(ground_o_i,2)+1:n+3*size(ground_o_i,2),:) = rot90(ground_spikeLoc(:,size(ground_spikeLoc,1)-n+1:end,:),1);
        spikeLoc_temp(n+size(o_i,1)+1:end,n+3*size(ground_o_i,1)+1:n+4*size(ground_o_i,1),:) = rot90(ground_spikeLoc(size(ground_spikeLoc,1)-n+1:end,:,:),2);

        % Pad with ceiling data
        o_i_temp(1:n,n+1:size(ceiling_o_i,1)+n,:) = fliplr(rot90(ceiling_o_i(:,size(ceiling_o_i,1)-n+1:end,:),1));
        o_i_temp(1:n,n+size(ceiling_o_i,1)+1:n+2*size(ceiling_o_i,1),:) = fliplr(ceiling_o_i(1:n,:,:));
        o_i_temp(1:n,n+2*size(ceiling_o_i,1)+1:n+3*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_o_i(:,1:n,:),-1));
        o_i_temp(1:n,n+3*size(ceiling_o_i,1)+1:n+4*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_o_i(size(ceiling_o_i,1)-n+1:end,:,:),2));
        spikeLoc_temp(1:n,n+1:size(ceiling_o_i,1)+n,:) = fliplr(rot90(ceiling_spikeLoc(:,size(ceiling_spikeLoc,1)-n+1:end,:),1));
        spikeLoc_temp(1:n,n+size(ceiling_o_i,1)+1:n+2*size(ceiling_o_i,1),:) = fliplr(ceiling_spikeLoc(1:n,:,:));
        spikeLoc_temp(1:n,n+2*size(ceiling_o_i,1)+1:n+3*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_spikeLoc(:,1:n,:),-1));
        spikeLoc_temp(1:n,n+3*size(ceiling_o_i,1)+1:n+4*size(ceiling_o_i,1),:) = fliplr(rot90(ceiling_spikeLoc(size(ceiling_spikeLoc,1)-n+1:end,:,:),2));

        % Pad with wall data on either end
        o_i_temp(n+1:n+size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
        o_i_temp(n+1:n+size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
        spikeLoc_temp(n+1:n+size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
        spikeLoc_temp(n+1:n+size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);

        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [n+1 n+size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;

    case 'Pillar1'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};

        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;

        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(25:32,25-n:24,:),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(25-n:24,25:32,:);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(25:32,33:32+n,:),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(33:32+n,25:32,:),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(25:32,25-n:24,:),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(25-n:24,25:32,:);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(25:32,33:32+n,:),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(33:32+n,25:32,:),2);

        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
        spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);

        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;

    case 'Pillar2'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};

        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;

        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(25:32,9-n:8,:),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(25-n:24,9:16,:);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(25:32,17:16+n,:),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(33:32+n,9:16,:),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(25:32,9-n:8,:),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(25-n:24,9:16,:);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(25:32,17:16+n,:),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(33:32+n,9:16,:),2);

        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
        spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);

        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;

    case 'Pillar3'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};

        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;

        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(9:16,25-n:24,:),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(9-n:8,25:32,:);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(9:16,33:32+n,:),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(17:16+n,25:32,:),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(9:16,25-n:24,:),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(9-n:8,25:32,:);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(9:16,33:32+n,:),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(17:16+n,25:32,:),2);

        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
        spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);

        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;

    case 'Pillar4'
        groundsection_ind = strcmp(gazeSections,'Ground');
        ground_o_i = grid_o_i{groundsection_ind};
        ground_spikeLoc = grid_spikeLoc{groundsection_ind};

        % Move original map to middle
        o_i_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        o_i_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = o_i;
        spikeLoc_temp = nan(size(o_i,1)+n,size(o_i,2)+2*n,size(o_i,3));
        spikeLoc_temp(1:size(o_i,1), n+1:n+size(o_i,2),:) = spikeLoc;

        % Pad with ground data
        o_i_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_o_i(9:16,9-n:8,:),-1);
        o_i_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_o_i(9-n:8,9:16,:);
        o_i_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_o_i(9:16,17:16+n,:),1);
        o_i_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_o_i(17:16+n,9:16,:),2);
        spikeLoc_temp(size(o_i,1)+1:end,n+1:(size(o_i,2)/4)+n,:) = rot90(ground_spikeLoc(9:16,9-n:8,:),-1);
        spikeLoc_temp(size(o_i,1)+1:end,n+(size(o_i,2)/4)+1:n+2*(size(o_i,2)/4),:) = ground_spikeLoc(9-n:8,9:16,:);
        spikeLoc_temp(size(o_i,1)+1:end,n+2*(size(o_i,2)/4)+1:n+3*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(9:16,17:16+n,:),1);
        spikeLoc_temp(size(o_i,1)+1:end,n+3*(size(o_i,2)/4)+1:n+4*(size(o_i,2)/4),:) = rot90(ground_spikeLoc(17:16+n,9:16,:),2);

        % Pad with pillar data on either end
        o_i_temp(1:size(o_i,1),1:n,:) = o_i(:,size(o_i,2)-n+1:end,:);
        o_i_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = o_i(:,1:n,:);
        spikeLoc_temp(1:size(o_i,1),1:n,:) = spikeLoc(:,size(o_i,2)-n+1:end,:);
        spikeLoc_temp(1:size(o_i,1),size(o_i_temp,2)-n+1:end,:) = spikeLoc(:,1:n,:);

        % Save indices of original grid [from_x to_x; from_y to_y]
        retrievemap = [1 size(o_i,1); ...
                       n+1 n+size(o_i,2)];
        % Send vars for adaptive smoothing
        o_i = o_i_temp;
        spikeLoc = spikeLoc_temp;

end

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
% 
% ISE=sv.data.ISE
% ISE_sh=sv.data.ISEsh;
% ISE_thr=prctile([ISE;sv.data.ISEsh],95);
% actual_image=sv.data.actual_image;
% shuffled_images=sv.data.shuffled_images;
% save('/Users/yuhsuan//Desktop/test_ise.mat');