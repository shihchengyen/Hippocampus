function placeandview_shuffle(total, chunk)

        poolobj = parpool(6);
        poolsize = poolobj.NumWorkers;
        disp(['number of workers: ' num2str(poolsize)]);
        
    % defining constants
    
        view_bin_count = 5122;
        shuffle_chunk = 50;
        shuffle_chunk = chunk;

    % temp args
    
        Args.ShuffleLimits = [0.1 0.9];
        Args.NumShuffles = 100;
        Args.NumShuffles = total;
    
    % loading required files
    
        pillars = load('/Volumes/Hippocampus/Data/picasso-misc/pillars.mat');
%         pillars = load('pillars.mat');
        pillars = pillars.pillar_bins;
        rp = rplparallel('auto');
        rp = rp.data;
        vms = vmsv('auto');
        vms = vms.data;
        um = umaze('auto');
        um = um.data;
        st = load('spiketrain.mat');
        st = st.timestamps;
        st = st./1000;
        
        % TEMP TEMP TEMP FOR OLD VMSV
%         a = load('bindepths.mat');
%         vms.binDepths = a.binDepths;

    % shuffling spike train
    
        maxTime = rp.timeStamps(end,3);
        tShifts = [0 ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
        full_arr = repmat(st, Args.NumShuffles+1, 1);
        full_arr = full_arr + tShifts';
        keepers = length(st) - sum(full_arr>maxTime, 2);
        for row = 2:size(full_arr,1)
            full_arr(row,:) = [full_arr(row,1+keepers(row):end)-maxTime-1 full_arr(row,1:keepers(row))];
        end    
    
    % getting place and view sessiontime (tracking changes in unique place- 
    % view pairs across time

        temp_view = vms.sessionTime_generated(:,1:2);
        temp_view = [temp_view zeros(size(temp_view,1),1) ones(size(temp_view,1),1)];
        temp_place = um.sessionTime(:,1:2);
        temp_place = [temp_place(:,1) zeros(size(temp_place,1),1) temp_place(:,2) 2.*ones(size(temp_place,1),1)];
        temp_comb = [temp_view; temp_place];
        temp_comb = sortrows(temp_comb, 1);

        for row = 2:size(temp_comb,1)
            if temp_comb(row,4) == 1
                temp_comb(row,3) = temp_comb(row-1,3);
            else
                temp_comb(row,2) = temp_comb(row-1,2);
            end
        end
    
    % binning spike trains into combined sessiontime to find view and place
    % at every spike

        flat_spiketimes = NaN(2,size(full_arr,1)*size(full_arr,2));
        temp = full_arr';
        flat_spiketimes(1,:) = temp(:);
        
        clear temp;
        
        flat_spiketimes(2,:) = repelem(1:size(full_arr,1), size(full_arr,2));
        edge_end = 0.5+size(full_arr,1);
        [N,~,~] = histcounts2(flat_spiketimes(1,:), flat_spiketimes(2,:), temp_comb(:,1), 0.5:1:edge_end);
        
        clear full_arr;
        clear flat_spiketimes;
        
        
        temp_comb(1:end-1,4) = temp_comb(2:end,1) - temp_comb(1:end-1,1);
        N = [N;zeros(1,size(N,2))];

        to_remove = (temp_comb(:,2)==0) | (temp_comb(:,3)==0) | (isnan(temp_comb(:,2))) | (isnan(temp_comb(:,3)));    
        N(to_remove,:) = [];
        temp_comb(to_remove,:) = [];

        spikes_actual_temp = [temp_comb(:,[2 3]) N(:,1)];
        
        occurs = [spikes_actual_temp; [5122 1600 0]];
        occurs = accumarray(occurs(:,1:2), occurs(:,3)>0);
        data.occurs = occurs;
        
        spikes_temp = [temp_comb(:,[2 3]) N(:,2:end)];
        dur_temp = temp_comb(:,2:4);
        dur_temp = [dur_temp; [view_bin_count 1600 0]]; % to get properly size array later
        dur_shuffle = accumarray(dur_temp(:,1:2),dur_temp(:,3),[]);
    
        size_n_2 = size(N,2);
        clear N;
        
    % to speed the griddifying up later on, we find the linear-grid 
    % mapping by applying slow function to bin indices
    % this is followed up with another linear-grid mapping that includes
    % padding
    
        index_tracker = 1:view_bin_count;

        % Restructure bins from linear to separate grids
        lin_grid_reference = cell(size(vms.binDepths,1),1);

        for jj = 1:size(vms.binDepths,1) % for each grid
            % Initialise empty matrices
            gridded_index = nan(vms.binDepths(jj,1),vms.binDepths(jj,2));

            % Assign linear bin to grid bin
            for mm = 1:vms.binDepths(jj,1)*vms.binDepths(jj,2) % For every point in linear map
                if mod(mm,vms.binDepths(jj,2)) == 0
                    y = vms.binDepths(jj,2);
                else
                    y = mod(mm,vms.binDepths(jj,2));
                end
                x = ceil(mm/vms.binDepths(jj,2));
                indbins_lin = mm + sum(vms.binDepths(1:jj-1,1).*vms.binDepths(1:jj-1,2));
                % Assign
                gridded_index(x,y) = index_tracker(indbins_lin);

            end
            % Collect output
            lin_grid_reference{jj} = gridded_index;
        end        
        
        padded_grid_reference = cell(size(vms.binDepths,1),1);
        cropper = cell(size(vms.binDepths,1),1);
        
        gazeSections = {'Cue' 'Hint' 'Ground' 'Ceiling' 'Walls' 'Pillar1' 'Pillar2' 'Pillar3' 'Pillar4'};
        for jj = 1:size(vms.binDepths,1) % for each grid
            if jj > 2 % exclude first 2 grids
                n = 5;
                [cropper{jj},~,padded_grid_reference{jj}] = padgrids(n,lin_grid_reference{jj},lin_grid_reference{jj},lin_grid_reference,lin_grid_reference,gazeSections,jj);
            end
        end        
        
    % initialize view x place table x subset of shuffles (in 500 chunks) -
    % view_bin_count x 1600 x 500 approx 30 gb, should reduce if running on lower ram
    % machine, followed by smoothing.
    
        spike_rate_shuffle = cell(1,ceil((size_n_2-1)/shuffle_chunk)); % dummy variable, no longer in use
        spike_rate_actual = accumarray([spikes_actual_temp(:,1:2);[view_bin_count 1600]],[spikes_actual_temp(:,3); 0],[])./dur_shuffle;
%         save('actual_sr_raw.mat','spike_rate_actual');
        data.sr_raw = spike_rate_actual;
        sr_actual_smoothed = smoothing(spike_rate_actual, lin_grid_reference, padded_grid_reference, cropper, pillars);        
%         save('actual_sr_smoothed.mat','sr_actual_smoothed');
        data.sr_smooth = sr_actual_smoothed;
        sr_actual_smoothed = reshape(sr_actual_smoothed, size(sr_actual_smoothed,1)*size(sr_actual_smoothed,2), 1);        
        
%         save('actual_time_spent_raw.mat','dur_shuffle');
        data.duration_raw = dur_shuffle;
        time_spent = smoothing(dur_shuffle, lin_grid_reference, padded_grid_reference, cropper, pillars);
        time_spent_smoothed = time_spent;
%         save('actual_time_spent.mat','time_spent_smoothed');
        data.duration_smooth = time_spent_smoothed;
        time_spent = reshape(time_spent,size(time_spent,1)*size(time_spent,2),1);
%         mean_sic = nan(length(spike_rate_shuffle),1);
        actual_sic = sic_batch(sr_actual_smoothed, time_spent);
        data.actual_sic = actual_sic;
%         save('actual_sic.mat','actual_sic');
        
%             spike_rate_shuffle_chunk = smoothing(spike_rate_shuffle_chunk, lin_grid_reference, padded_grid_reference, cropper, pillars);
%             spike_rate_shuffle_chunk = reshape(spike_rate_shuffle_chunk, size(spike_rate_shuffle_chunk,1)*size(spike_rate_shuffle_chunk,2), size(spike_rate_shuffle_chunk,3));
%             sic_chunk = sic_batch(spike_rate_shuffle_chunk, time_spent);        
        spikes_temp = sparse(spikes_temp);
        sic_shuffles = cell(1,length(spike_rate_shuffle));
        
        parfor chunk = 1:length(spike_rate_shuffle)
            columns_right = min([Args.NumShuffles (shuffle_chunk*chunk)])  +2;
            columns_left = shuffle_chunk*(chunk-1)+1  +2;
            disp([columns_left-2 columns_right-2]);
%             spikes_temp = full(spikes_temp);
            to_accum = [[repmat(spikes_temp(:,1:2),columns_right-columns_left+1,1) repelem((1:columns_right-columns_left+1)',size(spikes_temp,1),1) reshape(spikes_temp(:,columns_left:columns_right),[],1)]; [view_bin_count 1600 1 0]];
%             spikes_temp = sparse(spikes_temp);
            to_accum = full(to_accum);
            spike_rate_shuffle_chunk = accumarray(to_accum(:,1:3),to_accum(:,4),[])./repmat(dur_shuffle,1,1,columns_right-columns_left+1);
            to_accum = sparse(to_accum);
            % griddify, pad each grid, smooth, transform back to view_bin_countx1600xn
            spike_rate_shuffle_chunk = smoothing(spike_rate_shuffle_chunk, lin_grid_reference, padded_grid_reference, cropper, pillars);
            spike_rate_shuffle_chunk = reshape(spike_rate_shuffle_chunk, size(spike_rate_shuffle_chunk,1)*size(spike_rate_shuffle_chunk,2), size(spike_rate_shuffle_chunk,3));
            sic_chunk = sic_batch(spike_rate_shuffle_chunk, time_spent);
            sic_shuffles{chunk} = sic_chunk;
%             if chunk == 1
%                 sic_shuffles = sic_chunk';
%             else
%                 sic_shuffles = [sic_shuffles; sic_chunk'];
%             end
%             mean_sic(chunk) = mean(sic_shuffles);
%             clear spike_rate_shuffle_chunk;
%             save('temp.mat','sic_shuffles','mean_sic');
        end    
        data.mean_sic_shift = mean_sic;
        data.sic_shuffles = sic_shuffles;
        save('pnv_shuffle.mat', '-struct', 'data');

end

function [sic_array] = sic_batch(spike_rate_shuffle_chunk, time_spent)

    % spike_rate_shuffle_chunk is view_bin_countx1600 * shuffle_number, 2D
    % time_spent is view_bin_countx1600, 1D
    
    % converting duration to ratio
    time_spent = time_spent/sum(time_spent);
    time_spent = time_spent';
    
    % computing mean firing rate per shuffle
    average_per_shuffle = mean(spike_rate_shuffle_chunk,1);
    
    % removing bins with no firing rate
    
    second_half = spike_rate_shuffle_chunk.*log2(spike_rate_shuffle_chunk./repmat(average_per_shuffle,length(time_spent),1));
    second_half(isnan(second_half)) = 0;
    sic_array = time_spent*second_half;
    sic_array = sic_array./average_per_shuffle;

end

function [output_array] = smoothing(combined_array, lin_grid_reference, padded_grid_reference, cropper, pillars)

    output_array = nan(size(combined_array));
    output_array(1:2,:,:) = combined_array(1:2,:,:);
    
    % combined_array is view x place x shuffle

    % Restructure bins from linear view linear space to gridded view linear
    % space (view x space x shuffle) -> (padded-view x padded-view x space x shuffle)
    
        padded_array = cell(size(padded_grid_reference));
        combined_array = cat(1,combined_array, NaN(1,size(combined_array,2),size(combined_array,3))); % to be used 5-10 lines down

        for grid_ind = 1:length(padded_grid_reference)
            if grid_ind > 2 % exclude first 2 sections
                grid_ref = padded_grid_reference{grid_ind};
                grid_ref = grid_ref(:); % flatten to 1d
                spare = 5123; % to replace nan values temporarily
                grid_ref(isnan(grid_ref)) = spare;
                indexed_in = combined_array(grid_ref,:,:);
                padded_array{grid_ind} = reshape(indexed_in,size(padded_grid_reference{grid_ind},1),size(padded_grid_reference{grid_ind},2),size(indexed_in,2),size(combined_array,3));
            end
        end
        combined_array = [];
        
        
    % Restructure from (padded-view x padded-view x space x shuffle) to 
    % (padded-view x padded-view x space x space x shuffle)
    
        for grid_ind = 1:length(padded_grid_reference)
            if grid_ind > 2 % exclude first 2 sections
                padded_array{grid_ind} = reshape(padded_array{grid_ind}, size(padded_array{grid_ind},1), size(padded_array{grid_ind},2), size(lin_grid_reference{3},1), size(lin_grid_reference{3},1), size(padded_array{grid_ind},4));
            end
        end
    
    % Smoothing, then reshaping back into original (view x space x shuffle)
    
        disp('starting smooth process');
        for grid_ind = 1:length(padded_grid_reference)
            if grid_ind > 2 % exclude first 2 sections
                disp(grid_ind);
                disp('smoothing');
%                 temp = padded_array{grid_ind};
%                 padded_array{grid_ind} = [];
                padded_array{grid_ind}(isnan(padded_array{grid_ind})) = 0;
                % smoothing
                padded_array{grid_ind} = convn(padded_array{grid_ind},ones(5,5,5,5,1)./25,'same');
                disp('smoothing end');          
                
                % removing impossible smoothing (cutting into pillars)
                pillars_full = repmat(pillars, 1, 1, size(padded_array{grid_ind},1), size(padded_array{grid_ind},2), size(padded_array{grid_ind},5));
                pillars_full = permute(pillars_full, [3 4 1 2 5]);
                padded_array{grid_ind}(pillars_full) = 0;
                pillars_full = [];
                
                % removing excess bins from view dimensions
                crop = cropper{grid_ind};
                padded_array{grid_ind} = padded_array{grid_ind}(crop(1,1):crop(1,2),crop(2,1):crop(2,2),:,:,:);
                padded_array{grid_ind} = permute(padded_array{grid_ind}, [1 2 5 3 4]);
                padded_array{grid_ind} = reshape(padded_array{grid_ind}, size(padded_array{grid_ind},1), size(padded_array{grid_ind},2), size(padded_array{grid_ind},3), size(padded_array{grid_ind},4)*size(padded_array{grid_ind},5));
                padded_array{grid_ind} = permute(padded_array{grid_ind}, [1 2 4 3]);
                % squeezing view dimensions to one column, also do so for
                % reference indices
                padded_array{grid_ind} = reshape(padded_array{grid_ind}, size(padded_array{grid_ind},1)*size(padded_array{grid_ind},2), size(padded_array{grid_ind},3), size(padded_array{grid_ind},4));
                lgr = lin_grid_reference{grid_ind};
                lgr = lgr(:);
                % slotting everything back into the 3d array for output                
                output_array(lgr,:,:) = padded_array{grid_ind};
            end
        end    

    output_array(isnan(output_array)) = 0;
        
end


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

end



















