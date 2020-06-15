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
				'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4,'AdaptiveSmooth',1, 'FiltOption',1, 'ThresVel',0, 'UseAllTrials',1);
            
Args.flags = {'Auto','ArgsOnly','FRSIC','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps','NumShuffles','FiltOption','AdaptiveSmooth','ThresVel','UseAllTrials'};                           

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
    pv = vmpv('auto', varargin{:});
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

        if Args.ThresVel == 1
            if Args.FiltOption == 0
                conditions = conditions & (pv.data.thres_vel.place_good_rows > 0);
            elseif Args.FiltOption == 1
                conditions = conditions & (pv.data.thres_vel.place_good_rows > 0.5);
            else
                conditions = conditions & ismember(stc(:,2), pv.data.thres_vel.place_good_bins);
            end
        else
            if Args.FiltOption == 0
                conditions = conditions & (pv.data.all_vel.place_good_rows > 0);
            elseif Args.FiltOption == 1
                conditions = conditions & (pv.data.all_vel.place_good_rows > 0.5);
            else
                conditions = conditions & ismember(stc(:,2), pv.data.all_vel.place_good_bins);
            end    
        end


        disp('conditioning done');
        if repeat == 1
            dstc = diff(stc(:,1));
            stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
            stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
            stc_changing_ind = stc_changing_ind(1:end-1,:);
        end

        consol_arr = zeros(Args.GridSteps * Args.GridSteps,Args.NumShuffles + 1);

        interval = 1;

        for sp = 1:size(flat_spiketimes,1)

%             if rem(sp, 10000000) == 0
%                 disp(num2str(100*sp/size(flat_spiketimes,1)))
%             end

            while interval < size(stc_changing_ind,1)
                if flat_spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && flat_spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
                    break;
                end
                interval = interval + 1;
            end   

            bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),2);
            bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))));

            bins_hit(~(bins_hit>0)) = [];

            consol_arr(bins_hit,flat_spiketimes(sp,2)) = consol_arr(bins_hit,flat_spiketimes(sp,2)) + 1;

        end        
        
        firing_counts_full = consol_arr';
        stc_ss = stc(conditions,[2 4]);
        stc_ss(~(stc_ss(:,1) > 0),:) = [];
        stc_ss = [stc_ss; [1600 0]];
        gpdur = accumarray(stc_ss(:,1),stc_ss(:,2))';
        bins_sieved = 1:(Args.GridSteps * Args.GridSteps);

            if Args.AdaptiveSmooth

                firing_rates_full_raw = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);
                to_save = NaN(1,Args.GridSteps^2);
                to_save(bins_sieved) = firing_rates_full_raw(1,:);
                data.maps_raw = to_save;
                nan_track = isnan(to_save);

                alpha = 1e3;

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

                wip = ones(Args.NumShuffles+1,1);

                gpdur1 = zeros(1,Args.GridSteps^2);
                gpdur1(bins_sieved) = gpdur;

                preset_to_zeros = reshape(gpdur1, Args.GridSteps, Args.GridSteps); % will be set to nans afterwards, just swapped to zero to quickly cut 'done' shuffles
                preset_to_zeros(find(preset_to_zeros>0)) = 1;
                preset_to_zeros = ~preset_to_zeros;
                preset_to_zeros = repmat(preset_to_zeros, 1,1,Args.NumShuffles+1);

                gpdur1 = repmat(gpdur1,Args.NumShuffles + 1,1);
                gpdur1 = reshape(gpdur1, Args.NumShuffles + 1, Args.GridSteps,Args.GridSteps);
                gpdur1 = permute(gpdur1,[2,3,1]);

                firing_counts_full1 = zeros(Args.NumShuffles + 1, Args.GridSteps^2);
                firing_counts_full1(:,bins_sieved) = firing_counts_full;
                firing_counts_full1 = reshape(firing_counts_full1, Args.NumShuffles + 1, Args.GridSteps,Args.GridSteps);
                firing_counts_full1 = permute(firing_counts_full1,[2,3,1]);

                to_compute = 1:0.5:Args.GridSteps/2;
                possible = NaN(length(to_compute),2,Args.GridSteps,Args.GridSteps,Args.NumShuffles + 1);
                to_fill = NaN(size(possible,3), size(possible,4), size(possible,5));
                to_fill(preset_to_zeros) = 0;
                to_fill_time = NaN(size(possible,3), size(possible,4), size(possible,5));
                to_fill_time(preset_to_zeros) = 0; 
                to_fill_radius = NaN(size(possible,3), size(possible,4), size(possible,5));
                to_fill_radius(preset_to_zeros) = 0;

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

%                     if sum(sum(sum(isnan(to_fill(:,:,:))))) == 0
%                         disp('breaking');
%                         break;
%                     end

                end

                to_fill(isnan(to_fill)) = 0;
                to_fill = permute(to_fill, [3 1 2]);
                to_fill = reshape(to_fill, Args.NumShuffles + 1, Args.GridSteps^2);
                to_fill = to_fill(:,bins_sieved);

                to_fill_time(isnan(to_fill_time)) = 0;
                to_fill_time = permute(to_fill_time, [3 1 2]);
                to_fill_time = reshape(to_fill_time, Args.NumShuffles + 1, Args.GridSteps^2);
                to_fill_time = to_fill_time(:,bins_sieved);
                
                to_fill_radius = permute(to_fill_radius, [3 1 2]);
                to_fill_radius = reshape(to_fill_radius, Args.NumShuffles + 1, Args.GridSteps^2);

                firing_rates_full = to_fill;

                % smoothing part ends
                
                if repeat == 1
                    to_save = NaN(1,Args.GridSteps^2);
                    to_save(bins_sieved) = firing_rates_full(1,:);
                    to_save(nan_track) = nan;
                    data.maps_adsmooth = to_save;
                    to_save = NaN(size(firing_rates_full,1)-1,Args.GridSteps^2);
                    to_save(:,bins_sieved) = firing_rates_full(2:end,:);
                    data.maps_all = to_save;
                    to_save = NaN(size(firing_rates_full,1)-1,Args.GridSteps^2);
                    to_save(:,bins_sieved) = to_fill_time(2:end,:);
                    data.dur_map_all = to_save;
                    to_save = NaN(1,Args.GridSteps^2);
                    to_save(bins_sieved) = to_fill_time(1,:);
                    data.dur_map_actual = to_save;
                    data.radii = to_fill_radius;
                elseif repeat == 2
                    to_save = NaN(1,Args.GridSteps^2);
                    to_save(bins_sieved) = firing_rates_full(1,:);
                    to_save(nan_track) = nan;
                    data.maps_adsmooth1 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.GridSteps^2);
                    %to_save(:,bins_sieved) = firing_rates_full(2:end,:);
                    %data.maps_all1 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.GridSteps^2);
                    %to_save(:,bins_sieved) = to_fill_time(2:end,:);
                    %data.dur_map_all1 = to_save;
                    to_save = NaN(1,Args.GridSteps^2);
                    to_save(bins_sieved) = to_fill_time(1,:);
                    data.dur_map_actual1 = to_save;
                    data.radii1 = to_fill_radius;
                elseif repeat == 3
                    to_save = NaN(1,Args.GridSteps^2);
                    to_save(bins_sieved) = firing_rates_full(1,:);
                    to_save(nan_track) = nan;
                    data.maps_adsmooth2 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.GridSteps^2);
                    %to_save(:,bins_sieved) = firing_rates_full(2:end,:);
                    %data.maps_all2 = to_save;
                    %to_save = NaN(size(firing_rates_full,1)-1,Args.GridSteps^2);
                    %to_save(:,bins_sieved) = to_fill_time(2:end,:);
                    %data.dur_map_all2 = to_save;
                    to_save = NaN(1,Args.GridSteps^2);
                    to_save(bins_sieved) = to_fill_time(1,:);
                    data.dur_map_actual2 = to_save;
                    data.radii2 = to_fill_radius;
                end

            else
                firing_rates_full = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);
                to_fill_time = repmat(gpdur,size(firing_counts_full,1),1); % HM added

                to_save = NaN(1,Args.GridSteps^2);
                to_save(bins_sieved) = firing_rates_full(1,:);
                if repeat == 1
                    data.maps_raw = to_save;
                elseif repeat == 2
                    data.maps_raw1 = to_save;
                elseif repeat == 3
                    data.maps_raw2 = to_save;
                end
            end

                gpdur1 = zeros(Args.NumShuffles+1,Args.GridSteps^2);
                gpdur1(:,bins_sieved) = to_fill_time;
                Pi1 = gpdur1./sum(gpdur1,2); % consider nansum to play safe
        %         Pi1 = repmat(Pi1, Args.NumShuffles+1, 1);

                lambda_i = NaN(Args.NumShuffles+1,Args.GridSteps^2);
                lambda_i(:,bins_sieved) = firing_rates_full;
                lambda_i(isnan(lambda_i)) = 0;
                lambda_bar = sum(Pi1 .* lambda_i,2);
                % divide firing for each position by the overall mean
                FRratio = lambda_i./repmat(lambda_bar,1,Args.GridSteps^2);
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

        %     histogram(sic_out);

                % ISE part
                lambda_i = NaN(Args.NumShuffles+1,Args.GridSteps^2);
                lambda_i(:,bins_sieved) = firing_rates_full;                
                
                if repeat == 1
                    ise_out = ise(lambda_i(1,:), lambda_i(2:end,:), Args.GridSteps, Args.GridSteps);
                    data.ISE = ise_out(1);
                    data.ISEsh = ise_out(2:end);
                elseif repeat == 2
                    ise_out = ise(lambda_i, [], Args.GridSteps, Args.GridSteps);
                    data.ISE1 = ise_out;
                elseif repeat == 3
                    ise_out = ise(lambda_i, [], Args.GridSteps, Args.GridSteps);
                    data.ISE2 = ise_out;
                end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if repeat == 1
            data.SIC = sic_out(1);
            data.SICsh = sic_out';
        %     data.median_occ_firings = median_stats';
        %     data.variance_occ_firings = var_stats';
        %     data.perc_occ_firings = perc_stats';
        %     data.occ_data = occ_data;
        elseif repeat == 2
            data.SIC1 = sic_out;
        elseif repeat == 3
            data.SIC2 = sic_out;
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
