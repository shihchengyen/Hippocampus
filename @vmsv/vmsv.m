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
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',5000, ...
                'FRSIC',0, 'UseMedian',0, 'MinObs',5, ...
                'NumFRBins',4,'AdaptiveSmooth',1, 'FiltLowOcc',1);
            
Args.flags = {'Auto','ArgsOnly','FRSIC','UseAllTrials','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps','NumShuffles','FiltLowOcc','AdaptiveSmooth'};                           

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
     
    gaz1 = load('temp_gaze.mat');
    gaz.data.binGazeLin = gaz1.a1;
    gaz.data.timestamps = gaz1.a2;
    gaz.data.gpDurGaze = gaz1.a3;
    gaz.data.sessionTimeGaze = gaz1.a4;
    gaz2 = load('gaze1.mat');
    gaz.data.binGazeGrid = gaz2.a5;
    gaz.data.binDepths = gaz2.a6;
    gaz.data.binLocLin = gaz2.a7;
    gaz.data.binLocGrid = gaz2.a8;
    gaz3 = load('gaze2.mat');
    gaz.data.gazeSections = gaz3.a9;
    
    
	rp = load('../../../rplparallel.mat');
    rp = rp.rp.data;
%     rp = rplparallel('auto');
%     rp = rp.data;
    spiketrain = load(Args.RequiredFile);
%     gaz = gaze('auto',varargin{:});  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spiketimes = spiketrain.timestamps/1000; % now in seconds
    maxTime = rp.timeStamps(end,3);
    tShifts = [0 ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
    full_arr = repmat(spiketimes, Args.NumShuffles+1, 1);
    full_arr = full_arr + tShifts';
    keepers = length(spiketimes) - sum(full_arr>maxTime, 2);
    for row = 2:size(full_arr,1)
        full_arr(row,:) = [full_arr(row,1+keepers(row):end)-maxTime-1 full_arr(row,1:keepers(row))];
    end
    
    % calculating proportion of occupied time in each grid position across
    % entire session. trial based here, compared to occurrence based in
    % vmpc (should be standardized, needs time).
    
        diff1 = diff(gaz.data.binGazeLin);
        size(diff1)
        change_index = find(diff1~=0);
        size(change_index)

        testing = NaN(length(change_index),2);
        testing(:,2) = gaz.data.binGazeLin(change_index+1);
        testing(:,1) = gaz.data.timestamps(change_index+1);

        searching_for_first_trial = find(testing>rp.timeStamps(1,2));
        disp('test');
        testing = [0 0; rp.timeStamps(1,2) testing(searching_for_first_trial(1),2); testing(searching_for_first_trial(1):end,:)];

        for i = 1:size(rp.timeStamps,1)-1

            to_remove = testing(:,1) > rp.timeStamps(i,3);
            to_remove2 = testing(:,1) < rp.timeStamps(i+1,2);
            to_remove = and(to_remove, to_remove2);
            to_remove = find(to_remove);

            if ~isempty(to_remove)

                end_trial = [rp.timeStamps(i,3) 0];
                start_trial = [rp.timeStamps(i+1,2) testing(to_remove(length(to_remove)),2)];
                testing = [testing(1:to_remove(1)-1,:); end_trial; start_trial; testing(to_remove(length(to_remove))+1:end,:)];
            end
        end
        
        testing2 = zeros(size(testing,1),3);
        testing2(:,1:2) = testing;
        testing2(1:end-1,3) = diff(testing(:,1));
        
        gaz.data.sessionTimeGaze = testing2;
        zero_indices = find(gaz.data.sessionTimeGaze(:,2)==0);
    
    
    if Args.FiltLowOcc
        bins_sieved = find(sum(~isnan(gaz.data.gpDurGaze),2)>Args.MinObs);
        data.MinTrials = Args.MinObs;
    else
        bins_sieved = unique(gaz.data.sessionTimeGaze(2:end,2));
        bins_sieved(bins_sieved==0) = [];
        bins_sieved(isnan(bins_sieved)) = [];
        data.MinTrials = NaN;
    end
    
    gpdur = nansum(gaz.data.gpDurGaze, 2);
    gpdur_s = gpdur;
    gpdur = gpdur(bins_sieved)'; 
    Pi = gpdur/sum(gpdur);

    flat_spiketimes = NaN(2,size(full_arr,1)*size(full_arr,2));
    temp = full_arr';
    flat_spiketimes(1,:) = temp(:);
    flat_spiketimes(2,:) = repelem(1:size(full_arr,1), size(full_arr,2));
    edge_end = 0.5+size(full_arr,1);
    [N,Hedges,Vedges] = histcounts2(flat_spiketimes(1,:), flat_spiketimes(2,:), gaz.data.sessionTimeGaze(:,1), 0.5:1:edge_end);
    
    N(zero_indices(1:end-1),:) = [];
    N = N';

    location = gaz.data.sessionTimeGaze(:,2)';
    location(zero_indices) = [];
    duration1 = gaz.data.sessionTimeGaze(:,3)';
    duration1(zero_indices) = [];
    
    grid_numbers = bins_sieved;
    firing_counts_full = NaN(size(full_arr,1), length(grid_numbers));
    
%     median_stats = NaN(Args.GridSteps^2,1);
%     var_stats = NaN(Args.GridSteps^2,1);
%     perc_stats = NaN(Args.GridSteps^2,5);
    
    for grid_ind = 1:length(grid_numbers)
        tmp = N(:,location==grid_numbers(grid_ind));
        firing_counts_full(:,grid_ind) = sum(tmp,2);
    end   
    
    firing_counts_full1 = NaN(size(firing_counts_full,1),size(gaz.data.gpDurGaze,1));
    firing_counts_full1(:,grid_numbers) = firing_counts_full;
    lin_spikeLoc_Gaze = firing_counts_full1';
    lin_o_i_Gaze = repmat(gpdur_s,1,Args.NumShuffles+1);
    disp('yes');
    
    %%%
            % Restructure bins from linear to separate grids
        grid_o_i_Gaze = cell(size(gaz.data.binDepths,1),1);
        grid_spikeBin_Gaze = grid_o_i_Gaze;
        grid_smoothed_Gaze = grid_o_i_Gaze;

        for jj = 1:size(gaz.data.binDepths,1) % for each grid
            % Initialise empty matrices
            o_i = nan(gaz.data.binDepths(jj,1),gaz.data.binDepths(jj,2),Args.NumShuffles+1);
            spikeBin = o_i;
            map = o_i;
            % Assign linear bin to grid bin
            for mm = 1:gaz.data.binDepths(jj,1)*gaz.data.binDepths(jj,2) % For every point in linear map
                if mod(mm,gaz.data.binDepths(jj,2)) == 0
                    y = gaz.data.binDepths(jj,2);
                else
                    y = mod(mm,gaz.data.binDepths(jj,2));
                end
                x = ceil(mm/gaz.data.binDepths(jj,2));
                indbins_lin = mm + sum(gaz.data.binDepths(1:jj-1,1).*gaz.data.binDepths(1:jj-1,2));
                % Assign
                o_i(x,y,:) = lin_o_i_Gaze(indbins_lin,:);
                spikeBin(x,y,:) = lin_spikeLoc_Gaze(indbins_lin,:);

            end
            % Collect output 
            grid_o_i_Gaze{jj} = o_i;
            grid_spikeBin_Gaze{jj} = spikeBin;
        end
     
%%%%%%%%%%%%
            grid_o_i_Gaze_original = grid_o_i_Gaze;
            retrievemap = cell(size(grid_o_i_Gaze,1),1);
            for jj = 1:size(grid_o_i_Gaze,1) % For each separate grid
            
                if gaz.data.binDepths(jj,1)*gaz.data.binDepths(jj,2) > 2 % For non-cue/non-hint grids
                    % Pad each grid map with adjoining bins from other grids
                    % Pad with <<5>> extra bin rows
                    n = 5;
                    [retrievemap{jj},grid_o_i_Gaze{jj},grid_spikeBin_Gaze{jj}] = padgrids(n,grid_o_i_Gaze{jj},grid_spikeBin_Gaze{jj},grid_o_i_Gaze,grid_spikeBin_Gaze,gaz.data.gazeSections,jj);
    
                end
            end
            
%%%%%%%%%%%%                
                
    if Args.AdaptiveSmooth
        
        alpha = 1e2;
        grid_smoothed_Gaze{1} = grid_spikeBin_Gaze{1}./grid_o_i_Gaze{1};
        grid_smoothed_Gaze{2} = grid_spikeBin_Gaze{2}./grid_o_i_Gaze{2};
        
        for jj = 3:size(grid_o_i_Gaze,1) % for each grid
            
            disp('debug');
            wip = ones(Args.NumShuffles+1,1);
            gpdur1 = grid_o_i_Gaze{jj};
            preset_to_zeros = gpdur1(:,:,1);
            preset_to_zeros(find(preset_to_zeros>0)) = 1;
            preset_to_zeros(find(preset_to_zeros~=1)) = 0;
            preset_to_zeros = ~preset_to_zeros;
            preset_to_zeros = repmat(preset_to_zeros, [1,1,size(gpdur1,3)]);
            
            firing_counts_full1 = grid_spikeBin_Gaze{jj}; % kw_issue (need to set nans to zeros?)
            gpdur1(isnan(gpdur1)) = 0;
            firing_counts_full1(isnan(firing_counts_full1)) = 0;
            disp('debug');
            to_compute = 1:0.5:(max(size(gpdur1(:,:,1)))/4 + min(size(gpdur1(:,:,1)))/4);
            possible = NaN(length(to_compute),2,size(firing_counts_full1,1),size(firing_counts_full1,2),Args.NumShuffles + 1);
            to_fill = NaN(size(possible,3), size(possible,4), size(possible,5));
            to_fill(preset_to_zeros) = 0;
            
                    for idx = 1:length(to_compute)

                        f=fspecial('disk',to_compute(idx));
                        f(f>=(max(max(f))/3))=1;
                        f(f~=1)=0;

                        possible(idx,1,:,:,:) = repmat(imfilter(gpdur1(:,:,1), f), 1,1,Args.NumShuffles+1);   %./scaler;
                        possible(idx,2,:,:,find(wip)) = imfilter(firing_counts_full1(:,:,find(wip)), f);   %./scaler;

                        logic1 = squeeze(alpha./(possible(idx,1,:,:,:).*sqrt(possible(idx,2,:,:,:))) <= to_compute(idx));
                        slice1 = squeeze(possible(idx,1,:,:,:));
                        slice2 = squeeze(possible(idx,2,:,:,:));

                        to_fill(logic1 & isnan(to_fill)) = slice2(logic1 & isnan(to_fill))./slice1(logic1 & isnan(to_fill));

                        disp('smoothed with kernel size:');
                        disp(to_compute(idx));
                        disp('grids left');
                        disp(sum(sum(sum(isnan(to_fill(:,:,:))))));

                        check = squeeze(sum(sum(isnan(to_fill),2),1));
                        wip(check==0) = 0;

                        if sum(sum(sum(isnan(to_fill(:,:,:))))) == 0
                            disp('breaking');
                            break;
                        end

                    end       
                    
                    clear possible;
            to_fill(isnan(to_fill)) = 0;
            to_fill = to_fill(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2),:);
            grid_smoothed_Gaze{jj} = to_fill;
            
%             tmp10 = grid_o_i_Gaze{jj};
%             tmp10 = tmp10(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2),:);
%             grid_o_i_Gaze{jj} = tmp10;
%             tmp10 = grid_spikeBin_Gaze{jj};
%             tmp10 = tmp10(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2),:);
%             grid_spikeBin_Gaze{jj} = tmp10;
        end

        disp('checkpoint');
        
        % smoothing part ends
        data.maps_adsmooth = grid_smoothed_Gaze;
        
    else % kw_notes definitely broken
        firing_rates_full = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);

        to_save = NaN(1,Args.GridSteps^2);
        to_save(bins_sieved) = firing_rates_full(1,:);
        data.maps_raw = to_save;
    end
    
    
    
%%%%%%%%%%%% new SIC part
    
    total_grids = 0;
    for jj = 1:size(grid_smoothed_Gaze,1)
        total_grids = total_grids + size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2);
    end
    gpdur1 = zeros(1,total_grids);
    lambda_i = NaN(Args.NumShuffles+1,total_grids);
    filling_index = 0;
    for jj = 1:size(grid_o_i_Gaze,1)
        temp4 = reshape(grid_o_i_Gaze_original{jj}, [size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2) Args.NumShuffles+1]);
        temp4 = temp4(:,1);
        temp3 = reshape(grid_smoothed_Gaze{jj}, [size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2) Args.NumShuffles+1]);
        gpdur1(filling_index+1:filling_index+size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2)) = temp4(1,:);
        lambda_i(:,filling_index+1:filling_index+size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2)) = temp3';
        filling_index = filling_index + size(grid_smoothed_Gaze{jj},1)*size(grid_smoothed_Gaze{jj},2);
    end    
    
    Pi1 = gpdur1/sum(gpdur1);
    Pi1 = repmat(Pi1, Args.NumShuffles+1, 1);
    lambda_i(isnan(lambda_i)) = 0;
      
    
%%%%%%%%%%%%
        
        lambda_bar = sum(Pi1 .* lambda_i,2);
        % divide firing for each position by the overall mean
        FRratio = lambda_i./repmat(lambda_bar,1,total_grids);
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data.SIC = sic_out(1);
    data.SICsh = sic_out;
%     data.median_occ_firings = median_stats';
%     data.variance_occ_firings = var_stats';
%     data.perc_occ_firings = perc_stats';
    data.gridSteps = Args.GridSteps;
    
	% create nptdata so we can inherit from it    
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
