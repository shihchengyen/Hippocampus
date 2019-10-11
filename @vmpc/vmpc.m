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
Args.classname = 'vmpc';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vmpc';

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
    
	uma = umaze('auto',varargin{:});
%     uma = load('../../../unitymaze.mat');
%     uma = uma.um;
%     uma.data.zero_indices = find(uma.data.sessionTime(:,2)==0);
    
	rp = rplparallel('auto',varargin{:});
	spiketrain = load(Args.RequiredFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spiketimes = spiketrain.timestamps/1000; % now in seconds
    maxTime = rp.data.timeStamps(end,3);
    tShifts = [0 ((rand([1,Args.NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
    full_arr = repmat(spiketimes, Args.NumShuffles+1, 1);
    full_arr = full_arr + tShifts';
    keepers = length(spiketimes) - sum(full_arr>maxTime, 2);
    for row = 2:size(full_arr,1)

        full_arr(row,:) = [full_arr(row,1+keepers(row):end)-maxTime-1 full_arr(row,1:keepers(row))];

    end
    
    % calculating proportion of occupied time in each grid position across
    % entire session.
    
    if Args.FiltLowOcc
        bins_sieved = uma.data.sTPu;
        data.MinTrials = uma.data.Args.MinObs;
    else
        bins_sieved = uma.data.sortedGPindinfo(:,1);
        data.MinTrials = NaN;
    end
    
    gpdur = sum(uma.data.gpDurations, 2);
    gpdur = gpdur(bins_sieved)'; 
    Pi = gpdur/sum(gpdur);

    flat_spiketimes = NaN(2,size(full_arr,1)*size(full_arr,2));
    temp = full_arr';
    flat_spiketimes(1,:) = temp(:);
    flat_spiketimes(2,:) = repelem(1:size(full_arr,1), size(full_arr,2));
    edge_end = 0.5+size(full_arr,1);
    [N,Hedges,Vedges] = histcounts2(flat_spiketimes(1,:), flat_spiketimes(2,:), uma.data.sessionTime(:,1), 0.5:1:edge_end);
    size(full_arr)
    N(uma.data.zero_indices(1:end-1),:) = [];
    N = N';

    location = uma.data.sessionTime(:,2)';
    location(uma.data.zero_indices) = [];
    
    grid_numbers = bins_sieved;
    firing_counts_full = NaN(size(full_arr,1), length(grid_numbers));
    
    median_stats = NaN(Args.GridSteps^2,1);
    var_stats = NaN(Args.GridSteps^2,1);
    perc_stats = NaN(Args.GridSteps^2,5);
    
    for grid_ind = 1:length(grid_numbers)
        tmp = N(:,location==grid_numbers(grid_ind));
        var_stats(grid_numbers(grid_ind)) = var(tmp(1,:));
        median_stats(grid_numbers(grid_ind)) = median(tmp(1,:));
        perc_stats(grid_numbers(grid_ind),:) = prctile(tmp(1,:), [2.5 25 50 75 97.5]);
        firing_counts_full(:,grid_ind) = sum(tmp,2);
    end   
    
    if Args.AdaptiveSmooth
        
        firing_rates_full_raw = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);
        to_save = NaN(1,Args.GridSteps^2);
        to_save(bins_sieved) = firing_rates_full_raw(1,:);
        data.maps_raw = to_save;

        alpha = 1e4;
%         valid = zeros(1,Args.GridSteps^2);
%         valid(bins_sieved) = 1;
%         valid = reshape(valid, Args.GridSteps,Args.GridSteps);
        
        
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
          
%             if sum(sum(sum(isnan(to_fill(:,:,:))))) <= 0.05*(Args.GridSteps^2)*Args.NumShuffles % engage secondary mode of calculating
%                 to_fill = fill_remainder(to_fill, gpdur1, firing_counts_full1, to_compute(idx)+0.5, round(max(to_compute)), alpha);
%                 break;
%             end
                
        end
        
        to_fill(isnan(to_fill)) = 0;
        to_fill = permute(to_fill, [3 1 2]);
        to_fill = reshape(to_fill, Args.NumShuffles + 1, Args.GridSteps^2);
        to_fill = to_fill(:,bins_sieved);
        firing_rates_full = to_fill;

        % smoothing part ends
        to_save = NaN(1,Args.GridSteps^2);
        to_save(bins_sieved) = firing_rates_full(1,:);
        data.maps_adsmooth = to_save;
        
    else
        firing_rates_full = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);

        to_save = NaN(1,Args.GridSteps^2);
        to_save(bins_sieved) = firing_rates_full(1,:);
        data.maps_raw = to_save;
    end
    
        gpdur1 = zeros(1,Args.GridSteps^2);
        gpdur1(bins_sieved) = gpdur;
        Pi1 = gpdur1/sum(gpdur1);
        Pi1 = repmat(Pi1, Args.NumShuffles+1, 1);
        
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data.SIC = sic_out(1);
    data.SICsh = sic_out;
    data.median_occ_firings = median_stats';
    data.variance_occ_firings = var_stats';
    data.perc_occ_firings = perc_stats';
    
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

% function [to_fill] = fill_remainder(to_fill, gpdur1, firing_counts_full1, starting_size, max_to_consider, alpha)
% 
% disp('starting with this size');
% disp(starting_size);
% temp = padarray(to_fill, [max_to_consider, max_to_consider, 0], 'both');
% 
% targets = isnan(temp);
% [I, J, K] = find(targets);
% disp('total number of grids to fill');
% disp(length(J));
% 
% gpdur1 = padarray(gpdur1, [max_to_consider, max_to_consider, 0], 'both');
% firing_counts_full1 = padarray(firing_counts_full1, [max_to_consider, max_to_consider, 0], 'both');
% 
% for tind = 1:length(J)
%     
%     disp('doing grid number');
%     disp(tind);
%     size1 = starting_size;
%     
%     disp('with this size');
%     while size1 < max_to_consider
%         disp(size1);
%         f=fspecial('disk',size1);
%         f(f>=(max(max(f))/3))=1;
%         f(f~=1)=0;
%         
%         a1 = int32(I(tind))-int32(floor(size1));
%         a2 = int32(I(tind))+int32(floor(size1));
%         a3 = int32(J(tind))-int32(floor(size1));
%         a4 = int32(J(tind))+int32(floor(size1));
%         
%         d = sum(sum(gpdur1(a1:a2,a3:a4,K(tind)).*f));
%         s = sum(sum(firing_counts_full1(a1:a2,a3:a4,K(tind)).*f));
%         if size1 >= alpha/(d*sqrt(s))
%             temp(I(tind),J(tind),K(tind)) = s/d;
%             break;
%         end
%         
%         size1 = size1 + 0.5;
%     end
% end
% 
% to_fill = temp(max_to_consider+1:end-max_to_consider,max_to_consider+1:end-max_to_consider);



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
