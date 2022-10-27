function [obj, varargout] = vmhd(varargin)
%@vmhd Constructor function for vmhd class
%   OBJ = vmhd(varargin)
%
%   OBJ = vmhd('auto') attempts to create a vmhd object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmhd %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmhd('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell', 'RequiredFile','spiketrain.mat', ...
				'GridSteps',40, 'DirSteps',60,...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4,'SmoothType','Boxcar', 'SelectiveCriteria','RV', ...
                'UseMinObs',1, 'ThresVel',1, 'UseAllTrials',1, 'Alpha', 10000);
            
Args.flags = {'Auto','ArgsOnly','FRSIC','UseAllTrials','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps','DirSteps','NumShuffles','UseMinObs','AdaptiveSmooth','ThresVel','UseAllTrials', 'Alpha'};                           

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmhd';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vmd';

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
	uma = umaze('auto',varargin{:});
    ufile = unityfile('auto',varargin{:});
	rp = rplparallel('auto',varargin{:});
    
    %%%% PATCH
    cd ..; cd ..; cd ..;
    pv = load('1vmpv.mat');
    pv = pv.pv;
    %%%%%%%
    
    cd(ori);
    spiketrain = load(Args.RequiredFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nTrials = size(uma.data.dirDurations,2);
%     midTrial = ceil(nTrials/2);
%     sessionTimeInds = find(uma.data.sessionTimeDir(:,2) == 0); % Look for trial end markers
%     sessionTimeInds(1) = []; % Remove first marker which marks first cue-off, instead of trial end

    NumShuffles = Args.NumShuffles;
    for repeat = 1:3 % 1 = full trial, 2 = 1st half, 3 = 2nd half
        
        if repeat == 1
            disp('Full session:');
        elseif repeat == 2
            disp('First half:');
        elseif repeat == 3
            disp('Second half:');
        end
        
        if repeat > 1
            NumShuffles = 0;
        end
        
        if repeat == 1
            stc = pv.data.sessionTimeC;
        end
        
        % spike shuffling
        
        spiketimes = spiketrain.timestamps/1000; % now in seconds
%         maxTime = rp.data.timeStamps(end,3);
        maxTime = pv.data.rplmaxtime;
        tShifts = [0 ((rand([1,NumShuffles])*diff(Args.ShuffleLimits))+Args.ShuffleLimits(1))*maxTime];
        full_arr = repmat(spiketimes, NumShuffles+1, 1);
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
            disp('    Filtering...');
            stc(:,5) = [diff(stc(:,1)); 0];
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
            bins_sieved = 1:Args.DirSteps; % Currently don't have a threshold for min number of obs for head dir bins
%             bins_sieved = pv.data.place_good_bins; 
            conditions = conditions & (pv.data.pv_good_rows); % Make sure minobs take into account both place and view
        else
            bins_sieved = 1:Args.DirSteps;
        end
        
%         % calculating proportion of occupied time in each grid position across
%         % entire session.
% 
%         if Args.FiltLowOcc
%             bins_sieved = uma.data.sTDu;
%             data.MinTrials = uma.data.Args.MinObs;
%         else
%             bins_sieved = uma.data.sortedDirindinfo(:,1);
%             data.MinTrials = NaN;
%         end
        
        if repeat == 1 
            % Group into intervals those consecutive rows where same head dir bin is occupied
            dstc = diff(stc(:,1));
            stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
            stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
            stc_changing_ind = stc_changing_ind(1:end-1,:);
        end
        
        consol_arr = zeros(Args.DirSteps,NumShuffles + 1);
        interval = 1;
        if repeat == 1
            disp('    Assigning spikes to bins...');
        end
        for sp = 1:size(flat_spiketimes,1)

            while interval < size(stc_changing_ind,1)
                if flat_spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && ... % > start timestamp of this interval but < start timestamp of next
                        flat_spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
                    break;
                end
                interval = interval + 1; % didn't fall in this interval, search in the next interval
            end   

%             bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),2); % find the relevant place bin
%             bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2)))); % take out bins that don't satisfy filters
%             bins_hit(~(bins_hit>0)) = []; % take out bins where place bin = 0 
%             consol_arr(bins_hit,flat_spiketimes(sp,2)) = consol_arr(bins_hit,flat_spiketimes(sp,2)) + 1;

            
%             bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),[2 3]); % find the relevant place and view bin
            bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),[2 3 4]);
            bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))),:); % take out bins that don't satisfy filters
%             if ~isempty(bins_hit)
%                 disp('debug');
%             end
            if any(~(bins_hit(:,2)>0))
                bins_hit(~(bins_hit(:,2)>0),:) = []; % take out bins where HD bin = 0
%                 if any(~(bins_hit(:,1)>0)) || any(~(bins_hit(:,3)>0)) % If want to only consider samples with place and view data
%                     disp('problem place/view bins still exist');
%                     bins_hit(~(bins_hit(:,1)>0),:) = []; % take out bins where place bin = 0 
%                     bins_hit(~(bins_hit(:,3)>0),:) = []; % take out bins where view bin = nan
%                 end
            end
%             if isempty(bins_hit)
%                 continue;
%             end
            consol_arr(bins_hit(:,2),flat_spiketimes(sp,2)) = consol_arr(bins_hit(:,2),flat_spiketimes(sp,2)) + 1;

        end        
        
        firing_counts = consol_arr';

        % Remove non-place and non-view rows for duration
        stc_ss = stc(find(conditions==1),[2 3 4 5]); % [place hd view dur] with new vmpv that has head dir
%         stc_ss(~(stc_ss(:,1) > 0),:) = []; % remove place bin = 0
%         stc_ss(isnan(stc_ss(:,3)),:) = []; % remove NaN view bins
        stc_ss(~(stc_ss(:,2) > 0),:) = []; % remove hd bin = 0
        stc_ss(:,[1 3]) = []; % [place dur];
        stc_ss = [stc_ss; [60 0]];
        gpdurfull = accumarray(stc_ss(:,1),stc_ss(:,2))';
        
        % Remove low observation bins
        firing_counts_full = zeros(NumShuffles+1,Args.DirSteps);
        gpdur = zeros(1,Args.DirSteps);
        firing_counts_full(:,bins_sieved) = firing_counts(:,bins_sieved);
        gpdur(1,bins_sieved) = gpdurfull(1,bins_sieved);
        
        firing_rates_full_raw = firing_counts_full./repmat(gpdur,size(firing_counts_full,1),1);

        % Save raw maps
        to_save = firing_rates_full_raw(1,:);
        dur_raw = gpdur;
        spk_raw = firing_counts_full(1,:);
        if repeat == 1
            data.maps_raw = to_save;
            data.dur_raw = dur_raw;
            data.spk_raw = spk_raw;
        elseif repeat == 2
            data.maps_raw1 = to_save;
            data.dur_raw1 = dur_raw;
            data.spk_raw1 = spk_raw;
        elseif repeat == 3
            data.maps_raw2 = to_save;
            data.dur_raw2 = dur_raw;
            data.spk_raw2 = spk_raw;
        end

%             if Args.AdaptiveSmooth
                
                if repeat == 1
                    disp('Smoothing...');
                end
%                 nan_track = isnan(to_save);
                
                % Reshape shuffled maps
                to_smooth = firing_rates_full_raw;
                % Smooth with moving window average of n bins
                n = 5;
%                 [firing_rates_full]=smoothDirMap(to_smooth,n,Args.DirSteps); % legacy
                [firing_rates_full]=smoothdir(to_smooth,n,Args.DirSteps);
                % smoothing part ends
                
                if repeat == 1
%                     data.maps_adsm = firing_rates_full(1,:);
%                     data.maps_adsmsh = firing_rates_full(2:end,:);
%                     data.dur_adsm = to_fill_time(1,:);
%                     data.dur_adsmsh = to_fill_time(2:end,:);
%                     data.radii = to_fill_radius(1,:);
%                     data.radiish = to_fill_radius(2:end,:);
                    data.maps_sm = firing_rates_full(1,:);
                    data.maps_smsh = firing_rates_full(2:end,:);
%                     data.maps_dksm = maps_dksm(1,:);
%                     data.maps_dksmsh = maps_dksm(2:end,:);
                elseif repeat == 2
%                     data.maps_adsm1 = firing_rates_full(1,:);
%                     data.dur_adsm1 = to_fill_time(1,:); 
%                     data.radii1 = to_fill_radius(1,:);
                    data.maps_sm1 = firing_rates_full(1,:);
%                     data.maps_dksm1 = maps_dksm(1,:);
                elseif repeat == 3
%                     data.maps_adsm2 = firing_rates_full(1,:);
%                     data.dur_adsm2 = to_fill_time(1,:);
%                     data.radii2 = to_fill_radius(1,:);
                    data.maps_sm2 = firing_rates_full(1,:);
%                     data.maps_dksm2 = maps_dksm(1,:);
                end

%             else
%                 firing_rates_full = firing_rates_full_raw;
%                 to_fill_time = repmat(gpdur,NumShuffles+1,1); % HM added
% 
%             end

            % Rayleigh vector
            binSize=(pi*2)/length(firing_rates_full(1,:));
            binAngles=(0:binSize:( (359.5/360)*2*pi )) + binSize/2;
            binWeights=firing_rates_full./(max(firing_rates_full,[],2));
            S=sum( sin(binAngles).*binWeights , 2);
            C=sum( cos(binAngles).*binWeights , 2);
            R=sqrt(S.^2+C.^2);
            meanR=R./sum(binWeights,2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        if repeat == 1
%             data.SIC_adsm = sic_adsm(1);
%             data.SICsh_adsm = sic_adsm(2:end,1);
            data.crit_sm = meanR(1,1);
            data.critsh_sm = meanR(2:end,1);
%             data.SIC_dksm = sic_dksm(1);
%             data.SICsh_dksm = sic_dksm(2:end,1);
        %     data.median_occ_firings = median_stats';
        %     data.variance_occ_firings = var_stats';
        %     data.perc_occ_firings = perc_stats';
        %     data.occ_data = occ_data;
        elseif repeat == 2
%             data.SIC_adsm1 = sic_adsm;
            data.crit_sm1 = meanR;
%             data.SIC_dksm1 = sic_dksm;
        elseif repeat == 3
%             data.SIC_adsm2 = sic_adsm;
            data.crit_sm2 = meanR;
%             data.SIC_dksm2 = sic_dksm;
        end
            
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    end

            % create nptdata so we can inherit from it   
            data.DirBins = Args.DirSteps;
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

% function [map_sm]=smoothDirMap(map_raw,n,dirSteps)
% % Sliding window average of n bins
% % raw map input needs to be in column form. 
% % M dir bins by N shuffles
% dim1 = size(map_raw,1);
% dim2 = size(map_raw,2);
% flip = false;
% if dim1 ~= dirSteps
%     flip = true;
%     map_raw = map_raw';
% end
% % Smooth a dir map.
% if n==1; map_sm=map_raw; return; end
% p = (n-1)/2;                                               % Pad for circular smooth
% pad_map = [map_raw(end-p+1:end,:); map_raw; map_raw(1:p,:)];                    %  ..
% map_sm = mean( im2col(pad_map, [n 1], 'sliding') );
% % Reshape
% map_sm = reshape(map_sm,dirSteps,size(map_sm,2)/dirSteps);
% if flip
%     map_sm = map_sm';
% end
