function [obj, varargout] = vmcorr(varargin)
%@vmcorr Constructor function for vmcorr class
%   OBJ = vmcorr(varargin)
%
%   OBJ = vmcorr('auto') attempts to create a vmcorr object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmcorr %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = vmcorr('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell', 'RequiredFile','spiketrain.mat', ...
				'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',100000, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4, 'UseMinObs',0, 'ThresVel',1, 'UseAllTrials',1, 'StartOrig',0,'AlphaPlace',10000,'AlphaView',1000,...
                'ConvergeLim',0.001,'LLHIterLim',1000,'UseIterLim',0); 
            
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
Args.classname = 'vmcorr';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'vmcorr';

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
        % Temporary workaround so that we can use file named '1vmpv.mat'
        cd ..; cd ..; cd ..;
        pv = load('1vmpv.mat');
        pv = pv.pv;
    cd(ori);
    spiketrain = load(Args.RequiredFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    viewSections = {'Cue', 'Hint', 'Ground', 'Ceiling', 'Walls', 'Pillar1', 'Pillar2', 'Pillar3', 'Pillar4'};
    viewbinDepths = [1 1;
        1 1;
        40 40;
        40 40;
        8 160;
        5 32;
        5 32;
        5 32;
        5 32];
    placebinDepths = [pv.data.Args.GridSteps pv.data.Args.GridSteps];
    headdirectionbinDepths = [pv.data.headdirectionbins 1];

    % Numbers of pair permutations given input numbers of spatial variables
    spatialvars = {'place','view','headdirection'};
    spatialvarpairs = {{'place','view'}, {'place','headdirection'},{'headdirection','view'}}; % view must be in var2 position because of backfilling spikes later

    % Filter

    for repeat = 1:3 % 1 = full trial, 2 = 1st half, 3 = 2nd half
        
        for combi = 1:size(spatialvarpairs,2) % For each pair of spatial vars
            
            Var1 = spatialvarpairs{combi}{1};
            Var2 = spatialvarpairs{combi}{2};
            combiname = [lower(Var1(1)) lower(Var2(1))];
            
            % Debug
            if combi == 3
                disp('hd x v');
            end
        
            if repeat == 1
    %             stc = pv.data.sessionTimeC;
                disp('Full session:');
            elseif repeat == 2
                disp('1st half:');
            elseif repeat == 3
                disp('2nd half:');
            end

            stc = pv.data.sessionTimeC; % time, place, hd, view
            stc(:,5) = [diff(stc(:,1)); 0]; % column 5 as duration

            % spike binning
            disp('Binning spikes ...');
            spiketimes = spiketrain.timestamps/1000; % now in seconds
            spiketimes = spiketimes';
            spiketimes(spiketimes(:,1) < stc(1,1),:) = [];   
            binned = histcounts(spiketimes, stc(:,1))'; % Method 1 of 2 for binning spikes
            stc(:,6) = [binned; 0];

            % filtering rows from sessionTimeC
            disp('Filtering ...');
            conditions = ones(size(stc,1),1);
            if Args.UseAllTrials == 0
                conditions = conditions & pv.data.good_trial_markers;
            end
            if Args.ThresVel > 0
                conditions = conditions & get(pv,'SpeedLimit',Args.ThresVel);
            end
            if Args.UseMinObs %%%% FIX
                switch Var1
                    case 'place'
                        bins_sieved_var1 = pv.data.place_good_bins;
                        bins_removed_var1 = setdiff(1:size(pv.data.place_intervals_count,1),bins_sieved_var1);
                    case 'view'
                        bins_sieved_var1 = pv.data.view_good_bins;
                        bins_removed_var1 = setdiff(1:size(pv.data.view_intervals_count,1),bins_sieved_var1);
                    case 'headdirection' %%% FIX
                        
                end
                switch Var2
                    case 'place'
                        bins_sieved_var2 = pv.data.place_good_bins;
                        bins_removed_var2 = setdiff(1:size(pv.data.place_intervals_count,1),bins_sieved_var2);
                    case 'view'
                        bins_sieved_var2 = pv.data.view_good_bins;
                        bins_removed_var2 = setdiff(1:size(pv.data.view_intervals_count,1),bins_sieved_var2);
                    case 'headdirection' %%FIX
                        
                end
                    conditions = conditions & (pv.data.pv_good_rows); % Make sure maps take into account both place and view filters
            else
                bins_sieved_var1 = 1:pv.data.([Var1 'bins']);
                bins_removed_var1 = [];
                bins_sieved_var2 = 1:pv.data.([Var2 'bins']);
                bins_removed_var2 = [];
            end
            if repeat == 2
                conditions = conditions & (pv.data.halving_markers==1);
            elseif repeat == 3
                conditions = conditions & (pv.data.halving_markers==2);
            end

            %%
    %         if repeat == 1 % Full session

    %         % Method 2 of 2 for binning spikes - 
    %         % Used for checking for consistency against the main Method 1 which follows from vmpc/sv
    %         dstc = diff(stc(:,1));
    %         stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
    %         stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
    %         stc_changing_ind = stc_changing_ind(1:end-1,:);
    % 
    %         full_spikes2 = zeros(size(pv.data.view_intervals_count,1),Args.GridSteps * Args.GridSteps);
    %         place_spikes2 = zeros(1,Args.GridSteps * Args.GridSteps);
    %         view_spikes2 = zeros(1,size(pv.data.view_intervals_count,1));
    %         interval = 1;
    %         for sp = 1:size(spiketimes,1)
    % 
    %             while interval < size(stc_changing_ind,1)
    %                 if spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
    %                     break;
    %                 end
    %                 interval = interval + 1;
    %             end   
    % 
    %             bins_hit = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),[2 3]);
    %             bins_hit = bins_hit(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))),:);
    %             bins_hit(~(bins_hit(:,1)>0),:) = [];
    %             bins_hit(~(bins_hit(:,2)>0),:) = [];
    % 
    %             place_spikes2(1,bins_hit(:,1)) = place_spikes2(1,bins_hit(:,1)) + 1;
    %             view_spikes2(1,bins_hit(:,2)) = view_spikes2(1,bins_hit(:,2)) + 1;
    %             full_spikes2(bins_hit(:,2),bins_hit(:,1)) = full_spikes2(bins_hit(:,2),bins_hit(:,1)) + 1;
    % 
    %         end        
    % 
    %         % Backfill duration for view bins that occupy the same time bin
    %         stc_v = stc; 
    %         stc_v(stc_v(:,4)==0,4) = nan;
    %         stc_v(:,4) = fillmissing(stc_v(:,4), 'next'); % If > 1 view bin for 1 place bin, time is recorded with the last view bin
    % 
    %         % Remove non-place and non-view rows for duration
    %         stc_p = stc;
    %         stc_p = stc_p(find(conditions==1),[2 3 4 5]); % [place view dur spk]
    %         stc_p(~(stc_p(:,1)>0),:) = []; % remove place bin = 0
    %         stc_p(~(stc_p(:,2)>0),:) = []; % remove NaN view bins
    %         stc_p = [stc_p; [size(pv.data.place_intervals_count,1) size(pv.data.view_intervals_count,1) 0 0]];
    %         stc_v = stc_v(find(conditions==1),[2 3 4 5]); % [place view dur spk]
    %         stc_v(~(stc_v(:,1)>0),:) = []; % remove place bin = 0
    %         stc_v(~(stc_v(:,2)>0),:) = []; % remove NaN view bins
    %         stc_v = [stc_v; [size(pv.data.place_intervals_count,1) size(pv.data.view_intervals_count,1) 0 0]];
    %         gpdurplace = accumarray(stc_p(:,1),stc_p(:,3))';
    %         gpdurview = accumarray(stc_v(:,2),stc_v(:,3))';
    % 
    %         % Remove low observation bins
    %         place_spikes_full = zeros(1,Args.GridSteps*Args.GridSteps);
    %         gpdurp = zeros(1,Args.GridSteps*Args.GridSteps);
    %         place_spikes_full(:,bins_sieved_p) = place_spikes2(:,bins_sieved_p);
    %         gpdurp(1,bins_sieved_p) = gpdurplace(1,bins_sieved_p);
    %         view_spikes_full = zeros(1,size(pv.data.view_intervals_count,1));
    %         gpdursv = zeros(1,size(pv.data.view_intervals_count,1));
    %         view_spikes_full(:,bins_sieved_sv) = view_spikes2(:,bins_sieved_sv);
    %         gpdursv(1,bins_sieved_sv) = gpdurview(1,bins_sieved_sv);
    % 
    %         rates_place = place_spikes_full./gpdurp;
    %         rates_view = view_spikes_full./gpdursv;
    %         end

            % Method 1 of 2 for binning spikes - 
            % Uses a more efficient method than that used in vmpc/sv
            % Since it's different from previous objects, I've included Method
            % 2 for checking consistency with vmpc/sv objects 

            %% Consolidate pv array into view by place bin array

            % remove filtered rows
            stc_ssv = stc(conditions == 1,:);
            stc_colnames = {'time','place','headdirection','view','dur','spk'};
            keepcol = [find(ismember(stc_colnames,spatialvarpairs{combi}{1})) ...
                find(ismember(stc_colnames,spatialvarpairs{combi}{2}))];
%             stc_ssv = stc(conditions==1,[keepcol 5 6]);

            % Initialise variables
            var1_durations1 = nan(1,pv.data.([Var1 'bins']));
            var1_spikes1 = zeros(1,pv.data.([Var1 'bins']));
            full_durations1 = zeros(pv.data.([Var2 'bins']),pv.data.([Var1 'bins']));
            full_spikes1 = zeros(pv.data.([Var2 'bins']),pv.data.([Var1 'bins']));
            for ii = 1:pv.data.([Var1 'bins'])

                inds = stc_ssv(:,keepcol(1))==ii;
                subsample = stc_ssv(inds,:);
                % Consider only samples where both variables are sampled (nan for nonsampled view, 0 for nonsampled place, 0 for nonsampled hd)
                for cc = 1:length(keepcol)
                    if strcmp(spatialvarpairs{combi}(cc),'view')
                        subsample(isnan(subsample(:,keepcol(cc))),:) = [];
                    else 
                        subsample(subsample(:,keepcol(cc)) < 1) = [];
                    end
                end
                if ~isempty(subsample) 
                    % Get spikes and duration for place only
                    var1_durations1(1,ii) = sum(subsample(:,5));
                    var1_spikes1(1,ii) = sum(subsample(:,6));
                    % back-filling spikes for Var2
                    subsample(subsample(:,6)==0,6) = nan;
    %                 subsample(:,7) = circshift(subsample(:,5)~=0 ,-1); % Use only if spike is recorded in first bin of view set and time in last bin (like how it was before)
                    subsample(:,7) = subsample(:,5)~=0;
                    subsample(isnan(subsample(:,6)) & subsample(:,7), 6) = 0;
                    subsample(:,7) = [];
                    subsample(:,6) = fillmissing(subsample(:,6), 'next'); % In current version of pv, time and spikes are recorded in last repeat bin
                    % back-filling time for Var2
                    subsample(subsample(:,5)==0,5) = nan;
                    subsample(:,5) = fillmissing(subsample(:,5), 'next'); % If > 1 view bin for 1 place bin, time is recorded with the last view bin
                    % padding with max Var2 bin
                    if subsample(end,keepcol(2)) ~= pv.data.([Var2 'bins'])
    %                     subsample = [subsample; [NaN 5122 NaN]]; % Used to work, but now is value is nan, sum also becomes nan.
                        segment = [0 0 0 0 0 0];
                        segment(1,keepcol) = [ii pv.data.([Var2 'bins'])];
                        subsample = [subsample; segment];
                    end
                    % remove bad view spots
                    subsample(isnan(subsample(:,keepcol(ismember(spatialvarpairs{combi},'view')))),:) = [];
                    % sum durations
    %                 full_durations1(:,ii) = accumarray(subsample(:,2), subsample(:,3),[],[],NaN);
                    full_durations1(:,ii) = accumarray(subsample(:,keepcol(2)), subsample(:,5),[],[],0);
                    % sum spikes % Method 1 of 2 to bin spikes (counter check Method 2
                    full_spikes1(:,ii) = accumarray(subsample(:,keepcol(2)), subsample(:,6),[],[],0);
                end
            end

            % Filter out low-sampled bins
            var1_durations1(isnan(var1_durations1)) = 0; % Necessary because NaNs seem to mess up the smoothing. Can't do the same for full_durations because it is being used to compute llh and will erroneously give inf. 
    %         full_durations(isnan(full_durations)) = 0; 
            full_durations1(isnan(full_durations1)) = 0;
            var1_durations1(bins_removed_var1) = 0;
            full_durations1(bins_removed_var2,:) = 0;
            full_durations1(:,bins_removed_var1) = 0;
            var1_spikes1(bins_removed_var1) = 0;
            full_spikes1(bins_removed_var2,:) = 0;
            full_spikes1(:,bins_removed_var1) = 0;

            var1_array_orig = var1_spikes1./var1_durations1;
            var2_array_orig = sum(full_spikes1,2)./sum(full_durations1,2);

    %         % Consistency check between Methods 1 and 2 for binning spikes
    %         temp_p = rates_place - p_array_orig;
    %         temp_p(isnan(temp_p)) = 0;
    %         temp_v = rates_view - sv_array_orig';
    %         temp_v(1,1:2) = 0;
    %         temp_v(isnan(temp_v)) = 0;
    %         if repeat == 1 && (sum(temp_p)~=0 || sum(temp_v)~=0)
    %             disp(sum(temp_p));
    %             disp(sum(temp_v));
    %         end

            %% Compute covariance matrix
            if repeat == 1
                disp('Computing covariance matrix');
                var2mappervar1bin = nan(size(full_durations1));
                for ii = 1:size(var1_durations1,2)
                    var2mappervar1bin(:,ii) = full_spikes1(:,ii)./full_durations1(:,ii);
                end
                covmat = cov(var2mappervar1bin,'partialrows');
                covmat_norm = covmat./nanmax(nanmax(abs(covmat)));
                % Replace NaNs with zeros in covariance matrix for norm calculations
                covmat_nonan = covmat;
                covmat_nonan(isnan(covmat_nonan)) = 0;
                % Calculate norms
                l1norm = norm(covmat_nonan,1); % maximum of column sum
                l2norm = norm(covmat_nonan,2); % maximum single value
            end

            %% Distributive hypothesis testing

            % Null hypothesis: Under the assumption that firing of cell is
            % ideally specific to Var 1, predicted Var 2 firing can be 
            % calculated and should not be modulated by Var 2. 
            % (e.g. place cell firing is ideally location-specific, 
            % direction firing can be predicted from place map, and should
            % be non-selective for direction. Test for difference between
            % real HD map and predicted HD map. If difference is large,
            % with predicted deviating little from baseline,
            % then the cell is truly directional. If the difference is
            % small, and predicted has same peaks as observed, then there
            % is little directional selectivity. 

            % Removing hint and cue effects
            var2_array_orig_trun = var2_array_orig;
            full_durations1_trun = full_durations1;
            if strcmp(Var2,'view')
                 % disregard cue and hint bins
                var2_array_orig_trun(1:2) = nan;
                full_durations1_trun(1:2,:) = 0;
            end

            % Predicted rate as a function of var2
            topterm = nansum(var1_array_orig.*full_durations1_trun,2);
            bottomterm = nansum(full_durations1_trun,2);
            var2_array_pred = topterm./bottomterm; % raw map

            switch Var2
                case 'view'
                    % Adaptive smooth predicted view map
                    % transform from linear to grid
                    toptermG = lineartogrid(topterm,'view',viewbinDepths);
                    bottomtermG = lineartogrid(bottomterm,'view',viewbinDepths);
                    % Pad sv map with 5 extra rows
                    n = 5;
                    padpillar = false;
                    [emptyfloorref_pad,~] = padsvmap(n,bottomtermG,viewSections,padpillar);
                    padpillar = true;
                    [bottomtermGpad,retrievemap] = padsvmap(n,bottomtermG,viewSections,padpillar);
                    [toptermGpad,~] = padsvmap(n,toptermG,viewSections,padpillar);
                    sv_array_pred_adsmGpad = cell(size(viewbinDepths,1),1);
                    for gg = 1:size(viewbinDepths,1)
                        if gg == 1 || gg == 2
                            sv_array_pred_adsmGpad{gg} = toptermG{gg}/bottomtermG{gg};
                        else
                            sv_array_pred_adsmGpad{gg} = adsmooth(bottomtermGpad{gg},toptermGpad{gg},Args.AlphaView);
                        end
                    end
                    % Unpad
                    [sv_array_pred_adsmG] = unpadsvmap(sv_array_pred_adsmGpad,retrievemap,bottomtermG);
                    % Linearise
                    var2_array_pred_adsm = gridtolinear(sv_array_pred_adsmG,'view',viewbinDepths);
                case 'headdirection'
                    n = 5; % Boxcar window of 5
                    [var2_array_pred_adsm] = smoothdir(var2_array_pred,n,pv.data.([Var2 'bins']));
                case 'place'
                    toptermG = cell2mat(lineartogrid(topterm','place',[Args.GridSteps Args.GridSteps]));
                    bottomtermG = cell2mat(lineartogrid(bottomterm','place',[Args.GridSteps Args.GridSteps]));
                    [p_array_pred_adsmG,~,~] = adsmooth(bottomtermG,toptermG,Args.AlphaPlace);
                    var2_array_pred_adsm = gridtolinear({p_array_pred_adsmG},'place',[Args.GridSteps Args.GridSteps])';

            end

            % If place cell firing is only mod by location, and view influence is attributable only to inhomogeous sampling, 
            % then dr = 0. If dr is significant, place cell is view-modulated.
            ratio = log((1+var2_array_orig_trun)./(1+var2_array_pred)); % Muller 1994 
    %         ratio = log(1+sv_array_orig_trun)./(1+sv_array_pred); % Cacucci 2004
            dr_var2 = nansum(abs(ratio))/sum(bottomterm>0); 
            % Correlation between orig and predicted view maps
            vis = ~isnan(var2_array_orig_trun);
            dcorr_var2 = corr2(var2_array_orig_trun(vis),var2_array_pred(vis));

            % Predicted rate as a function of var1
            topterm = nansum(var2_array_orig_trun.*full_durations1_trun,1);
            bottomterm = nansum(full_durations1_trun,1);
            var1_array_pred = topterm./bottomterm; % raw map

            % Adaptive smooth predicted var1 map
            switch Var1
                case 'place' 
                    toptermG = cell2mat(lineartogrid(topterm','place',[Args.GridSteps Args.GridSteps]));
                    bottomtermG = cell2mat(lineartogrid(bottomterm','place',[Args.GridSteps Args.GridSteps]));
                    [p_array_pred_adsmG,~,~] = adsmooth(bottomtermG,toptermG,Args.AlphaPlace);
                    var1_array_pred_adsm = gridtolinear({p_array_pred_adsmG},'place',[Args.GridSteps Args.GridSteps])';
                case 'headdirection'
                    n = 5; % Boxcar window of 5
                    [var1_array_pred_adsm] = smoothdir(var1_array_pred,n,pv.data.([Var1 'bins']));
                case 'view'
                    % Adaptive smooth predicted view map
                    % transform from linear to grid
                    toptermG = lineartogrid(topterm,'view',viewbinDepths);
                    bottomtermG = lineartogrid(bottomterm,'view',viewbinDepths);
                    % Pad sv map with 5 extra rows
                    n = 5;
                    padpillar = false;
                    [emptyfloorref_pad,~] = padsvmap(n,bottomtermG,viewSections,padpillar);
                    padpillar = true;
                    [bottomtermGpad,retrievemap] = padsvmap(n,bottomtermG,viewSections,padpillar);
                    [toptermGpad,~] = padsvmap(n,toptermG,viewSections,padpillar);
                    sv_array_pred_adsmGpad = cell(size(viewbinDepths,1),1);
                    for gg = 1:size(viewbinDepths,1)
                        if gg == 1 || gg == 2
                            sv_array_pred_adsmGpad{gg} = toptermG{gg}/bottomtermG{gg};
                        else
                            sv_array_pred_adsmGpad{gg} = adsmooth(bottomtermGpad{gg},toptermGpad{gg},Args.AlphaView);
                        end
                    end
                    % Unpad
                    [sv_array_pred_adsmG] = unpadsvmap(sv_array_pred_adsmGpad,retrievemap,bottomtermG);
                    % Linearise
                    var1_array_pred_adsm = gridtolinear(sv_array_pred_adsmG,'view',viewbinDepths);
            end

            % If view cell firing is only mod by view location, and place influence is attributable only to inhomogeous sampling, 
            % then dr = 0. If dr is significant, view cell is place-modulated.
            ratio = log((1+var1_array_orig)./(1+var1_array_pred)); % Muller 1994
    %         ratio = log(1+p_array_orig)./(1+p_array_pred); % Cacucci 2004 
            dr_var1 = nansum(abs(ratio))/sum(bottomterm>0);
            % Correlation between orig and predicted place maps
            vis = ~isnan(var1_array_orig);
            dcorr_var1 = corr2(var1_array_orig(vis)',var1_array_pred(vis)');


            %% Establish independence of place and view maps
            disp('Maximising llh ...');
            mlm = true;
            initialwith = Var2;
            convergewithvar2 = false;
            convergewithvar1 = false;
            while mlm == true

                llh_vec = [];
                reinitialise = false;

                var2_spikes_temp = full_spikes1;
                var2_spikes_temp(isnan(full_spikes1)) = 0; % to prevent error on factorial

                % Maximum-likelihood maps. 
                if Args.StartOrig
                    % Start with Original maps
                    var1_array = var1_array_orig;  
                    var2_array = var2_array_orig;
                    var1_spk = var1_spikes1;
                    var2_spk = nansum(full_spikes1,2);
                else
                    % Start with uniform array
                    var1_array = ones(1,size(var1_spikes1,2));
                    var1_array(isnan(var1_array_orig)) = nan;
                    var2_array = ones(size(full_durations1,1),1);
                    var2_array(isnan(var2_array_orig)) = nan;
                    var1_spk = var1_array.*nansum(full_durations1,1);
                    var1_spk(isnan(var1_spk)) = 0;
                    var2_spk = var2_array.*nansum(full_durations1,2);
                    var2_spk(isnan(var2_spk)) = 0;
                end
                % Initialise storage. % Pad with original raw map first in
                % both cases so we can smooth them too for comparison with
                % pc/sv smoothed maps
                var1_array_set = [var1_array_orig; var1_array];
                var2_array_set = [var2_array_orig, var2_array];
                var1_spk_set = [var1_spikes1; var1_spk];
                var2_spk_set = [nansum(full_spikes1,2), var2_spk];

                % Starting llh ( upper limit of factorial 170, anything beyond that becomes Inf )
    %             full_durations1(full_durations1==0) = nan; % zeros mess up llh calculation
                for ii = 1:size(var1_array_set,1)
                    llh = sum( nansum(full_spikes1.*log(var1_array_set(ii,:).*full_durations1.*var2_array_set(:,ii))) - nansum(var1_array_set(ii,:).*full_durations1.*var2_array_set(:,ii)) - nansum(log(factorial(var2_spikes_temp))) );
                    llh_vec(end+1,1) = llh;
                end
                if isinf(llh) % Evaluate the 2nd llh in the set since that is the starting point of the iterations
                    pick = 2;
                    if max(max(var2_spikes_temp)) > 170
                        picklabel = ['spikecount' num2str(max(max(var2_spikes_temp)))];
                        disp(['llh inf from start: spike count ' num2str(max(max(var2_spikes_temp)))]);
                    else
                        picklabel = 'inf1';
                        disp('llh inf from start: other reasons');
                    end
                    % end iterations
                    mlm = false;
                    break;
                end

                while strcmp(initialwith,Var2) && ~convergewithvar2 % Try to find max llh first with flat sv array. If fail to converge, try with flat p array

                    dur_adj = nan(size(full_durations1,1),size(full_durations1,2));
                    for ii = 1:size(full_durations1,2) % var1 (place)
                        for jj = 1:size(full_durations1,1) % var2
                            if ~isnan(full_durations1(jj,ii))
                                if isnan(var2_array(jj))
    %                                 dur_adj(jj,ii) = nan;
                                    dur_adj(jj,ii) = 0;
                                elseif var2_array(jj) == 0
                                    dur_adj(jj,ii) = full_durations1(jj,ii); 
                                else
                                    dur_adj(jj,ii) = full_durations1(jj,ii)*var2_array(jj);
                                end
                            end
                        end
                    end
                    var1_array = nansum(full_spikes1,1)./nansum(dur_adj,1);
                    spk = var1_array.*nansum(full_durations1,1);
                    spk(isnan(spk)) = 0;
                    var1_array_set = [var1_array_set; var1_array];
                    var1_spk_set = [var1_spk_set; spk];

                    dur_adj = nan(size(full_durations1,1),size(full_durations1,2));
                    for jj = 1:size(full_durations1,1) % var2
                        for ii = 1:size(full_durations1,2) % var1 (place)
                            if ~isnan(full_durations1(jj,ii))
                                if isnan(var1_array(ii))
                                    dur_adj(jj,ii) = 0;
    %                                 dur_adj(jj,ii) = nan;
                                elseif var1_array(ii) == 0
                                    dur_adj(jj,ii) = full_durations1(jj,ii); 
                                else
                                    dur_adj(jj,ii) = full_durations1(jj,ii)*var1_array(ii);
                                end
                            end
                        end
                    end
                    var2_array = nansum(full_spikes1,2)./nansum(dur_adj,2);
                    spk = var2_array.*nansum(full_durations1,2);
                    spk(isnan(spk)) = 0;
                    var2_array_set = [var2_array_set var2_array];
                    var2_spk_set = [var2_spk_set spk];

                    llh = sum( nansum(full_spikes1.*log(var1_array.*full_durations1.*var2_array)) - nansum(var1_array.*full_durations1.*var2_array) - nansum(log(factorial(var2_spikes_temp))) );
                    llh_vec(end+1,1) = llh;
    %                 llh = sum( nansum(full_spikes1.*log(p_array_set(ii,:).*full_durations1.*sv_array_set(:,ii))) - nansum(p_array_set(ii,:).*full_durations1.*sv_array_set(:,ii)) - nansum(log(factorial(view_spikes_temp))) );

                    if Args.UseIterLim % Converge if can find max within 1000 iterations

                        % Run iterations in sets of 100. 
                        % Stop if maximum of llh lands in second-to-last set
                        % For every set of 100 iterations
                        if size(llh_vec,1) >= 101 && rem(size(llh_vec,1),100) == 1

                            % Find max of last and next-to-last sets
                            maxlast = max(llh_vec(size(llh_vec,1)-100+1:end));
                            disp(['iter ' num2str(size(llh_vec,1)) ' : currlocalmaxllh = ' num2str(maxlast)]);

                            if size(llh_vec,1) >= 201 
                                maxnextlast = max(llh_vec(size(llh_vec,1)-200+1:size(llh_vec,1)-100));
                                % If llh doesn't converge
                                if size(llh_vec,1) > Args.LLHIterLim
                                    % Switch to using flat view array to start iterating
                                    initialwith = Var1;
                                    reinitialise = true;
                                    disp(['LLH did not converge with initial ' Var2 ' array. Try initialising with ' ...
                                        Var1 ' array.']);
                                    break;
                                end
                                % Test if global max is found
                                if maxnextlast > maxlast 
                                    convergewithvar2 = true;
                                    % Record global max
                                    pick = find(llh_vec == max(llh_vec(2:end)));
                                    if size(pick,1) > 1
                                        pick = pick(1);
                                        picklabel = ['undiff'];
                                    elseif isinf(pick)
                                        % Skip on to using differnt initial array
                                        initialwith = Var1;
                                        reinitialise = true;
                                        disp(['Failed to converge to non-infinite values with initial ' Var2  ...
                                            'array. Try initialising with ' Var1 ' array']);
                                        break;
                                    else
                                        picklabel = num2str(pick);
                                        % display results
                                        disp(['Converged with initial view array, pick = ' num2str(pick)]);
                                        % break out of loop
                                        mlm = false;
                                        break;
                                    end
                                end
                            end
                        end

                    else % Converge if subsequent iterations < 0.1% of delta between iterations 1 and 2

                        % Report every set of 100 iterations
                        if size(llh_vec,1) >= 101 && rem(size(llh_vec,1),100) == 1
                            disp(['Iter ' num2str(size(llh_vec,1))]);
                        end
                        % Check for convergence
                        if abs((llh-llh_vec(end-1))/(llh_vec(2)-llh_vec(3))) < Args.ConvergeLim && size(llh_vec,1) < Args.LLHIterLim % Converged
                            convergewithvar2 = true;
                            % Record global max
                            pick = find(llh_vec == max(llh_vec(2:end)));
                            if size(pick,1) > 1
                                pick = pick(1);
                                picklabel = 'undiff';
                            elseif isinf(pick)
                                % Skip to using different initial array
                                initialwith = Var1;
                                reinitialise = true;
                                disp(['Failed to converge to non-infinite values with initial ' ...
                                    Var2 'array. Try initialising with ' Var1 ' array']);
                                break;
                            else
                                picklabel = num2str(pick);
                                % display results
                                disp(['Converged with initial ' Var2 ' array, pick = ' ...
                                    num2str(pick) ' out of ' num2str(size(llh_vec,1))]);
                                % break out of loop
                                mlm = false;
                                break;
                            end
                        elseif isinf(llh)
                            % Skip to using different initial array
                            initialwith = Var1;
                            reinitialise = true;
                            disp(['Failed to converge to non-infinite values with initial ' ...
                                Var2 ' array. Try initialising with ' Var1 ' array']);
                            break;
                        elseif size(llh_vec,1) > Args.LLHIterLim
                            % No convergence. Switch to using flat view array to start iterating
                            initialwith = Var1;
                            reinitialise = true;
                            disp(['LLH did not converge with initial ' ...
                                Var2 ' array. Try initialising with ' Var1 ' array.']);
                            break;
                        end
                    end

                end
                if reinitialise
                    continue;
                end

                while strcmp(initialwith,Var1) && ~convergewithvar1 % If first try with flat sv array fails to converge, try with flat p array

                    dur_adj = nan(size(full_durations1,1),size(full_durations1,2));
                    for jj = 1:size(full_durations1,1) % view
                        for ii = 1:size(full_durations1,2) % place
                            if ~isnan(full_durations1(jj,ii))
                                if isnan(var1_array(ii))
                                    dur_adj(jj,ii) = 0;
                                elseif var1_array(ii) == 0
                                    dur_adj(jj,ii) = full_durations1(jj,ii); 
                                else
                                    dur_adj(jj,ii) = full_durations1(jj,ii)*var1_array(ii);
                                end
                            end
                        end
                    end
                    var2_array = nansum(full_spikes1,2)./nansum(dur_adj,2);
                    spk = var2_array.*nansum(full_durations1,2);
                    spk(isnan(spk)) = 0;
                    var2_array_set = [var2_array_set var2_array];
                    var2_spk_set = [var2_spk_set spk];

                    dur_adj = nan(size(full_durations1,1),size(full_durations1,2));
                    for ii = 1:size(full_durations1,2) % var1 (place)
                        for jj = 1:size(full_durations1,1) % var2
                            if ~isnan(full_durations1(jj,ii))
                                if isnan(var2_array(jj))
                                    dur_adj(jj,ii) = 0;
                                elseif var2_array(jj) == 0
                                    dur_adj(jj,ii) = full_durations1(jj,ii); 
                                else
                                    dur_adj(jj,ii) = full_durations1(jj,ii)*var2_array(jj);
                                end
                            end
                        end
                    end
                    var1_array = nansum(full_spikes1,1)./nansum(dur_adj,1);
                    spk = var1_array.*nansum(full_durations1,1);
                    spk(isnan(spk)) = 0;
                    var1_array_set = [var1_array_set; var1_array];
                    var1_spk_set = [var1_spk_set; spk];

                    llh = sum( nansum(full_spikes1.*log(var1_array.*full_durations1.*var2_array)) - nansum(var1_array.*full_durations1.*var2_array) - nansum(log(factorial(var2_spikes_temp))) );
                    llh_vec(end+1,1) = llh;

                    if Args.UseIterLim % Converge if can find max within 1000 iterations

                        % Run iterations in sets of 100. 
                        % Stop if maximum of llh lands in second-to-last set
                        % For every set of 100 iterations
                        if size(llh_vec,1) >= 101 && rem(size(llh_vec,1),100) == 1

                            % Find max of last and next-to-last sets
                            maxlast = max(llh_vec(size(llh_vec,1)-100+1:end));
                            disp(['iter ' num2str(size(llh_vec,1)) ' : currlocalmaxllh = ' num2str(maxlast)]);

                            if size(llh_vec,1) >= 201
                                maxnextlast = max(llh_vec(size(llh_vec,1)-200+1:size(llh_vec,1)-100));
                                % If llh doesn't converge
                                if size(llh_vec,1) > Args.LLHIterLim
                                    pick = size(llh_vec,1);
                                    picklabel = 'nonconverged';
                                    disp('LLH did not converge with either initial arrays. Abandon.');
                                    mlm = false;
                                    break;
                                end
                                % Test if global max is found
                                if maxnextlast > maxlast % Yes, global max is found
                                    % Record global max
                                    convergewithvar1 = true;
                                    pick = find(llh_vec == max(llh_vec(2:end)));
                                    if size(pick,1) > 1
                                        pick = pick(1);
                                        picklabel = ['undiff'];
                                    elseif isinf(pick)
                                        picklabel = 'inf';
                                        disp('Failed to converge to non-infinite values with either initial arrays. Abandon');
                                    else
                                        picklabel = num2str(pick);
                                        % display results
                                        disp(['Converged with initial ' Var1 ' array, pick = ' num2str(pick)]);
                                    end
                                    % break out of loop
                                    mlm = false;
                                    break;
                                end
                            end
                        end

                    else % Converge if subsequent iterations < 0.1% of delta between iterations 1 and 2

                        % Report every set of 100 iterations
                        if size(llh_vec,1) >= 101 && rem(size(llh_vec,1),100) == 1
                            disp(['Iter ' num2str(size(llh_vec,1))]);
                        end
                        % Check for convergence
                        if abs((llh-llh_vec(end-1))/(llh_vec(2)-llh_vec(3))) < Args.ConvergeLim && size(llh_vec,1) < Args.LLHIterLim % Converged
                            convergewithvar1 = true;
                            % Record global max
                            pick = find(llh_vec == max(llh_vec(2:end)));
                            if size(pick,1) > 1
                                pick = pick(1);
                                picklabel = 'undiff';
                            elseif isinf(pick)
                                % Terminate
                                picklabel = 'inf';
                                disp('Failed to converge to non-infinite values with either array. Terminate');
                                mlm = false;
                            else
                                picklabel = num2str(pick);
                                % display results
                                disp(['Converged with initial ' ...
                                    Var1 ' array, pick = ' num2str(pick) ' out of ' num2str(size(llh_vec,1))]);
                                % break out of loop
                                mlm = false;
                            end
                        elseif isinf(llh)
                            % Terminate
                            pick = size(llh_vec,1);
                            picklabel = 'inf';
                            disp('Failed to converge to non-infinite values with either array. Terminate');
                            mlm = false;
                        elseif size(llh_vec,1) > Args.LLHIterLim % No convergence
                            pick = size(llh_vec,1);
                            picklabel = 'nonconverged';
                            disp('LLH did not converge with either initial arrays. Terminate.');
                            mlm = false;
                        end
                    end

                end

            end

            % Scale each map to match total observed number of spikes
            var1fac = repmat(nansum(nansum(full_spikes1)),size(var1_spk_set,1),1) ./ nansum(var1_spk_set,2);
            var2fac = repmat(nansum(nansum(full_spikes1)),1,size(var2_spk_set,2))./ nansum(var2_spk_set,1);
            var1fac(1) = 1; % Original place map
            var2fac(1) = 1; % Original view map
            if Args.StartOrig 
                var1fac(2) = 1;
                var2fac(2) = 1;
            end
            var1_array_set = repmat(var1fac,1,size(var1_array_set,2)) .* var1_array_set;
            var1_spk_set = repmat(var1fac,1,size(var1_spk_set,2)) .* var1_spk_set;
            var2_array_set = repmat(var2fac,size(var2_array_set,1),1) .* var2_array_set;
            var2_spk_set = repmat(var2fac,size(var2_spk_set,1),1) .* var2_spk_set;

            % Remove low obs bins
            var1_array_set(:,bins_removed_var1) = nan;
            var1_spk_set(:,bins_removed_var1) = nan;
            var2_array_set(bins_removed_var2,:) = nan;
            var2_spk_set(bins_removed_var2,:) = nan;

            % Limit output map set so as not to overwhelm file size
            numIterLlh = size(var1_array_set,1);
            maxNumSmooth = 10;
            if pick <= maxNumSmooth
                smoothset = 1:pick;
            else
                smoothset = [1 pick-maxNumSmooth+2:pick];
            end
            numIterSm = length(smoothset);
            smoothpick = numIterSm;

            var1_array_smoothset = var1_array_set(smoothset,:);
            var1_spk_smoothset = var1_spk_set(smoothset,:);
            var2_array_smoothset = var2_array_set(:,smoothset);
            var2_spk_smoothset = var2_spk_set(:,smoothset);

            %%%%%%%%%%%%%%%



            % NOTE: Currently unable to make Kian Wei's batch smoothing code
            % work. 
            % So for now just looping

            disp('Smoothing ...');

            %% Smoothing
            
            for cc = 1:2 % for each var
                
                if cc == 1
                    arr = var1_array_smoothset;
                    spk = var1_spk_smoothset;
                    dur = var1_durations1;
                    if size(arr,1) < size(arr,2)
                        arr = arr'; spk = spk'; dur = dur';
                    end
                else
                    arr = var2_array_smoothset;
                    spk = var2_spk_smoothset;
                    dur = nansum(full_durations1,2);
                    if size(arr,1) < size(arr,2)
                        arr = arr'; spk = spk';
                    end
                end
                
                switch spatialvarpairs{combi}{cc}
                    
                    case 'place'
                        placesmooth = 'adaptive';
                        switch placesmooth

                            case 'adaptive'

                                % Adaptive Smoothing
                                maps_sm = nan(size(arr,2),size(arr,1));
                                dur_sm = nan(size(arr,2),size(arr,1));
                                durG = lineartogrid(dur,'place',[Args.GridSteps Args.GridSteps]);
                                durG = durG{1};
                                for ii = 1:numIterSm
                                    spkG = lineartogrid(spk(:,ii),'place',[Args.GridSteps Args.GridSteps]);
                                    spkG = spkG{1}; 
                                    [maps_adsmG,~,dur_adsmG] = adsmooth(durG,spkG,Args.AlphaPlace);
                                    dur_adsmG(isnan(dur_adsmG)) = 0;
                                    maps_sm(ii,:) = gridtolinear({maps_adsmG},'place',[Args.GridSteps Args.GridSteps]);
                                    dur_sm(ii,:) = gridtolinear({dur_adsmG},'place',[Args.GridSteps Args.GridSteps]);
                                end
                                % Adaptive SIC
                                crit_out = skaggs_sic(maps_sm',dur_sm');

                            case 'boxcar'

                                arrG = cell2mat(lineartogrid(arr,'place',[Args.GridSteps Args.GridSteps]));
                                durG = repmat(durG,1,1,size(arr,1)); %% ???
                                unvis = ~(durG > 0);

                                % Boxcar smoothing
                                maps_bcsmG = smooth(arrG,5,unvis,'boxcar');
                                dur_bcsmG = smooth(durG,5,unvis,'boxcar');
                                maps_sm = (gridtolinear({maps_bcsmG},'place',[Args.GridSteps Args.GridSteps]))';
                                dur_bcsm = (gridtolinear({dur_bcsmG},'place',[Args.GridSteps Args.GridSteps]))';
                                dur_bcsm(isnan(dur_bcsm)) = 0;
                                % Boxcar SIC
                                crit_out = skaggs_sic(maps_sm',dur_bcsm');

                            case 'disk'

                                % Disk smoothing
                                maps_dksmG = smooth(arrG,5,unvis,'disk');
                                dur_dksmG = smooth(durG,5,unvis,'disk');
                                maps_sm = (gridtolinear({maps_dksmG},'place',[Args.GridSteps Args.GridSteps]))';
                                dur_dksm = (gridtolinear({dur_dksmG},'place',[Args.GridSteps Args.GridSteps]))';
                                dur_dksm(isnan(dur_dksm)) = 0;
                                % Disk SIC
                                crit_out = skaggs_sic(maps_sm',dur_dksm');
                        end
                        %         % ISE 
                        %         lambda_i = var1_maps_adsm; 
                        %         var1_ise_out = ise(lambda_i(1,:), lambda_i(2:end,:), Args.GridSteps, Args.GridSteps);
                    
                    case 'view'
                        viewsmooth = 'adaptive';

                        maps_sm = nan(numIterSm,size(arr,1));
                        dur_sm = nan(numIterSm,size(arr,1));
            %             var2_maps_bcsm = nan(numIterSm,size(var2_array_smoothset,1));
            %             var2_dur_bcsm = nan(numIterSm,size(var2_array_smoothset,1));
            %             var2_maps_dksm = nan(numIterSm,size(var2_array_smoothset,1));
            %             var2_dur_dksm = nan(numIterSm,size(var2_array_smoothset,1));

                        for ii = 1:numIterSm % For each iteration

                            % Assign linear bin to grid bin - left to right, bottom to top
                            durG = lineartogrid(dur,'view',viewbinDepths);
                            spkG = lineartogrid(spk(:,ii),'view',viewbinDepths);
                            ratesG = lineartogrid(arr(:,ii),'view',viewbinDepths);

                            % Pad sv map with 5 extra rows
                            n = 5;
                            padpillar = false;
                            [emptyfloorref_pad,~] = padsvmap(n,durG,viewSections,padpillar);
                            padpillar = true;
                            [durGpad,retrievemap] = padsvmap(n,durG,viewSections,padpillar);
                            [spkGpad,~] = padsvmap(n,spkG,viewSections,padpillar);
                            [ratesGpad,~] = padsvmap(n,ratesG,viewSections,padpillar);

                            if strcmp(viewsmooth,'adaptive')

                                % Adaptive smooth
                                maps_adsmGpad = cell(size(durGpad));
                                dur_adsmGpad = cell(size(durGpad));
                                for jj = 1:size(viewbinDepths,1)
                                    if jj == 1 || jj == 2
                                        maps_adsmGpad{jj} = spkGpad{jj}/durGpad{jj};
                                        dur_adsmGpad{jj} = durGpad{jj};
                                    else
                                        [maps_adsmGpad{jj},spk_adsmG,dur_adsmGpad{jj}] = adsmooth(durGpad{jj},spkGpad{jj},Args.AlphaView);
                                    end
                                end
                                % Unpad smoothed map
                                maps_adsmG = unpadsvmap(maps_adsmGpad,retrievemap,durG);
                                dur_adsmG = unpadsvmap(dur_adsmGpad,retrievemap,durG);
                                % Convert grid map back to linear sv map
                                maps_sm(ii,:) = gridtolinear(maps_adsmG,'view',viewbinDepths);
                                dur_sm(ii,:) = gridtolinear(dur_adsmG,'view',viewbinDepths);
                                dur_sm(isnan(dur_sm)) = 0;
                                % Adaptive SIC 
                                crit_out = skaggs_sic(maps_sm',dur_sm');
                            else 
                                % Boxcar/disk smooth
                                maps_bcsmGpad = cell(size(durGpad));
                                dur_bcsmGpad = cell(size(durGpad));
                                maps_dksmGpad = cell(size(durGpad));
                                dur_dksmGpad = cell(size(durGpad));
                                for jj = 1:size(viewbinDepths,1)
                                    if jj == 1 || jj == 2
                                        maps_bcsmGpad{jj} = ratesGpad{jj};
                                        maps_dksmGpad{jj} = ratesGpad{jj};
                                        dur_bcsmGpad{jj} = durGpad{jj};
                                        dur_dksmGpad{jj} = durGpad{jj};
                                    else
                                        unvis = ~(emptyfloorref_pad{jj}>0) | isnan(ratesGpad{jj});
                                        if strcmp(viewsmooth,'boxcar')
                                            % Boxcar smoothing
                                            maps_bcsmGpad{jj} = smooth(ratesGpad{jj},5,unvis,'boxcar');
                                            dur_bcsmGpad{jj} = smooth(durGpad{jj},5,unvis,'boxcar');
                                        elseif strcmp(viewsmooth,'disk')
                                            % Disk smoothing
                                            maps_dksmGpad{jj} = smooth(ratesGpad{jj},5,unvis,'disk');
                                            dur_dksmGpad{jj} = smooth(durGpad{jj},5,unvis,'disk');
                                        end
                                    end
                                end

                                if strcmp(viewsmooth,'boxcar')
                                    % Unpad smoothed map
                                    maps_bcsmG = unpadsvmap(maps_bcsmGpad,retrievemap,durG);
                                    dur_bcsmG = unpadsvmap(dur_bcsmGpad,retrievemap,durG);
                                    % Convert grid map back to linear sv map
                                    maps_sm(ii,:) = gridtolinear(maps_bcsmG,'view',viewbinDepths);
                                    dur_sm(ii,:) = gridtolinear(dur_bcsmG,'view',viewbinDepths);
                                    dur_sm(isnan(var2_dur_bcsm)) = 0;
                                    % Boxcar SIC
                                    crit_out = skaggs_sic(var2_maps_bcsm',dur_sm');
                                elseif strcmp(viewsmooth,'disk')
                                    % Unpad smoothed map
                                    maps_dksmG = unpadsvmap(maps_dksmGpad,retrievemap,durG);
                                    dur_dksmG = unpadsvmap(dur_dksmGpad,retrievemap,durG);
                                    % Convert grid map back to linear sv map
                                    maps_sm(ii,:) = gridtolinear(maps_dksmG,'view',viewbinDepths);
                                    dur_sm(ii,:) = gridtolinear(dur_dksmG,'view',viewbinDepths);
                                    dur_sm(isnan(var2_dur_dksm)) = 0;
                                    % Disk SIC
                                    crit_out = skaggs_sic(var2_maps_dksm',dur_sm');
                                end
                            end

                        end
                        
                    case 'headdirection'
                        
                        n = 5;
                        headdirectionsmooth = 'boxcar';
                        % Smooth
                        [maps_sm] = smoothdir(arr,n,pv.data.headdirectionbins);
                        [dur_sm] = smoothdir(repmat(dur,1,numIterSm),n,pv.data.headdirectionbins);
                        maps_sm = maps_sm';
                        dur_sm = dur_sm';

                        % Rayleigh vector
                        binSize=(pi*2)/length(maps_sm(1,:));
                        binAngles=(0:binSize:( (359.5/360)*2*pi )) + binSize/2;
                        binWeights=maps_sm./(max(maps_sm,[],2));
                        S=nansum( sin(binAngles).*binWeights , 2);
                        C=nansum( cos(binAngles).*binWeights , 2);
                        R=sqrt(S.^2+C.^2);
                        meanR=R./nansum(binWeights,2);
                        crit_out = meanR';
                end
            
                if cc == 1
                    var1_maps_sm = maps_sm;
                    var1_dur_sm = dur_sm;
                    var1_crit_out = crit_out;
                else
                    var2_maps_sm = maps_sm;
                    var2_dur_sm = dur_sm;
                    var2_crit_out = crit_out;
                end
                    
            end

            % Store output data
            if repeat == 1
                
                data.(combiname).var1 = Var1;
                data.(combiname).var2 = Var2;
                data.(combiname).([lower(Var1) 'smooth']) = eval([lower(Var1) 'smooth']);
                data.(combiname).([lower(Var2) 'smooth']) = eval([lower(Var2) 'smooth']);
                data.(combiname).([lower(Var1) 'binDepths']) = eval([Var1 'binDepths']);
                data.(combiname).([lower(Var2) 'binDepths']) = eval([Var2 'binDepths']);

                data.(combiname).(['maps_dist_' lower(Var1(1))]) = var1_array_pred;
                data.(combiname).(['maps_dist_' lower(Var1(1)) '_adsm']) = var1_array_pred_adsm;
                data.(combiname).(['maps_dist_' lower(Var2(1))]) = var2_array_pred';
                data.(combiname).(['maps_dist_' lower(Var2(1)) '_adsm']) = var2_array_pred_adsm';
                data.(combiname).(['distratio_' lower(Var1(1))]) = dr_var1;
                data.(combiname).(['distratio_' lower(Var2(1))]) = dr_var2;
                data.(combiname).(['distcorr_' lower(Var1(1))]) = dcorr_var1;
                data.(combiname).(['distcorr_' lower(Var2(1))]) = dcorr_var2;

                data.(combiname).(['convergewith' lower(Var2(1))]) = convergewithvar2;
                data.(combiname).(['convergewith' lower(Var1(1))]) = convergewithvar1;
                data.(combiname).NumIterLlh = numIterLlh;
                data.(combiname).smoothpick = smoothpick;
                data.(combiname).llh = llh_vec;
                data.(combiname).llhpick = pick;
                data.(combiname).llhpicklabel = picklabel;
                data.(combiname).(['maps_raw_corr' lower(Var1(1))]) = var1_array_set(pick,:);
                data.(combiname).(['maps_raw_corr' lower(Var2(1))]) = var2_array_set(:,pick)';
                data.(combiname).(['maps_raw_corr' lower(Var1(1)) 'set']) = var1_array_smoothset;
                data.(combiname).(['maps_raw_corr' lower(Var2(1)) 'set']) = var2_array_smoothset';
    %             data.ratesplaceMethod2 = rates_place;
    %             data.ratesviewMethod2 = rates_view;

                data.(combiname).(['maps_sm_corr' lower(Var1(1))]) = var1_maps_sm(smoothpick,:);
                data.(combiname).(['maps_sm_corr' lower(Var1(1)) 'set']) = var1_maps_sm;
                data.(combiname).(['dur_sm_corr' lower(Var1(1))]) = var1_dur_sm(smoothpick,:);
                data.(combiname).(['crit_sm_corr' lower(Var1(1))]) = var1_crit_out(smoothpick);
                data.(combiname).(['crit_sm_corr' lower(Var1(1)) 'set']) = var1_crit_out';

                data.(combiname).(['maps_sm_corr' lower(Var2(1))]) = var2_maps_sm(smoothpick,:);
                data.(combiname).(['maps_sm_corr' lower(Var2(1)) 'set']) = var2_maps_sm;
                data.(combiname).(['dur_sm_corr' lower(Var2(1))]) = var2_dur_sm(smoothpick,:);
                data.(combiname).(['crit_sm_corr' lower(Var2(1))]) = var2_crit_out(smoothpick);
                data.(combiname).(['crit_sm_corr' lower(Var2(1)) 'set']) = var2_crit_out';

                data.(combiname).covmat = covmat;
                data.(combiname).covmat_norm = covmat_norm;
                data.(combiname).l1norm = l1norm;
                data.(combiname).l2norm = l2norm;
            elseif repeat == 2

                data.(combiname).(['maps_dist_' lower(Var1(1)) '1']) = var1_array_pred;
                data.(combiname).(['maps_dist_' lower(Var1(1)) '_adsm' '1']) = var1_array_pred_adsm;
                data.(combiname).(['maps_dist_' lower(Var2(1)) '1']) = var2_array_pred';
                data.(combiname).(['maps_dist_' lower(Var2(1)) '_adsm' '1']) = var2_array_pred_adsm';
                data.(combiname).(['distratio_' lower(Var1(1)) '1']) = dr_var1;
                data.(combiname).(['distratio_' lower(Var2(1)) '1']) = dr_var2;
                data.(combiname).(['distcorr_' lower(Var1(1)) '1']) = dcorr_var1;
                data.(combiname).(['distcorr_' lower(Var2(1)) '1']) = dcorr_var2;

                data.(combiname).(['maps_sm_corr' lower(Var1(1)) '1']) = var1_maps_sm(smoothpick,:);
                data.(combiname).(['maps_sm_corr' lower(Var1(1)) 'set' '1']) = var1_maps_sm;
                data.(combiname).(['crit_sm_corr' lower(Var1(1)) '1']) = var1_crit_out(smoothpick);
                data.(combiname).(['crit_sm_corr' lower(Var1(1)) 'set' '1']) = var1_crit_out';

                data.(combiname).(['maps_sm_corr' lower(Var2(1)) '1']) = var2_maps_sm(smoothpick,:);
                data.(combiname).(['maps_sm_corr' lower(Var2(1)) 'set' '1']) = var2_maps_sm;
                data.(combiname).(['crit_sm_corr' lower(Var2(1)) '1']) = var2_crit_out(smoothpick);
                data.(combiname).(['crit_sm_corr' lower(Var2(1)) 'set' '1']) = var2_crit_out';
                
            elseif repeat == 3

                data.(combiname).(['maps_dist_' lower(Var1(1)) '2']) = var1_array_pred;
                data.(combiname).(['maps_dist_' lower(Var1(1)) '_adsm' '2']) = var1_array_pred_adsm;
                data.(combiname).(['maps_dist_' lower(Var2(1)) '2']) = var2_array_pred';
                data.(combiname).(['maps_dist_' lower(Var2(1)) '_adsm' '2']) = var2_array_pred_adsm';
                data.(combiname).(['distratio_' lower(Var1(1)) '2']) = dr_var1;
                data.(combiname).(['distratio_' lower(Var2(1)) '2']) = dr_var2;
                data.(combiname).(['distcorr_' lower(Var1(1)) '2']) = dcorr_var1;
                data.(combiname).(['distcorr_' lower(Var2(1)) '2']) = dcorr_var2;

                data.(combiname).(['maps_sm_corr' lower(Var1(1)) '2']) = var1_maps_sm(smoothpick,:);
                data.(combiname).(['maps_sm_corr' lower(Var1(1)) 'set' '2']) = var1_maps_sm;
                data.(combiname).(['crit_sm_corr' lower(Var1(1)) '2']) = var1_crit_out(smoothpick);
                data.(combiname).(['crit_sm_corr' lower(Var1(1)) 'set' '2']) = var1_crit_out';

                data.(combiname).(['maps_sm_corr' lower(Var2(1)) '2']) = var2_maps_sm(smoothpick,:);
                data.(combiname).(['maps_sm_corr' lower(Var2(1)) 'set' '2']) = var2_maps_sm;
                data.(combiname).(['crit_sm_corr' lower(Var2(1)) '2']) = var2_crit_out(smoothpick);
                data.(combiname).(['crit_sm_corr' lower(Var2(1)) 'set' '2']) = var2_crit_out';
                
            end
        end
        
    end
        
    %% SAVE   
    
    % create nptdata so we can inherit from it    
    data.spatialvars = spatialvars;
    data.spatialvarpairs = spatialvarpairs;
    data.gridSteps = Args.GridSteps;
%     data.binDepths = viewbinDepths;
    data.placebins = pv.data.placebins;
    data.viewbins = pv.data.viewbins;
    data.headdirectionbins = pv.data.headdirectionbins;
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


