function [obj, varargout] = vmcorr(varargin)

%%%% Method 1


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
				'GridSteps',40, 'pix',1,...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
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

    
    % Numbers of pair permutations given input numbers of spatial variables
    spatialvars = {'place','view','headdirection'};
    spatialvarpairs = {{'place','view'}, {'place','headdirection'},{'headdirection','view'}}; % view must be in var2 position because of backfilling spikes later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ori = pwd;
    data.origin = {pwd}; 
%     pv = vmpv('auto', varargin{:});
        % Temporary workaround so that we can use file named '1vmpv.mat'
        cd ..; cd ..; cd ..;
        pv = load([num2str(Args.pix) 'vmpv.mat']);
        pv = pv.pv;
    cd(ori);

    % Set up view binning for smoothing
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


    spiketrain = load(Args.RequiredFile);
    stc = pv.data.sessionTimeC; % time, place, hd, view

    % Bin spikes % For Method 1
    spiketimes = spiketrain.timestamps/1000; % now in seconds
    maxTime = pv.data.rplmaxtime;
    spiketimes(spiketimes>maxTime) = [];
    spiketimes(spiketimes(:,1) < stc(1,1),:) = []; 
    spiketimes = spiketimes';
    stc(:,5) = [diff(stc(:,1)); 0]; % column 5 as duration
    binned = histcounts(spiketimes, stc(:,1))'; % Method 1 of 2 for binning spikes
    stc(:,6) = [binned; 0];

    % For view only, backfill time and spikes
    stc_view = stc;
    % back-filling spikes for view bins that occupy the same time bin
    stc_view(stc_view(:,6)==0,6) = nan;
    stc_view(:,7) = stc_view(:,5)~=0;
    stc_view(isnan(stc_view(:,6)) & stc_view(:,7), 6) = 0;
    stc_view(:,7) = [];
    stc_view(:,6) = fillmissing(stc_view(:,6), 'next');
    stc_view(isnan(stc_view(:,6)),6) = 0;
    % back-filling duration for view bins that occupy the same time bin
    stc_lasttime = stc_view(end,5); % Make sure if last duration sample is zero, it remains zero, not nan
    stc_view(stc_view(:,5)==0,5) = nan;
    stc_view(end,5) = stc_lasttime;
    stc_view(:,5) = fillmissing(stc_view(:,5), 'next'); % [timestamp place hd view dur spk]

    % Run for full session and each half
    for repeat = 1:3 % 1 = full trial, 2 = 1st half, 3 = 2nd half

        if repeat == 1
            disp('Full session:');
        elseif repeat == 2
            disp('1st half:');
        elseif repeat == 3
            disp('2nd half:');
        end
        
        for combi = 1:size(spatialvarpairs,2) % For each pair of spatial vars

            Var1 = spatialvarpairs{combi}{1};
            Var2 = spatialvarpairs{combi}{2};
            combiname = [lower(Var1(1)) lower(Var2(1))];
            disp([' ' combiname ' - ']);

            % filtering rows from sessionTimeC
            disp('    Filtering ...');
            conditions = ones(size(stc,1),1);
            if repeat == 2
                conditions = conditions & (pv.data.halving_markers==1);
            elseif repeat == 3
                conditions = conditions & (pv.data.halving_markers==2);
            end
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
                        bins_sieved = 1:Args.DirSteps; % Currently don't have a threshold for min number of obs for head dir bins
            %             bins_sieved = pv.data.place_good_bins; 
                        conditions = conditions & (pv.data.pv_good_rows); % Make sure minobs take into account both place and view
                end
                switch Var2
                    case 'place'
                        bins_sieved_var2 = pv.data.place_good_bins;
                        bins_removed_var2 = setdiff(1:size(pv.data.place_intervals_count,1),bins_sieved_var2);
                    case 'view'
                        bins_sieved_var2 = pv.data.view_good_bins;
                        bins_removed_var2 = setdiff(1:size(pv.data.view_intervals_count,1),bins_sieved_var2);
                    case 'headdirection' %%FIX
                        bins_sieved = 1:Args.DirSteps; % Currently don't have a threshold for min number of obs for head dir bins
            %             bins_sieved = pv.data.place_good_bins; 
                        conditions = conditions & (pv.data.pv_good_rows); % Make sure minobs take into account both place and view
                end
                    conditions = conditions & (pv.data.pv_good_rows); % Make sure maps take into account both place and view filters
            else
                bins_sieved_var1 = 1:pv.data.([Var1 'bins']);
                bins_sieved_var2 = 1:pv.data.([Var2 'bins']);
                bins_removed_var1 = [];
                bins_removed_var2 = [];
            end
            
            % remove filtered rows
            stc_filt = stc(conditions == 1,:);
            stc_view_filt = stc_view(conditions == 1,:);
            % Consider only samples where all spatial variables are
            % sample (in line with the way vmpc/sv/hd are filtered)
            % (nan for nonsampled view, 0 for nonsampled place, 0 for nonsampled hd)
            reject = ~(stc_filt(:,2)>0) | ~(stc_filt(:,3)>0) | isnan(stc_filt(:,4));
            stc_filt(reject,:) = [];
            % checkind = find(conditions == 1);
            % checkind(reject) = [];
            reject_view = ~(stc_view_filt(:,2)>0) | ~(stc_view_filt(:,3)>0) | isnan(stc_view_filt(:,4));
            stc_view_filt(reject_view,:) = [];

            %% Consolidate pv array into view by place bin array

            stc_colnames = {'time','place','headdirection','view','dur','spk'};
            keepcol = [find(ismember(stc_colnames,spatialvarpairs{combi}{1})) ...
                find(ismember(stc_colnames,spatialvarpairs{combi}{2}))];

            % Initialise variables
            disp('    Binning ...');

            % Pad last row so that accumarray can work
            stc_filt_append = [stc_filt;[0 1600 60 5122 0 0]];
            stc_view_filt_append = [stc_view_filt;[0 1600 60 5122 0 0]];

            % Consolidate durations and spikes
            var1_durations1 = accumarray(stc_filt_append(:,keepcol(1)),stc_filt_append(:,5),[],[],0);
            var1_spikes1 = accumarray(stc_filt_append(:,keepcol(1)),stc_filt_append(:,6),[],[],0);
            var1_durations1 = reshape(var1_durations1,1,length(var1_durations1));
            var1_spikes1 = reshape(var1_spikes1,1,length(var1_spikes1));
            var2_durations1 = accumarray(stc_filt_append(:,keepcol(2)),stc_filt_append(:,5),[],[],0);
            var2_spikes1 = accumarray(stc_filt_append(:,keepcol(2)),stc_filt_append(:,6),[],[],0);
            var2_durations1 = reshape(var2_durations1,length(var2_durations1),1);
            var2_spikes1 = reshape(var2_spikes1,length(var2_spikes1),1);
            if strcmp(Var2,'view')
                full_durations1 = accumarray(stc_view_filt_append(:,[keepcol(2) keepcol(1)]), stc_view_filt_append(:,5), [],[],0);
                full_spikes1 = accumarray(stc_view_filt_append(:,[keepcol(2) keepcol(1)]), stc_view_filt_append(:,6), [],[],0);
            else
                full_durations1 = accumarray(stc_filt_append(:,[keepcol(2) keepcol(1)]), stc_filt_append(:,5), [],[],0);
                full_spikes1 = accumarray(stc_filt_append(:,[keepcol(2) keepcol(1)]), stc_filt_append(:,6), [],[],0);
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
            if strcmp(Var2,'view')
                var2_array_orig = sum(full_spikes1,2)./sum(full_durations1,2);
            else
                var2_array_orig = var2_spikes1./var2_durations1;
            end

            % Debug
            if repeat == 1
                cd(['FiltVel/' num2str(Args.pix) 'px']);
                switch Var1
                    case 'place'
                        Var1obj = load('vmpc.mat');
                        Var1obj = Var1obj.vmp;
                    case 'view'
                        Var1obj = load('vmsv.mat');
                        Var1obj = Var1obj.vms;
                    case 'headdirection'
                        Var1obj = load('vmhd.mat');
                        Var1obj = Var1obj.vmd;
                        % cd(ori);
                        % Var1obj = vmhd('auto','redo','NumShuffles',2);
                        % % Var1obj = Var1obj.data;
                end
                switch Var2
                    case 'place'
                        Var2obj = load('vmpc.mat');
                        Var2obj = Var2obj.vmp;
                    case 'view'
                        Var2obj = load('vmsv.mat');
                        Var2obj = Var2obj.vms;
                    case 'headdirection'
                        Var2obj = load('vmhd.mat');
                        Var2obj = Var2obj.vmd;
                        % cd(ori);
                        % Var2obj = vmhd('auto','redo','NumShuffles',2);
                        % Var2obj = Var2obj.data;
                end
                cd(ori);
                diffs1 = reshape(var1_array_orig,length(var1_array_orig),1) - reshape(Var1obj.data.maps_raw,length(Var1obj.data.maps_raw),1);
                diffs2 = reshape(var2_array_orig,length(var2_array_orig),1) - reshape(Var2obj.data.maps_raw,length(Var2obj.data.maps_raw),1);
                errorcode = {};
                if (sum(diffs1,[],'omitnan') ~= 0 || sum(diffs2,[],'omitnan') ~= 0) 
                    disp('error in binning spikes');
                    if sum(diffs1,[],'omitnan') ~= 0
                        Var1spk = var1_spikes1';
                        Var1spk(:,2) = Var1obj.data.spk_raw;
                        Var1spk(:,3) = Var1spk(:,1) - Var1spk(:,2);
                        Var1dur = var1_durations1';
                        Var1dur(:,2) = Var1obj.data.dur_raw;
                        Var1dur(:,3) = Var1dur(:,1) - Var1dur(:,2);
                        if sum(abs(Var1spk(:,3))) > 0
                            disp('Var1 spk error = '); 
                            disp(Var1spk(Var1spk(:,3)~=0,3));
                            errorcode{end+1,1} = ['Var1spk ' num2str(sum(Var1spk(:,3)~=0)) 'pt ' num2str(sum(abs(Var1spk(:,3))))];
                        end
                        if sum(abs(Var1dur(:,3))) > 0
                            disp('Var1 dur error = ')
                            disp(Var1dur(Var1dur(:,3)~=0,3));
                            errorcode{end+1,1} = ['Var1dur ' num2str(sum(Var1dur(:,3)~=0)) 'pt ' num2str(sum(abs(Var1dur(:,3))))];
                        end
                    end
                    if sum(diffs2,[],'omitnan') ~= 0
                        if strcmp(Var2,'view)')
                            Var2spk = sum(full_spikes1,2);
                            Var2dur = sum(full_durations1,2);
                        else 
                            Var2spk = var2_spikes1;
                            Var2dur = var2_durations1;
                        end
                        Var2spk(:,2) = Var2obj.data.spk_raw;
                        Var2spk(:,3) = Var2spk(:,1) - Var2spk(:,2);
                        Var2dur(:,2) = Var2obj.data.dur_raw;
                        Var2dur(:,3) = Var2dur(:,1) - Var2dur(:,2);
                        if sum(abs(Var2spk(:,3))) > 0
                            disp('Var2 spk error = '); 
                            % disp(Var2spk(Var2spk(:,3)~=0,3));
                            errorcode{end+1,1} = ['Var2spk ' num2str(sum(Var2spk(:,3)~=0)) 'pt ' num2str(sum(abs(Var2spk(:,3))))];
                        end
                        if sum(abs(Var2dur(:,3))) > 0
                            disp('Var2 dur error = ')
                            % disp(Var2dur(Var2dur(:,3)~=0,3));
                            errorcode{end+1,1} = ['Var2dur ' num2str(sum(Var2dur(:,3)~=0)) 'pt ' num2str(sum(abs(Var2dur(:,3))))];
                        end
                    end
                end
                % % debug stc
                % % Are there any instances of 
                % stc_test = stc;
                % stc_test(:,7) = false; % mark change of place/hd bin
                % stc_test(:,8) = false; % mark suspicious rows
                % for tt = 2:size(stc_test,1)
                %     % Mark rows where there is a change of place/hd bin
                %     if stc_test(tt,2) ~= stc_test(tt-1,2) || stc_test(tt,3) ~= stc_test(tt-1,3)
                %         stc_test(tt,7) = true;
                %         if stc_test(tt,5) == 0
                %             % find next time stamp
                %             count = tt;
                %             if count < size(stc_test,1)
                %                 while stc_test(count+1,1) == stc_test(count,1) && count < size(stc_test,1)
                %                     count = count + 1;
                %                 end
                %             end
                %             % mark rows where there is no duration allocated to
                %             % this pair of place/hd bin
                %             if stc_test(count,5) == 0
                %                 stc_test(tt,8) = true;
                %             end
                %         end
                %     end
                % end
            end
            
            %% Compute covariance matrix
            if repeat == 1
                disp('    Computing covariance matrix ...');
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

            disp('    Distributive hypothesis testing ...');
            % Removing hint and cue effects
            var2_array_orig_trun = var2_array_orig;
            full_durations1_trun = full_durations1;
            % if strcmp(Var2,'view')
            %      % disregard cue and hint bins
            %     var2_array_orig_trun(1:2) = nan;
            %     full_durations1_trun(1:2,:) = 0;
            % end

            % Predicted rate as a function of var2
            % If predicted map of var2 based on behavior/sampling conforms 
            % closely to observed map of var2, then observed var2
            % selectivity is artefactual. 
            topterm = sum(var1_array_orig.*full_durations1_trun,2,'omitnan');
            bottomterm = sum(full_durations1_trun,2,'omitnan');
            var2_array_pred = topterm./bottomterm; % raw map

            % Smoothing predicted maps
            switch Var2
                case 'view'
                    % Adaptive smooth 
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
                    dur_pred_adsmGpad = cell(size(viewbinDepths,1),1);
                    for gg = 1:size(viewbinDepths,1)
                        if gg == 1 || gg == 2
                            sv_array_pred_adsmGpad{gg} = toptermG{gg}/bottomtermG{gg};
                            dur_pred_adsmGpad{gg} = bottomtermG{gg};
                        else
                            [sv_array_pred_adsmGpad{gg},~,dur_pred_adsmGpad{gg}] = adsmooth(bottomtermGpad{gg},toptermGpad{gg},Args.AlphaView);
                        end
                    end
                    % Unpad
                    [sv_array_pred_adsmG] = unpadsvmap(sv_array_pred_adsmGpad,retrievemap);
                    [dur2_pred_adsmG] = unpadsvmap(dur_pred_adsmGpad,retrievemap);
                    % Linearise
                    var2_array_pred_adsm = gridtolinear(sv_array_pred_adsmG,'view',viewbinDepths);
                    dur2_pred_adsm = gridtolinear(dur2_pred_adsmG,'view',viewbinDepths);
                    dur2_pred_adsm(isnan(dur2_pred_adsm)) = 0;
                case 'headdirection'
                    n = 5; % Boxcar window of 5
                    [var2_array_pred_adsm] = smoothdir(var2_array_pred,n,pv.data.([Var2 'bins']));
                case 'place'
                    % Adaptive smooth
                    toptermG = cell2mat(lineartogrid(topterm','place',[Args.GridSteps Args.GridSteps]));
                    bottomtermG = cell2mat(lineartogrid(bottomterm','place',[Args.GridSteps Args.GridSteps]));
                    [p_array_pred_adsmG,~,dur2_pred_adsmG] = adsmooth(bottomtermG,toptermG,Args.AlphaPlace);
                    var2_array_pred_adsm = gridtolinear({p_array_pred_adsmG},'place',[Args.GridSteps Args.GridSteps])';
                    dur2_pred_adsm = gridtolinear({dur2_pred_adsmG},'place',[Args.GridSteps Args.GridSteps]);
                    dur2_pred_adsm(isnan(dur2_pred_adsm)) = 0;
            end
            
            % SIC for predicted/distributed maps
            if ~strcmp(Var2,'headdirection')
                sic2_pred_adsm = skaggs_sic(var2_array_pred_adsm,dur2_pred_adsm);
            else
                % Rayleigh vector
                binSize=(pi*2)/length(var2_array_pred_adsm(1,:));
                binAngles=(0:binSize:( (359.5/360)*2*pi )) + binSize/2;
                binWeights=var2_array_pred_adsm./(max(var2_array_pred_adsm,[],2));
                S=sum( sin(binAngles).*binWeights , 2);
                C=sum( cos(binAngles).*binWeights , 2);
                R=sqrt(S.^2+C.^2);
                sic2_pred_adsm=R./sum(binWeights,2);
            end

            % If place cell firing is only mod by location, and view influence is attributable only to inhomogeous sampling, 
            % then dr = 0. If dr is significant, place cell is view-modulated.
            % ratio measures goodness of fit between obs and pred, 0 =
            % perfect prediction. 
            ratio = log((1+var2_array_orig_trun)./(1+var2_array_pred)); % Muller 1994 
    %         ratio = log(1+sv_array_orig_trun)./(1+sv_array_pred); % Cacucci 2004
            dr_var2 = sum(abs(ratio),[],'omitnan')/sum(bottomterm>0); % influence of var2 on var1
            % Correlation between orig and predicted view maps
            vis = ~isnan(var2_array_orig_trun);
            dcorr_var2 = corr2(var2_array_orig_trun(vis),var2_array_pred(vis));

            % Predicted rate as a function of var1
            topterm = sum(var2_array_orig_trun.*full_durations1_trun,1,'omitnan');
            bottomterm = sum(full_durations1_trun,1,'omitnan');
            var1_array_pred = topterm./bottomterm; % raw map

            % Adaptive smooth predicted var1 map
            switch Var1
                case 'place' 
                    toptermG = cell2mat(lineartogrid(topterm','place',[Args.GridSteps Args.GridSteps]));
                    bottomtermG = cell2mat(lineartogrid(bottomterm','place',[Args.GridSteps Args.GridSteps]));
                    [p_array_pred_adsmG,~,dur1_pred_adsmG] = adsmooth(bottomtermG,toptermG,Args.AlphaPlace);
                    var1_array_pred_adsm = gridtolinear({p_array_pred_adsmG},'place',[Args.GridSteps Args.GridSteps])';
                    dur1_pred_adsm = gridtolinear({dur1_pred_adsmG},'place',[Args.GridSteps Args.GridSteps]);
                    dur1_pred_adsm(isnan(dur1_pred_adsm)) = 0;
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
                    dur_pred_adsmGpad = cell(size(viewbinDepths,1),1);
                    for gg = 1:size(viewbinDepths,1)
                        if gg == 1 || gg == 2
                            sv_array_pred_adsmGpad{gg} = toptermG{gg}/bottomtermG{gg};
                            dur_pred_adsmGpad{gg} = bottomtermG{gg};
                        else
                            [sv_array_pred_adsmGpad{gg},~,dur_pred_adsmGpad{gg}] = adsmooth(bottomtermGpad{gg},toptermGpad{gg},Args.AlphaView);
                        end
                    end
                    % Unpad
                    [sv_array_pred_adsmG] = unpadsvmap(sv_array_pred_adsmGpad,retrievemap);
                    [dur1_pred_adsmG] = unpadsvmap(dur_pred_adsmGpad,retrievemap);
                    % Linearise
                    var1_array_pred_adsm = gridtolinear(sv_array_pred_adsmG,'view',viewbinDepths);
                    dur1_pred_adsm = gridtolinear(dur1_pred_adsmG,'view',viewbinDepths);
                    dur1_pred_adsm(isnan(dur1_pred_adsm)) = 0;
            end
            % SIC for predicted/distributed maps
            if ~strcmp(Var1,'headdirection')
                sic1_pred_adsm = skaggs_sic(var1_array_pred_adsm,dur1_pred_adsm);
            else
                % Rayleigh vector
                binSize=(pi*2)/length(var1_array_pred_adsm(1,:));
                binAngles=(0:binSize:( (359.5/360)*2*pi )) + binSize/2;
                binWeights=var1_array_pred_adsm./(max(var1_array_pred_adsm,[],2));
                S=sum( sin(binAngles).*binWeights , 2);
                C=sum( cos(binAngles).*binWeights , 2);
                R=sqrt(S.^2+C.^2);
                sic1_pred_adsm=R./sum(binWeights,2);
            end

            % If view cell firing is only mod by view location, and place influence is attributable only to inhomogeous sampling, 
            % then dr = 0. If dr is significant, view cell is place-modulated.
            ratio = log((1+var1_array_orig)./(1+var1_array_pred)); % Muller 1994
    %         ratio = log(1+p_array_orig)./(1+p_array_pred); % Cacucci 2004 
            dr_var1 = sum(abs(ratio),[],'omitnan')/sum(bottomterm>0); % influence of var1 on var2
            % Correlation between orig and predicted place maps
            vis = ~isnan(var1_array_orig);
            dcorr_var1 = corr2(var1_array_orig(vis)',var1_array_pred(vis)');


            %% Establish independence of place and view maps
            disp('    Maximising llh ...');
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
                    pidjtij = var1_array_set(ii,:).*full_durations1.*var2_array_set(:,ii);
                    term1 = full_spikes1.*log(pidjtij);
                    term3 = log(factorial(var2_spikes_temp));
                    midterm = term1 - pidjtij - term3;
                    llh = sum(sum(midterm,1,'omitnan'),2,'omitnan');
                    llh_vec(end+1,1) = llh;
                end
                if isinf(llh) % Evaluate the 2nd llh in the set since that is the starting point of the iterations
                    pick = 2;
                    if max(max(var2_spikes_temp)) > 170
                        picklabel = ['spikecount' num2str(max(max(var2_spikes_temp)))];
                        disp(['        llh inf from start: spike count ' num2str(max(max(var2_spikes_temp)))]);
                    else
                        picklabel = 'inf1';
                        disp('        llh inf from start: other reasons');
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

                    pidjtij = var1_array.*full_durations1.*var2_array;
                    term1 = full_spikes1.*log(pidjtij);
                    term3 = log(factorial(var2_spikes_temp));
                    midterm = term1 - pidjtij - term3;
                    llh = sum(sum(midterm,1,'omitnan'),2,'omitnan');
                    % llh = sum( nansum(full_spikes1.*log(var1_array.*full_durations1.*var2_array)) - nansum(var1_array.*full_durations1.*var2_array) - nansum(log(factorial(var2_spikes_temp))) );
                    llh_vec(end+1,1) = llh;
    %                 llh = sum( nansum(full_spikes1.*log(p_array_set(ii,:).*full_durations1.*sv_array_set(:,ii))) - nansum(p_array_set(ii,:).*full_durations1.*sv_array_set(:,ii)) - nansum(log(factorial(view_spikes_temp))) );

                    if Args.UseIterLim % Converge if can find max within 1000 iterations

                        % Run iterations in sets of 100. 
                        % Stop if maximum of llh lands in second-to-last set
                        % For every set of 100 iterations
                        if size(llh_vec,1) >= 101 && rem(size(llh_vec,1),100) == 1

                            % Find max of last and next-to-last sets
                            maxlast = max(llh_vec(size(llh_vec,1)-100+1:end));
                            disp(['        iter ' num2str(size(llh_vec,1)) ' : currlocalmaxllh = ' num2str(maxlast)]);

                            if size(llh_vec,1) >= 201 
                                maxnextlast = max(llh_vec(size(llh_vec,1)-200+1:size(llh_vec,1)-100));
                                % If llh doesn't converge
                                if size(llh_vec,1) > Args.LLHIterLim
                                    % Switch to using flat view array to start iterating
                                    initialwith = Var1;
                                    reinitialise = true;
                                    disp(['        LLH did not converge with initial ' Var2 ' array. Try initialising with ' ...
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
                                        disp(['        Failed to converge to non-infinite values with initial ' Var2  ...
                                            'array. Try initialising with ' Var1 ' array']);
                                        break;
                                    else
                                        picklabel = num2str(pick);
                                        % display results
                                        disp(['        Converged with initial view array, pick = ' num2str(pick)]);
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
                            disp([        'Iter ' num2str(size(llh_vec,1))]);
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
                                disp(['        Failed to converge to non-infinite values with initial ' ...
                                    Var2 'array. Try initialising with ' Var1 ' array']);
                                break;
                            else
                                picklabel = num2str(pick);
                                % display results
                                disp(['        Converged with initial ' Var2 ' array, pick = ' ...
                                    num2str(pick) ' out of ' num2str(size(llh_vec,1))]);
                                % break out of loop
                                mlm = false;
                                break;
                            end
                        elseif isinf(llh)
                            % Skip to using different initial array
                            initialwith = Var1;
                            reinitialise = true;
                            disp(['        Failed to converge to non-infinite values with initial ' ...
                                Var2 ' array. Try initialising with ' Var1 ' array']);
                            break;
                        elseif size(llh_vec,1) > Args.LLHIterLim
                            % No convergence. Switch to using flat view array to start iterating
                            initialwith = Var1;
                            reinitialise = true;
                            disp(['        LLH did not converge with initial ' ...
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

                    pidjtij = var1_array.*full_durations1.*var2_array;
                    term1 = full_spikes1.*log(pidjtij);
                    term3 = log(factorial(var2_spikes_temp));
                    midterm = term1 - pidjtij - term3;
                    llh = sum(sum(midterm,1,'omitnan'),2,'omitnan');
                    % llh = sum( nansum(full_spikes1.*log(var1_array.*full_durations1.*var2_array)) - nansum(var1_array.*full_durations1.*var2_array) - nansum(log(factorial(var2_spikes_temp))) );
                    llh_vec(end+1,1) = llh;

                    if Args.UseIterLim % Converge if can find max within 1000 iterations

                        % Run iterations in sets of 100. 
                        % Stop if maximum of llh lands in second-to-last set
                        % For every set of 100 iterations
                        if size(llh_vec,1) >= 101 && rem(size(llh_vec,1),100) == 1

                            % Find max of last and next-to-last sets
                            maxlast = max(llh_vec(size(llh_vec,1)-100+1:end));
                            disp(['        iter ' num2str(size(llh_vec,1)) ' : currlocalmaxllh = ' num2str(maxlast)]);

                            if size(llh_vec,1) >= 201
                                maxnextlast = max(llh_vec(size(llh_vec,1)-200+1:size(llh_vec,1)-100));
                                % If llh doesn't converge
                                if size(llh_vec,1) > Args.LLHIterLim
                                    pick = size(llh_vec,1);
                                    picklabel = 'nonconverged';
                                    disp('        LLH did not converge with either initial arrays. Abandon.');
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
                                        disp('        Failed to converge to non-infinite values with either initial arrays. Abandon');
                                    else
                                        picklabel = num2str(pick);
                                        % display results
                                        disp(['        Converged with initial ' Var1 ' array, pick = ' num2str(pick)]);
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
                            disp(['        Iter ' num2str(size(llh_vec,1))]);
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
                                disp('        Failed to converge to non-infinite values with either array. Terminate');
                                mlm = false;
                            else
                                picklabel = num2str(pick);
                                % display results
                                disp(['        Converged with initial ' ...
                                    Var1 ' array, pick = ' num2str(pick) ' out of ' num2str(size(llh_vec,1))]);
                                % break out of loop
                                mlm = false;
                            end
                        elseif isinf(llh)
                            % Terminate
                            pick = size(llh_vec,1);
                            picklabel = 'inf';
                            disp('        Failed to converge to non-infinite values with either array. Terminate');
                            mlm = false;
                        elseif size(llh_vec,1) > Args.LLHIterLim % No convergence
                            pick = size(llh_vec,1);
                            picklabel = 'nonconverged';
                            disp('        LLH did not converge with either initial arrays. Terminate.');
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

            disp('    Smoothing ...');

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
                                maps_adsmG = unpadsvmap(maps_adsmGpad,retrievemap);
                                dur_adsmG = unpadsvmap(dur_adsmGpad,retrievemap);
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
                                    maps_bcsmG = unpadsvmap(maps_bcsmGpad,retrievemap);
                                    dur_bcsmG = unpadsvmap(dur_bcsmGpad,retrievemap);
                                    % Convert grid map back to linear sv map
                                    maps_sm(ii,:) = gridtolinear(maps_bcsmG,'view',viewbinDepths);
                                    dur_sm(ii,:) = gridtolinear(dur_bcsmG,'view',viewbinDepths);
                                    dur_sm(isnan(var2_dur_bcsm)) = 0;
                                    % Boxcar SIC
                                    crit_out = skaggs_sic(var2_maps_bcsm',dur_sm');
                                elseif strcmp(viewsmooth,'disk')
                                    % Unpad smoothed map
                                    maps_dksmG = unpadsvmap(maps_dksmGpad,retrievemap);
                                    dur_dksmG = unpadsvmap(dur_dksmGpad,retrievemap);
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
                data.(combiname).errorcode = errorcode;

                data.(combiname).(['maps_dist_' lower(Var1(1))]) = var1_array_pred; % predicted rate as function of var1, raw
                data.(combiname).(['maps_dist_' lower(Var1(1)) '_adsm']) = var1_array_pred_adsm; % predicted rate as function of var1, adaptive-smoothed
                data.(combiname).(['maps_dist_' lower(Var2(1))]) = var2_array_pred';
                data.(combiname).(['maps_dist_' lower(Var2(1)) '_adsm']) = var2_array_pred_adsm';
                data.(combiname).(['crit_sm_dist' lower(Var1(1))]) = sic1_pred_adsm;
                data.(combiname).(['crit_sm_dist' lower(Var2(1))]) = sic2_pred_adsm;
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


