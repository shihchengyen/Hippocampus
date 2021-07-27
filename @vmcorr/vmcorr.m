function [obj, varargout] = vmcorr(varargin)
%@vmpc Constructor function for vmpc class
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
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',0, ...
                'FRSIC',0, 'UseMedian',0, ...
                'NumFRBins',4,'AdaptiveSmooth',1, 'UseMinObs',0, 'ThresVel',1, 'UseAllTrials',1, 'StartOrig',0,'AlphaPlace',10000,'AlphaView',1000,...
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

    Args.NumShuffles = 0;
%     for repeat = 1:1 % 1 = full session, 2 = 1st half, 3 = 2nd half

%% Consolidate pv array into view by place bin array

%         if repeat == 1
            stc = pv.data.sessionTimeC;
%         end
        
        % spike shuffling

        spiketimes = spiketrain.timestamps/1000; % now in seconds
        spiketimes = spiketimes';

        spiketimes(spiketimes(:,1) < stc(1,1),:) = [];      
        
        % selecting rows from sessionTimeC
        stc(:,4) = [diff(stc(:,1)); 0];
        
        conditions = ones(size(stc,1),1);

        if Args.UseAllTrials == 0
            conditions = conditions & pv.data.good_trial_markers;
        end

        if Args.ThresVel > 0
            conditions = conditions & get(pv,'SpeedLimit',Args.ThresVel);
        end
        
        if Args.UseMinObs
            bins_sieved_p = pv.data.place_good_bins;
            bins_removed_p = setdiff(1:size(pv.data.place_intervals_count,1),bins_sieved_p);
            bins_sieved_sv = pv.data.view_good_bins;
            bins_removed_sv = setdiff(1:size(pv.data.view_intervals_count,1),bins_sieved_sv);
            conditions = conditions & (pv.data.pv_good_rows); % Make sure maps take into account both place and view filters
        else
            bins_sieved_p = 1:(Args.GridSteps * Args.GridSteps);
            bins_removed_p = [];
            bins_sieved_sv = 1:size(pv.data.view_intervals_count,1);
            bins_removed_sv = [];
        end

        disp('conditioning done');

        dstc = diff(stc(:,1));
        stc_changing_ind = [1; find(dstc>0)+1; size(stc,1)];
        stc_changing_ind(:,2) = [stc_changing_ind(2:end)-1; nan];
        stc_changing_ind = stc_changing_ind(1:end-1,:);

        % Initialise storage
        view_spikes = zeros(size(pv.data.view_intervals_count,1),Args.GridSteps * Args.GridSteps);
        place_spikes = zeros(1,Args.GridSteps * Args.GridSteps);
        interval = 1;
        for sp = 1:size(spiketimes,1)

            if rem(sp, 10000000) == 0
                disp(num2str(100*sp/size(spiketimes,1)))
            end

            while interval < size(stc_changing_ind,1)
                if spiketimes(sp,1) >= stc(stc_changing_ind(interval,1),1) && spiketimes(sp,1) < stc(stc_changing_ind(interval+1,1),1)
                    break;
                end
                interval = interval + 1;
            end   

            bins_hit_place = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),2);
            bins_hit_place = bins_hit_place(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))));
            bins_hit_place(~(bins_hit_place>0)) = [];

            bins_hit_view = stc(stc_changing_ind(interval,1):stc_changing_ind(interval,2),3);
            bins_hit_view = bins_hit_view(logical(conditions(stc_changing_ind(interval,1):stc_changing_ind(interval,2))));
            bins_hit_view(~(bins_hit_view>0)) = [];
            
            if ~isempty(bins_hit_view) && ~isempty(bins_hit_place) % Only spikes where both place and view are sampled are considered
                place_spikes(1,bins_hit_place) = place_spikes(1,bins_hit_place) + 1;
                view_spikes(bins_hit_view,bins_hit_place) = view_spikes(bins_hit_view,bins_hit_place) + 1;
            end
         
        end        
        
        stc_ssv = stc(find(conditions==1),[2 3 4]); % [place view dur]
        stc_ssv(~(stc_ssv(:,1) > 0),:) = [];
        stc_ssv = [stc_ssv; [0 size(pv.data.view_intervals_count,1) 0]];
        
        place_durations = nan(size(place_spikes,1), size(place_spikes,2));
        view_durations = nan(size(pv.data.view_intervals_count,1),size(place_spikes,2));
        for ii = 1:Args.GridSteps*Args.GridSteps
            
            inds = stc_ssv(:,1)==ii;
            subsample = [stc_ssv(inds,:)]; % [place view dur]
            % Consider only samples where both place and view are sampled
            subsample(isnan(subsample(:,2)),:) = [];
            if ~isempty(subsample) 
                % Get spikes and duration for place only
                place_durations(1,ii) = sum(subsample(:,3));
                % back-filling time for view
                subsample(subsample(:,3)==0,3) = nan;
                subsample(:,3) = fillmissing(subsample(:,3), 'previous');
                % padding with 5122 bin
                if subsample(end,1) ~= 5122
%                     subsample = [subsample; [NaN 5122 NaN]]; % Used to work, but now is value is nan, sum also becomes nan.
                    subsample = [subsample; [ii 5122 0]];
                end
                % remove bad view spots
                subsample(isnan(subsample(:,2)),:) = [];
                % sum durations
                view_durations(:,ii) = accumarray(subsample(:,2), subsample(:,3),[],[],NaN);
            end
        end
        % Filter out low-sampled bins
        place_durations(isnan(place_durations)) = 0; % Necessary because NaNs seem to mess up the smoothing. Can't do the same for view_durations because it is being used to compute llh and will erroneously give inf. 
%         view_durations(isnan(view_durations)) = 0; 
        place_durations(bins_removed_p) = 0;
        view_durations(bins_removed_sv,:) = 0;
        view_durations(:,bins_removed_p) = 0;
        place_spikes(bins_removed_p) = 0;
        view_spikes(bins_removed_sv,:) = 0;
        view_spikes(:,bins_removed_p) = 0;
        
        p_array_orig = place_spikes./place_durations;
        sv_array_orig = nansum(view_spikes,2)./nansum(view_durations,2);
        
        %% Compute covariance matrix
        disp('Computing covariance matrix');
        viewmapperplacebin = nan(size(view_durations));
        for ii = 1:size(place_durations,2)
            viewmapperplacebin(:,ii) = view_spikes(:,ii)./view_durations(:,ii);
        end
        covmat = cov(viewmapperplacebin,'partialrows');
        covmat_norm = covmat./nanmax(nanmax(abs(covmat)));
        % Replace NaNs with zeros in covariance matrix for norm calculations
        covmat_nonan = covmat;
        covmat_nonan(isnan(covmat_nonan)) = 0;
        % Calculate norms
        norml1 = norm(covmat_nonan,1); % maximum of column sum
        norml2 = norm(covmat_nonan,2); % maximum single value
            
        %% Establish independence of place and view maps
        mlm = true;
        initialwith = 'spatialview';
        convergewithsv = false;
        convergewithp = false;
        llh_vec = [];
        while mlm == true
        
            reinitialise = false;
            
            view_spikes_temp = view_spikes;
            view_spikes_temp(isnan(view_spikes)) = 0; % to prevent error on factorial

            % Maximum-likelihood maps. 
            if Args.StartOrig
                % Start with Original maps
                p_array = p_array_orig;  
                sv_array = sv_array_orig;
                p_spk = place_spikes;
                sv_spk = nansum(view_spikes,2);
            else
                % Start with uniform array
                p_array = ones(1,size(place_spikes,2));
                p_array(isnan(p_array_orig)) = nan;
                sv_array = ones(size(view_durations,1),1);
                sv_array(isnan(sv_array_orig)) = nan;
                p_spk = p_array.*nansum(view_durations,1);
                p_spk(isnan(p_spk)) = 0;
                sv_spk = sv_array.*nansum(view_durations,2);
                sv_spk(isnan(sv_spk)) = 0;
            end
            % Initialise storage. % Pad with original raw map first in
            % both cases so we can smooth them too for comparison with
            % pc/sv smoothed maps
            p_array_set = [p_array_orig; p_array];
            sv_array_set = [sv_array_orig, sv_array];
            p_spk_set = [place_spikes; p_spk];
            sv_spk_set = [nansum(view_spikes,2), sv_spk];
            
            % Starting llh ( upper limit of factorial 170, anything beyond that becomes Inf )
            for ii = 1:size(p_array_set,1)
                llh = sum( nansum(view_spikes.*log(p_array_set(ii,:).*view_durations.*sv_array_set(:,ii))) - nansum(p_array_set(ii,:).*view_durations.*sv_array_set(:,ii)) - nansum(log(factorial(view_spikes_temp))) );
                llh_vec(end+1,1) = llh;
            end
            if isinf(llh) % Evaluate the 2nd llh in the set since that is the starting point of the iterations
                pick = 2;
                if max(max(view_spikes_temp)) > 170
                    picklabel = ['spikecount' num2str(max(max(view_spikes_temp)))];
                    disp(['llh inf from start: spike count ' num2str(max(max(view_spikes_temp)))]);
                else
                    picklabel = 'inf1';
                    disp('llh inf from start: other reasons');
                end
                % end iterations
                mlm = false;
                break;
            end
            disp('maximising llh');

            while strcmp(initialwith,'spatialview') && ~convergewithsv % Try to find max llh first with flat sv array. If fail to converge, try with flat p array
                
                dur_adj = nan(size(view_durations,1),size(view_durations,2));
                for ii = 1:size(view_durations,2) % place
                    for jj = 1:size(view_durations,1) % view
                        if ~isnan(view_durations(jj,ii))
                            if isnan(sv_array(jj))
                                dur_adj(jj,ii) = 0;
                            elseif sv_array(jj) == 0
                                dur_adj(jj,ii) = view_durations(jj,ii); 
                            else
                                dur_adj(jj,ii) = view_durations(jj,ii)*sv_array(jj);
                            end
                        end
                    end
                end
                p_array = nansum(view_spikes,1)./nansum(dur_adj,1);
                spk = p_array.*nansum(view_durations,1);
                spk(isnan(spk)) = 0;
                p_array_set = [p_array_set; p_array];
                p_spk_set = [p_spk_set; spk];

                dur_adj = nan(size(view_durations,1),size(view_durations,2));
                for jj = 1:size(view_durations,1) % view
                    for ii = 1:size(view_durations,2) % place
                        if ~isnan(view_durations(jj,ii))
                            if isnan(p_array(ii))
                                dur_adj(jj,ii) = 0;
                            elseif p_array(ii) == 0
                                dur_adj(jj,ii) = view_durations(jj,ii); 
                            else
                                dur_adj(jj,ii) = view_durations(jj,ii)*p_array(ii);
                            end
                        end
                    end
                end
                sv_array = nansum(view_spikes,2)./nansum(dur_adj,2);
                spk = sv_array.*nansum(view_durations,2);
                spk(isnan(spk)) = 0;
                sv_array_set = [sv_array_set sv_array];
                sv_spk_set = [sv_spk_set spk];

                llh = sum( nansum(view_spikes.*log(p_array.*view_durations.*sv_array)) - nansum(p_array.*view_durations.*sv_array) - nansum(log(factorial(view_spikes_temp))) );
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
                                % Switch to using flat view array to start iterating
                                initialwith = 'place';
                                reinitialise = true;
                                disp('LLH did not converge with initial view array. Try initialising with place array.');
                                break;
                            end
                            % Test if global max is found
                            if maxnextlast > maxlast 
                                convergewithsv = true;
                                % Record global max
                                pick = find(llh_vec == max(llh_vec(2:end)));
                                if size(pick,1) > 1
                                    pick = pick(1);
                                    picklabel = ['undiff'];
                                elseif isinf(pick)
                                    % Skip on to using differnt initial array
                                    initialwith = 'place';
                                    reinitialise = true;
                                    disp('Failed to converge to non-infinite values with initial view array. Try initialising with place array');
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
                        convergewithsv = true;
                        % Record global max
                        pick = find(llh_vec == max(llh_vec(2:end)));
                        if size(pick,1) > 1
                            pick = pick(1);
                            picklabel = 'undiff';
                        elseif isinf(pick)
                            % Skip to using different initial array
                            initialwith = 'place';
                            reinitialise = true;
                            disp('Failed to converge to non-infinite values with initial view array. Try initialising with place array');
                            break;
                        else
                            picklabel = num2str(pick);
                            % display results
                            disp(['Converged with initial view array, pick = ' num2str(pick) ' out of ' num2str(size(llh_vec,1))]);
                            % break out of loop
                            mlm = false;
                            break;
                        end
                    elseif isinf(llh)
                        % Skip to using different initial array
                        initialwith = 'place';
                        reinitialise = true;
                        disp('Failed to converge to non-infinite values with initial view array. Try initialising with place array');
                        break;
                    elseif size(llh_vec,1) > Args.LLHIterLim
                        % No convergence. Switch to using flat view array to start iterating
                        initialwith = 'place';
                        reinitialise = true;
                        disp('LLH did not converge with initial view array. Try initialising with place array.');
                        break;
                    end
                end

            end
            if reinitialise
                continue;
            end
            
            while strcmp(initialwith,'place') && ~convergewithsv % If first try with flat sv array fails to converge, try with flat p array

                dur_adj = nan(size(view_durations,1),size(view_durations,2));
                for jj = 1:size(view_durations,1) % view
                    for ii = 1:size(view_durations,2) % place
                        if ~isnan(view_durations(jj,ii))
                            if isnan(p_array(ii))
                                dur_adj(jj,ii) = 0;
                            elseif p_array(ii) == 0
                                dur_adj(jj,ii) = view_durations(jj,ii); 
                            else
                                dur_adj(jj,ii) = view_durations(jj,ii)*p_array(ii);
                            end
                        end
                    end
                end
                sv_array = nansum(view_spikes,2)./nansum(dur_adj,2);
                spk = sv_array.*nansum(view_durations,2);
                spk(isnan(spk)) = 0;
                sv_array_set = [sv_array_set sv_array];
                sv_spk_set = [sv_spk_set spk];
                
                dur_adj = nan(size(view_durations,1),size(view_durations,2));
                for ii = 1:size(view_durations,2) % place
                    for jj = 1:size(view_durations,1) % view
                        if ~isnan(view_durations(jj,ii))
                            if isnan(sv_array(jj))
                                dur_adj(jj,ii) = 0;
                            elseif sv_array(jj) == 0
                                dur_adj(jj,ii) = view_durations(jj,ii); 
                            else
                                dur_adj(jj,ii) = view_durations(jj,ii)*sv_array(jj);
                            end
                        end
                    end
                end
                p_array = nansum(view_spikes,1)./nansum(dur_adj,1);
                spk = p_array.*nansum(view_durations,1);
                spk(isnan(spk)) = 0;
                p_array_set = [p_array_set; p_array];
                p_spk_set = [p_spk_set; spk];

                llh = sum( nansum(view_spikes.*log(p_array.*view_durations.*sv_array)) - nansum(p_array.*view_durations.*sv_array) - nansum(log(factorial(view_spikes_temp))) );
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
                                convergewithp = true;
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
                                    disp(['Converged with initial place array, pick = ' num2str(pick)]);
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
                        convergewithp = true;
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
                            disp(['Converged with initial place array, pick = ' num2str(pick) ' out of ' num2str(size(llh_vec,1))]);
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
        pfac = repmat(nansum(nansum(view_spikes)),size(p_spk_set,1),1) ./ nansum(p_spk_set,2);
        svfac = repmat(nansum(nansum(view_spikes)),1,size(sv_spk_set,2))./ nansum(sv_spk_set,1);
        pfac(1) = 1; % Original place map
        svfac(1) = 1; % Original view map
        if Args.StartOrig 
            pfac(2) = 1;
            svfac(2) = 1;
        end
        p_array_set = repmat(pfac,1,size(p_array_set,2)) .* p_array_set;
        p_spk_set = repmat(pfac,1,size(p_spk_set,2)) .* p_spk_set;
        sv_array_set = repmat(svfac,size(sv_array_set,1),1) .* sv_array_set;
        sv_spk_set = repmat(svfac,size(sv_spk_set,1),1) .* sv_spk_set;
        
        % Remove low obs bins
        p_array_set(:,bins_removed_p) = nan;
        p_spk_set(:,bins_removed_p) = nan;
        sv_array_set(bins_removed_sv,:) = nan;
        sv_spk_set(bins_removed_sv,:) = nan;
        
        % Limit output map set so as not to overwhelm file size
        numIterLlh = size(p_array_set,1);
        maxNumSmooth = 10;
        if pick <= maxNumSmooth
            smoothset = 1:pick;
        else
            smoothset = [1 pick-maxNumSmooth+2:pick];
        end
        numIterSm = length(smoothset);
        smoothpick = numIterSm;
        
        p_array_smoothset = p_array_set(smoothset,:);
        p_spk_smoothset = p_spk_set(smoothset,:);
        sv_array_smoothset = sv_array_set(:,smoothset);
        sv_spk_smoothset = sv_spk_set(:,smoothset);
        
        %%%%%%%%%%%%%%%
        
        % Store llh output
        data.convergewithsv = convergewithsv;
        data.convergewithp = convergewithp;
        data.llh = llh_vec;
        data.llhpick = pick;
        data.llhpicklabel = picklabel;
        data.maps_raw_corrp = p_array_set(pick,:);
        data.maps_raw_corrsv = sv_array_set(:,pick)';
        data.maps_raw_corrpset = p_array_smoothset;
        data.maps_raw_corrsvset = sv_array_smoothset';

        
        %% PLACE Smoothing
        
        % NOTE: Currently unable to make Kian Wei's batch smoothing code
        % work. 
        % So for now just looping

        % Smooth
        p_maps_adsm = nan(size(p_array_smoothset,1),size(p_array_smoothset,2));
        place_durationsG = reshape(place_durations,Args.GridSteps,Args.GridSteps);
        for ii = 1:numIterSm
            disp(['Smoothing place map ' num2str(smoothset(ii))]);
            place_spikesG = reshape(p_spk_smoothset(ii,:),Args.GridSteps,Args.GridSteps);
            [maps_adsmG,spk_adsmG,dur_adsmG] = adsmooth(place_durationsG,place_spikesG,Args.AlphaPlace);
            p_maps_adsm(ii,:) = reshape(maps_adsmG,1,Args.GridSteps*Args.GridSteps);
        end
        
        % SIC
        pdur1 = repmat(place_durations,numIterSm,1);
        Pi1 = pdur1./sum(pdur1,2); % consider nansum to play safe
        lambda_i = p_maps_adsm;
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

        % ISE 
        lambda_i = p_maps_adsm; 
        ise_out = ise(lambda_i(1,:), lambda_i(2:end,:), Args.GridSteps, Args.GridSteps);
        
        % Store output data
        data.maps_adsm_corrp = p_maps_adsm(smoothpick,:);
        data.maps_adsm_corrpset = p_maps_adsm;
        data.SIC_corrp = sic_out(smoothpick);
        data.SIC_corrpset = sic_out;
        data.ISE_corrp = ise_out(smoothpick);
        data.ISE_corrpset = ise_out;
        
        %% View Smoothing

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
        sv_maps_adsm = nan(size(sv_array_smoothset,1),size(sv_array_smoothset,2));
        durG = cell(size(binDepths,1),1);
        spkG = cell(size(binDepths,1),1);
        mapG = cell(size(binDepths,1),1);
        svdurlin = nansum(view_durations,2);
        for ii = 1:numIterSm % For each iteration
            disp(['Smoothing view map ' num2str(smoothset(ii))]);
            for jj = 1:9 % For each grid
                % Initialize empty grids
                dur = nan(binDepths(jj,1),binDepths(jj,2));
                spk = dur;
                % Assign linear bin to grid bin. Note that this goes left
                % to right, top to bottom
                for mm = 1:binDepths(jj,1)*binDepths(jj,2) % For every point in linear map
                    if mod(mm,binDepths(jj,2)) == 0
                        y = binDepths(jj,2);
                    else
                        y = mod(mm,binDepths(jj,2));
                    end
                    x = ceil(mm/binDepths(jj,2));
                    indbins_lin = mm + sum(binDepths(1:jj-1,1).*binDepths(1:jj-1,2));
                    % Assign
                    dur(x,y) = svdurlin(indbins_lin,1);
                    spk(x,y) = sv_spk_smoothset(indbins_lin,ii);
                end
                % Collect output
                durG{jj} = dur;
                spkG{jj} = spk;
            end
            retrievemap = cell(size(binDepths,1),1);
            durG_temp = durG; % keep original for reference during padding
            spkG_temp = spkG;
            for jj = 3:size(binDepths,1) % For each separate grid, minus the cue and hint
                % Pad each grid map with adjoining bins from other grids
                % Pad with <<5>> extra bin rows
                n = 5;
                [retrievemap{jj},durG{jj},spkG{jj}] = padgrids(n,durG_temp{jj},spkG_temp{jj},durG_temp,spkG_temp,gazeSections,jj);
%                 [retrievemap{jj},durG{jj},spkG{jj}] = padgrids(n,durG{jj},spkG{jj},durG,spkG,gazeSections,jj);
                % Smooth
                [maps_adsmG,spk_adsmG,dur_adsmG] = adsmooth(durG{jj},spkG{jj},Args.AlphaView);
                % Remove padding 
                unpad_maps_adsmG = maps_adsmG(retrievemap{jj}(1,1):retrievemap{jj}(1,2),retrievemap{jj}(2,1):retrievemap{jj}(2,2));
                % Put grid map back into linear map
                set = reshape(flipud(rot90(unpad_maps_adsmG)),size(unpad_maps_adsmG,1)*size(unpad_maps_adsmG,2),1);
                lin_inds = sum(binDepths(1:jj-1,1).*binDepths(1:jj-1,2))+1:sum(binDepths(1:jj,1).*binDepths(1:jj,2));
                sv_maps_adsm(lin_inds,ii) = reshape(set,1,binDepths(jj,1)*binDepths(jj,2));
            end
        end
        
        % SIC 
        svdur1 = repmat(svdurlin,1,numIterSm);
        Pi1 = svdur1./sum(svdur1,1);
        lambda_i = sv_maps_adsm;
        lambda_bar = nansum(Pi1 .* lambda_i,1);
        % divide firing for each position by the overall mean
        FRratio = lambda_i./repmat(lambda_bar,5122,1);
        % compute first term in SIC
        SIC1 = Pi1 .* lambda_i; 
        SIC2 = log2(FRratio);
        zeros_placing = SIC1==0;  

        bits_per_sec = SIC1 .* SIC2 ./ lambda_bar;
        bits_per_sec(zeros_placing) = NaN;
        lambda_bar_ok = lambda_bar>0;
        lambda_bar_bad = ~lambda_bar_ok;
        sic_out = nansum(bits_per_sec, 1);
        sic_out(lambda_bar_bad) = NaN;
        sic_out = sic_out';

        % ISE portion - taken directly from vmsv
        % create overall map and insert padded portions in, to account for
        % cross-portion pairs
        tic;
        canvas = nan(51, 161, numIterSm);
        firing_rates = sv_maps_adsm';

        % flooring
        floor_padded = nan(42,42,numIterSm);
        floor_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);
        floor_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,3203:3203+39),[2 1]), 40, 1, numIterSm),1);
        floor_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,3243:3243+39),[2 1]), 1, 40, numIterSm);
        floor_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,3283:3283+39),[2 1]), 40, 1, numIterSm);
        floor_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,3323:3323+39),[2 1]), 1, 40, numIterSm), 2);
        canvas(10:end,1:42,:) = floor_padded;

        % ceiling
        ceiling_padded = nan(42,42,numIterSm);
        ceiling_padded(2:end-1, 2:end-1, :) = flip(permute(reshape(firing_rates(:,1603:3202),size(firing_rates,1),40,40), [3 2 1]), 1);
        ceiling_padded(2:end-1,1,:) = flip(reshape(permute(firing_rates(:,4323:4323+39),[2 1]), 40, 1, numIterSm),1);
        ceiling_padded(1,2:end-1,:) = reshape(permute(firing_rates(:,4363:4363+39),[2 1]), 1, 40, numIterSm);
        ceiling_padded(2:end-1,end,:) = reshape(permute(firing_rates(:,4403:4403+39),[2 1]), 40, 1, numIterSm);
        ceiling_padded(end,2:end-1,:) = flip(reshape(permute(firing_rates(:,4443:4443+39),[2 1]), 1, 40, numIterSm), 2);
        canvas(10:end,44:85,:) = ceiling_padded;

        % walls
        walls_padded = nan(8,161,numIterSm);
        walls_padded(:,1:end-1,:) = flip(permute(reshape(firing_rates(:,3203:3203+1280-1), numIterSm, 40*4, 8),[3 2 1]), 1);
        walls_padded(:,end,:) = walls_padded(:,1,:);
        canvas(1:8,:,:) = walls_padded;

        % used to pad pillar base more easily
        floor_base = flip(permute(reshape(firing_rates(:,3:1602),size(firing_rates,1),40,40), [3 2 1]), 1);

        % pillars
        PTL_padded = nan(6,33,numIterSm);
        PTL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4963:4963+160-1), numIterSm, 8*4, 5),[3 2 1]), 1);
        % small diagonal issue here, diagonal floor bins at the corners are put
        % side by side, only 16 such occurrences in total, neglected for now.
        PTL_padded(end,1:8,:) = flip(permute(floor_base(9:16,8,:),[2 1 3]),2);
        PTL_padded(end,9:16,:) = floor_base(8,9:16,:);
        PTL_padded(end,17:24,:) = permute(floor_base(9:16,17,:),[2 1 3]);
        PTL_padded(end,25:32,:) = flip(floor_base(17,9:16,:),2);
        PTL_padded(:,end,:) = PTL_padded(:,1,:);
        canvas(10:10+6-1,87:87+32,:) = PTL_padded;

        PTR_padded = nan(6,33,numIterSm);
        PTR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4803:4803+160-1), numIterSm, 8*4, 5),[3 2 1]), 1);
        PTR_padded(end,1:8,:) = flip(permute(floor_base(9:16,24,:),[2 1 3]),2);
        PTR_padded(end,9:16,:) = floor_base(8,25:32,:);
        PTR_padded(end,17:24,:) = permute(floor_base(9:16,33,:),[2 1 3]);
        PTR_padded(end,25:32,:) = flip(floor_base(17,25:32,:),2);
        PTR_padded(:,end,:) = PTR_padded(:,1,:);
        canvas(10:10+6-1,121:121+32,:) = PTR_padded;

        PBL_padded = nan(6,33,numIterSm);
        PBL_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4643:4643+160-1), numIterSm, 8*4, 5),[3 2 1]), 1);
        PBL_padded(end,1:8,:) = flip(permute(floor_base(25:32,8,:),[2 1 3]),2);
        PBL_padded(end,9:16,:) = floor_base(24,9:16,:);
        PBL_padded(end,17:24,:) = permute(floor_base(25:32,17,:),[2 1 3]);
        PBL_padded(end,25:32,:) = flip(floor_base(33,9:16,:),2);
        PBL_padded(:,end,:) = PBL_padded(:,1,:);
        canvas(17:17+6-1,87:87+32,:) = PBL_padded;

        PBR_padded = nan(6,33,numIterSm);
        PBR_padded(1:end-1,1:end-1,:) = flip(permute(reshape(firing_rates(:,4483:4483+160-1), numIterSm, 8*4, 5),[3 2 1]), 1);
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

        disp(['time taken to pad map for ISE: ' num2str(toc)]);
        tic;

        ise_out = ise(actual_image, shuffled_images, 51, 161);
        disp(['time taken to compute ISE: ' num2str(toc)]);

        % Store output data
        sv_maps_adsm = sv_maps_adsm';
        data.maps_adsm_corrsv = sv_maps_adsm(smoothpick,:);
        data.maps_adsm_corrsvset = sv_maps_adsm;
        data.SIC_corrsv = sic_out(smoothpick);
        data.SIC_corrsvset = sic_out;
        data.ISE_corrsv = ise_out(smoothpick);
        data.ISE_corrsvset = ise_out;
        data.covmat = covmat;
        data.covmat_norm = covmat_norm;
        data.norml1 = norml1;
        data.norml2 = norml2;
        
    %% SAVE   
    
    % create nptdata so we can inherit from it    
    data.gridSteps = Args.GridSteps;
    data.NumIterLlh = numIterLlh;
    data.smoothpick = smoothpick;
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

function [smoothedRate,smoothedSpk,smoothedDur]=adsmooth(dur,spk,alpha)
% Adaptive smoothing of rate maps.
%
%       [smoothedRate,smoothedSpk,smoothedPos]=rates_adaptivesmooth(posMap,spkMap,alpha)
%
% Each bin is smoothed using a flat, circular kernal. The circle radius 
% is set for each bin, indivdually, such that 
%
%   radius => alpha ./ ( sqrt(nSpike) .* nDwell )
%
% where nSpike and nDwell are the number of spikes, and the amount of dwell time (in s) within the kernel.
%
% smoothedRate, smoothedSpk, smoothedPos are the smoothed maps (spike and pos maps are smoothed 
% with the same kernal group as for the rate map.

% Check for empty spk maps %
if sum(sum(spk))==0
    smoothedDur=dur;    smoothedDur(dur==0)=nan;
    smoothedSpk=spk;    smoothedSpk(dur==0)=nan;
    smoothedRate=spk;   smoothedRate(dur==0)=nan;
    return
end
% Pre-assign output %
smoothedDur=zeros(size(dur));
smoothedSpk=zeros(size(dur));
% Visited env template: use this to get numbers of visited bins in filter at edge of environemnt %
vis=zeros(size(dur));
vis(dur>0)=1;
% Pre-assign map which records which bins have passed %
smoothedCheck=false(size(dur));
smoothedCheck(dur==0)=true; % Disregard unvisited - mark as already done.
% Pre-assign list of radii used (this is for reporting purposes, not used for making maps) %
radiiUsedList=nan(1,sum(sum(dur>0)));
radiiUsedCount=1;

%%% Run increasing radius iterations %%%
r=1; % Circle radius
boundary=0; % IMFILTER boundary condition
while any(any(~smoothedCheck))
    % Check radius isn't getting too big (if >map/2, stop running) %
    if r>max(size(dur))/2
%     if r>20
        smoothedSpk(~smoothedCheck)=nan;
        smoothedDur(~smoothedCheck)=nan;
        break
    end
    % Construct filter kernel ...
    % Place: Flat disk, where r>=distance to bin centre %
    f=fspecial('disk',r); 
    f(f>=(max(max(f))/3))=1;
    f(f~=1)=0;   
    % Filter maps (get N spikes and pos sum within kernel) %
    fSpk=imfilter(spk,f,boundary);
    fDur=imfilter(dur,f,boundary);
    fVis=imfilter(vis,f,boundary);
    % Which bins pass criteria at this radius? %
    warning('off', 'MATLAB:divideByZero');
    binsPassed=alpha./(sqrt(fSpk).*fDur) <= r;
    warning('on', 'MATLAB:divideByZero');
    binsPassed=binsPassed & ~smoothedCheck; % Only get the bins that have passed in this iteration.
    % Add these to list of radii used %
    nBins=sum(binsPassed(:));
    radiiUsedList(radiiUsedCount:radiiUsedCount+nBins-1)=r;
    radiiUsedCount=radiiUsedCount+nBins;
    % Assign values to smoothed maps %
    smoothedSpk(binsPassed)=fSpk(binsPassed)./fVis(binsPassed);
    smoothedDur(binsPassed)=fDur(binsPassed)./fVis(binsPassed);
    % Record which bins were smoothed this iteration %
    smoothedCheck(binsPassed)=true;
    % Increase circle radius (half-bin steps) %
    r=r+0.5; % Increase radius in 0.5 bin steps.
end

% Assign Output %
warning('off', 'MATLAB:divideByZero');
smoothedRate=smoothedSpk./smoothedDur;
warning('on', 'MATLAB:divideByZero');
smoothedRate(dur==0)=nan;
smoothedDur(dur==0)=nan;
smoothedSpk(dur==0)=nan;

% Report radii sizes %
if 0
    hAllFigs = get(0, 'children');
    hFig = findobj(hAllFigs, 'flat', 'tag', 'adaptiveSmoothPlotWindow');
    if isempty(hFig);
        hFig=figure;
        set(hFig,'tag','adaptiveSmoothPlotWindow');
    else
        figure(hFig);
    end
    hist(radiiUsedList,1:10);
    uiwait(hFig,1.5);
end



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
% Patch: Remove NaNs because adaptive smoothing seems to have trouble with
% boundaries of NaNs
o_i(isnan(o_i)) = 0;
spikeLoc(isnan(spikeLoc)) = 0;
