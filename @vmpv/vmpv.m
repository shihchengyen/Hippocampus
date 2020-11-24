function [obj, varargout] = vmpv(varargin)
%@vmpv Constructor function for vmpv class
%   OBJ = vmpv(varargin)
%
%   OBJ = vmpv('auto') attempts to create a vmpv object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmpv %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%example [as, Args] = vmpv('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Session','RequiredFile','binData.hdf', ...
				'GridSteps',40, 'overallGridSize',25, ...
                'MinObsPlace',5,'MinObsView',5,'MinDurPlace',0.05,'MinDurView',0.01);
            
Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps', 'MinObsPlace','MinObsView','overallGridSize','MinDurPlace','MinDurView'};                         

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'vmpv';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'pv';

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
	rp = rplparallel('auto',varargin{:});
    viewdata = h5read('binData.hdf','/data');
    size(viewdata)
    cd(ori);

    % major section 1 - getting combined sessiontime
    
    pst = uma.data.sessionTime;
    vst = viewdata';
    vst(:,1) = vst(:,1)/1000;
    cst_1 = [pst(:,1:2) nan(size(pst,1),1) ones(size(pst,1),1)];
    cst_2 = [vst(:,1) nan(size(vst,1),1) vst(:,2) 2.*ones(size(vst,1),1)];
    cst = [cst_2; cst_1];
    
    % hybrid insertion sort instead of default sortrows for speed
    % does insertion all at one go, at the end
    % can possibly be improved through some sort of linked-list?
    
    disp('sorting in');
    tic;
    index_in = nan(1,size(cst_1,1));
    idx1 = 1;
    history1 = 1;
    for r = size(cst_2,1)+1:size(cst,1) % for each place time
        if idx1 > 889108590
            disp('debug');
        end
        for i = history1:size(cst_2,1) % for each view time
            if cst(r,1) == cst(i,1) || cst(r,1) < cst(i,1) % if place time less than view time
                history1 = i;
                index_in(idx1) = i + idx1 - 1; % where to insert the place row in between view rows
                break;
            end
        end            
        if isnan(index_in(idx1))
            index_in(idx1) = idx1 + size(cst_2,1);
            history1 = size(cst_2,1);
        end
        if rem(idx1, 1000) == 0
            disp([num2str(idx1) '/' num2str(size(cst_1,1))]);
        end
        idx1 = idx1 + 1;
    end
    to_insert = cst_1;
    cst(setdiff(1:size(cst,1),index_in),:) = cst(1:size(cst_2,1),:); % setdiff one length more than expected here
    cst(index_in,:) = to_insert;
    disp(['sorting in done, time elapsed: ' num2str(toc)]);
 
    % cst = sortrows(cst, [1 4]); % sort by time, then put place in front
    disp('filling (no progress marker)');
    tic;
    cst(:,2) = fillmissing(cst(:,2),'previous'); % fills in view rows with place it current is at
    place_rows = find(cst(:,4) == 1);
    disp(['filling done, time elapsed: ' num2str(toc)]);
    disp('place pruning start (no progress marker)');
    tic;
    for row = length(place_rows):-1:1 
        % in the event that place changed at the same time as view's new sample, 
        % we want to remove the place row, as the place information has
        % already been transferred down to the subsequent view rows.
        % note that visual inspection might see it as the same when it
        % isn't, due to rounding in GUI.
        if place_rows(row) ~= size(cst,1)
            if cst(place_rows(row),1) ~= cst(place_rows(row)+1,1)
                place_rows(row) = [];
            end
        end
    end
    cst(place_rows,:) = [];
    disp(['place pruning end, time elapsed: ' num2str(toc)]);
    
    % now we need to replace surviving place rows, with duplicates of
    % itself, corresponding to the number of views in the previous time
    % bin
    place_rows = find(cst(:,4) == 1);
    first_non_nan = find(isnan(cst(:,3))==0);
    first_non_nan = first_non_nan(1);
    place_rows(place_rows < first_non_nan) = [];
    place_bins = cst(place_rows,2);
    
    reference_rows = place_rows - 1; 
    reference_times = cst(reference_rows,1); % timestamps from which to pull view bins from
    actual_times = cst(place_rows,1);
        
    % iterating and inserting rows takes too much time, below method to
    % make it feasible. first we count the number of duplicated place rows
    % to be added, then preallocate the final-sized array. 
    disp('counting max array size');
    tic;
    membership = ismember(cst(:,1),reference_times);
    cst_full = nan(size(cst,1)+sum(membership)-length(place_rows),4);
    searching_portion = cst(membership,:);
    % size of gaps needed is calculated here
    insertion_gaps = nan(length(reference_times),1);
    for idx = 1:length(insertion_gaps)
        if rem(idx, 1000) == 0
            disp([num2str(idx) '/' num2str(length(insertion_gaps))]);
        end
        insertion_gaps(idx) = sum(searching_portion(:,1)==reference_times(idx));
    end
    disp(['counting max done, time elapsed: ' num2str(toc)]);
    % slotting in original cst into enlarged one
    full_start = 1;
    original_start = 1;
    disp('slotting');
    tic;
    for idx = 1:length(reference_rows)
        if rem(idx, 1000) == 0
            disp([num2str(idx) '/' num2str(length(reference_rows))]);
        end
        chunk_to_insert = cst(original_start:reference_rows(idx),:);
        cst_full(full_start:full_start-1+size(chunk_to_insert,1),:) = chunk_to_insert;
        original_start = reference_rows(idx)+2;
        full_start = full_start-1+size(chunk_to_insert,1)+insertion_gaps(idx)+1;
    end
    cst_full(full_start:end,:) = cst(reference_rows(end)+2:end,:);
    disp(['slotting end, time elapsed: ' num2str(toc)]);
    % constructing the new portion
    disp('inserting (no progress marker)');
    tic;
    inserting_portion = searching_portion;
    inserting_portion(:,1) = repelem(actual_times,insertion_gaps);
    inserting_portion(:,2) = repelem(place_bins,insertion_gaps);
    inserting_portion(:,4) = 4;
    cst_full(isnan(cst_full(:,1)),:) = inserting_portion;
    disp(['inserting end, time elapsed: ' num2str(toc)]);
    if ~isempty(find(diff(cst_full(:,1))<0))
        error('combined sessiontime misaligned!');
    end
    cst_full = cst_full(:,1:3);
    disp('guaranteeing unique and ascending');
    tic;
    dti = find(diff(cst_full(:,1)));
    dti(:,2) = [1; 1 + dti(1:end-1,1)];
    dti = dti(:,[2 1]);
    dti = [dti; [dti(end,2)+1 size(cst_full,1)]];
    to_remove = [];
    for chunk = 1:size(dti, 1)
        cst_full(dti(chunk, 1): dti(chunk,2),:) = sortrows(cst_full(dti(chunk, 1): dti(chunk,2),:), [1 3]);
        identify_dup = diff(cst_full(dti(chunk, 1): dti(chunk,2),:));
        if dti(chunk, 1) - dti(chunk, 2) ~= 0
            if ~isempty(find(sum(identify_dup,2)==0))
                to_remove = [to_remove; find(sum(identify_dup,2)==0)+dti(chunk,1)];
            end
        end
    end
    cst_full(to_remove,:) = [];
    disp(['duplicates found: ' num2str(length(to_remove))]);
    disp(['guaranteeing unique and ascending done, time elapsed: ' num2str(toc)]);
    % cst_full = unique(cst_full,'rows');
    % cst_full = sortrows(cst_full,[1 3]);
    
    % before we save the combined sessiontime, try to reduce size and
    % optimize searching later on by targeting consecutive timestamps with
    % exactly the same view spots.
    
    % first, chunk by number of view bins per timestamp
    time_transitions = find(diff(cst_full(:,1))~=0) + 1;
    time_repeats = diff(time_transitions);
%     unique_times = unique(cst_full(:,1));
    time_repeats = [1; time_repeats; size(cst_full,1)+1-time_transitions(end)];
    time_repeats = [cst_full([1; time_transitions],1) time_repeats];
    %time_repeats = [unique(cst_full(:,1)) [1; time_repeats; NaN]];
    
    % third (for now) column will track candidate redundancies
    time_repeats(:,3) = [0; diff(time_repeats(:,2))==0];
    time_repeats(:,4) = 0;
    time_repeats = [[1; time_transitions] time_repeats];
    % at this point, interpret as below:
    % 1st col - row where timestamp first occurred in cst_full
    % 2nd col - new timestamp
    % 3rd col - number of occurrences for this timestamp
    % 4th col - marks out timestamps with same number of occurrences as
    % predecessor
    % 5th col - timestamps with exact same view bins as predecessor (to be
    % filled now)
    possible = find(time_repeats(:,4)==1);
    cst_full(isnan(cst_full(:,3)),3) = -10; % convert nan views to -10 for now
    for idx = 1:length(possible)
        if rem(idx,100000)==0
            disp(['view ' num2str(100*idx/length(possible))]); % percentage done output
        end
        % long formula basically checks for perfectly identical view bins
        checker = ismember(cst_full(time_repeats(possible(idx),1):time_repeats(possible(idx),1)+time_repeats(possible(idx),3)-1,3), cst_full(time_repeats(possible(idx)-1,1):time_repeats(possible(idx)-1,1)+time_repeats(possible(idx),3)-1,3));
        if sum(checker) == length(checker)
            time_repeats(possible(idx),5) = 1;
        end
    end
    % last column - 6th filters for same place bin
    time_repeats = [time_repeats zeros(size(time_repeats,1),1)];
    possible = find(time_repeats(:,5)==1);
    for idx = 1:length(possible)
        if rem(idx,100000)==0
            disp(['place ' num2str(100*idx/length(possible))]); % percentage done output
        end
        if cst_full(time_repeats(possible(idx),1),2) == cst_full(time_repeats(possible(idx)-1,1),2)
            time_repeats(possible(idx),6) = 1;
        end        
    end
    % actually remove timings, from combined sessiontime
    timing_to_remove = time_repeats(time_repeats(:,6)==1,2);
    cst_indices = find(ismember(cst_full(:,1), timing_to_remove));
    disp(['removing ' num2str(length(cst_indices)) ' rows']);
    cst_full(cst_indices,:) = [];   

    cst_full(cst_full(:,3)==-10,3) = NaN; % back-conversion
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data.sessionTimeC = cst_full;
    stc = cst_full;
    
    % major section 2 - getting view intervals, accounting and not
    % accounting for speed thresholding
                
        max_allo = 5000;
        view_split1 = nan(5122,max_allo,2);
        view_split2 = ones(5122,1);
        view_split3 = zeros(5122,1); % tracks starting time of this bin's interval if active, 0 if not.
        idx_tracker = nan(5122,max_allo);
        %stacked = 0;                
        last_processed_time = [0 0];
        for row = 1:size(stc,1)
%             if stc(row,1) > 18
%                 disp('stopper');
%             end
            if ~isnan(stc(row,3))
                if stc(row,1) ~= last_processed_time(2)
                    last_processed_time = [last_processed_time(2) stc(row,1)];
                end
                target = stc(row,3);
                if view_split3(target) == 0 % && stc(row,2) ~= -1
                    if stc(row,2) > 0 || stc(row,2) == -1 %~= 0 % mod line (was > 0)
                        view_split1(target,view_split2(target),1) = stc(row,1);
                        view_split3(target) = stc(row,1);
                        idx_tracker(target,view_split2(target)) = row;
                    end
                elseif (view_split3(target)~=last_processed_time(1) && view_split3(target)~=last_processed_time(2)) || stc(row,2) == 0 %== 0 % mod line (was < 1) % should be < 1

                        for subr = idx_tracker(target,view_split2(target)):size(stc,1)
                            if stc(subr,1) > stc(idx_tracker(target,view_split2(target)),1)
                                view_split1(target,view_split2(target),2) = stc(subr,1);
                                break;
                            end
                        end
            %             view_split1(target,view_split2(target),2) = view_split3(target);
            %             idx_tracker(target,view_split2(target)) = row;
                        if view_split2(target) == size(view_split1,2) - 1
                            view_split1 = cat(2, view_split1, nan(5122,max_allo,2));
                            idx_tracker = cat(2, idx_tracker, nan(5122,max_allo));
                        end
                        view_split2(target) = view_split2(target) + 1;
                        if stc(row,2) > 0 || stc(row,2) == -1 %~= 0 % mod line (was > 0)
                            view_split3(target) = stc(row,1);
                            view_split1(target, view_split2(target), 1) = stc(row,1);
                            idx_tracker(target,view_split2(target)) = row;
                        else
                            view_split3(target) = 0;
                        end

                else
                    view_split3(target) = stc(row,1);
                    idx_tracker(target,view_split2(target)) = row;
                end  
            end
        end

        for view = 1:5122
            disp(['view shifting: ' num2str(view)]);
            if sum(isnan(view_split1(view, view_split2(view),:)))==1
                found = 0;
                for subr = idx_tracker(view,view_split2(view)):size(stc,1)
                    if stc(subr,1) > stc(idx_tracker(view,view_split2(view)),1)
                        view_split1(view,view_split2(view),2) = stc(subr,1);
                        found = 1;
                        break;
                    end
                end       
                if found == 0
                    view_split1(view, view_split2(view), 2) = stc(end,1);
                end
%                 view_split1(view, view_split2(view), 2) = stc(end,1);
            elseif sum(isnan(view_split1(view, view_split2(view),:)))==2
                view_split2(view) = view_split2(view) - 1;
            end
            
        end
                view_split1 = view_split1(:,1:min([size(view_split1,2) max(view_split2)+2]),:);
                data.view_intervals = view_split1; % upper bound for intervals is the first timestamp where it is no longer seen
                data.view_intervals_count = view_split2;                
            
    % major section 3 - getting place intervals, with and without
    % thresholding speed
    
    sessionTimeC = stc;
    
    preallo = 150;
    
    place_ignore_speed_intervals = nan(1600, preallo, 2);
    place_ignore_speed_intervals_count = ones(1600,1);

    curr_place = 0;
    curr_start_time = 0;

    for row = 1:size(sessionTimeC,1)

        if sessionTimeC(row,2) ~= curr_place && sessionTimeC(row,2) ~= -1
            interval = [curr_start_time; sessionTimeC(row,1)];
            if curr_place ~= 0 && curr_place ~= -1
                place_ignore_speed_intervals(curr_place, place_ignore_speed_intervals_count(curr_place), 1) = interval(1);
                place_ignore_speed_intervals(curr_place, place_ignore_speed_intervals_count(curr_place), 2) = interval(2);
                place_ignore_speed_intervals_count(curr_place) = place_ignore_speed_intervals_count(curr_place) + 1;
                if place_ignore_speed_intervals_count(curr_place) == size(place_ignore_speed_intervals,1) - 5
                    place_ignore_speed_intervals = cat(2,place_ignore_speed_intervals,nan(1600, preallo, 2));
                end
            end
            curr_place = sessionTimeC(row,2);
            curr_start_time = sessionTimeC(row,1);
            if row == size(sessionTimeC,1) % in the rare case that the place changes on the last stc row
                place_ignore_speed_intervals(curr_place,place_ignore_speed_intervals_count(curr_place),1) = curr_start_time;
                place_ignore_speed_intervals(curr_place,place_ignore_speed_intervals_count(curr_place),2) = uma.data.sessionTime(end,1);
                place_ignore_speed_intervals_count(curr_place) = place_ignore_speed_intervals_count(curr_place) + 1;
            end
        elseif row == size(sessionTimeC,1)
            if curr_place ~= 0 && curr_place ~= -1
                place_ignore_speed_intervals(curr_place,place_ignore_speed_intervals_count(curr_place),1) = curr_start_time;
                place_ignore_speed_intervals(curr_place,place_ignore_speed_intervals_count(curr_place),2) = sessionTimeC(end,1);
                place_ignore_speed_intervals_count(curr_place) = place_ignore_speed_intervals_count(curr_place) + 1;
            end
        end

    end

    place_ignore_speed_intervals_count = place_ignore_speed_intervals_count - 1;    
    place_ignore_speed_intervals = place_ignore_speed_intervals(:,1:max(place_ignore_speed_intervals_count),:);
    data.place_intervals = place_ignore_speed_intervals;
    data.place_intervals_count = place_ignore_speed_intervals_count;
    
    % major section 4 - labeling instance by instance validity for view
    % domain
    
        view_ignore_speed_intervals = data.view_intervals;
        view_ignore_speed_intervals_count = data.view_intervals_count;

        long_interval_cutoff = Args.MinDurView;

        view_good_ignore_speed = zeros(size(sessionTimeC,1),1);

        curr_int_ignore_speed = ones(5122,1); % track the last used interval per bin

        for strow = 1:size(sessionTimeC,1)

            timestamp = sessionTimeC(strow,1);
            view_bin = sessionTimeC(strow,3);

            if isnan(view_bin)
                continue;
            end

            % ignoring speed requirement
            while timestamp > view_ignore_speed_intervals(view_bin,curr_int_ignore_speed(view_bin),2) && curr_int_ignore_speed(view_bin) <= view_ignore_speed_intervals_count(view_bin) % pan through intervals
                curr_int_ignore_speed(view_bin) = curr_int_ignore_speed(view_bin) + 1;
            end
            if curr_int_ignore_speed(view_bin) <= view_ignore_speed_intervals_count(view_bin) % still have possible intervals
                if timestamp >= view_ignore_speed_intervals(view_bin,curr_int_ignore_speed(view_bin),1) && timestamp < view_ignore_speed_intervals(view_bin,curr_int_ignore_speed(view_bin),2) % within current interval
                    if view_ignore_speed_intervals(view_bin,curr_int_ignore_speed(view_bin),2) - view_ignore_speed_intervals(view_bin,curr_int_ignore_speed(view_bin),1) >= long_interval_cutoff
                        view_good_ignore_speed(strow) = 1; % resides in long interval
                    else
                        view_good_ignore_speed(strow) = 0.5; % in short interval
                    end
                end
            end

            if rem(strow,10000000) == 0
                disp(['view instance: ' num2str(strow)]);
            end

        end
        
        data.view_good_rows = view_good_ignore_speed;

        
    % major section 5 - row by row filtering for intervals in place domain

        long_interval_cutoff = Args.MinDurPlace;

        place_good_ignore_speed = zeros(size(sessionTimeC,1),1);

        curr_int_ignore_speed = ones(1600,1); % track the last used interval per bin
        
        hist_bin = 0;
        for strow = 1:size(sessionTimeC,1)
        
            timestamp = sessionTimeC(strow,1);
            place_bin = sessionTimeC(strow,2);

            if place_bin == -1
                place_bin = hist_bin;
            end
            
            if isnan(place_bin) || place_bin == 0
                continue;
            end

            
            % ignoring speed requirement
            while timestamp > place_ignore_speed_intervals(place_bin,curr_int_ignore_speed(place_bin),2) && curr_int_ignore_speed(place_bin) <= place_ignore_speed_intervals_count(place_bin) % pan through intervals
                curr_int_ignore_speed(place_bin) = curr_int_ignore_speed(place_bin) + 1;
            end
            if curr_int_ignore_speed(place_bin) <= place_ignore_speed_intervals_count(place_bin) % still have possible intervals
                if timestamp >= place_ignore_speed_intervals(place_bin,curr_int_ignore_speed(place_bin),1) % within current interval
                    if place_ignore_speed_intervals(place_bin,curr_int_ignore_speed(place_bin),2) - place_ignore_speed_intervals(place_bin,curr_int_ignore_speed(place_bin),1) >= long_interval_cutoff
                        place_good_ignore_speed(strow) = 1; % resides in long interval
                    else
                        place_good_ignore_speed(strow) = 0.5; % in short interval
                    end
                end
            end
            
            hist_bin = place_bin;
            
            if rem(strow,10000000) == 0
                disp(['place instance ignore_speed: ' num2str(strow)]);
            end

        end    
    
        data.place_good_rows = place_good_ignore_speed;
        
    % Consolidate good rows across place and view
    
        data.pv_good_rows = view_good_ignore_speed & place_good_ignore_speed;
        
    % minor section 1 - thresholding minobs long duration    
        

    dur_pi = data.place_intervals(:,:,2) - data.place_intervals(:,:,1);
    dur_vi = data.view_intervals(:,:,2) - data.view_intervals(:,:,1);
    
    count_pi = sum(dur_pi > Args.MinDurPlace,2);
    count_vi = sum(dur_vi > Args.MinDurView,2);
    
    data.place_good_bins = find(count_pi > Args.MinObsPlace);
    data.view_good_bins = find(count_vi > Args.MinObsView);
        
    % minor section 2 - marking out first half second half of session

    rpts = rp.data.timeStamps;
    last_row = round(size(rpts,1)/2);
    last_half_timing = rpts(last_row,3);
    halving_markers = double(sessionTimeC(:,1) > last_half_timing) + 1;
    data.halving_markers = halving_markers;
    data.unityData = uma.data.unityData;
    data.unityTime = uma.data.unityTime;
    
    % minor section 3 - marking out good trial rows
    
    good_trials = uma.data.processTrials;
    good_trial_markers = zeros(size(sessionTimeC,1),1);

    possible_int = 1;
    for row = 1:length(good_trial_markers)    
        for idx = possible_int:length(good_trials)
            if sessionTimeC(row,1) >= rpts(good_trials(idx),1)
                if sessionTimeC(row,1) < rpts(good_trials(idx),3)
                    good_trial_markers(row) = 1;
                    possible_int = idx;
                    break;
                end
            end
        end
        if rem(row,10000000)==0
            disp(['good trials: ' num2str(row)]);
        end
    end    
    data.good_trial_markers = good_trial_markers;
    
    % misc section, for shuffling and stuff
    
    data.rplmaxtime = rp.data.timeStamps(end,3);
    
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
