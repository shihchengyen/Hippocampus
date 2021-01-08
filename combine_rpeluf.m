function combine_rpeluf(day_dir)
% combine_rpeluf(day_dir)
%   day_dir example: '/Volumes/Hippocampus/Data/picasso-misc/20180628'
%   
%   For a single day directory, searches all the sessions (session0*) for
%   an object file and merges them into a single file. This is separately
%   carried out for rplparallel.mat, eyelink.mat and unityfile.mat each.
%
%   Between consecutive sessions, data falling in between the unity
%   session end (e.g. session n end) to start times (session n+1 start) is removed.
%   
%   Run aligning_objects.m after combine_rpeluf.
%

cd(day_dir)
sessions = dir('session0*');
if size(sessions,1) == 1
    error('There is only 1 session folder.')
end

disp(['day_dir: ' day_dir])



%%%%%%%%%%%%%%%%%%%%%%%
%%% eyelink section %%%
%%%%%%%%%%%%%%%%%%%%%%%

cd([day_dir '/session01'])
r = eyelink('auto');

r.data.time_removed_all = [0];
r.data.ind_removed_all = [0];

for s = 2:size(sessions,1)
    cd([day_dir '/session0' num2str(s)])
    p = r;
    q = eyelink('auto');
    if isempty(q)
        error(['Error in ' pwd ': empty eyelink'])
    end
    
    % useful fields for most objects
    r.data.numSets = p.data.numSets + q.data.numSets;
    r.data.noOfTrials = p.data.noOfTrials + q.data.noOfTrials;

    % object specific fields
    % This eyelink plus.m is written for combining same-day sessions

    % First, remove all information from the last trial of p to the
    % session start of q
    temp_trial_timestamps = p.data.trial_timestamps(~isnan(p.data.trial_timestamps));
    p_last_trial = temp_trial_timestamps(end);
    scanning_ind = p_last_trial;
    q_start_ind = 0;
    while q_start_ind == 0
        if q.data.timestamps(scanning_ind) >= q.data.session_start
            q_start_ind = scanning_ind - sum(p.data.ind_removed_all); % this gives q start in terms of p indices
        else
            scanning_ind = scanning_ind + 1;
        end
    end

    time_removed = p.data.timestamps(q_start_ind) - p.data.timestamps(p_last_trial);
    ind_removed = q_start_ind - p_last_trial;
    if time_removed <= 0 || ind_removed <= 0
        error(['eyelink session0' num2str(s) ': time_removed or ind_removed is negative'])
    end
    r.data.eye_pos(p_last_trial+1:q_start_ind-1,:) = [];
    r.data.timestamps(p_last_trial+1:q_start_ind-1) = [];
    r.data.timestamps(p_last_trial+1:end) = r.data.timestamps(p_last_trial+1:end) - time_removed;
    r.data.session_start = [p.data.session_start; r.data.timestamps(p_last_trial+1)];

    r.data.trial_timestamps = concat (p.data.trial_timestamps, (q.data.trial_timestamps-sum(r.data.ind_removed_all)-ind_removed), 'Rowise');
    r.data.trial_codes = concat (p.data.trial_codes, q.data.trial_codes, 'Rowise');

    % We must adjust the appropriate fix_times and fix_event entries.
    for setno = r.data.numSets:r.data.noOfSessions
        new_fix_times = r.data.fix_times(r.data.fix_times(:,3*(setno-1)+1)~=0, 3*(setno-1)+1:3*(setno-1)+2);
        new_fix_times = new_fix_times - double(time_removed);
        r.data.fix_times(1:size(new_fix_times,1),3*(setno-1)+1:3*(setno-1)+2) = new_fix_times;
    end
    % Remove the fix_times and fix_event(s) which occur during the "cut
    % out" time
    remember_sacc = 0;
    for row = 1:size(p.data.fix_times,1)
        if p.data.fix_times(row,3*(r.data.numSets-2)+2) >= p.data.timestamps(p_last_trial)
            r.data.fix_times(row,3*(r.data.numSets-2)+1:3*(r.data.numSets-2)+3) = 0;
            r.data.fix_event(row,r.data.numSets-1) = 0;
            r.data.sacc_event(row,r.data.numSets-1) = 0;
            remember_sacc = row;
        end
    end
    if remember_sacc ~= 0
        r.data.sacc_event(remember_sacc:end,r.data.numSets-1) = 0;
    end
    % When last session is reached, concatenate row wise for raycasting
    if r.data.numSets == r.data.noOfSessions
        final_fix_times = r.data.fix_times(:,1:3);
        final_fix_event = r.data.fix_event(:,1);
        final_sacc_event = r.data.sacc_event(:,1);
        for row = size(final_sacc_event,1):-1:1
            if final_sacc_event(row) ~= 0
                sacc_last_row = row;
                break
            end
        end
        final_sacc_event = final_sacc_event(1:sacc_last_row); % to handle the case where the entries are 0 but not at the end
        
        for setno = 2:r.data.noOfSessions
            final_fix_times = concatenate(final_fix_times,r.data.fix_times(:,3*(setno-1)+1:3*(setno-1)+3),'Rowise');
            final_fix_event = concatenate(final_fix_event,r.data.fix_event(:,setno),'Rowise');
            final_sacc_event = concatenate(final_sacc_event,r.data.sacc_event(:,setno),'Rowise');
            for row = size(final_sacc_event,1):-1:1
                if final_sacc_event(row) ~= 0
                    sacc_last_row = row;
                    break
                end
            end
            final_sacc_event = final_sacc_event(1:sacc_last_row);
        end
        fixt_remove = find(final_fix_times(:,1)==0);
        final_fix_times(fixt_remove,:) = [];
        final_fix_event = final_fix_event(final_fix_event ~= 0);

        r.data.fix_times = final_fix_times;
        r.data.fix_event = final_fix_event;
        r.data.sacc_event = final_sacc_event;
    end
    r.data.time_removed_all = [p.data.time_removed_all time_removed];
    r.data.ind_removed_all = [p.data.ind_removed_all ind_removed];
    
    %contains all the eye positions for all the sessions 
    %r.data.eye_pos = concat(p.data.eye_pos, q.data.eye_pos);  
    %OOF
    r.data.timeouts = concat(p.data.timeouts, q.data.timeouts, 'Rowise');

end

cd([day_dir '/session01'])
copyfile eyelink.mat eyelink_original.mat % backup original
el = r;
save('eyelink.mat','el')
disp('Merged eyelink saved!')



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% rplparallel section %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd([day_dir '/session01'])
r = rplparallel('auto');

for s = 2:size(sessions,1)
    cd([day_dir '/session0' num2str(s)])
    p = r;
    q = rplparallel('auto');
    if isempty(q)
        error(['Error in ' pwd ': empty rplparallel'])
    end
    if q.data.Args.Data.markers(1) ~= 84
        error(['Error in rplparallel session0' num2str(s) ': first marker not 84'])
    end
    
    % useful fields for most objects
    r.data.numSets = p.data.numSets + q.data.numSets;

    % object specific fields
    r.data.markers = [p.data.markers; q.data.markers];
    % "cut out" the duration before the 84 marker of q
    q_session_start = q.data.Args.Data.timeStamps(1);
    r.data.timeStamps = [p.data.timeStamps; (q.data.timeStamps - q_session_start + p.data.Args.Data.timeStamps(end))];
    r.data.Args.Data.markers = [p.data.Args.Data.markers q.data.Args.Data.markers];
    r.data.Args.Data.timeStamps = [p.data.Args.Data.timeStamps (q.data.Args.Data.timeStamps - q_session_start + p.data.Args.Data.timeStamps(end))];
    % save the "cut out" duration for later (for e.g., used to combine
    % spiketrain.mat)
    if ~isfield(p.data,'session_start_all') || length(p.data.session_start_all) == 1
        r.data.session_start_all = [p.data.Args.Data.timeStamps(1) q_session_start];
    else
        r.data.session_start_all = [p.data.session_start_all q_session_start];
    end
    % recalculate trialIndices
    r.data.trialIndices = floor(r.data.timeStamps*r.data.SampleRate);

end

cd([day_dir '/session01'])
copyfile rplparallel.mat rplparallel_original.mat % backup original
rp = r;
save('rplparallel.mat','rp')
disp('Merged rplparallel saved!')



%%%%%%%%%%%%%%%%%%%%%%%%%
%%% unityfile section %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

cd([day_dir '/session01'])
r = unityfile('auto');

for s = 2:size(sessions,1)
    cd([day_dir '/session0' num2str(s)])
    p = r;
    q = unityfile('auto');
    if isempty(q)
        error(['Error in ' pwd ': empty unityfile'])
    end
    
    % useful fields for most objects
    r.data.numSets = p.data.numSets + q.data.numSets;


    % object specific fields
    r.data.unityData = [p.data.unityData; q.data.unityData];
    r.data.unityTriggers = [p.data.unityTriggers; (q.data.unityTriggers + size(p.data.unityData,1))];
    r.data.unityTime = [p.data.unityTime; (q.data.unityTime(2:end) - q.data.unityTime(1) + p.data.unityTime(end))];

    pnewTrialTime = nan(max([size(p.data.unityTrialTime,1) size(q.data.unityTrialTime,1)]), size(p.data.unityTrialTime,2));
    pnewTrialTime(1:size(p.data.unityTrialTime,1),1:size(p.data.unityTrialTime,2)) = p.data.unityTrialTime;
    qnewTrialTime = nan(max([size(p.data.unityTrialTime,1) size(q.data.unityTrialTime,1)]), size(q.data.unityTrialTime,2));
    qnewTrialTime(1:size(q.data.unityTrialTime,1),1:size(q.data.unityTrialTime,2)) = q.data.unityTrialTime;
    r.data.unityTrialTime = cat(2, pnewTrialTime, qnewTrialTime);

end

cd([day_dir '/session01'])
copyfile unityfile.mat unityfile_original.mat % backup original
uf = r;
save('unityfile.mat','uf')
disp('Merged unityfile saved!')


end

