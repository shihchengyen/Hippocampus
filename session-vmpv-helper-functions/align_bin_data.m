function aligned_bin_data = align_bin_data(vargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align_bin_data() This is part of a set of helper functions used to construct VMPV objects
% at the session level on days with multiple sessions
% Designed flow of functions:
% sessionInfo() --> splitBinData('1') --> align_bin_data()

% Description on how to use the helper functions %

% Workflow:

% sessionInfo() is run at /Session01 but can be modified to run at
% /2018xxxx instead

% sessionInfo() - creates a sessionInfo.mat ('sesh') object containing details on
%                   number of trials per session, based on raw session*.txt
%                   files found in each /Session0_ folder

% data structure within sesh obj:

% sessionCount: Integer value of number of sessions for the day
% sessionTrialUTime: N x 2 matrix , N = No. of sessions
%        Each Column Represents:
%
%        | Calculated No. of Trials| Number of Samples recorded | 

% After day-level binData has been constructed from raycasting,
% splitBinData('1') is run at Session01 to break 1binData_new.csv
% (day-level raycasted bin data) into session level and store the
% respective session-levels binData files at the correct session folders

% align_bin_data() is run at /Session02 onwards to make sure bin Data
% starts from the same timestamp as 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Args = struct('pix', 1, 'filename', 'binData.csv');
    bin_file = [num2str(Args.pix) Args.filename];
    
    %check if rplparallel is present and load rplparallel
    rp = rplparallel('auto');
    true_timestamps = rp.data.timeStamps';
    true_timestamps = true_timestamps(:) * 1000; % in ms

    %check if _binData.mat is present and loads
    if exist(bin_file, 'file') == 2
        fprintf('"%s file found". Loading existing object ... \n',bin_file);
        %bin_data = load(bin_file);
        bin_data = readmatrix(bin_file);
        %bin_data=bin_data.viewdata;
    else
        fprintf('"%s" not found. \n', bin_file);
    end

   
    %process binData
    bin_data_timestamps = bin_data(:,1);
    bin_data_timestamps = bin_data_timestamps - (bin_data_timestamps(1)-true_timestamps(1));
    bin_data_bins = bin_data(:,2);
    bin_data_markers = bin_data(:,3);
    bin_data_trial_timestamps = find(bin_data(:,3) ~= 0);
    bin_data_trial_timestamps_flat = bin_data_trial_timestamps';
    sorted_timestamps = sort(bin_data_timestamps);
    is_duplicate_in_sorted = diff(sorted_timestamps) == 0;
 
% 3. Extract the values at these positions and get unique ones
% The 'unique' ensures you only get each duplicate value once
duplicate_values = unique(sorted_timestamps(is_duplicate_in_sorted));
    
    for i = 1:length(true_timestamps)-1

        true_start = true_timestamps(i);
        true_end = true_timestamps(i+1);
        true_diff = true_end - true_start;
        current_start = bin_data_trial_timestamps_flat(i);
        current_end = bin_data_trial_timestamps_flat(i+1);
        current_chunk = double(bin_data_timestamps(current_start:current_end));
        current_diff = double(current_chunk(end) - current_chunk(1));
        current_start_time = current_chunk(1);
        current_end_time = current_chunk(end);
        current_chunk = (current_chunk - current_start_time)* true_diff/current_diff; % now scaled to rpl timing  
        current_chunk = current_chunk + current_start_time;
        shifting_needed = current_chunk(end) - current_end_time;
        bin_data_timestamps(current_start:current_end) = uint32(current_chunk);
        bin_data_timestamps(current_end+1:length(bin_data_timestamps))=bin_data_timestamps(current_end+1:length(bin_data_timestamps))+shifting_needed;
    end
    

    
    fprintf("updated 1bindata.csv \n");
    viewdata = [bin_data_timestamps, bin_data_bins, bin_data_markers];
    viewdata = viewdata(~isnan(viewdata(:,2)),:);
    disp(bin_file);
    writematrix(viewdata, bin_file);

end
