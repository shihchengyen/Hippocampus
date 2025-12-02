function splitBinData(pix)
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
    %pix is a string, reflects pixels used
    %e.g. '1' not 1
    if isnumeric(pix)
        pix = num2str(pix);
    end

    sesh = sessionInfo();
    session_dir = pwd; %starting directory to return to at the end
    cd ../
    day_dir_list = dir();
    mask = ((strcmp({day_dir_list.name},'Session01')) | (strcmp({day_dir_list.name},'session01')));
    cd (day_dir_list(mask).name) 
    %navigates to session01 / Session01 and look for _binData.hdf
    fileName = [pix 'binData_new.csv'];
    %viewdata = h5read(fileName,'/data');
    viewdata = readmatrix(fileName);
    %viewdata = viewdata';
 
    viewdata_copy = viewdata;
    raycasted_csv=readtable("unityfile_eyelink_new.csv");
   
    %slower approach

    %{
    row = 1;
    for i=1:sesh.data.sessionCount
        next_dir = "Session0" + string(i);
        cd ../
        cd (next_dir)
        marker_count = 0;
        session_start = 0;
        session_end = 0;
        markersCurrentSession = sesh.data.sessionTrialUTime(i,1)*2;       
        while marker_count ~= markersCurrentSession
            if row == height(raycasted_csv)
                disp("error")
            end
            event_marker = raycasted_csv{row,18};
            if ~strcmp(event_marker,"") && ~strcmp(event_marker,"F")
                %loop for timestamp with trial markers
                numeric_value = event_marker{1};
                string_value = num2str(numeric_value);
                first_char = string_value(1);
                if strcmp(first_char,'3') | strcmp(first_char, '4') | strcmp(first_char, '1')
                    %look for trial start/ trial end markers
                    marker_count = marker_count + 1;
                    if (marker_count == 1)
                        session_start = raycasted_csv{row,2};
                    end
                end
            end
            
            if marker_count == markersCurrentSession
                session_end = raycasted_csv{row,2};
            end
            row = row + 1;
        end
        keep_rows = (viewdata_copy(:, 1) <= session_end) & (viewdata_copy(:,1)>=session_start);
        viewdata = viewdata_copy(keep_rows,:);
        saveFileName = [pix 'binData.mat'];
        save (saveFileName, 'viewdata')         
    end
    %}


% Convert relevant columns to arrays/cell arrays once for efficiency
event_marker_col_raw = raycasted_csv{:, 18}; % Get the entire column 18 - event markers
timestamp_col = [raycasted_csv{:, 2}];      % Get the entire column 2 - timestamps

% Identify markers that are NOT empty and NOT "F"
valid_marker_mask = ~cellfun('isempty', event_marker_col_raw) & ~strcmp(event_marker_col_raw, "F");
timestamp_col = timestamp_col(valid_marker_mask);
event_marker_col_raw=event_marker_col_raw(valid_marker_mask);
event_marker_col_raw = cellfun(@(x) floor(str2double(x) / 10), event_marker_col_raw, 'UniformOutput', false);
raycasted_csv = [timestamp_col, [event_marker_col_raw{:}]'];
%raycasted_csv=raycasted_csv(raycasted_csv(:,2)~=2,:);
%number of rows of raycasted_csv = 3*total number of trials
row = 1;
for i=1:sesh.data.sessionCount
    next_dir = "session0" + string(i);
    cd ../
    cd (next_dir)
    markersCurrentSession = sesh.data.sessionTrialUTime(i,1)*3;
    session_start = raycasted_csv(row,1);
    session_end = raycasted_csv(row+markersCurrentSession-1,1);
    keep_rows = (viewdata_copy(:, 1) <= session_end) & (viewdata_copy(:,1)>=session_start);
    viewdata = viewdata_copy(keep_rows,:);
    viewdata(:,3) = 0;
    [lia, locb] = ismember(viewdata(:,1), raycasted_csv(:,1));
    viewdata(lia, 3) = raycasted_csv(locb(lia), 2);
    saveFileName = [pix 'binData.csv']; % Or .xlsx, etc.
    writematrix(viewdata, saveFileName);
    cd (session_dir) %splitBinData always returns to the starting directory
    row = row + markersCurrentSession;
end
