function obj = sessionInfo()
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
    args.classname = 'sessionInfo';
    args.matname = [args.classname '.mat'];
    args.matvarname = 'sesh';
% To decide the method to create or load the object
    start_dir=pwd;
    cd ../
    objectFileName = args.matname;
    objectVarName = args.matvarname;
    if exist(objectFileName, 'file') == 2
        % If the file exists, load the object
        fprintf('"%s" file found. Loading existing object...\n', objectFileName);
        % Load the file. 
        obj_0=load(objectFileName); 
        obj=obj_0.sesh;
        fprintf('Object "%s" loaded successfully from "%s".\n', objectVarName, objectFileName);
        cd (start_dir)

    else
    % If the file does not exist, create the object via your process
        cd (start_dir)
        fprintf('"%s" not found. Creating new object...\n', objectFileName);
        obj = createObject(args);
    end
end

function obj = createObject(args)
    % example object
    start_dir = pwd;
    cd ../
    day_dir = pwd;
    
    %--------- Counts number of session using session folders ------ %
    pattern = '^session0.*$';
    
    % Get directory listing
    listing = nptDir();
    
    % Get only directories
    is_folder = [listing.isdir];
    folders = listing(is_folder);
    
    % Get folder names as a cell array of strings
    folder_names = {folders.name};
    
    % Use cellfun with regexprep to find matches
    matching_logical_idx = ~cellfun('isempty', regexp(folder_names, pattern, 'match', 'ignorecase'));
    
    % Count the number of true values in the logical index
    sessionCount = sum(matching_logical_idx);
    sessionFolders = folder_names(matching_logical_idx);
    data.sessionCount = sessionCount;
    data.sessionTrialUTime = zeros(sessionCount, 2);
    %sessionTrialTime will consist, for each session (row):
    % col 1: number of trials | col 2: number of
    % timestamps recorded in unity raw data
    for i=1:length(sessionFolders)
        cd (sessionFolders{i})
        rd = dir('RawData*');
        if(~isempty(rd))
            cd(rd(1).name)
            dlist = nptDir('session_1*txt');
            % get entries in directory
            dnum = size(dlist,1);
            % check if the right conditions were met to create object
            if(dnum>0)
                % dnum - number of session raw files
                % this is a valid object
                % these are fields that are useful for most objects
                % these are object specific fields
                 %count of individual records within session
                trial_count = 0;%counts of trials within session
                row_count = 0; %total number of rows in all session.txt files
                for k = 1:dnum % go through all 'session_' files and extract data   
                    rawdata = dlmread(dlist(k).name,'',15,0); % start reading at row 15 (skip logged parameters)
                    row_count = row_count + height(rawdata);
                end
                for k=1:height(rawdata)
                    if (rawdata(k,1)~=0 && rawdata(k,1)~=84)
                        trial_count = trial_count+1;
                    end
                end
                data.sessionTrialUTime(i,1) = floorDiv(trial_count,3);
                data.sessionTrialUTime(i,2) = row_count;
            else
                disp ("No valid raw session files")
                obj = createEmptyObject(args);
            end
        else
            disp("No Raw txt file found in day/session")
            obj = createEmptyObject(args);
        end
        cd (day_dir)
    end
    obj = struct();
    obj.data = data;
    obj.classname = args.classname;
    eval([args.matvarname ' = obj;']);
    save(args.matname, args.matvarname);
    cd (start_dir)

end
function obj = createEmptyObject(args)
    % Create empty object with default fields
    data.sessionTrialUTime = [];
    data.sessionCount = 0;
    % Create the object as a structure
    obj = struct();
    obj.data = data;
    obj.classname = args.classname;
end