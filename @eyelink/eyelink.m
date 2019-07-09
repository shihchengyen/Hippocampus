%this is an object template file
%obj == et which was the object returned from the eyestarget function

function [obj, varargout] = eyelink(varargin)
%@eyelink Constructor function for eyelink class
%   OBJ = eyelink(varargin)
%
%   OBJ = eyelink('auto') attempts to create a eyelink object by first
%   creating an edf file for the day, and then extracting the eye
%   positions. It also has fields that store the trials and sessions in a
%   day.
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on eyelink %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = eyelink('save','redo')
%
%dependencies:

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0,'Navigation', 0, 'Calibration',0);
Args.flags = {'Auto','ArgsOnly', 'Navigation', 'Calibration'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {};

%varargin contains the arguments entered into the function. Args is the
%array in which the parameter-value pairs are stored
[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'eyelink';
Args.matname = [Args.classname '.mat']; %this will get the file which stores the object, if the object was saved
Args.matvarname = 'df';

% To decide the method to create or load the object
[command,robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyEyelink(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyEyelink(Args);
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
dnum = size(dlist,1); %this returns the number of rows. i.e. the number of directories

% check if the right conditions were met to create object
if(dnum>0) %you can create the object now

    if(Args.Navigation)
        %get to the day directory and create the object
        cwd = pwd; %save the current working directory
        str = getDataOrder('day');
        cd (str);

        %look for the edf file
        dlist = nptDir ('.edf');
        %convert the edf file into a MATLAB accessible format
        edfdata = edfmex (dlist(1).name);

        %get the experiment start time
        %expTime = edfdata.FEVENT(1).sttime;

        %%Storing the eye positions: extract the x ad y positions of the eye
        %%and remove the eye movements where the monkey was not looking at the
        %%screen, i.e. the eye positions were outside the screen bounds
        x = (edfdata.FSAMPLE.gx(1,:))';
        x(x>=1921)= NaN;
        %x(x<=-1)= NaN;

        y = (edfdata.FSAMPLE.gy(1,:))';
        y(y>=1081)= NaN;
        %y(y<=-1)= NaN;

        eyePos = horzcat(x,y);

        %Create a matrix that has all the possible messages and a vector that
        %contains the index that all of them happen at.- this is done
        %to extract the start and end of the experiment and all the trials
        type = {edfdata.FEVENT(:).type};
        type = cell2mat(type);
        type = type-24; %the index for a 'Start/Cue/End' message is 24
        messageEvent = find (~type); %stores the indices where 'Start/Cue/End' messages were generated.
        messageEvent = messageEvent';
        m = {edfdata.FEVENT(messageEvent(:)).message}'; %this vector stores all the messages as cell chars

        %clearing the first few messages till Trigger Version #
        m(1:7) = [];
        messageEvent (1:7) = [];

        %Trigger Version 84 signals the start of the session, and so we shall
        %make different time stamps for the different sessions
        s = {'Trigger Version 84'};
        sessionIndex = find(strcmp(m, s)); %sessionIndex has the index inside messageEvent where a new session starts
        noOfSessions = size (sessionIndex, 1); %this stores the number of sessions in the day
        fprintf ('No of Sessions ');
        disp(noOfSessions);

        %In this for loop, we interatively call completeData function to check
        %on the various
        trialTimestamps = zeros(size(m,1), 3*noOfSessions); 
        noOfTrials = zeros(1, noOfSessions);
        missingData = [];

        for i=1:noOfSessions 

            sessionName = dir('session*');
            if(contains(sessionName(i).name, num2str(i)) == 1) 
                fprintf(strcat('Session Name: ', sessionName(i).name, '\n'));
                idx = sessionIndex(i,1);

                if (i==noOfSessions)  
                    [corrected_times,tempMissing] = completeData(edfdata, m(idx:end, 1), messageEvent(idx:end,1), sessionName(i).name);
                else
                    idx2 = sessionIndex(i+1,1);
                    [corrected_times,tempMissing] = completeData(edfdata, m(idx:idx2, 1), messageEvent(idx:idx2,1), sessionName(i).name);
                end

                l = 1 + (i-1)*3;
                u = 3 + (i-1)*3;
                row = size(corrected_times,1);
                trialTimestamps (1:row, l:u) = corrected_times;
                noOfTrials (1,i) = size(corrected_times, 1);
                missingData = vertcat(missingData, tempMissing);
            end 
        end

         %edit the size of the array and remove all zero rows '
         trialTimestamps = trialTimestamps(any(trialTimestamps,2),:);

         if(size(missingData,1) ~= 0) %if the table is not empty 
             str = strcat('missingData_', (dlist(1).name), '.csv');
             writetable(missingData, (str));
         end 

         %%Make a matrix with all the timeouts in all the trials in the session
         %%which we can check when we are graphing lines for the end trial (refer plot.m)
         c = {'Timeout'};
         timeouts = contains (m,c);
         timeouts = {edfdata.FEVENT(messageEvent(find(timeouts))).sttime}';
         timeouts = cell2mat (timeouts); %stores all the timeouts in the session

         %%Now, get the time for all the eye events that occur in all
         %%experiments
         %startTime = edfdata.FEVENT(1).sttime; % get time when recording was started (time when 'MODE RECORD' message was sent)

         %Now, we are going to try and the eye Fix and Sacc events according to
         %sessions and trial, as well as their durations. They are going to be
         %stored in separate vectors, sacc and fix

         events = {edfdata.FEVENT(:).codestring}';
         indexFix = find(strcmp(events,'ENDFIX')); % get index of fixations (end fixations)
         indexSacc = find(strcmp(events,'ENDSACC')); % get index of saccades (end saccades)


         for j=1:noOfSessions

             if (j==noOfSessions)
                 idx2 = indexSacc(indexSacc > messageEvent(sessionIndex(j)));
                 idx1 = indexFix(indexFix > messageEvent(sessionIndex(j)));
             else
                 idx2 = indexSacc(indexSacc > messageEvent(sessionIndex(j,1)) & indexSacc < messageEvent(sessionIndex(j+1,1)));
                 idx1 = indexFix(indexFix > messageEvent(sessionIndex(j,1)) & indexFix< messageEvent(sessionIndex(j+1,1))); %extract the relevant indices
             end

            fixEvents = cell (size(idx1,1),2);
            fixEvents (:,1) =  {edfdata.FEVENT(idx1).sttime}';
            fixEvents (:,2)=  {edfdata.FEVENT(idx1).entime}';
            fixEvents = cell2mat(fixEvents);
            fixEvents(:,3)  = fixEvents(:,2) - fixEvents(:,1); %get the duration
            fixEvents (:, 1:2) = [];
            fix(1:size(idx1,1), j) = fixEvents;

            saccEvents = cell (size(idx2,1),2);
            saccEvents (:,1) =  {edfdata.FEVENT(idx2).sttime}';
            saccEvents (:,2)=  {edfdata.FEVENT(idx2).entime}';
            saccEvents = cell2mat(saccEvents);
            saccEvents(:,3)  = saccEvents(:,2) - saccEvents(:,1); %get the duration
            saccEvents (:, 1:2) = [];
            sacc(1:size(idx2,1),j) = saccEvents;

         end

         %remove all the excess 0 rows
         fix = fix(any(fix,2),:);
         sacc = sacc(any(sacc,2), :);

         %Assign values to the data object structure
         noOfSessions = size(trialTimestamps,2)/3;

         for idx=1:noOfSessions

             %this selects the session director and cds into it
             strName = sessionName(idx).name;
             cd (strName);

             l = 1+(idx-1)*3;
             u = l+2;
             data.trial_timestamps = trialTimestamps(:, l:u); %contains all start, cue and and times for all the trials
             data.trial_timestamps = data.trial_timestamps(any(data.trial_timestamps,2),:);

             data.sacc_event = sacc;
             data.fix_event = fix;
             data.timestamps = (edfdata.FSAMPLE.time)'; %all the time that the experiments were run for
             data.eye_pos = eyePos;  %contains all the eye positions in all the sessions
             data.noOfSessions = noOfSessions;
             data.timeouts = timeouts;
             data.noOfTrials= noOfTrials(1, idx);

             % create nptdata so we can inherit from it
             data.numSets = 1;    %eyelink is a session object = each session has only object 
             data.Args = Args;
             n = nptdata(data.numSets,0,pwd);
             d.data = data;
             obj = class(d,Args.classname,n);
             saveObject(obj,'ArgsC',Args);

             %after saving object, we cd back to the parent directory 
             cd ..
         end 
         
%This creates an object containing the data from the calibration tests
%performed in the day
    elseif (Args.Calibration) 
        %get to the day directory and create the object
        cwd = pwd; %save the current working directory
        str = getDataOrder('day');
        cd (str);

        dlist = dir ('*_*.edf')
        edfdata = edfmex (dlist(1).name); %get the name of the calibration edf
        
        startTime = edfdata.FSAMPLE.time(1);
        type = {edfdata.FEVENT(:).type};
        type = cell2mat(type); 
        type = type - 24;
        messageEvent = find(~type); 
        messages = {edfdata.FEVENT(messageEvent(:)).message}'; %stores a digit indicating if it's start/reward/end
        
        eventTimes = {edfdata.FEVENT(messageEvent(:)).sttime}'; %stores the time the events took place
        sz = ceil(size(eventTimes,1)/3); %each trial has three events in it 
        trialTimestamps = zeros (sz,4);
        idx = 1; 
        %This loop goes through all the messages and extracts all the major
        %events in the session - start trial, start of reward, end trial or
        %failed trial.
        for i = 1:size(eventTimes,1)
            
            if(ismember(messages(i,1), {'1  0  0  0  0  0  0  0'}) == 1) %start of the trial session 
                trialTimestamps (idx,1) = cell2mat(eventTimes(i,1));
            elseif(ismember(messages(i,1),{'0  0  0  0  0  1  1  0'})==1) %start of reward 
                trialTimestamps(idx,2) = cell2mat(eventTimes(i,1));
            elseif(ismember(messages(i,1), {'0  0  1  0  0  0  0  0'})==1) %end of trial session
                trialTimestamps(idx,3) = cell2mat(eventTimes(i,1));
                idx = idx+1;
            elseif(ismember(messages(i,1),{'0  0  0  0  0  1  1  1'}) == 1)%failed trial
                trialTimestamps(idx,4) = cell2mat(eventTimes(i,1));
            end 
        end
        %get rid of the zeros - rows
        trialTimestamps = trialTimestamps(any(trialTimestamps,2),:);
        if (size(any(trialTimestamps(:,4)) == 0))
            trialTimestamps(:,4) = [];
        end
        
        %USe trial timestamps to index into the eye positions so as to be
        %able to draw them in the plot function
        indices = trialTimestamps-double(startTime);
        indices (indices<0) = NaN;
        
        x = (edfdata.FSAMPLE.gx(1,:))';
        x(x>=1921)= NaN;
        %x(x<=-1)= NaN;

        y = (edfdata.FSAMPLE.gy(1,:))';
        y(y>=1081)= NaN;
        %y(y<=-1)= NaN;

        eyePos = horzcat(x,y);

        data.trial_timestamps = trialTimestamps; 
        data.indices = indices;
        data.eyePos = eyePos;
        data.noOfSessions = 1; 
    
        % create nptdata so we can inherit from it
        data.numSets = 1;    %eyelink is a session object = each session has only object 
        data.Args = Args;
        n = nptdata(data.numSets,0,pwd);
        d.data = data;
        obj = class(d,Args.classname,n);
        saveObject(obj,'ArgsC',Args);
        cd (cwd);

    end
             
else
	 % create empty object
	 obj = createEmptyEyelink(Args);
end

%---------------------------------------
%This function creates an empty object for storing eyelink info
function obj = createEmptyEyelink(Args)

% these are object specific fields
data.trial_timestamps = []; %contains all start, cue and and times for all the trials
data.sacc_event = [];
data.fix_event = [];
data.timestamps = []; %all the time that the experiments were run for
data.eye_pos = [];  %contains all the eye positions in all the sessions
%data.eyeEvent_time = eventTimestamps; %contains all the indices of all the events
data.noOfSessions = 0;
%data.exp_start_time = double(startTime);
data.timeouts = [];
data.noOfTrials= 0;


% create nptdata so we can inherit from it
% useful fields for most objects
data.numSets = 0;
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);


        
