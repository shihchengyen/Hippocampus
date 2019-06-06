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

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0);
Args.flags = {'Auto','ArgsOnly'};
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

    %look for the edf file
    dlist = nptDir ('.edf');
    %get entires in the directory list
    dnum = size (dlist,1);
    %convert the edf file into a MATLAB accessible format 
	edfdata = edfmex (dlist(1).name);

    %get the experiment start time
    expTime = edfdata.FEVENT(1).sttime;

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

    %%Find the number of trials in total, and use this to construct
    %%sessions
    %%trialTimestamps is the vector storing all the sessions. Each session
    %%has 3 columns for Start/Cue/End and all the sessions are arranged
    %%columnwise 
    
    n = size (m, 1);
    n = ceil((n-2)/3); %ceil is to ensure that any decimal numbers may be rounded
    trialTimestamps = zeros(n, 3*noOfSessions);
    noOfTrials = zeros(1, noOfSessions);

    start = NaN (n,1);
    endtimes = NaN (n,1);
    cue = NaN (n,1);
    
     if (noOfSessions==1) %if there is only one session
            
         idx = 1; %counter of indices 
         wasCue = 0;
         i = 1;
         n = size (m,1); %to parse messages completely 
         
         %This while loop is to ensure that we can deal with missing cue
         %offset messages 
         while i < n %incrementing counter to parse through m
            
             %disp (idx);
             if (contains(m(i,1), {'Start Trial'}) == 1)
                start(idx,1) = edfdata.FEVENT(messageEvent(i,1)).sttime;
             elseif (contains(m(i,1), {'End Trial'}) == 1 || contains(m(i,1), {'Timeout'}) == 1)
                 endtimes (idx, 1) = edfdata.FEVENT(messageEvent (i,1)).sttime;
                 if (wasCue == 0)
                     cue(idx,1) = NaN; 
                 end
                 wasCue = 0;
                 idx = idx+1;
             elseif (contains (m (i,1), {'Cue Offset'}) == 1)
                 cue (idx,1) = edfdata.FEVENT(messageEvent(i,1)).sttime;
                 wasCue = 1;
             end 
             
             i = i+1;  
         end
         
 
        %now concactenate these into one matric
        noOfTrials = size (start, 1); %number of trials in the session
        trialTimestamps = horzcat(start, cue, endtimes);
         

     else %multiple sessions
         for i = 1:noOfSessions

            if ( i == noOfSessions) %this is to account for the last few indices to access from messageEvent 
                start = messageEvent (sessionIndex(i,1)+1:3:end-1);
                endtimes = messageEvent (sessionIndex(i,1)+2:3:end-1);
                cuetimes = messageEvent (sessionIndex(i,1)+3:3:end-1); %turn them into vertical matrices

                start = {edfdata.FEVENT(start(:)).sttime}';
                endtimes= {edfdata.FEVENT(endtimes(:)).sttime}';
                cuetimes = {edfdata.FEVENT(cuetimes(:)).sttime}';

                start = cell2mat (start);
                endtimes = cell2mat(endtimes);
                cuetimes= cell2mat (cuetimes);

                s1 = size(start,1);
                s2 = size (trialTimestamps,1);
                noOfTrials(1, i) = s1;

                if (s1 < s2) %ensure that all the matrices being concatenated are of the same length 
                    display('Concact');
                    start (s1+1:s2,:) = 0;
                    endtimes (s1+1:s2,:) = 0;
                    cuetimes (s1+1:s2,:) = 0;
                end

                %now concactenate these into one matric
                l = 1 + (i-1)*3;
                u = 3 + (i-1)*3;
                trialTimestamps (:, l:u) = horzcat(start, endtimes, cuetimes);

            else

                start = messageEvent (sessionIndex(i,1)+1:3:sessionIndex(i+1)-1);
                endtimes = messageEvent (sessionIndex(i,1)+2:3:sessionIndex(i+1)-1);
                cuetimes = messageEvent (sessionIndex(i,1)+3:3:sessionIndex(i+1)-1); %turn them into vertical matrices

                start = {edfdata.FEVENT(start(:)).sttime}';
                endtimes= {edfdata.FEVENT(endtimes(:)).sttime}';
                cuetimes = {edfdata.FEVENT(cuetimes(:)).sttime}';

                start = cell2mat (start);
                endtimes = cell2mat(endtimes);
                cuetimes= cell2mat (cuetimes);

                s1 = size(start,1) ;
                s2 = size (trialTimestamps,1);
                noOfTrials(1, i) = s1;

                if (s1 < s2)
                    display('Concact');
                    start (s1+1:s2, :) = 0;
                    endtimes (s1+1:s2, :) = 0;
                    cuetimes (s1+1:s2, :) = 0;
                end

                %now concactenate these into one matrix
                l = 1 + (i-1)*3;
                u = 3 + (i-1)*3;
                trialTimestamps (:,l:u)= horzcat(start, endtimes, cuetimes);
            end

         end

     end

     %edit the size of the array and remove all zero rows '
     if(size(trialTimestamps,2) > 3*noOfSessions)
       trialTimestamps(:,1:3*noOfSessions) = [];
     end

     trialTimestamps = trialTimestamps(any(trialTimestamps,2),:);

     %%Make a matrix with all the timeouts in all the trials in the session
     %%which we can check when we are graphing lines for the end trial (refer plot.m) 
     c = {'Timeout'};
     timeouts = contains (m,c);
     timeouts = {edfdata.FEVENT(messageEvent(find(timeouts))).sttime}';
     timeouts = cell2mat (timeouts); %stores all the timeouts in the session

     %%Now, get the time for all the eye events that occur in all
     %%experiments
     startTime = edfdata.FEVENT(1).sttime; % get time when recording was started (time when 'MODE RECORD' message was sent)
     
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
     data.trial_timestamps = trialTimestamps - double (startTime); %contains all start, cue and and times for all the trials
     data.sacc_event = sacc;
     data.fix_event = fix;
     data.timestamps = (edfdata.FSAMPLE.time)'; %all the time that the experiments were run for
     data.eye_pos = eyePos;  %contains all the eye positions in all the sessions
     %data.eyeEvent_time = eventTimestamps; %contains all the indices of all the events
     data.noOfSessions = noOfSessions;
     %data.exp_start_time = double(startTime);
     data.timeouts = timeouts;
     data.noOfTrials= noOfTrials;
   

     % create nptdata so we can inherit from it
	 data.numSets = noOfSessions;    %this needs to contain the number of sessions for the InspectGUI
     data.Args = Args;
	 n = nptdata(data.numSets,0,pwd);
	 d.data = data;
	 obj = class(d,Args.classname,n);
	 saveObject(obj,'ArgsC',Args);
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
