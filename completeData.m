%%This function fills in the missing timestamps and messages for incomplete
%%edf files. The function returns elTrials, which stores the normalised
%%time for the session, while missingData returns a table containing the
%%missing messages and timestamps. This is turned into a .csv file
%%containing all the missing events timestamps and messages across all
%%sessions in the day
function [elTrials, missingData] = completeData(edfdata, m, messageEvent, sessionName)

    %WE extract the start time of the session in case the start time of the 
    %1st trial is missing (the 1st index of messageEvent is the session
    %index
    sessionStart = edfdata.FEVENT(messageEvent(1,1)).sttime; 

    %the start of the experiment is taken to normalise the data 
    expTime = edfdata.FEVENT(1).sttime;
    expTime = uint32(expTime);
    
    %stores the starting time of all the events 
    eltimes = {edfdata.FEVENT(messageEvent(2:end)).sttime}'; 
    eltimes = cell2mat(eltimes); 
    eltimes = uint32(eltimes);
    
    %%Correct the data that has been extracted from edfdata for ONE session
    
    %create a new matrix that contains all trial messages only
    messages = m(2:end);
    
    %Make the ripple object to access its fields
    cd (sessionName)
    rpl = rplparallel('auto');
    cd ..
    
    %markers will store all the event numbers in the trial, as taken 
    %from the ripple object. This will be used to mark which events are missing
    %in th eyelink object. (1-start, 2-cue, 3/4 - end/timeout 
    markers = rpl.data.markers;
    n = size(markers, 1);
    noOfmessages = size(messages, 1); %stores the number of messages recorded by eyelink in the session
    rpltimeStamps = rpl.data.timeStamps; %stores the event time for ripple
    
     if (n*3 > noOfmessages) %if there are missing messages
         missing = zeros (n, size(markers,2)); %stores the event that is missing
         rpldurations = zeros (n, size(markers,2)); %stores the rpl durations for filling missing eyelink timestamps
         elTrials = zeros (n, size(markers,2)); %stores the eyelink timestamps for the events in a trial
     end
   
    %calculate the durations
    %rpldurations(1,1) is assumed to be the time difference between the
    %start of the session and the trial start time
    rpldurations (:,1) = [rpltimeStamps(1,1); rpltimeStamps(2:end, 1) - rpltimeStamps(1:end-1, 3)]; %end of last trial - start of this trial
    rpldurations (:,2) = rpltimeStamps(:,2) - rpltimeStamps(:,1); %cue time - start time
    rpldurations (:,3) = rpltimeStamps (:,3) - rpltimeStamps(:,2); %end - cue time
    
    idx = 1;
    n = n*3;
    newMessages = cell(n,1); %stores all the missing messages 
    
    for i = 1:n
        %convert a linear index into a 2D index 
        r = floor(i/3) + 1; %row 
        c = rem(i,3); %column
        
        if (c == 0)
            r = r-1;
            c = 3;
        end 
        
        if (contains(messages(idx,1), {num2str(markers(r,c))}) == 0) %if the trial is missing
            missing(r,c) = markers(r,c); %mark what event is missing (start(1), cue(2), or end(3/4)) 
             if (c==1)
                 if (r~=1) %if it is not the first trial and start time of the session
                     elTrials(r,c) = elTrials (r-1, 3) + 1000*rpldurations(r, c);
                 else %the start time of the first trial in the session is missing 
                     elTrials(r,c) = sessionStart + 1000*rpldurations(1,1); %the start time of trial + diff bet session and start time
                 end
                  newMessages(i,1) = {strcat('Start Trial', {' '}, num2str (markers(r,c)))};
             elseif (c==2) %if the cue offset time is missing
                 elTrials(r,c) = elTrials (r, 1) + 1000*rpldurations (r, c);
                 newMessages(i,1) = {strcat('Cue Offset', {' '}, num2str (markers(r,c)))};
             elseif (c==3) %if the end trial time is missing 
                 elTrials(r,c) = elTrials (r, 2) + 1000*rpldurations (r, c);
                 
                 if(markers(r,c) == 31)
                    newMessages(i,1) = {strcat('End Trial', {' '}, num2str (markers(r,c)))};
                 elseif(markers(r,c)==41)
                    newMessages(i,1) = {strcat('Timeout', {' '}, num2str (markers(r,c)))}; 
                 end
            end
            
        else %if the trial is found 
            missing (r,c) = 0;     
            elTrials (r,c) = eltimes(idx,1);
            idx = idx+1;
        end
    end
    
    %Now, turn the time into a csv file
    disp('Missing messages');   
    n= size(find(missing),1)
    elTrials = uint32(elTrials); %convert to uint to get rid of floating points
 
    if (n ~= 0) %if there are missing messages in this session, make the missingData matrix to add to the .csv file 
        type = cell(n,1);
        type(:,:) = {24};
        correctedTimes = num2cell(sortrows(elTrials(find(missing))));
        newMessages = newMessages(~cellfun('isempty',newMessages));
        missingData = horzcat(type, correctedTimes, newMessages); %this needs to be turned into a csv file 
        missingData = cell2table(missingData, 'VariableNames', {'Type', 'Timestamps', 'Messages'}); 
        clear correctedTimes;
    else %if there are not any missing messages 
        disp('No missing messages');
        missingData = cell(n,3);
        missingData (:,:) = {0};
        missingData = cell2table(missingData, 'VariableNames', {'Type', 'Timestamps', 'Messages'});
    end 
    
    elTrials = elTrials - expTime;
    
    %%We now plot a histogram of the discrepancies in the start-cue,
    %%cue-end and end-start durations for ripple and eyelink objects for
    %%the same session to ensure that the data recorded is consistent 
    eldurations (:,1) = [0; elTrials(2:end, 1) - elTrials(1:end-1, 3)]; %end of last trial - start of this trial
    eldurations (:,2) = elTrials(:,2) - elTrials(:,1); %cue time - start time
    eldurations (:,3) = elTrials (:,3) - elTrials(:,2); %end - cue time
    
    eldurations = eldurations/1000; %conversion to seconds
    eldurations = double(eldurations);
    discrepancies = abs(rpldurations - eldurations); %stores the discrepancy between the two times in seconds 
    
    %%Plot the distributions of the durations in ms
    figure (1)
    %edges = [lower:0.0001:upper];
    histogram (discrepancies, 50);
    hold on;
    title('Histogram for rplparallel and eyelink durations');
    xlabel ('s');
    ylabel ('occurence');
    hold off 
    
    