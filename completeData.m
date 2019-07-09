%%This function fills in the missing timestamps and messages for incomplete
%%edf files. The function returns elTrials, which stores the normalised
%%time for the session, while missingData returns a table containing the
%%missing messages and timestamps. This is turned into a .csv file
%%containing all the missing events timestamps and messages across all
%%sessions in the day
function [elTrials, missingData, flag] = completeData(edfdata, m, messageEvent, sessionName, moreSessionsFlag)

    %WE extract the start time of the session in case the start time of the 
    %1st trial is missing (the 1st index of messageEvent is the session
    %index
    sessionStart = edfdata.FEVENT(messageEvent(1,1)).sttime; 

    %the start of the experiment is taken to normalise the data 
    expTime = edfdata.FEVENT(1).sttime;
    expTime = uint32(expTime);
   
    %%Correct the data that has been extracted from edfdata for ONE sessiom
    %create a new matrix that contains all trial messages only
    messageEvent(contains(m, 'Trigger')) = [];
    m(contains(m, 'Trigger')) = [];
    
    messageEvent(contains(m, 'end')) = [];
    m(contains(m, 'end')) = [];
    
    messageEvent(contains(m, 'ERROR')) = [];
    m(contains(m, 'ERROR')) = [];
    messages = m;
    
    %stores the starting time of all the events 
    eltimes = {edfdata.FEVENT(messageEvent(1:end)).sttime}'; 
    eltimes = cell2mat(eltimes); 
    eltimes = uint32(eltimes);
    
    %Make the ripple object to access its fields
    cd (sessionName)
    rpl = rplparallel('auto');
   
        
    if (rpl.data.numSets ~= 0 && isempty(rpl.data.timeStamps)~=1) %no missing rplparallel.mat
        
        %markers will store all the event numbers in the trial, as taken 
        %from the ripple object. This will be used to mark which events are missing
        %in th eyelink object. (1-start, 2-cue, 3/4 - end/timeout)
        markers = rpl.data.markers;
        rpltimeStamps = rpl.data.timeStamps; %stores the event time for ripple
        n = size(markers, 1);
        %CHECK IF THE RPLPARALLEL OBJECT IS FORMATTED CORRECTLY OR IS
        %MISSING INFORMATION
        if(n==1) %if the formatting is 1xSIZE
            df = rpl; 
            save('rplparallel0.mat', 'df');
            markers(1) = [];
            rpltimeStamps(1) = [];
            rpltimeStamps(find(~markers)) = [];
            markers(markers==0) = [];
            n = size(markers,2)/3;
            if (rem(size(markers),3) ~= 0)
               markers = rpl.data.markers;
               rpltimeStamps = rpl.data.timeStamps;
               [markers,rpltimeStamps] = callEyelink(markers, messages, eltimes-expTime, rpltimeStamps);
               
            else 
                markers = reshape (markers, [3 n]);
                rpltimeStamps = reshape(rpltimeStamps, [3 n]);
                markers = markers';
                rpltimeStamps = rpltimeStamps';
            end 
             n = size(markers, 1);
             
             df.data.markers = markers;
             df.data.timeStamps = rpltimeStamps;
             save('rplparallel.mat', 'df');
             
        elseif (n*3 < size(m,1)) %If rplparallel obj is missing data   
            df = rpl;  
            [markers,rpltimeStamps] = callEyelink(markers, messages, eltimes-expTime, rpltimeStamps);
            
            n = size(markers, 1);
            df.data.markers = markers;
            df.data.timeStamps = rpltimeStamps;
            save('rplparallel.mat', 'df');
        end
        %save object and then return
        cd ..;
        
        noOfmessages = size(messages, 1); %stores the number of messages recorded by eyelink in the session 
         
        missing = zeros (n, size(markers,2)); %stores the event that is missing
        rpldurations = zeros (n, size(markers,2)); %stores the rpl durations for filling missing eyelink timestamps
        elTrials = zeros (n, size(markers,2)); %stores the eyelink timestamps for the events in a trial
         
        %To check if there are more sessions than must be - put a session
        %NUMBER check as well to ensure that you aren't nerfing normal
        %files
        
        if(moreSessionsFlag ~= 0)
           if ((n*3) - size(messageEvent,1)-2 >= 100)
               flag = 1;
               missingData = cell(n,3);
               elTrials = cell (n,3);
               missingData = cell2table(missingData, 'VariableNames', {'Type', 'Timestamps', 'Messages'});
               return; 
           else 
               flag = 0;
           end 
        end 
        flag = 0; 
        %calculate the durations
        %rpldurations(1,1) is assumed to be the time difference between the
        %start of the session and the trial start time
        rpldurations (:,1) = [rpltimeStamps(1,1); rpltimeStamps(2:end, 1) - rpltimeStamps(1:end-1, 3)]; %end of last trial - start of this trial
        rpldurations (:,2) = rpltimeStamps(:,2) - rpltimeStamps(:,1); %cue time - start time
        rpldurations (:,3) = rpltimeStamps (:,3) - rpltimeStamps(:,2); %end - cue time
        
        idx = 1;
        n = n*3; %size of the markers 
        newMessages = cell(n,1); %stores all the missing messages 

        for i = 1:n
            %convert a linear index into a 2D index 
            r = floor(i/3) + 1; %row 
            c = rem(i,3); %column

            if (c == 0)
                r = r-1;
                c = 3;
            end
                
   
            if ((idx == size(messages,1) && i~=n) || (contains(messages(idx,1), {num2str(markers(r,c))}) == 0 && markers(r,c)~=0 )) %if the trial is missing
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
                     if(floor(markers(r,c)/10) == 3)
                        newMessages(i,1) = {strcat('End Trial', {' '}, num2str (markers(r,c)))};
                     elseif(floor(markers(r,c)/10)==4)
                        newMessages(i,1) = {strcat('Timeout', {' '}, num2str (markers(r,c)))}; 
                     end
                end
            elseif (markers(r,c) == 0)
                continue; 
            else %if the trial is found 
                missing (r,c) = 0;     
                elTrials (r,c) = eltimes(idx,1);
                idx = idx+1;
            end
            
        end 
        
        %SAVING THE OBJECT  
        %Now, turn the time into a csv file
        [row ~] = find(~elTrials);
        elTrials(row, :) = []; 
        disp('Missing messages');   
        n = size(find(missing),1)
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
        eldurations = double(eldurations);
        eldurations = eldurations/1000; %conversion to seconds
        
        discrepancies = abs(rpldurations - eldurations); %stores the discrepancy between the two times in seconds

        %%Plot the distributions of the durations in ms
        figure (1)
        histogram (discrepancies, 50);
        hold on;
        title('Histogram for rplparallel and eyelink durations');
        xlabel ('s');
        ylabel ('occurence');
        hold off 
        
    else %missing rplparallel  
        %you are going to assume that there are no missing messages
        cd ..; 
        fprintf ('Empty Object. Just fill up time array\n');
        n = size(messages,1);
        elTrials = zeros (n, 3); %stores the eyelink timestamps for the events in a trial
        missing = zeros(n,3); 
        
        for i = 1:n
        
            r = floor(i/3) + 1; %row 
            c = rem(i,3); %column

            if (c == 0)
                r = r-1;
                c = 3;
            end 
           
            missing (r,c) = 0;  
            if (contains(messages(i,1), {'Start Trial'})==1)
                elTrials (r,1) = eltimes(i,1) - expTime;
            elseif (contains(messages(i,1), {'End Trial'}) == 1)
                elTrials (r,3) = eltimes(i,1) - expTime;
            elseif (contains(messages(i,1), {'Timeout'}) == 1)
                elTrials (r,3) = eltimes(i,1) - expTime;
            elseif (contains(messages(i,1), {'Cue Offset'}) == 1)
                elTrials (r,2) = eltimes(i,1) - expTime;
            end
            
        end
        
        disp('No missing messages');
        missingData = cell(n,3);
        missingData (:,:) = {0};
        missingData = cell2table(missingData, 'VariableNames', {'Type', 'Timestamps', 'Messages'});
    

    end 

