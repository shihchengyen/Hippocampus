%callEyelink function:
%(1)Checks if rplparallel object is missing information, and if it is,
%use the eyelink data to complete the missing rpl markers. 
%(2)If eyelink is missing markers, and rplparallel has them, it ignores it
%(3)Common trials that are missing - they are deleted 
function [markersNew, timesNew] = callEyelink(markersRaw, messages, eltimes, rpltimeStamps ) %(markersRaw)

%stores the timestamps from eyelink 
eldurations = [0; eltimes(2:end)-eltimes(1:end-1)];
eldurations = double(eldurations);
eldurations = eldurations/1000; %conversion to seconds
eltimes = eldurations;

%Get rid of the 0s that separate the markers and the first 84 markers 
%Do the same for the rplparallel timestamps to reconstruct them.
if (markersRaw(1) == 84)
    markersRaw(1) = [];
    rpltimeStamps(1) = [];
end
rpltimeStamps(find(~markersRaw)) = [];
markersRaw(markersRaw==0) = [];

disp(size(markersRaw)); 
n = size(markersRaw,2); %stores all the trial elementsin it; 

if(n < size(messages,1)) 

    m = size(messages, 1);
    
    %first, check if edf file is missing something too
    remainder = rem(m, 3); 
    if (remainder ~=0)
        fprintf('Edf file incomplete\n ');
    end 
    markersNew = zeros (m, 1);
    timesNew = zeros(m,3); 
    idx = 1;
    idx2 = 1;
    sz = m+n;
    count = 1;
    
    for i=1:sz
       
        %Convert to row-column indexing for markersNew 
        r = floor(i/3) + 1; %row 
        c = rem(i,3); %column

        if (c == 0)
            r = r-1;
            c = 3;
        end
        
        if (idx2<=m)%prevent accessing more than the array size 
            %Convert the string to a number for easier manipulation
            message = (cell2str(messages(idx2,1)));
            message = str2num(message(end-4:end-2));
        end 
        
        if ((floor(message/10) == c) || (floor(message/10)==4 && c==3)) %ensures that messages itself isn't missing a trial
            if (idx<=n && message == markersRaw(1, idx)) %if the marker exists in rplparallel
                markersNew(r,c) = markersRaw(1, idx);
                timesNew(r,c) = rpltimeStamps(1, idx); 
                idx = idx+1;
            else %rplparallel is missing a marker
                fprintf('Missing in rplparallel but found in messages\n');
                if (c==1 && r~=1)
                    timesNew (r,c) = timesNew(r-1, 3)+eltimes(idx2,1);
                elseif(c==1 && r==1) %assume the first existing time is cue offset 
                    timesNew (r,c) = rpltimeStamps(1,idx) - eltimes(idx2+1,1); 
                else
                    timesNew(r,c) =timesNew(r,c-1) + eltimes(idx2,1);
                end 
                
                markersNew(r,c) = message;
            end 
            idx2 = idx2+1;
        else %check if markersRaw has it instead 
            if(floor(markersRaw(1, idx)/10) == 3 || floor(markersRaw(1, idx)/10) == 4 && c ==3) 
                if (c~= 1 && rem(markersRaw(1,idx), 10) == rem(markersRaw(1,idx-1),10))
                    markersNew(r,c) = markersRaw(1, idx); 
                    timesNew(r,c) = rpltimeStamps(1, idx); 
                    fprintf ('Missing Data from messages. but found in rplparallel\n');
                    disp([r,c])
                end 
            else 
                markersNew(r,c) = 0; %unity maze??? 
                timesNew(r,c) = 0;
                fprintf ('Use unitymaze\n');
                count = count+1;
            end
            idx = idx+1;
        end   
        
        %once done checking 
        if(idx>n && idx2>m) 
          break;
        end 
    end %end of for loop
    
%RPLPARALLEL MARKERS ARE MORE THAN EYELINK
else 
    
    m = size(messages,1);
    markersNew = zeros(n+m, 3);
    timesNew = zeros(n+m,3); 
    idx = 1; %parsing through markers 
    idx2 = 1; %parsing through messages 
    count = 0;
    sz = n+m;
    
    for i=1:sz
        
        %index for markersNew
        r = floor(i/3) + 1; %row 
        c = rem(i,3); %column

        if (c == 0)
            r = r-1;
            c = 3;
        end
        
        if(idx2<=m)
            message = (cell2str(messages(idx2,1)));
            message = str2num(message(end-4:end-2));
        end
        
        if(floor(markersRaw(1,idx)/10) == c) %if it is going in the correct column
           if (idx2<= m && message == markersRaw(1, idx)) %if messages has it we can move ahead 
               idx2 = idx2+1;
           else 
               count = count+1; %used in debugging to make sure the edf errors agree with those found here.
           end 
            markersNew (r,c) = markersRaw(1,idx);
            timesNew (r,c) = rpltimeStamps(1, idx);
            idx = idx+1;
        elseif (floor(markersRaw(1,idx)/10) == 4 && c == 3) %Timeout condition
            if (idx2<= m && message == markersRaw(1, idx)) %if messages has it we can move ahead 
               idx2 = idx2+1;
            end 
            markersNew (r,c) = markersRaw(1,idx);
            timesNew(r,c) = rpltimeStamps(1,idx); 
            idx = idx+1;
        else %it is missing from rplparallel 
            %check if messages has it 
            if(floor(message/10) == c || (floor(message/10) == 4 && c == 3)) %message has it 
                fprintf('Missing in rplparallel. But found in messages\n'); 
                markersNew(r,c) = message;
                if (c==1)
                    timesNew (r,c) = timesNew(r-1, 3)+eltimes(idx2,1);
                else
                    timesNew(r,c) =timesNew(r,c-1) + eltimes(idx2,1);
                end
                idx2 = idx2+1;
            else %even messages doesnt have it 
                fprintf('Delete Trial\n');
                markersNew(r,c) = 0;
                timesNew(r,c) = 0;
                count= count+1;
            end
            
        end
        
        %to ensure that we can break from the loop and don't need to waste
        %execution time
        if(idx > n && idx> m) 
          break;
        end 
    end
end

markersNew = markersNew(any(markersNew,2),:);  
timesNew = timesNew(any(timesNew,2),:);  
disp(count); 

%check if the last trial was missed by both 
if (~isempty(find(~markersNew)))
   fprintf('Last event missing. Use Unity Maze\n');
   try 
       um = unitymaze('auto');
   catch
       fprintf('Could not load unity maze\n');
       [row col] = find(~markersNew);
       markersNew(row, :) = []; 
       timesNew(row, :) = []; 
   end 
end 





