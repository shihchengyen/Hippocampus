function [unityData] = synchronise(unityData, markers)

 try 
    el = load ('eyelink.mat');
 catch
     fprintf('No eyelink\n');
 end 
 
trialTimestamps = el.el.data.trial_timestamps;

%find the durations between the different 
eldurations(:,1) = [0;trialTimestamps(2:end,1)-trialTimestamps(1:end-1,3)];
eldurations(:,2) = trialTimestamps(1:end,2) - trialTimestamps(1:end,1);
eldurations(:,3) = [trialTimestamps(1:end,3) - trialTimestamps(1:end,2)];

%convert to seconds
eldurations = eldurations/1000;

umdurations = zeros(size(markers));
sz = size(markers, 1);
uD = zeros(size(unityData,1),size(unityData,2));

for i = 1:sz
   
    %For one trial
    
    %Find the number of elements, and total duration between start and cue 
    noElements2 = size(unityData(markers(i,1):markers(i,2),2),1);
    umdurations(i,2) = sum(unityData(markers(i,1):markers(i,2),2));
    
    %Find the number of elements, and total duration between cue and end 
    noElements3 = size(unityData(markers(i,2):markers(i,3),2),1);
    umdurations(i,3) = sum(unityData(markers(i,2):markers(i,3),2));
    
    %Find the number of elements and total duration between prev end and
    %this start 
    if(i==1)
        umdurations(i,1) = 0;
        uD (1:markers(i,1),2)=unityData(1:markers(i,1),2);
    else
        noElements1 = size(unityData(markers(i-1,3):markers(i,1),2),1);
        umdurations(i,1) = sum(unityData(markers(i-1,3):markers(i,1),2));
    end 
    
    %Find the discrepancies between eyelink and unitymaze durations
     umdurations(i,:) = eldurations(i,:) - umdurations(i,:);
     
     %Divide the discrepancies with the found numer of Elements for each
     %part of the trial, and add this to the unityData data to correct it
     uD(markers(i,1):markers(i,2),2) = unityData(markers(i,1):markers(i,2),2) + (umdurations(i,2))/noElements2;
     uD(markers(i,2):markers(i,3),2) = unityData(markers(i,2):markers(i,3),2) + (umdurations(i,3))/noElements3; 
     if(i~=1)
         uD(markers(i-1,3):markers(i,1),2) =  unityData(markers(i-1,3):markers(i,1),2) + (umdurations(i,1))/noElements1; 
     end 
     
end 

uD(markers(end,end):end,2)=unityData(markers(end,end):end,2);
uD(:,1) = unityData(:,1);
uD (:,3:5)= unityData(:,3:5);

unityData = uD;
%correct unityData
% unityData (markers(end,end)+1:end,:) = [];
% unityData (1:markers(1,1)-1,:) = [];



