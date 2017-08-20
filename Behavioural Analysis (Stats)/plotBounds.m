%% plots robot trajectory
% boundaries work for Tee Left/Tee Left Block/Tee Right/Tee Right Block
% solid lines indicate completed trial; dotted lines indicate timeout
% x marks indicate rewards

clear all
close all
clc

xBound = [-2.56 -2.56 -12.72 -12.72 12.28 12.28 2.23 2.23 -2.56];
zBound = [-12.37 4.95 4.95 10.09 10.09 4.95 4.95 -12.37 -12.37];

[files, path] = uigetfile('*.txt','select data file','MultiSelect','on');

if iscell(files)
    iter = length(files);
else
    cellstr(files);
end

totalTO = 0;
totalTrials = 0;

totalTData = [];
totalTDataTO = [];

for j = 1:iter
    file = files{j};
    data = tblread([path file]);
    
    trials = find(data(:,1) == 1);
    trialEnd = find(data(:,1) == 3);
    timeOut = find(data(:,1) == 4);
    reward = find(data(:,1) == 2);
    tData = [];
    tDataTO = [];
    
    figure
    hold all
    plot(xBound,zBound,'k','LineWidth',2)
    
    if isempty(timeOut)
        
        for i = 1:length(trials)
            plot(data(trials(i):trialEnd(i),3),data(trials(i):trialEnd(i),4),'LineWidth',1.5)
            plot(data(trialEnd(i),3),data(trialEnd(i),4),'gx','LineWidth',1.5)
            
            tData{i} = data(trials(i):trialEnd(i),:);
        end
    else
        
        for i = 1:length(trialEnd)
            trialStart = find(data(1:trialEnd(i),1) == 1,1,'last');
            plot(data(trialStart:trialEnd(i),3),data(trialStart:trialEnd(i),4),'LineWidth',1.5)
            plot(data(trialEnd(i),3),data(trialEnd(i),4),'gx','LineWidth',1.5)
            
            tData{i} = data(trialStart:trialEnd(i),:);
        end
        
        for i = 1:length(timeOut)
            trialStart = find(data(1:timeOut(i),1) == 1,1,'last');
            plot(data(trialStart:timeOut(i),3),data(trialStart:timeOut(i),4),':','LineWidth',1.5)
            
            tDataTO{i} = data(trialStart:timeOut(i),:);
        end
    end
    
    for i = 1:length(reward)
        plot(data(reward(i),3),data(reward(i),4),'gx','LineWidth',1.5)
    end
    
    title(file)
    
    totalTO = totalTO + length(timeOut);
    totalTrials = totalTrials + length(trialEnd) + length(timeOut);
    
    totalTData{j} = tData;
    totalTDataTO{j} = tDataTO;
    
    fid = fopen([path file]);
    sessionType = fgets(fid);
    sessionTypeAll{j} = sessionType;
end

f = figure;
hold all
plot(xBound,zBound,'k','LineWidth',2)

for i = 1:length(totalTData)
    for j = 1:length(totalTData{i})
        plot(totalTData{i}{j}(:,3),totalTData{i}{j}(:,4),'LineWidth',1.5)
    end
end

for i = 1:length(totalTDataTO)
    for j = 1:length(totalTDataTO{i})
        plot(totalTDataTO{i}{j}(:,3),totalTDataTO{i}{j}(:,4),':','LineWidth',1.5)
    end
end
title('All Trials')

success = 1 - totalTO/totalTrials


uisave({'files','sessionTypeAll','totalTData','totalTDataTO','totalTO','totalTrials','xBound','zBound'})
[fileName, filePath] = uiputfile('*.mat');
saveas(f,[filePath fileName(1:end-4)],'fig')
