clear all

global tdist;
global tturn;
global ttime; % take note: time stamps in data files are inconsistent
global perf;

%Extract data folders 
currdir = dir(pwd);
dirFlags = [currdir.isdir]; %get a logical vector that tells which is a directory
subfolder = currdir(dirFlags); %extract those that are directories
subfolder(1:2) = []; %remove ','and '..' folders which refer to folders containing subfolders (current directory and up one directory)

acrossDays = cell(size(subfolder,1),1); %preallocate array to store summary data across days
sum = zeros(size(subfolder,1),7);
sumL = zeros(size(subfolder,1),7);
sumR = zeros(size(subfolder,1),7);
sumCL = zeros(size(subfolder,1),7);
sumCR = zeros(size(subfolder,1),7);
sumWL = zeros(size(subfolder,1),7);
sumWR = zeros(size(subfolder,1),7);

%Loop through data folders, get summary data for each day
for i = 1:size(subfolder,1);
    folder = strcat('C:\Users\Tabitha\Desktop\BehaviouralData Analysis (Stats)\', subfolder(i).name);
    cd(folder);
    dlist = dir('session*');
    
    withinDay = zeros(size(dlist,1),5); %preallocate array for DVs

    %Loop through session files, get summary data for each session
    for j = 1:size(dlist,1);
        
        data = dlmread(dlist(j).name,'',10,0); % start reading at row 11 (skip logged parameters)
        
        % get start and end row
        startrow = find(data(:,1) == 1); % look for trigger '1' (start trial) in first column
        endrow = find(data(:,1) > 2); % look for trigger '3' (end trial) or '4' (timeout) in last column
        data(endrow + 1:size(data,1),:) = [];% delete rows after end of trial
        data(1:startrow-1,:) = [];% delete rows before start of trial
        
        Distance(data);
        Angle(data);
        TimeTaken(data);
        Performance(data);

        %Get maze direction and type from text file header
        fileID = fopen(dlist(j).name,'r');
        header = textscan(fileID,'%s', 2, 'Delimiter','\n'); % read first two headerlines (to accommodate change in output log)
        
            if strfind(char(header{1}{1}),'Left') 
                mazedir = 0;
            elseif strfind(char(header{1}{2}),'Left') 
                mazedir = 0;
            elseif strfind(char(header{1}{1}),'Right') 
                mazedir = 1;
            elseif strfind(char(header{1}{2}),'Right') 
                mazedir = 1;
            end
            
        fclose(fileID);

        acrossDays(i,1) = cellstr(subfolder(i).name);
        withinDay(j,:) = [tdist tturn ttime perf mazedir]; % parameters: date, session info, total dist, total turn, total time taken, performance, maze direction
    end
    
    % Overall analysis
        avDist = mean(withinDay(:,1)); sdDist = std(withinDay(:,1));
        avTurn = mean(withinDay(:,2)); sdTurn = std(withinDay(:,2));
        avTime = mean(withinDay(:,3)); sdTime = std(withinDay(:,3));
        corrTrials = find(withinDay(:,4)==3); avPerf = size(corrTrials,1)/size(withinDay,1);
        
        sum(i,:) = [avDist sdDist avTurn sdTurn avTime sdTime avPerf];
        
    % Separate analysis for maze directions
        % mazedir = 0 (left)
        rowL = find(withinDay(:,5)==0);
        avDistL = mean(withinDay(rowL,1)); sdDistL = std(withinDay(rowL,1));
        avTurnL = mean(withinDay(rowL,2)); sdTurnL = std(withinDay(rowL,2));
        avTimeL = mean(withinDay(rowL,3)); sdTimeL = std(withinDay(rowL,3));
        corrTrialsL =find(withinDay(rowL,4)==3); avPerfL = size(corrTrialsL,1)/size(rowL,1);

        sumL(i,:) = [avDistL sdDistL avTurnL sdTurnL avTimeL sdTimeL avPerfL];
        
        % mazedir = 1 (right)
        rowR = find(withinDay(:,5)==1);
        avDistR = mean(withinDay(rowR,1)); sdDistR = std(withinDay(rowR,1));
        avTurnR = mean(withinDay(rowR,2)); sdTurnR = std(withinDay(rowR,2));
        avTimeR = mean(withinDay(rowR,3)); sdTimeR = std(withinDay(rowR,3));
        corrTrialsR =(find(withinDay(rowR,4)==3)); avPerfR = size(corrTrialsR,1)/size(rowR,1);

        sumR(i,:) = [avDistR sdDistR avTurnR sdTurnR avTimeR sdTimeR avPerfR];
    
    % Separate analysis for performance
        % correct trials = 3, left
        rowCL = find(withinDay(:,4)==3 & withinDay(:,5)==0);
        avDistCL = mean(withinDay(rowCL,1)); sdDistCL = std(withinDay(rowCL,1));
        avTurnCL = mean(withinDay(rowCL,2)); sdTurnCL = std(withinDay(rowCL,2));
        avTimeCL = mean(withinDay(rowCL,3)); sdTimeCL = std(withinDay(rowCL,3));
        corrTrialsCL =find(withinDay(rowCL,4)==3); avPerfCL = size(corrTrialsCL,1)/size(rowCL,1);
        
        sumCL(i,:) = [avDistCL sdDistCL avTurnCL sdTurnCL avTimeCL sdTimeCL avPerfCL];
        
        % correct trials = 3, right
        rowCR = find(withinDay(:,4)==3 & withinDay(:,5)==1);
        avDistCR = mean(withinDay(rowCR,1)); sdDistCR = std(withinDay(rowCR,1));
        avTurnCR = mean(withinDay(rowCR,2)); sdTurnCR = std(withinDay(rowCR,2));
        avTimeCR = mean(withinDay(rowCR,3)); sdTimeCR = std(withinDay(rowCR,3));
        corrTrialsCR =find(withinDay(rowCR,4)==3); avPerfCR = size(corrTrialsCR,1)/size(rowCR,1);
        
        sumCR(i,:) = [avDistCR sdDistCR avTurnCR sdTurnCR avTimeCR sdTimeCR avPerfCR];
        
        % wrong trials = 4, left
        rowWL = find(withinDay(:,4)==4 & withinDay(:,5)==0);
        avDistWL = mean(withinDay(rowWL,1)); sdDistWL = std(withinDay(rowWL,1));
        avTurnWL = mean(withinDay(rowWL,2)); sdTurnWL = std(withinDay(rowWL,2));
        avTimeWL = mean(withinDay(rowWL,3)); sdTimeWL = std(withinDay(rowWL,3));
        corrTrialsWL =find(withinDay(rowWL,4)==4); avPerfWL = size(corrTrialsWL,1)/size(rowWL,1);
        
        sumWL(i,:) = [avDistWL sdDistWL avTurnWL sdTurnWL avTimeWL sdTimeWL avPerfWL];
        
        % wrong trials = 4, right
        rowWR = find(withinDay(:,4)==4 & withinDay(:,5)==1);
        avDistWR = mean(withinDay(rowWR,1)); sdDistWR = std(withinDay(rowWR,1));
        avTurnWR = mean(withinDay(rowWR,2)); sdTurnWR = std(withinDay(rowWR,2));
        avTimeWR = mean(withinDay(rowWR,3)); sdTimeWR = std(withinDay(rowWR,3));
        corrTrialsWR =find(withinDay(rowWR,4)==4); avPerfWR = size(corrTrialsWR,1)/size(rowWR,1);
        
        sumWR(i,:) = [avDistWR sdDistWR avTurnWR sdTurnWR avTimeWR sdTimeWR avPerfWR];
        
end

%% Plot performance (bar graphs)
% x = 1:size(sumL,1);
% performanceL = sumL(:,7);
% w1 = 0.5;
% bar(x,distanceL,w1,'FaceColor',[0.2 0.2 0.5]);
% hold on
% performanceR = sumR(:,7);
% w2 = 0.25;
% bar(x,distanceR,w2,'FaceColor',[0 0.7 0.7]);

% grid on
% legend({'Left','Right'},'Location','northwest');
% ylabel('% Correct trials (reward obtained)');

% ax = gca;
% ax.Title.String = 'Performance (Forced alternating T-maze)';
% ax.Title.String = 'Performance(Cued alternating T-maze)';
% ax.Title.String = 'Performance (Cued alternating T-maze, stem reward removed)';
% ax.XTick = 1:size(sumL,1);
% ax.XTickLabels = acrossDays(:,1);
% ax.XTickLabelRotation = 90;
% ylim(ax,[0 1]);

%% Error function not working properly in Matlab, use excel to plot distance and angular graphs
