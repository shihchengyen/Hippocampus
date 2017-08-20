global tdist
global tturn
global ttime

% Get user input of day to process
    dayInfo = int2str(input('Enter day to analyze: '));
    folder = strcat('C:\Users\Tabitha\Desktop\BehaviouralData Analysis (Trajectory)\', dayInfo);
    cd(folder);

dlist = dir('session*');

for i = 1:size(dlist,1);
    
    % read file but skip 1st row
    data = dlmread(dlist(i).name,'',1,0);

    % get start and end row
    startrow = find(data(:,1) == 1); % look for trigger '1' (start trial) in first column
    endrow = find(data(:,1) > 2); % look for trigger '3' (end trial) or '4' (timeout) in last column
    data(endrow + 1:size(data,1),:) = [];% delete rows after end of trial
    data(1:startrow-1,:) = [];% delete rows before start of trial

Distance(data); 
Angle(data); 
TimeTaken(data); 

dist = strcat('Total distance: ', num2str(tdist)); disp(dist);
turn = strcat('Total angular turn: ', num2str(tturn)); disp(turn);
time = strcat('Total time taken: ', num2str(ttime)); disp(time);

ReadUnity(data);
pause;
hold off; %refresh plot

end