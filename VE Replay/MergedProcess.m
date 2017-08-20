%% Process frames
% Read output data file
data = dlmread('session01_T05_2010201611247.csv',',',10,1);

startrow = find(data(:,1) == 1); % get trial start row
endrow = (find(data(:,1) > 2) + 1); % get trial end row + one row buffer

dupLimit = round((data(startrow:endrow,2))*1000); % get number of times needed to duplicate each screen shot by

for i = 1:length(dupLimit)-1;
    dupLimit(i+1,2) = sum(dupLimit(1:i,1));
end

% Crop and duplicate frames
list = dir('*.png');
for j = 1:length(list);
    img = imread(list(j).name);
    imgC = imcrop(img,[1915 0 1913 1002]); % crop dimensions: [x-start y-start x-length y-length]
        
        for k = 1:dupLimit(j,1); % duplicate frames and save to array (array name: frame) in order
            frame{k+dupLimit(j,2)} = imgC;
        end
end

%% Overlay eye movements
% Load edf data
edfdata = edfmex('161020.edf');

% Get x-coord (replace values out of screen)
x = (edfdata.FSAMPLE.gx(1,:))';
x(x>=1913)= NaN;
x(x<=-1)= NaN;

% Get y-coord (replace values out of screen)
y = (edfdata.FSAMPLE.gy(1,:))';
y(y>=1032)= NaN;
y(y<=29)= NaN;
y = y-29;
%y = abs(y-1080); %do not invert y-axis if superimposing onto image

% Store in two column vector
eyePos = horzcat(x,y);

% Get row indices of 'ENDFIX' and 'ENDSACC' in FEVENT
events = {edfdata.FEVENT(:).codestring}';
indexFix = find(strcmp(events,'ENDFIX')); % get index of fixations (end fixations)
indexSacc = find(strcmp(events,'ENDSACC')); % get index of saccades (end saccades)
indexEvents = sortrows(vertcat(indexFix,indexSacc)); % merged index of all events

% Get row indices of eye events in FSAMPLE
startTime = edfdata.FEVENT(1).sttime; % get time when recording was started (time when 'MODE RECORD' message was sent)

timestamps = cell(size(indexEvents,1),2); % preallocate array to store start and end index of events
timestamps(:,1) = {edfdata.FEVENT(indexEvents).sttime}'; % extract start timestamps from FEVENT
timestamps(:,2) = {edfdata.FEVENT(indexEvents).entime}'; % extract end timestamps from FEVENT
timestamps = cell2mat(timestamps); 
timestamps = timestamps - startTime + 1; % get start and end index in FSAMPLE
timestamps(:,3) = timestamps(:,2) - timestamps(:,1) + 1; % get duration

% Mark out start and end events in eyePos 
eyePos(timestamps(:,1),3) = 1; % start 
eyePos(timestamps(:,2),3) = 2; % end
eyePos(timestamps(:,1),4) = timestamps(:,3); % duration 

% Select one trial only// JUST FOR THIS DEMO
% eyePos = eyePos(11391:38848,:);
timestamps = timestamps(654:854,:);

%For loop
%Get very first eye event in current trial
start = find(eyePos(:,3) == 1);

% Set up writerobj
writerObj = VideoWriter('Test.avi');
writerObj.FrameRate=100;
open(writerObj);

try
    for i = 1:size(timestamps,1);   
        disp(i);

        for j = 1:timestamps(i,3);

            img = frame{timestamps(i,1)-timestamps(1,1)+(j-1)+start(1)}; % read in image// REVISIT HOW START AND END OF TRIAL IS DEFINED
            ax1 = subplot(1,1,1);
            imshow(img,'Border','tight');
            hold(ax1,'on');

            plot(eyePos(timestamps(i,1)+(j-1),1),eyePos(timestamps(i,1)+(j-1),2),'g.','MarkerSize',70); % plot current eye position
            hold on;
            plot(eyePos(timestamps(i,1):timestamps(i,1)+(j-2),1),eyePos(timestamps(i,1):timestamps(i,1)+(j-2),2),'g.','MarkerSize',35); % plot all (earlier) eye positions of current event
            xlim([0 1913]);
            ylim([0 1002]);  
            pause(.0001);
            
            newImg = getframe(gcf);
            writeVideo(writerObj, newImg);

            hold off
        end
    end
catch
    close(writerObj)
end
close(writerObj);


