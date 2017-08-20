% Load edf data
edfdata = edfmex('160920.edf');

% Get x-coord (replace values out of screen)
x = (edfdata.FSAMPLE.gx(1,:))';
x(x>=1921)= NaN;
x(x<=-1)= NaN;

% Get y-coord (replace values out of screen)
y = (edfdata.FSAMPLE.gy(1,:))';
y(y>=1081)= NaN;
y(y<=-1)= NaN;
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
timestamps = timestamps(52:251,:);

% Read in image (.png) list
dlist = dir('frame*');

%For loop
%Get very first eye event in current trial
start = find(eyePos(:,3) == 1);

for i = 1:size(timestamps,1);   
    
    for j = 1:timestamps(i,3);
        
        img = imread(dlist(timestamps(i,1)-timestamps(1,1)+(j-1)+start(1)).name); % read in image// REVISIT HOW START AND END OF TRIAL IS DEFINED
        ax1 = subplot(1,1,1);
        imshow(img,'Border','tight');
        hold(ax1,'on');
        
        plot(eyePos(timestamps(i,1)+(j-1),1),eyePos(timestamps(i,1)+(j-1),2),'g.','MarkerSize',70); % plot current eye position
        hold on;
        plot(eyePos(timestamps(i,1):timestamps(i,1)+(j-2),1),eyePos(timestamps(i,1):timestamps(i,1)+(j-2),2),'g.','MarkerSize',35); % plot all (earlier) eye positions of current event
        xlim([0 1920]);
        ylim([0 1080]);  
        pause(.0001);
        
        % Save screenshot
        screen = num2str(timestamps(i,1)-timestamps(1,1)+(j-1)+start(1), '%05d');
        newFilename = strcat('screen', screen);
        set(gcf,'PaperPositionMode','auto');
        print(newFilename,'-dpng','-r0'); % save as screenshot frame
        
        hold off
    end
   
%     else
%         plot(eyePos(i,1),eyePos(i,2),'g.','MarkerSize',80); % plot eye position
%         xlim([0 1920]);
%         ylim([0 1080]);  
%         pause(.0001);
%         hold on;
    
end

% % % Vectorized version
% % plot(eyePos(:,1),eyePos(:,2),'r.','MarkerSize',10);
% % xlim([0 1920]);
% % ylim([0 1080]);


