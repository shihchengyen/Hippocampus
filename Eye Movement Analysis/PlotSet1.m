% Read edf file
edfdata = edfmex('P11_24.edf'); % SPECIFY FILE TO READ

% Get time when recording was started
startTime = edfdata.FEVENT(52).sttime; %recording starts at row 52 

% Get row index of reward messages ('0 0 0 0 0 1 1 0')
messages = {edfdata.FEVENT(:).message}';
indexReward = find(strcmp(messages,'0  0  0  0  0  1  1  0'));

for i = 1:size(indexReward)
    timeReward(i,2) = edfdata.FEVENT(indexReward(i)).sttime; % time at reward
    timeReward(i,1) = timeReward(i,2) - 600 + 1; % time at start of fixation (600ms prior to reward)
end

% Get gaze x-coord (replace values out of screen)
x = (edfdata.FSAMPLE.gx(1,:))';
x(x>=100000000)= NaN;

% Get gaze y-coord (replace values out of screen)
y = (edfdata.FSAMPLE.gy(1,:))';
y(y>=100000000)= NaN;

%Plot 9 fixation windows/areas
for i = 1:9,
    
    if i == 1, j = 0; k = 0;
    elseif i == 2, j = 1; k = 0;
    elseif i == 3, j = 2; k = 0;
    elseif i == 4, j = 0; k = 1;
    elseif i == 5, j = 1; k = 1;
    elseif i == 6, j = 2; k = 1;
    elseif i == 7, j = 0; k = 2;
    elseif i == 8, j = 1; k = 2;
    elseif i == 9, j = 2; k = 2;
    end 
    
    outX = [125+760*j 275+760*j 275+760*j 125+760*j 125+760*j]; outY = [125+340*k 125+340*k 275+340*k 275+340*k 125+340*k]; % x- and y-coordinates of outer square
    outer = fill(outX,outY,[0.5,0.5,0.5]); set(outer,'EdgeColor','None'); 
    hold on    
    innerX = [175+760*j 225+760*j 225+760*j 175+760*j 175+760*j]; innerY = [175+340*k 175+340*k 225+340*k 225+340*k 175+340*k]; % x- and y-coordinates of inner square
    inner = fill(innerX,innerY,[0,0,0]); set(inner,'EdgeColor','None'); 
    hold on
end

% Specify plot colours for fixation clusters (Set1 consists of 30 trials)
colour(1,:) = [0.7,0.1,0.2]; % cherry cobbler
colour(2,:) = [0,0.6,0.4]; % glorious green
colour(3,:) = [0.1,0.3,0.4]; % midnight muse
colour(4,:) = [0.7,0.3,0.1]; % cajun craze
colour(5,:) = [0.2,0.2,0.4]; % concord crush
colour(6,:)= [1,0.5,0]; % pumpkin pie
colour(7,:) = [0.3,0.2,0.2]; % early espresso
colour(8,:) = [0.9,0.2,0.2]; % poppy parade
colour(9,:) = [0.6,0.3,0.4]; % rich razelberry
colour(10,:) = [0.2,0.3,0.3]; % handsome hunter
colour(11,:) = [0.7,0.6,0.4]; % baked brown sugar
colour(12,:) = [0,0.5,0.7]; % pacific point
colour(13,:) = [0.7,0.3,0.3]; % bravo burgundy
colour(14,:) = [1,0.8,0.2]; % summer sun
colour(15,:) = [0.4,0.5,0.2]; % gumball green
colour(16,:) = [0.4,0.6,0.7]; % marina mist
colour(17,:) = [0.5,0.6,0.3]; % old olive
colour(18,:) = [0.9,0.2,0.5]; % melon mambo
colour(19,:) = [0.5,0.5,0.7]; % lovely lilac
colour(20,:) = [1,0.9,0]; % yoyo yellow
colour(21,:) = [0.2,0.4,0.5]; % not quite navy
colour(22,:) = [0.8,0.4,0.3]; % dusty durango
colour(23,:) = [0,0.5,0.5]; % island indigo
colour(24,:) = [0.5,0.4,0.3]; % soft suede
colour(25,:) = [0.8,0.6,0.5]; % creamy caramel
colour(26,:) = [0.8,0.6,0.8]; % orchid opulence
colour(27,:) = [1,0.6,0.3]; %peach parfait
colour(28,:) = [0.6,0.1,0.1]; % raspberry ripple
colour(29,:) = [0.2,0.2,0.4]; %concord crush
colour(30,:) = [0.7,0.7,0.3]; % kiwi kiss 

% Plot fixation clusters
eventStart= cell2mat({edfdata.FEVENT(:).sttime}');

for i = 1:30,
    lowBound = find(eventStart(1:indexReward(i)) <= timeReward(i,1)); lowBoundInd(i,1) = max(lowBound); % get row index of the event containing the start of 600ms fixation
    
    for j = lowBoundInd(i,1):indexReward(i,1)-1, % number of eye events in current trial
        
        if indexReward(i,1)-lowBoundInd(i,1) < 2 % only one eye event in the 600ms prior to reward 
        	if strcmp(edfdata.FEVENT(j).codestring,'STARTFIX') % check if it is a fixation or saccade
            	plot(x(timeReward(i,2)-600-startTime:timeReward(i,2)-startTime,1),y(timeReward(i,2)-600-startTime:timeReward(i,2)-startTime,1),'Color', colour(i,:),'MarkerSize',20); % only plot if it is a fixation event
                xlim([0 1920]);
                ylim([0 1080]);
                set(gca, 'Ydir', 'reverse') % reverse axis (note: eyelink starts with 0,0 at the top left corner)
                hold on;   
            else
            end
            
        else
            
            pointerStart = edfdata.FEVENT(j).sttime - startTime; % get index to start plotting in edfdata.FSAMPLE
            pointerEnd = edfdata.FEVENT(j).entime - startTime; % get index to stop plotting in edfdata.FSAMPLE
        
            if j == lowBoundInd(i,1) % look for first fixation event
                if strcmp(edfdata.FEVENT(j).codestring,'ENDFIX')
                    plot(x(timeReward(i,1)-startTime:pointerEnd,1),y(timeReward(i,1)-startTime:pointerEnd,1),'Color', colour(i,:),'MarkerSize',20);
                    xlim([0 1920]);
                    ylim([0 1080]);
                    set(gca, 'Ydir', 'reverse') % reverse axis (note: eyelink starts with 0,0 at the top left corner)
                    hold on;           
                else
                end

            elseif j > lowBoundInd(i,1) && j < indexReward(i,1)-1 % look for subsequent fixation events
                if strcmp(edfdata.FEVENT(j).codestring,'ENDFIX')
                    plot(x(pointerStart:pointerEnd,1),y(pointerStart:pointerEnd,1),'Color', colour(i,:),'MarkerSize',20);
                    xlim([0 1920]);
                    ylim([0 1080]);
                    set(gca, 'Ydir', 'reverse') % reverse axis (note: eyelink starts with 0,0 at the top left corner)
                    hold on;           
                else
                end

            elseif j == indexReward(i,1)-1 % plot if reward is received in the middle of a fixation ('STARTFIX' right before reward message) with preceding fixation events
                if strcmp(edfdata.FEVENT(j).codestring,'STARTFIX') 
                    plot(x(pointerStart:timeReward(i,2)-startTime,1),y(pointerStart:timeReward(i,2)-startTime,1),'Color', colour(i,:),'MarkerSize',20);
                    hold on;
                else
                end
            end
        end
    end
end

