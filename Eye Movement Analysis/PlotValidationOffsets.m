% Read edf file
edfdata = edfmex('P11_28.edf'); % SPECIFY FILENAME

% Get validation offsets
for i = 31:39, %CHECK ROWS IN EDFDATA.FEVENT
    val = edfdata.FEVENT(i).message; % get validation messages
    val = regexp(val,'([^ ,:]*)','tokens'); % split validation messages by delimiters
    offsetX = str2double(char(val{1,12})); offsetY = str2double(char(val{1,13})); % get x and y offsets
    offsets(i-30,:) = [offsetX offsetY]; % concatenate offsets to single matrix
end

% Create static plot
static(:,1) = [216 960 1703 115 960 1804 216 960 1703]; % static x-coordinates
static(:,2) = [145 92 145 540 540 540 934 987 934]; % static y-coordinates
%static(:,2) = abs(static(:,2) - 1080);

for i = 1:size(static, 1); % plot static and measured coordinates
    plot(static(i,1), static(i,2),'r.','MarkerSize',10);
    hold on
    plot(static(i,1)-offsets(i,1), static(i,2)-offsets(i,2),'g.','MarkerSize',10);
    xlim([0 1920]);
    ylim([0 1080]);
    set(gca, 'Ydir', 'reverse') % reverse axis (note: eyelink starts with 0,0 at the top left corner)
    hold on
end