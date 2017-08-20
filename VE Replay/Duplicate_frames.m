data = dlmread('160920_Session1.csv',',',1,0);

startrow = find(data(:,2) == 1);
endrow = (find(data(:,2) > 2)) + 3; % give three extra frames buffer

frames = round((data(startrow:endrow,3))*1000); % get number of times needed to multiply each screen shot by

for i = 0:1074;
    num = int2str(i);
    filename = strcat(num,'.png');
    screenIm = imread(filename); % read in screen image file
    dupLimit = frames(i+1); % get number of times to duplicate the screen image file
    
    % Duplicate each image for the specified number of times
    for j = 1:dupLimit;
        num2 = num2str(i, '%04d');
        copy = num2str(j, '%02d');
        newFilename = strcat('frame', num2, '_', copy,'.png');
        imwrite(screenIm, newFilename);      
    end
end

