function TimeTaken(data)

global ttime

%data = dlmread('session04_T07_1172016145727.txt','',1,0);
lastrow = size(data(:,2),1);

seconds = mod(data(:,2),100); % get seconds
hrmin = floor(data(:,2)/100); minutes = mod(hrmin,100);% get minutes
hours = floor(data(:,2)/10000); % get hours

% Set dummy values for year, month, date, hours and minutes
t1 = [2016 12 31 hours(1) minutes(1) seconds(1)];
t2 = [2016 12 31 hours(lastrow) minutes(lastrow) seconds(lastrow)];
ttime = etime(t2,t1);

%e = etime(data(size(data,1),2), data(1,2));

%ttime = data(size(data,1),2)-data(1,2); %sum delta time in column 2
