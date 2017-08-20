function Distance2(data)

global tdist

%% Remove rows with duplicate x-, y-coordinates and angle (i.e. animal is stationary in the virtual environment)
% Vectorized version
% extract position data
posdata = data(:,3:5);
% take difference between rows in each column
diffdata = diff(posdata);

% sum differences across cols, so if there were no changes in any of the
% columns, the sum will be zero
dds = sum(diffdata,2); % second argument '2' returns column vector, '1' returns row vector
% start with 1st row, then find the rows with non-zero sums (i.e. some
% change in at least one of the columns, and then add 1 to the row number
pdata = posdata([1; find(dds)+1],:);

%% Calculator magnitude of each vector
% Vectorized version
% compute dx and dy
dv = diff(pdata(:,1:2));
% take the square, i.e. dx^2 and dy^2
dv2 = dv.^2;
% take the sum, i.e. dx^2 + dy^2
dpf = dv2 * [1; 1];
% take the sqrt
dpfs = sqrt(dpf);
tdist = sum(dpfs);

%disp(tdist);
