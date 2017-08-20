function Angle(data)

global tturn

%data = dlmread('session01_T04_1172016143727.txt','',1,0);

DiffArray = diff(data(:,5)); % create array containing difference values between consecutive rows for column 5 
DiffAngle = zeros(size(DiffArray,1),1); %preallocate array for DVs

for i = 1: size(DiffArray); %Note diff function calculates [X(2)-X(1), X(3)-X(2)...]
    %IMPORTANT: clarify threshold, currently set to +/-5
    if DiffArray(i)>10; % suggests counter-clockwise crossing of 0-360 boundary 
        DiffAngle(i)= (360-data(i+1,5)) + data(i,5); 
    elseif DiffArray(i)<-10; % suggests clockwise crossing of 0-360 boundary
        DiffAngle(i)= data(i+1,5) + (360-data(i,5)); 
    else
        DiffAngle(i)= DiffArray(i);
    end
end

tturn = sum(abs(DiffAngle));
%disp(tturn);