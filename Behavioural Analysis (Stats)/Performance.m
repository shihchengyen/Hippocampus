function Performance(data)

global perf

perf = data(size(data,1),1); % get performance (reward vs timeout)
%disp(perf)