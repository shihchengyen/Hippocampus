function [sort_index, gpinfo] = groupdata(data)
%groupdata Sorts and groups data
%   [SORT_INDEX, INDEX_INFO] = groupdata(A) sorts the data in A, and returns the indices
%   in SORT_INDEX so that A(SORT_INDEX) returns sort(A). INDEX_INFO is a matrix with the 
%   unique values in A sorted in increasing order in the first column, the first index in 
%   SORT_INDEX that contains that unique value in the second column, the last index in
%   SORT_INDEX that contains that unique value in the third column, and the number of 
%   indices that contain that unique value in the fourth column.
%   
%   Example: 
%      >> r = [1 4 3 1 2 4 1 4 2 5 2 4 1 2 1 3 4 3 3 4];
%      >> [si,gi] = groupdata(r');
%      >> gi =
%
%      1     1     5     5
%      2     6     9     4
%      3    10    13     4
%      4    14    19     6
%      5    20    20     1
%
%      >> r(si(gi(1,2):gi(1,3)))
%      ans =
%
%      1     1     1     1     1      

% get length of data
[data_length,data_width] = size(data);

% return error if the data is not a column vector
if(data_width>1)
	display('This function only works with single-column data')
	gpinfo = [];
	sort_index = [];
else
	% first sort the data
	[sorteddata,sort_index] = sort(data);
	% find the indices where the sorted values change, add 1 to account for the change in
	% indices when using diff, but pre-pend 1 to indicate the start of the first unique value
	group_start_index = [1; find(diff(sorteddata)~=0) + 1];
	% fill in the end index for each unique value
	start_end_indices = [ group_start_index [ [group_start_index(2:end)-1]; data_length]];
	% add in the number of indices in the 3rd column
	group_numbers = diff(start_end_indices,1,2) + 1;
	% add in the unique values in the 1st column
	gpinfo = [sorteddata(group_start_index) start_end_indices group_numbers];
end
