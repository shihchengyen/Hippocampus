function data = arrangeMarkers(data)
% Marker formats
% A - 1 (start)
%   - 3 (end trial correct)
%   - 4 (end trial time-out)
% e.g. [1     3     1     3     1     3     1     3     1     3]
% B - 1 (start)
%   - 3 (end trial correct)
%   - 4 (end trial time-out)
%   - 0 (found between separate markers)
% e.g. [1     0     4     0     1     0     4     0     1     0]
% C - 1 (start/cue on)
%   - 2 (cue off)
%   - 3 (end trial correct)
%   - 4 (end trial time-out)
%   - 0 (found between separate markers)
% e.g. [1     0     2     0     3     0     1     0     2     0]
% D - 1x (start/cue on, with x (i.e. 1, 2, 3, or 4) as the target)
%   - 2x (cue off, with x as the target)
%   - 3x (end trial correct, with x as the target)
%   - 4x (end trial time-out, with x as the target)
%   - 0 (found between separate markers)
% e.g. [11     0    21     0    31     0    13     0    23     0]
% 204 - 1x (start/cue on, with x (i.e. 1, 2, 3, 4, 5, or 6) as the target)
%   - 2x (cue off, with x as the target)
%   - 3x (end trial correct, with x as the target)
%   - 4x (end trial time-out, with x as the target)
%   - 0 (found between separate markers)
% e.g. [204	0	11     0    21     0    31     0    13     0    23     0]
% 84 - 1x (start/cue on, with x (i.e. 1, 2, 3, 4, 5, or 6) as the target)
%   - 2x (cue off, with x as the target)
%   - 3x (end trial correct, with x as the target)
%   - 4x (end trial time-out, with x as the target)
%   - 0 (found between separate markers)
% e.g. [84	0	11     0    21     0    31     0    13     0    23     0]

rawMarkers = data.markers;
% reshape into 2 rows
rm1 = reshape(rawMarkers,2,[]);
% check if there are 0's between markers
if(rm1(2,1)==0)
	% not Format A, but still need to determine which format the markers are in
	if(rm1(1,1)>1)
		switch(rm1(1,1))
			case(204)
				% Format 204
				% remove 1st column before reshaping
				data.markers = reshape(rm1(1,2:end),3,[])';
				% get start time for each trial
                % skip the first two, since the first indicates version
                % number, while the second is a 0 to clear the previous
                % marker
				rtime = data.timeStamps(3:end);
				rt1 = reshape(rtime,6,[]);
				data.timeStamps = rt1([1 3 5],:)';    
				data.trialIndices = floor(data.timeStamps*data.SampleRate);				
			case(84)
				% Format 204
				% remove 1st column before reshaping
				data.markers = reshape(rm1(1,2:end),3,[])';
				% get start time for each trial
                % skip the first two, since the first indicates version
                % number, while the second is a 0 to clear the previous
                % marker
				rtime = data.timeStamps(3:end);
				rt1 = reshape(rtime,6,[]);
				data.timeStamps = rt1([1 3 5],:)';    
				data.trialIndices = floor(data.timeStamps*data.SampleRate);				
			otherwise
				% Format D
				% format is: 1x for cue onset/start of trial; 2x for cue offset; 
				% and 3x or 4x for reward or error/timeout, where x is the correct
				% stimulus for that trial
				% so we will reshape markers into 3 columns
				data.markers = reshape(rm1(1,:),3,[])';
				% get start time for each trial
				rtime = data.timeStamps;
				rt1 = reshape(rtime,6,[]);
				data.timeStamps = rt1([1 3 5],:)';    
				data.trialIndices = floor(data.timeStamps*data.SampleRate);				
			end
	else
		% check if 2nd markers is cue off (2)
		if(rm1(1,2)==2)
			% Format C
			% format is: 1 for cue onset/start of trial; 2 for cue offset; 
			% and 3 or 4 for reward or error/timeout
			% so we will reshape markers into 3 columns
			data.markers = reshape(rm1(1,:),3,[])';
			% get start time for each trial
			rtime = data.timeStamps;
			rt1 = reshape(rtime,6,[]);
			data.timeStamps = rt1([1 3 5],:)';    
			data.trialIndices = floor(data.timeStamps*data.SampleRate);				
		else
			% Format B
			% format is 1 marker for start of trial, and another for either
			% correct or incorrect, so reshape into 2 columns
			data.markers = reshape(rm1(1,:),2,[])';
			% get start time for each trial
			rtime = data.timeStamps;
			% reshape to remove timestamps for 0's
			rt1 = reshape(rtime,2,[]);
			% select 1st row to remove 0's, and the reshape timestamps into 
			% 2 columns
			data.timeStamps = reshape(rt1(1,:),2,[])';
			% compute the data indices corresponding to the marker timestamps
			data.trialIndices = floor(data.timeStamps*data.SampleRate);				
		end
	end
else
	% Format A, so transpose and assign to data structure
	data.markers = rm1';
    data.timeStamps = reshape(data.timeStamps,2,[])';
    % get start time for each trial
    data.trialIndices = floor(reshape(data.timeStamps*data.SampleRate,2,[])');
end
