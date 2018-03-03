function submitSort(varargin)
	% submitSort Submits channels for spike sorting
	%   This function checks a spreadsheet to submit channels for spike sorting. This function should be called from
	%   session directory, or can be used in the following manner to submit several days of data for spike sorting:
	%       ProcessLevel(nptdata,'Levels','Days','Exclude',{'analog'},'nptLevelCmd',{'Session','submitSort'})
	%   This function uses the following optional arguments:
	%       CellLog	- 'Cell activity log.xlsx'
	%       ChannelsPerArray - 32
	%       SortCmd - 'source ~/.bash_profile; cp ~/Dropbox/Work/Matlab/hmmsort/hmmsort5.dag .; condor_submit_dag -maxpre 10
	%                  hmmsort5.dag'

Args = struct('CellLog','Cell activity log.xlsx','ChannelsPerArray',32, ...
	'SortCmd','source ~/.bash_profile; cp ~/Dropbox/Work/Matlab/hmmsort/hmmsort5.dag .; condor_submit_dag -maxpre 10 hmmsort5.dag');
Args.flags = {};

[Args,modvarargin] = getOptArgs(varargin,Args);

% change to the days directory to read the spreadsheet
[p,cwd] = getDataOrder('days','CDNow');

% read spreadsheet indicating channels with possible single units
num = xlsread(Args.CellLog);

% return to original directory
cd(cwd)
% get day of current directory
[p,dayname,e] = fileparts(getDataOrder('day'));
% get index of current directory
dayidx = find(num(1,:)==str2num(dayname));

% find rows with single units
ai =find(num(:,dayidx)==1);
% get the channel numbers
ch_nums = num(ai,1);

for chidx = 1:size(ch_nums,1)
	chan_num = ch_nums(chidx);
	% convert channel number to array number
	array_num = floor((chan_num-1)/Args.ChannelsPerArray)+1;
	array_dir = sprintf('array%02d',array_num);
	chan_dir = sprintf('channel%03d',chan_num);
	display([array_dir filesep chan_dir])
	cd(array_dir);
	cd(chan_dir);
	system(Args.SortCmd);
	cd(cwd)
end
