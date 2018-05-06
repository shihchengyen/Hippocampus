function [obj, varargout] = plot(obj,varargin)
%@hidata/plot Plot function for hidata object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		'Array',0, 'Session',0, 'Day',0, 'UseGMR',0, ...
		'Objects',{''}, 'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','Array','Session','Day', 'UseGMR'};
[Args,varargin2] = getOptArgs(varargin,Args,'removenumargs');

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

UseGMRlayout = 0;

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	if(Args.Array)
		% plot waveforms array by array
		% get starting row in ChannelIndex
		chnstart = obj.data.ArrayIndex(n);
		% get ending row in ChannelIndex
		chnend = obj.data.ArrayIndex(n+1);
		if(Args.UseGMR)
			UseGMRlayout = 1;
		end
	elseif(Args.Session)
		% plot waveforms session by session
		% get starting row in ChannelIndex
		chnstart = obj.data.SessionIndex(n);
		% get ending row in ChannelIndex
		chnend = obj.data.SessionIndex(n+1);
		if(Args.UseGMR)
			UseGMRlayout = 1;
		end
	elseif(Args.Day)
		% plot waveforms session by session
		% get starting row in ChannelIndex
		chnstart = obj.data.DayIndex(n);
		% get ending row in ChannelIndex
		chnend = obj.data.DayIndex(n+1);
	else
		chnstart = n - 1;
		chnend = n;
	end  % if(Args.Array)

	% get directories
	sdstr = get(obj,'SessionDirs');
	if(UseGMRlayout)
		% get number of objects
		nobjs = size(Args.Objects,1);
		% grab current directory
		cwd = pwd;	
		for pidx = (chnstart+1):chnend
			% get directory
			cdir = sdstr{pidx};
			% select the appropriate subplot
			subplotGMR(cdir,Args)
			% change to the appropriate directory
			cd(cdir)
			% load objects
			for i = 1:nobjs
				try
		            % get the number of cols for this object. Check here so the
		            % user avoid having to type empty cell arrays when it is not
		            % needed
		            [oRows,oCols] = size(Args.Objects);
		            % check if there is a 3rd column in this object
		            if(oCols>2)
		                % instantiate object with arguments in 3rd column of Objects
						thisObj = feval(Args.Objects{i,1},'auto',Args.Objects{i,3}{:});
		            else
		                % instantiate object with just the 'auto' argument
						thisObj = feval(Args.Objects{i,1},'auto');
		            end
		            if(~isempty(thisObj))
						if(oCols>1)
							plot(thisObj,Args.Objects{i,2}{:});
						else
							plot(thisObj);
						end
		            end  % if(~isempty(thisObj))
					clear thisObj
				catch
					% get last error
					lm = lasterr;
					fprintf('Warning: Problem plotting %s object!\n', ...
						Args.Objects{i,1});
					fprintf('%s\n',lm);
				end  % catch
			end  % for i = 1:nobjs
			cd(cwd)
		end	 % for pidx = (chnstart+1):chnend	
	else  % if(UseGMRlayout)
		% create a new nptdata object using the relevant directories
		nd = nptdata('SessionDirs',sdstr( (chnstart+1):chnend ));
		plot(nd,varargin2{:});
	end  % if(UseGMRlayout)
else  % if(~isempty(Args.NumericArguments))
	% plot all data
	sdstr = get(obj,'SessionDirs');
	% create a new nptdata object using the relevant directories
	nd = nptdata('SessionDirs',sdstr( (chnstart+1):chnend ));
	plot(nd,varargin2{:});
end  % if(~isempty(Args.NumericArguments))

% add code for plot options here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% @hidata/PLOT takes 'LabelsOff' as an example
if(~Args.LabelsOff)
	xlabel('Data Points')
	ylabel('Voltage (uV)')
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RR = eval('Args.ReturnVars');
lRR = length(RR);
if(lRR>0)
    for i=1:lRR
        RR1{i}=eval(RR{i});
    end 
    varargout = getReturnVal(Args.ReturnVars, RR1);
else
    varargout = {};
end

% layout for hippocampus channels
function subplotGMR(dirname,Args)	
	% allocate memory for layout
	alayout = zeros(124,1);
	slayout = alayout;
	
	% set up channel mapping for arrays
	a = reshape(1:40,8,5);
	b = rot90(a,-1);
	c = fliplr(b);
	clayout = c(2:4,8:-1:1)';
	clayout2 = c(5,[8 7]);
	clayout3 = c(1,6:-1:1);
	% array 1
	alayout(1:6) = c(1,7:-1:2);
	alayout(7:30) = clayout;
	alayout(31:32) = clayout2;
	% array 2
	alayout(33:38) = clayout3;
	alayout(39:62) = clayout;
	alayout(63:64) = clayout2;
	% array 3
	alayout(65:96) = alayout(33:64);
	% array 4
	alayout(97:102) = clayout3;
	alayout(103:118) = clayout(:,1:2);
	alayout(119:124) = c(4,7:-1:2);
	
	% array 1
	% alayout(1:6) = 7:-1:2;
	% alayout(7:14) = 16:-1:9;
	% alayout(15:22) = 24:-1:17;
	% alayout(23:30) = 32:-1:25;
	% alayout(31:32) = 40:-1:39;
	% array 2
	% alayout(33:38) = 6:-1:1;
	% alayout(39:46) = 16:-1:9;
	% alayout(47:54) = 24:-1:17;
	% alayout(55:62) = 32:-1:25;
	% alayout(63:64) = 40:-1:39;
	% array 3
	% alayout(65:70) = 6:-1:1;
	% alayout(71:78) = 16:-1:9;
	% alayout(79:86) = 24:-1:17;
	% alayout(87:94) = 32:-1:25;
	% alayout(95:96) = 40:-1:39;
	% array 4
	% alayout(97:102) = 6:-1:1;
	% alayout(103:110) = 16:-1:9;
	% alayout(111:118) = 24:-1:17;
	% alayout(119:124) = 31:-1:26;

	% set up channel mapping for session
	a = reshape(1:128,16,8);
	b = rot90(a,-1);
	c = fliplr(b);
	slayout(1:6) = c(2:7,1);
	slayout(7:118) = c(:,2:15);
	slayout(119:124) = c(2:7,16);
	
	% get channel str
	[fdir,cstr] = fileparts(dirname);
	% get channel number
	chnnum = sscanf(cstr,'channel%d');
	if(Args.Array)
		subplot(5,8,alayout(chnnum))
	elseif(Args.Session)
		subplot(8,16,slayout(chnnum))
	end
