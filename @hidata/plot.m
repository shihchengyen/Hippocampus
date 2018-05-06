function [obj, varargout] = plot(obj,varargin)
%@hidata/plot Plot function for hidata object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		'Array',0, 'Session',0, 'Day',0, ...
		'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','Array','Session','Day'};
[Args,varargin2] = getOptArgs(varargin,Args,'removenumargs');

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	% if(Args.Array || Args.Session || Args.Day)
		if(Args.Array)
			% plot waveforms array by array
			% get starting row in ChannelIndex
			chnstart = obj.data.ArrayIndex(n);
			% get ending row in ChannelIndex
			chnend = obj.data.ArrayIndex(n+1);
		elseif(Args.Session)
			% plot waveforms session by session
			% get starting row in ChannelIndex
			chnstart = obj.data.SessionIndex(n);
			% get ending row in ChannelIndex
			chnend = obj.data.SessionIndex(n+1);
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
		% create a new nptdata object using the relevant directories
		nd = nptdata('SessionDirs',sdstr( (chnstart+1):chnend ));
		plot(nd,varargin2{:});
else
	% plot all data
	sdstr = get(obj,'SessionDirs');
	% create a new nptdata object using the relevant directories
	nd = nptdata('SessionDirs',sdstr( (chnstart+1):chnend ));
	plot(nd,varargin2{:});
end

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
