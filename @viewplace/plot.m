function [obj, varargout] = plot(obj,varargin)
%@viewplace/plot Plot function for viewplace object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		'Array',0, 'Session',0, 'Day',0, ...
		'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','Array','Session','Day'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	if(Args.Array || Args.Session || Args.Day)
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
		end  % if(Args.Array)
		% get number of channels in this array
		numSets = chnend - chnstart;
	    for index = 1:numSets
	        if(numSets>1)
	            nptSubplot(numSets,index);
	        end
			chnindex = chnstart + index;
			xind = (obj.data.ChannelIndex(chnindex)+1):obj.data.ChannelIndex(chnindex+1);
			plot((obj.data.spikeForms(xind,:))','.-')
			sdstr = get(obj,'SessionDirs');
			title(getDataOrder('ShortName','DirString',sdstr{chnindex}))
		end  % for index = 1:numSets
	else  % if(Args.Array || Args.Session || Args.Day)
		% display png file for each unit
		% first get the directory
		sidx = find(obj.data.ChannelIndex>=n);
		sdstr = get(obj,'SessionDirs');
		cwd = pwd;
		sidx1 = sidx(1) - 1;
		cd(sdstr{sidx1})
		% need to figure which picture inside this directory we need to load
		pnum = n - obj.data.ChannelIndex(sidx1);
		% load the picture
		a = imread([sprintf('unit%02d',pnum) '_placefield.png']);
		imshow(a)
		cd(cwd)
	end  % if(Args.Array || Args.Session || Args.Day)
else
	% plot all data
	stem(obj.data.locSI)
	ci = obj.data.ChannelIndex(2:end);
	nci = size(ci,1);
	cix = repmat((ci+0.5)',2,1);
	ciy = repmat(ylim',1,nci);
	line(cix,ciy)
	if(~Args.LabelsOff)
		xlabel('Cells')
		ylabel('Place Selectivity Index')
	end
end

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
