function [obj, varargout] = plot(obj,varargin)
%@rplhighpass/plot Plot function for rplhighpass object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0, 'GroupPlots',1, 'GroupPlotIndex',1, 'Color','b', ...
			'SpikeData',[], 'SpikeTriggerIndex',26, 'SpikeHeight',100, ...
		    'FFT',0, 'XLims',[0 10000], 'LoadSort',0, 'Cmds','', ...
		    'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','NormalizeTrial','LoadSort','FFT'};
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
	if(Args.FFT)
		PlotFFT(obj.data.analogData,obj.data.analogInfo.SampleRate)
		if(~Args.LabelsOff)
			xlabel('Freq (Hz)')
			ylabel('Magnitude')
		end
		xlim(Args.XLims)
	else
		plot(obj.data.analogTime * 1000,obj.data.analogData,'.-')
		if(Args.LoadSort)
			l = load('hmmsort.mat');
			Args.SpikeData = l.mlseq;
		end
		if(~isempty(Args.SpikeData))
			hold on
			% get SpikeData
			mlseq = Args.SpikeData;
			spidx = Args.SpikeTriggerIndex;
			ncells = size(mlseq,1);
			clist = nptDefaultColors(1:ncells);
			% add spike trains
			for spi = 1:ncells
				st1 = find(mlseq(spi,:)==spidx);
				% add stem plot
				stem( obj.data.analogTime(st1) * 1000, repmat(Args.SpikeHeight,[size(st1),1]), 'Color', clist(spi,:))
			end
			hold off
		end
		if(~Args.LabelsOff)
			xlabel('Time (ms)')
			ylabel('Voltage (uv)')
		end
	end
else
	% plot all data
	if(Args.FFT)
		PlotFFT(obj.data.analogData,obj.data.analogInfo.SampleRate)
		if(~Args.LabelsOff)
			xlabel('Freq (Hz)')
			ylabel('Magnitude')
		end
		xlim(Args.XLims)
	else
		plot(obj.data.analogTime * 1000,obj.data.analogData,'.-')
		if(~Args.LabelsOff)
			xlabel('Time (ms)')
			ylabel('Voltage (uv)')
		end
	end
end

sdstr = get(obj,'SessionDirs');
title(getDataOrder('ShortName','DirString',sdstr{1}))

if(~isempty(Args.Cmds))
    % save the current figure in case Args.Cmds switches to another figure
    h = gcf;
    eval(Args.Cmds)
    % switch back to previous figure
    figure(h);
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
