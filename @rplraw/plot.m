function [obj, varargout] = plot(obj,varargin)
%@rplraw/plot Plot function for rplraw object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'FFT',0, 'XLims',[0 10000], 'Cmds','', ...
		  'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','NormalizeTrial','FFT'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

% analogTime needs to be in double-precision in order to avoid rounding
% error, but we don't really want to save it on disk as it takes up too
% much space and makes loading the object slower. So we will create it the
% first time we need it
% check if the field analogTime exists OR if the analogTime field exists
% and is not a double precision variable
% we will create the variable here.
if( ~isfield(obj.data,'analogTime') || ~isa(obj.data.analogTime,'double') )
	obj.data.analogTime = (0:(obj.data.analogInfo.NumberSamples-1))' ./ obj.data.analogInfo.SampleRate;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	if(Args.FFT)
		PlotFFT(obj.data.analogData,obj.data.analogInfo.SampleRate);
		if(~Args.LabelsOff)
			xlabel('Freq (Hz)')
			ylabel('Magnitude')
	    end
		xlim(Args.XLims)
	else
		plot(obj.data.analogTime *1000,obj.data.analogData,'.-')
		if(~Args.LabelsOff)
			xlabel('Time (ms)')
			ylabel('Voltage (uV)')
		end
	end
else
	% plot all data
	if(Args.FFT)
		PlotFFT(obj.data.analogData,obj.data.analogInfo.SampleRate);
		if(~Args.LabelsOff)
			xlabel('Freq (Hz)')
			ylabel('Magnitude')
	    end
		xlim(Args.XLims)
	else
		plot(obj.data.analogTime *1000,obj.data.analogData,'.-')
		if(~Args.LabelsOff)
			xlabel('Time (ms)')
			ylabel('Voltage (uV)')
		end
	end
end

set(gca,'TickDir','out')
sdstr = get(obj,'SessionDirs');
title(getDataOrder('ShortName','DirString',sdstr{1}))

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
