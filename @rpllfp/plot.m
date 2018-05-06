function [obj, varargout] = plot(obj,varargin)
%@rpllfp/plot Plot function for rpllfp object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'FFT',0, 'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','NormalizeTrial','FFT'};
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
	else
		plot(obj.data.analogTime * 1000,obj.data.analogData,'.-')
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
