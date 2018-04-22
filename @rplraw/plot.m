function [obj, varargout] = plot(obj,varargin)
%@rpllfp/plot Plot function for rplraw object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'FreqPlot',0, 'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','NormalizeTrial','FreqPlot'};
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
	if(Args.FreqPlot)
		PlotFFT(obj.data.analogData,obj.data.analogInfo.SampleRate);
		set(gca,'TickDir','out')
	else
		plot(obj.data.analogTime,obj.data.analogData,'.-')
	end
	if(~Args.LabelsOff)
		xlabel('Freq (Hz)')
		ylabel('Magnitude')
    end
    
    sdstr = get(obj,'SessionDirs');
	title(getDataOrder('ShortName','DirString',sdstr{n}))
else
	% plot all data
	plot(obj.data.analogTime,obj.data.analogData,'.-')
	if(~Args.LabelsOff)
		xlabel('Time (ms)')
		ylabel('Voltage (uV)')
	end
end

RR = eval('Args.ReturnVars');
for i=1:length(RR) RR1{i}=eval(RR{i}); end 
varargout = getReturnVal(Args.ReturnVars, RR1);
