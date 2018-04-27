function [obj, varargout] = plot(obj,varargin)
%@rplparallel/plot Plot function for rplparallel object.
%   OBJ = plot(OBJ) creates a plot of the Ripple markers

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
			'M1', 20, 'M2', 30, ...
		  'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly'};
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
	stem(obj.data.markers(3:2:end))
	line(xlim,repmat(Args.M1,1,2))
	line(xlim,repmat(Args.M2,1,2),'Color','r')
else
	% plot all data
	panGUI
	stem(obj.data.markers(3:2:end))
	line(xlim,repmat(Args.M1,1,2))
	line(xlim,repmat(Args.M2,1,2),'Color','r')
end

% add code for plot options here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% @rplparallel/PLOT takes 'LabelsOff' as an example
if(~Args.LabelsOff)
	xlabel('Marker Number')
	ylabel('Marker Value')
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