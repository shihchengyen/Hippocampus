function [obj, varargout] = plot(obj,varargin)
%@rplview/plot Plot function for rplview object.
%   OBJ = plot(OBJ,N) plots the Nth entity in the data file.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
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
	% plot one entity at a time
	n = Args.NumericArguments{1};
	[ns_status,hFile] = ns_OpenFile(obj.data.rawfname); 
	[ns_status,nsEI] = ns_GetEntityInfo(hFile, n);
	numSamples = nsEI.ItemCount;
	switch(nsEI.EntityType)
		case 1
			% Event entity
			ddata = NaN(1, numSamples); timeStamps = NaN(1,numSamples);
			for i = 1:numSamples
				[~, timeStamps(i), ddata(i)] = ns_GetEventData(hFile, n, i);
			end 
			stem(timeStamps,ddata)
			if(~Args.LabelsOff)
				xlabel('Time (s)')
				ylabel('Marker')
				title(nsEI.EntityLabel)
			end			
		case 2
			% Analog entity
			[ns_RESULT,~,analogData] = ns_GetAnalogData(hFile,n,1,numSamples);
			plot(analogData)
			% display analog info
			[ns_RESULT,analogInfo] = ns_GetAnalogInfo(hFile,n);
			analogInfo
			if(~Args.LabelsOff)
				xlabel('Time (s)')
				ylabel('Voltage (uV)')
				title(nsEI.EntityLabel)
			end
		case 3		
			% Segment entity
			ddata = NaN(1, numSamples); timeStamps = NaN(1,numSamples);
			for i = 1:numSamples
				[~, timeStamps(i), ddata(i)] = ns_GetEventData(hFile, ni, i);
			end 
			stem(timeStamps,ddata)
			if(~Args.LabelsOff)
				xlabel('Time (s)')
				ylabel('Marker')
				title(nsEI.EntityLabel)
			end			
		case 4
			% Neural event entity
		otherwise
	end
else
	% plot all data
end

% add code for plot options here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% @rplview/PLOT takes 'LabelsOff' as an example
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RR = eval('Args.ReturnVars');
for i=1:length(RR) RR1{i}=eval(RR{i}); end 
varargout = getReturnVal(Args.ReturnVars, RR1);
