function [obj, varargout] = plot(obj,varargin)
%@rplview/plot Plot function for rplview object.
%   OBJ = plot(OBJ,N) plots the Nth entity in the data file.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'PlotEntity',0, 'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','PlotEntity'};
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
	% display entity info
	nsEI
	numSamples = nsEI.ItemCount;
	switch(nsEI.EntityType)
		case 1
			if(Args.PlotEntity)
				% Event entity
				ddata = NaN(1, numSamples); 
				timeStamps = NaN(1,numSamples);
				for i = 1:numSamples
					[~, timeStamps(i), ddata(i)] = ns_GetEventData(hFile, n, i);
				end 
				stem(timeStamps,ddata)
				if(~Args.LabelsOff)
					xlabel('Time (s)')
					ylabel('Marker')
					title(nsEI.EntityLabel)
				end
			end % if(Args.PlotEntity)
			% display event info
			[ns_status,eventInfo] = ns_GetEventInfo(hFile,n);
			eventInfo
		case 2
			if(Args.PlotEntity)
				% Analog entity
				[ns_RESULT,~,analogData] = ns_GetAnalogData(hFile,n,1,numSamples);
				plot(analogData)
				if(~Args.LabelsOff)
					xlabel('Time (s)')
					ylabel('Voltage (uV)')
					title(nsEI.EntityLabel)
				end
			end % if(Args.PlotEntity)
			% display analog info
			[ns_RESULT,analogInfo] = ns_GetAnalogInfo(hFile,n);
			analogInfo
		case 3		
			if(Args.PlotEntity)
				% Segment entity
				% read first segment to get size
				[ns_RESULT, ~, ~, sample_count] = ns_GetSegmentData(hFile, n, 1);
				sTime = zeros(1,numSamples);
				sData = zeros(sample_count,numSamples);
				sCount = zeros(1,numSamples);
				unit_id = zeros(1,numSamples);
				for i = 1:numSamples
					[ns_RESULT, sTime(i), sData(:,i), sCount(i), unit_id(i)] ...
						= ns_GetSegmentData(hFile, n, i);
				end
				% Sorted spikes
				uniqueUnits = unique(unit_id(unit_id~=0));
				for i = 1:length(uniqueUnits)
					SortedSpikeData{i} = sData(:,unit_id==uniqueUnits(i));
				end
		
				if isempty(uniqueUnits)
					hold on;
					for i = 1:numSamples 
						plot(sData(:,i),'b'); 
					end
					hold off
				else
					ul = length(uniqueUnits);
					for i = 1:ul
						if(ul>1)
							subplot(ul,1,i)
						end
						plot(SortedSpikeData{i})
					end
				end            
				if(~Args.LabelsOff)
					xlabel('Sample'); 
					ylabel('Spike Data (\muV)');
					title(nsEI.EntityLabel)
				end
			end % if(Args.PlotEntity)
			% display segment info
			[ns_status,segmentInfo] = ns_GetSegmentInfo(hFile,n);
			segmentInfo
			[ns_status,segmentSourceInfo] = ns_GetSegmentSourceInfo(hFile,n,segmentInfo.SourceCount);
			segmentSourceInfo
		case 4
			% Neural event entity
		otherwise
	end
else
	% plot all data
end

%%
RR = eval('Args.ReturnVars');
for i=1:length(RR) RR1{i}=eval(RR{i}); end 
varargout = getReturnVal(Args.ReturnVars, RR1);
