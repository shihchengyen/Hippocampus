function [obj, varargout] = plot(obj,varargin)
%@viewsort/plot Plot function for viewsort object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		'Array',0, 'Session',0, 'Day',0,  ...
		'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','Array','Session','Day'};
[Args,varargin2] = getOptArgs(varargin,Args);

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
	if(Args.Array || Args.Session || Args.Day)
		if(Args.Array)
            % plot waveforms array by array
            % get starting row in ChannelIndex
            chnstart = obj.data.ArrayIndex(n);
            % get ending row in ChannelIndex
            chnend = obj.data.ArrayIndex(n+1);
%             if(Args.UseGMR)
%                 UseGMRlayout = 1;
%             end
		elseif(Args.Session)
			% plot waveforms session by session
            % get starting row in ChannelIndex
            chnstart = obj.data.SessionIndex(n);
            % get ending row in ChannelIndex
            chnend = obj.data.SessionIndex(n+1);
%             if(Args.UseGMR)
%                 UseGMRlayout = 1;
%             end
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
            
            % plot noise on top of waveform
            hold on
            line(repmat(xlim',1,2),repmat([-obj.data.Noise(index) obj.data.Noise(index)],2,1),'Color','r')
%             line(repmat(xlim',1,2),repmat([-obj.data.Noise],2,1));
%            plot(obj.data.Noise(xind(end),:),'r')
%            plot(0-(obj.data.Noise(xind(end),:)),'r')
            hold off            
		end  % for index = 1:numSets
	else  % if(Args.Array || Args.Session || Args.Day)
		xind = (obj.data.ChannelIndex(n)+1):obj.data.ChannelIndex(n+1);
		plot((obj.data.spikeForms(xind,:))','.-')
		sdstr = get(obj,'SessionDirs');
		title(getDataOrder('ShortName','DirString',sdstr{n}))
        % plot noise on top of waveform
        hold on
        line(repmat(xlim',1,2),repmat([-obj.data.Noise(n) obj.data.Noise(n)],2,1),'Color','r')
        hold off
	end  % if(Args.Array || Args.Session || Args.Day)
else
	% plot all data
	numSets = obj.data.numSets;
    for index = 1:numSets
        if(numSets>1)
            nptSubplot(numSets,index);
        end
		xind = (obj.data.ChannelIndex(index)+1):obj.data.ChannelIndex(index+1);
		plot((obj.data.spikeForms(xind,:))','.-')
		sdstr = get(obj,'SessionDirs');
		title(getDataOrder('ShortName','DirString',sdstr{index}))
        
        % add code for plot options here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % @viewsort/PLOT takes 'LabelsOff' as an example
        if(~Args.LabelsOff)
            xlabel('Data Points')
            ylabel('Voltage (uV)')
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
end

% % add code for plot options here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % @viewsort/PLOT takes 'LabelsOff' as an example
% if(~Args.LabelsOff)
% 	xlabel('Data Points')
% 	ylabel('Voltage (uV)')
% end
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
