function [obj, varargout] = plot(obj,varargin)
%@viewsort/plot Plot function for viewsort object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		'Array',0, 'Session',0, 'Day',0,  ...
		'Cmds', '', 'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','Array','Session','Day'};
[Args,varargin2] = getOptArgs(varargin,Args,'removenumargs');

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

delete(findall(gcf,'type','annotation'))

% plot waveform
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
			if(obj.data.Args.Saved==0)
                plot((obj.data.spikeForms(xind,:))','.-')
            else
                % plot saved spikeforms in bold
                plot((obj.data.spikeForms(xind,:))','.-', 'LineWidth', 2)
            end
			sdstr = get(obj,'SessionDirs');
			title(getDataOrder('ShortName','DirString',sdstr{chnindex}))
            
            % plot noise on top of waveform
            hold on
            line(repmat(xlim',1,2),repmat([-obj.data.Noise(index) obj.data.Noise(index)],2,1),'Color','r')
            % plot(obj.data.Noise(xind(end),:),'r')
            % plot(0-(obj.data.Noise(xind(end),:)),'r')
            hold off
            
            % add axis labels
            if(~Args.LabelsOff)
				xlabel('Data Points')
				ylabel('Voltage (uV)')
            end
		end  % for index = 1:numSets
	else  % if(Args.Array || Args.Session || Args.Day)
        xind = (obj.data.ChannelIndex(n)+1):obj.data.ChannelIndex(n+1);
        if(obj.data.Args.Saved==0)
            plot((obj.data.spikeForms(xind,:))','.-')
        else
            % plot saved spikeforms in bold
            plot((obj.data.spikeForms(xind,:))','.-', 'LineWidth', 2)
        end
        sdstr = get(obj,'SessionDirs');
        title(getDataOrder('ShortName','DirString',sdstr{n}))
        
        % plot noise on top of waveform
        hold on
        line(repmat(xlim',1,2),repmat([-obj.data.Noise(n) obj.data.Noise(n)],2,1),'Color','r')
        % plot(obj.data.Noise(xind(end),:),'r')
        % plot(0-(obj.data.Noise(xind(end),:)),'r')
        hold off
        
        % add axis labels
        if(~Args.LabelsOff)
            xlabel('Data Points')
            ylabel('Voltage (uV)')
        end
        
        if(obj.data.Args.Saved==0)
            % show ISI coefficients of variation and mean firing rate
            legendLabels = cell(1,size(xind,2));
            for k = 1:size(xind,2)
                if isnan(obj.data.coeffV_ISI(xind(k)))
                    legendLabels{k} = 'N/A';
                else
                    cvlabel = round(obj.data.coeffV_ISI(xind(k)),2);
                    frlabel = 1000./obj.data.meanISI(xind(k));
                    legendLabels{k} = sprintf('%u: %.2f %.2f', k-1, cvlabel,frlabel);
                end
            end
            lgd = legend(legendLabels, 'FontSize', 12);
            title(lgd,{'ISI Distribution','Coefficient of Variation:','Firing Rate:'})
            % subplot(round((numUnits+1)/2),2,d+1);
            % histogram(spike_ISI,1000,'facecolor',plotColour(d),'edgecolor',plotColour(d)); hold on
            % title(strcat('(',num2str(d),')',' ISI Distribution, Coefficient of Variation: ', num2str(round(coeffV_ISI,2))));
            % xlabel('ISI(ms)','FontSize',10); ylabel('Count','FontSize',10);
            
            % show spike similarities (dot product)
            spikeind = (obj.data.spikesimIndex(n)+1):obj.data.spikesimIndex(n+1);
            if (~isnan(obj.data.spikesim(spikeind)))
                perms = nchoosek(1:size(xind,2),2);
                % add 1 to the correlation coefficients as we are going to subtract 1 from
                % all the values in perms to switch from 1 index to 0 index
                perms(:,3:4) = [obj.data.spikesim(spikeind) obj.data.spikesp2pdiffs(spikeind)] + 1;
                dim = [0.2, 0.2, 0.1, 0.1];
                annotation('textbox',dim,'String', num2str(perms-1));
            end
        end %if(obj.data.Args.Saved==0)
    end  % if(Args.Array || Args.Session || Args.Day)
else
	% plot all data
	% numSets = obj.data.numSets;
    numSets = get(obj,'number');
    for index = 1:numSets
        if(numSets>1)
            nptSubplot(numSets,index);
        end
        xind = (obj.data.ChannelIndex(index)+1):obj.data.ChannelIndex(index+1);
        plot((obj.data.spikeForms(xind,:))','.-')
        sdstr = get(obj,'SessionDirs');
        title(getDataOrder('ShortName','DirString',sdstr{index}))
        % plot noise on top of waveform
        hold on
        line(repmat(xlim',1,2),repmat([-obj.data.Noise(index) obj.data.Noise(index)],2,1),'Color','r')
        % plot(obj.data.Noise(xind(end),:),'r')
        % plot(0-(obj.data.Noise(xind(end),:)),'r')
        hold off
        % add axis labels
        if(~Args.LabelsOff)
            xlabel('Data Points')
            ylabel('Voltage (uV)')
        end
    end
end

% evaluate 'Cmds'
if(~isempty(Args.Cmds))
    % save the current figure in case Args.Cmds switches to another figure
    h = gcf;
    % save current directory
    cwd = pwd;
    % change to corresponding session directory
    cd(sdstr{n})
    % run command
    eval(Args.Cmds)
    % switch back to previous figure
    figure(h);
    % switch back to previous directory
    cd(cwd);
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
