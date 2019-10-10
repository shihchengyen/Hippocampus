function [obj, varargout] = plot(obj,varargin)
%@vmspiketrain/plot Plot function for vmspiketrain object.
%   OBJ = plot(OBJ) creates a form of raster plot for locations in the virtual maze with 
%   spikes. 
%
%   example InspectGUI(vs,'GridSteps',5) 
%
%   Dependencies: unitymaze.m

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', 'DataFile','spiketrain.mat', ...
		  'Linear',0);
Args.flags = {'LabelsOff','ArgsOnly','Linear'};
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

	clist = get(gca,'ColorOrder');
	clistsize = size(clist,1);

	spike_xy = obj.data.spike_xy{n};
	sorted_sgpi = obj.data.sorted_sgpi{n};
	sorted_sgpi_info = obj.data.sorted_sgpi_info{n};

    % for-loop over grid positions
    for i=1:size(sorted_sgpi_info,1)
        % create subplot
        % index into spike_xy to get xy position
        indices = sorted_sgpi_info(i,2):sorted_sgpi_info(i,3);
        % get xy positions for this grid position
        sxy = spike_xy(sorted_sgpi(indices),:);
        if(Args.Linear)
            reps = sorted_sgpi_info(i,5);
            repi = obj.data.rep_num(sorted_sgpi(indices));
            
        else
        	plot(sxy(:,1),sxy(:,2),'.','Color',clist(mod(i-1,clistsize)+1,:))
			if(i==1)
				hold on
			end
		end
	end
    
    horGridBound = obj.data.horGridBound;
    vertGridBound = obj.data.vertGridBound;
    
    % plot grid
    % line(repmat(um.data.horGridBound,2,1),ylim,'k:')
    % line(repmat(um.data.vertGridBound,2,1),'k:')
    line( [repmat(horGridBound,2,1) repmat(horGridBound([1 end])',1,6)], ...
        [repmat(vertGridBound([1 end])',1,6) repmat(vertGridBound,2,1)], ...
        'Color','k','LineStyle',':')

    hold off
	
	% label the axis
	xlabel('X Dimension')
	ylabel('Y Dimension')

	% add an appropriate title
	sdstr = get(obj,'SessionDirs');
	title(getDataOrder('ShortName','DirString',sdstr{n}))
else
	% plot all data
end

% The following code allows any commands to be executed as part of each plot
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
