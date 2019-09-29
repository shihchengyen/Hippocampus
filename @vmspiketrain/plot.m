function [obj, varargout] = plot(obj,varargin)
%@vmspiketrain/plot Plot function for vmspiketrain object.
%   OBJ = plot(OBJ) creates a form of raster plot for locations in the virtual maze with 
%   spikes. 
%
%   example InspectGUI(vs,'GridSteps',5) 
%
%   Dependencies: unitymaze.m, 

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', 'DataFile','spiketrain.mat');
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

	sdstr = get(obj,'SessionDirs');
	cwd = pwd;
	cd(sdstr{n})

	um = unitymaze('auto',varargin2{:});
	% load file
	st = load(Args.DataFile);

	cd(cwd)

	% these are object specific fields
	ts = (st.timestamps/1000)';
	% perform histogram on spike times using unity time bins in order to get the closest xy
	% position
	[~,~,ubin] = histcounts(ts,um.data.unityTime);
	% since ts may contain spike times after the unity program has ended, we need to look only
	% at the non-zero values in ubin
	spike_xy = um.data.unityData(ubin(ubin>0),3:4);

	% perform histogram on spike times to see which grid position the above information should
	% go to
	% [n,tsbin] = histc(ts,sTime(:,1));
	[~,~,~,binH,binV] = histcounts2(spike_xy(:,1),spike_xy(:,2),um.data.horGridBound,um.data.vertGridBound);
	% compute grid position number
	spike_gridposition = binH + ((binV - 1) * um.data.gridSteps);

	% get the grid position for each spike
	% spike_gridposition = sTime(tsbin(tsbin>0),2);

	% sort according to grid position
	[sorted_spike_gridposition,sorted_sgpi] = sort(spike_gridposition);

	% find the points where grid position changes
	% first set of values in sorted_spike_gridposition is for 0, which is where there is no 
	% grid position, e.g. during cue presentation, etc., so the first change is where we want
	% to start. 
	change_ssgpindex = find(diff(sorted_spike_gridposition));
	change_ssgpindex(end+1) = size(spike_xy,1);

	clist = get(gca,'ColorOrder');
	clistsize = size(clist,1);

	% for-loop over grid positions
	for i=1:(size(change_ssgpindex,1)-1)
		% create subplot
		% index into spike_xy to get xy position
		indices = (change_ssgpindex(i)+1):change_ssgpindex(i+1);
		% get xy positions for this grid position
		sxy = spike_xy(sorted_sgpi(indices),:);
		plot(sxy(:,1),sxy(:,2),'.','Color',clist(mod(i-1,clistsize)+1,:))
		if(i==1)
			hold on
		end
    end
    
    % plot grid
    % line(repmat(um.data.horGridBound,2,1),ylim,'k:')
    % line(repmat(um.data.vertGridBound,2,1),'k:')
    line( [repmat(um.data.horGridBound,2,1) repmat(um.data.horGridBound([1 end])',1,6)], ...
        [repmat(um.data.vertGridBound([1 end])',1,6) repmat(um.data.vertGridBound,2,1)], ...
        'Color','k','LineStyle',':')

    hold off
	
	% label the axis
	xlabel('X Dimension')
	ylabel('Y Dimension')

	% add an appropriate title
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
