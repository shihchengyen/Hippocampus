function [obj, varargout] = plot(obj,varargin)
%@unitymaze/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Trial',0, 'FrameIntervals',0, ...
		  'FrameIntervalTriggers',[2 3], 'DurationDiff',0);
Args.flags = {'LabelsOff','ArgsOnly','Trial','FrameIntervals','DurationDiff'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

% Vertices coordinates: 
vertices = [-10 10; -5 10; 0 10; 5 10; 10 10; -10 5; 0 5; 10 5; -10 0; -5 0; 0 0; 5 0; 10 0; -10 -5; 0 -5; 10 -5; -10 -10; -5 -10; 0 -10; 5 -10; 10 -10;];

% Poster coordinates
posterpos = [-5 -7.55; -7.55 5; 7.55 -5; 5 7.55; 5 2.45; -5 -2.45]; % Posters 1 to 6 respectively

% Plot boundaries
xBound = [-12.5,12.5,12.5,-12.5,-12.5]; % (clockwise from top-left corner) outer maze wall
zBound = [12.5,12.5,-12.5,-12.5,12.5];
x1Bound =[-7.5,-2.5,-2.5,-7.5,-7.5]; % yellow pillar 
z1Bound =[7.5,7.5,2.5,2.5,7.5];
x2Bound =[2.5,7.5,7.5,2.5,2.5]; % red pillar
z2Bound =[7.5,7.5,2.5,2.5,7.5];
x3Bound =[-7.5,-2.5,-2.5,-7.5,-7.5]; % blue pillar
z3Bound =[-2.5,-2.5,-7.5,-7.5,-2.5];
x4Bound =[2.5,7.5,7.5,2.5,2.5]; % green pillar
z4Bound =[-2.5,-2.5,-7.5,-7.5,-2.5];

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	if(Args.Trial)
		if(Args.FrameIntervals)
			indices = obj.data.unityTriggers(n,Args.FrameIntervalTriggers);
            uData = obj.data.unityData((indices(1)+1):indices(2),2);
			stem(uData)
			xlabel('Frames')
			ylabel('Interval (s)')
			% use n and setIndex to find which session this trial belongs to
			sessionnum = find(obj.data.setIndex<=n);
			sdstr = get(obj,'SessionDirs');
			sessionstr = getDataOrder('ShortName','DirString',sdstr{sessionnum(end)});
            % in order to display the duration difference, we will need to
            % load the rplparallel object to get the Ripple timestamps
            % first check if the rplparallel object has been loaded before
            % and stored in UserData
            ud = get(gcf,'UserData');
            if(isfield(ud,'rp'))
            	rp = ud.rp;
            else
	            rp = rplparallel('auto',varargin2{:});
	            ud.rp = rp;
	            set(gcf,'UserData',ud);
	        end
            rpTrialDur = diff(rp.data.timeStamps(n,2:3),1,2);
            uet = cumsum(uData);
			title([sessionstr ' Trial ' num2str(n) ' Duration disparity: ' num2str(1000*(uet(end)-rpTrialDur)) ' ms'])
		else  % if(Args.FrameIntervals)
			plot(xBound,zBound,'k','LineWidth',1.5);
			hold on
			plot(x1Bound,z1Bound,'k','LineWidth',1);
			plot(x2Bound,z2Bound,'k','LineWidth',1);
			plot(x3Bound,z3Bound,'k','LineWidth',1);
			plot(x4Bound,z4Bound,'k','LineWidth',1);
		
			% plot current trial trajectory
			plot(obj.data.unityData( obj.data.unityTriggers(n,2):(obj.data.unityTriggers(n,3)-1) ,3), ...
				obj.data.unityData( obj.data.unityTriggers(n,2):(obj.data.unityTriggers(n,3)-1) ,4),'b','LineWidth',1); 
			% plot end point identifier
			plot(obj.data.unityData(obj.data.unityTriggers(n,3),3), ...
				obj.data.unityData(obj.data.unityTriggers(n,3),4),'k.','MarkerSize',20); 
			hold off
		
			rstr = num2str(obj.data.sumCost(n,2));
			sstr = num2str(obj.data.sumCost(n,1));
			ratiostr = num2str(obj.data.sumCost(n,2)/obj.data.sumCost(n,1));
			% use n and setIndex to find which session this trial belongs to
			sessionnum = find(obj.data.setIndex<=n);
			sdstr = get(obj,'SessionDirs');
			sessionstr = getDataOrder('ShortName','DirString',sdstr{sessionnum(end)});
				
			title([sessionstr ' T: ' num2str(n-obj.data.setIndex(sessionnum(end))) ' Route: ' rstr ' Shortest: ' sstr ' Ratio: ' ratiostr])
		end  % if(Args.FrameIntervals)
	elseif(Args.DurationDiff)  % if(Args.Trial)
		uTrigs = obj.data.unityTriggers((obj.data.setIndex(n)+1):obj.data.setIndex(n+1),:);
		uTime = obj.data.unityTime(:,n);
		% add 1 to start index since we want the duration between triggers
		startind = uTrigs(:,1)+1;
		endind = uTrigs(:,3);
		starttime = uTime(startind);
		endtime = uTime(endind);
		trialDurations = diff([starttime endtime],1,2)
		% in order to display the duration difference, we will need to
		% load the rplparallel object to get the Ripple timestamps
		% first save the current directory
		cwd = pwd;
		% use n and setIndex to find which session this trial belongs to
		sessionnum = find(obj.data.setIndex<=n);
		sdstr = get(obj,'SessionDirs');
		sessionDir = sdstr{sessionnum(end)};
		% change directory to find the appropriate rplparallel object
		cd(sessionDir)
		rp = rplparallel('auto',varargin2{:});
		% change back to the previous directory
		cd(cwd)
		rpTrialDur = diff(rp.data.timeStamps(:,[1 3]),1,2);
		% multiply by 1000 to convert to ms
		histogram((trialDurations - rpTrialDur)*1000);
		xlabel('Time (ms)')
		ylabel('Frequency')
		sessionstr = getDataOrder('ShortName','DirString',sessionDir);
		title(['UnityFile trial duration - Ripple trial duration ' sessionstr])
		set(gca,'YScale','log','YGrid','on')
		cylim = ylim;
		ylim([0.5 cylim(2)])
	else  % if(Args.Trial)
		xind = (obj.data.setIndex(n)+1):obj.data.setIndex(n+1);
		[h,l1,l2] = plotyy(xind,obj.data.sumCost(xind,1:2),xind,obj.data.sumCost(xind,2)./obj.data.sumCost(xind,1),'bar','stem');
		h(2).YColor = 'm';
		yticks = h(2).YTick;
		h(2).YTick = 0:1:yticks(end);
		h(2).YGrid = 'on';
		l1(1).FaceColor = 'y';
		l1(2).FaceColor = 'c';
		l2.Color = 'm';
		
		legend('Shortest','Route','Ratio')
		
		sdstr = get(obj,'SessionDirs');
		title(getDataOrder('ShortName','DirString',sdstr{n}))
	end
else
	% plot all data
	n = 1;
end

% add code for plot options here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% read file but skip 1st row
% dfilename = obj.data.dlist(n).name;
% data = dlmread(dfilename,'',1,0);

% col 3 & 4 are x & y

% col 5 is angle in degrees
% but to display correctly, we need to multiply by -1 and add 90 degrees
% then we will convert to radians as the sin and cos functions expect radians
% dangle = deg2rad(-data(:,5) + 90);

% plot using quiver function
% quiver(data(:,3),data(:,4),cos(dangle),sin(dangle))
% title(dfilename)
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
