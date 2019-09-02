function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0,'Trial',0,'Cmds', '');
Args.flags = {'LabelsOff','ArgsOnly', 'Trial'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

%draw the structure of the um layout 
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

%P1(-5,1.5,-7.55) P2(-7.55,1.5,5) P3(7.55,1.5,-5) P4(5,1.5,7.55) P5(5,1.5,2.45) P6(-5,1.5,-2.45)
catx = [-6.12; -3.88];
caty = [-7.55; -7.55];
camelx = [-7.55; -7.55];
camely = [3.88;6.12 ];
rabbitx = [7.55; 7.55];
rabbity = [-6.12; -3.88];
donkeyx = [3.88; 6.12];
donkeyy = [7.55; 7.55];
crocx = [3.88; 6.12];
crocy = [2.45; 2.45];
pigx = [-6.12; -3.88];
pigy = [-2.45; -2.45];

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	if(Args.Trial)
		% code to plot 1 kind of plot
		plot(xBound,zBound,'k','LineWidth',1.5);
		hold on
		plot(x1Bound,z1Bound,'k','LineWidth',1);
		plot(x2Bound,z2Bound,'k','LineWidth',1);
		plot(x3Bound,z3Bound,'k','LineWidth',1);
		plot(x4Bound,z4Bound,'k','LineWidth',1);
        
        %Plot the poster positions
        plot (catx, caty, 'r', 'LineWidth', 1.5);
        plot (camelx, camely, 'r', 'LineWidth', 1.5);
        plot (rabbitx, rabbity, 'r', 'LineWidth', 1.5);
        plot (donkeyx, donkeyy, 'r', 'LineWidth', 1.5);
        plot (crocx, crocy, 'r', 'LineWidth', 1.5);
        plot (pigx, pigy, 'r', 'LineWidth', 1.5);
        
        %Extract
        disp(n); 
        indices = obj.data.index(n,:); %stores indices with the start - cue - end event for nth trial
        %indices = indices(1,1)+1; %add the one because indices has indices to start-cue-end event message
        fixIndices = obj.data.fixIndex;
        
        %Plotting the locations looked at between cue offset and end Trial 
        cueEnd = fixIndices(fixIndices>=indices(1,2) & fixIndices<=indices(1,3));
        monkeyLoc = horzcat(obj.data.playerLocation(cueEnd,1), obj.data.playerLocation(cueEnd,3)); %stores x and z position of monkey
        gazeEnd = horzcat(obj.data.playerGazeLocation(cueEnd,1), obj.data.playerGazeLocation(cueEnd,3)); %stores x and z position of gaze
        objLookedAt = obj.data.fixatedObj(cueEnd);
        x = [monkeyLoc(:,1)'; gazeEnd(:,1)'];
        y = [monkeyLoc(:,2)'; gazeEnd(:,2)'];

        sz = size(monkeyLoc(:,1),1);
        g = linspace(0.7, 0, sz);
        g=g';
        scatter(monkeyLoc(:,1), monkeyLoc(:,2), 100, [g g g], 'filled'); 
        %Plot the gaze lines of the monkey, according to what it is looking
        %at 
        for i=1:sz
            if(contains(objLookedAt(i,1), {'Ground'})==1)
                plot (x(:,i), y(:,i), '-*','Color', [0 1.0 1.0], 'LineWidth',1);
            elseif(contains(objLookedAt(i,1), {'Ceiling'})==1)
                plot (x(:,i), y(:,i), '-*','Color', [1.0 0.6 0.2], 'LineWidth',1);
            elseif(contains(objLookedAt(i,1), {'wall'})==1) 
                plot (x(:,i), y(:,i), '-*','Color', [0 0.5 0], 'LineWidth',1);
            elseif (contains(objLookedAt(i,1), {'Poster'})==1) 
                plot (x(:,i), y(:,i),'-*', 'Color', [1.0 0 1.0], 'LineWidth',1);
            end
        end
       
        %Plot the Last Location of the monkey as well as the last location
        %that it looked at 
        message = obj.data.fixatedObj(indices(1,3));
        lpx = [obj.data.playerLocation(indices(1,3)-1,1) ; obj.data.playerGazeLocation(indices(1,3)-1,1)];
        lpy = [obj.data.playerLocation(indices(1,3)-1,3) ; obj.data.playerGazeLocation(indices(1,3)-1,3)];
        
        if(contains(message,{'End Trial'})==1)
           scatter(obj.data.playerLocation(indices(1,3)-1,1), obj.data.playerLocation(indices(1,3)-1,3), 120, [0 1.0 0], 'filled'); %monkey loc
           plot (lpx, lpy, 'b*-', 'LineWidth',1);
        else
           scatter(obj.data.playerLocation(index(1,1)-1,1), obj.data.playerLocation(index(1,3)-1,3), 120, [1.0 0 0], 'filled'); %monkey loc
           plot (lpx, lpy, 'r*-', 'LineWidth',1);
        end
        
        % label the axis
        title(['Trial: ',num2str(n)]);
		xlabel('X Axis')
		ylabel('Y Axis')
        hold off
	elseif(Args.Type2)
		% code to plot another kind of plot
		
		% label the axis
		xlabel('X Axis')
		ylabel('Y Axis')
	else
		% code to plot yet another kind of plot
		
		% label the axis
		xlabel('X Axis')
		ylabel('Y Axis')
    end	
    
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
