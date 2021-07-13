function [obj, varargout] = plot(obj,varargin)
    %@dirfiles/plot Plot function for dirfiles object.
    %   OBJ = plot(OBJ) creates a raster plot of the neuronal
    %   response.
    Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',10, ...
              'ReturnVars',{''},'Magnify',0,'Trial',0,'ArgsOnly',0, 'Cmds','');
    % flag: different ways to plot data--shortcut--true/false
    Args.flags = {'LabelsOff','ArgsOnly','Magnify'};
    [Args,varargin2] = getOptArgs(varargin,Args);
    % if user select 'ArgsOnly', return only Args structure for an empty object
    if Args.ArgsOnly
        Args = rmfield (Args, 'ArgsOnly');
        varargout{1} = {'Args',Args};
        return;
    end
     

    if(~isempty(Args.NumericArguments))
        eots = obj.data.eots;
        heads = obj.data.heads;
        %% Plotting path of each trial on map
        res = 0.025;

        % Read and display map image
        mapimage = imread(obj.data.mapimage);
        map = binaryOccupancyMap(mapimage, 1/res);
        map.GridOriginInLocal = [-13.80,-21.80];

        % show(map);
        % for n = 1:(length(eots)-1)
        % graph each map one by one
        
%       Conditional statement that selects trial number. Implemented
%       according to InspectFn.m: When there is only 1 numerical input, the
%       argument is Trial number; When there are two numerical inputs, the
%       first one is "Event Number", the second one is trial number.
        if(~mod(length(varargin),2) && isnumeric(varargin{2}))
            n = Args.NumericArguments{2};
        else
            n = Args.NumericArguments{1};
        end
        
        % figure(); %Only turn on when debugging plot.m
        hold on
        show(map);
        % The reward region of each cued poster is marked as yellow circle
        [cat_x, cat_y] = circle(-1.3727587461471558, 0.85040283203125, 1);
        [camel_x, camel_y] = circle(-1.6286569833755493, -2.674863815307617, 1);
        [croc_x, croc_y] = circle(0.5865049362182617, -4.530182838439941, 1);
        [rabbit_x, rabbit_y] = circle(-1.717581868171692, -5.609162330627441, 1);
        [pig_x, pig_y] = circle(-0.980493426322937, -8.913317680358887, 1);
        [donkey_x, donkey_y] = circle(-2.202519416809082, -11.747871398925781, 1);
        
        plot(cat_x, cat_y,'Color','yellow'); 
        plot(camel_x, camel_y,'Color','yellow');   
        plot(rabbit_x, rabbit_y,'Color','yellow');
        plot(donkey_x, donkey_y,'Color','yellow');
        plot(croc_x, croc_y,'Color','yellow');  
        plot(pig_x, pig_y,'Color','yellow'); 
        % case poster_ID --> highlight the reward region
        
        % Todo: mark the correct and incorrect trial stop
        cue = obj.data.animal_cues(n);
        switch cue
            case 1
                plot(cat_x, cat_y,'Color','red');   
            case 2
                plot(camel_x, camel_y,'Color','red');     
            case 3
                plot(rabbit_x, rabbit_y,'Color','red');
            case 4
                plot(donkey_x, donkey_y,'Color','red');
            case 5
                plot(croc_x, croc_y,'Color','red');        
            case 6
                plot(pig_x, pig_y,'Color','red'); 
            otherwise
                display('No matching posters!\n');
        end       
        
        % Show posters' locations
        text(-2.1,1,'Cat','Color','red','FontSize',10);
        text(-3.0,-2.7,'Camel','Color','red','FontSize',10);
        text(-3.2,-5.5,'Rabbit','Color','red','FontSize',10);
        text(0.4,-4.5,'Croc','Color','red','FontSize',10);
        text(-1.0,-8.7,'Pig','Color','red','FontSize',10);
        text(-2.8,-12.0,'Donkey','Color','red','FontSize',10);
        
        % Define the extra points at the end of each trial
        if n==length(heads)&& n<length(eots)
            extra = obj.data.poseStructs((eots(n)+1):eots(n+1));
        elseif n==length(heads) || eots(n)==heads(n+1)
            extra = {};
        else
            extra = obj.data.poseStructs((eots(n)+1):heads(n+1));
        end                
        
        % Define each trial
        rosT = obj.data.poseStructs(heads(n):eots(n));
        xPoints = cellfun(@(m) double(m.Pose.Pose.Position.X),rosT);
        yPoints = cellfun(@(m) double(m.Pose.Pose.Position.Y),rosT);
        % scatter(xPoints(2:end-1),yPoints(2:end-1),'o','filled');
        
         % Find Orientation
        orienT = obj.data.poseStructs(heads(n):eots(n));
        z = cellfun(@(m) double(m.Pose.Pose.Orientation.Z),orienT);

        yaw = 2*asin(z); % https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Intuition
        nx = cos(yaw);
        ny = sin(yaw);
        %% Begin plotting quiver plot
        %figure();
        %hold on;
        quiver_scale = 2;
        % norm = sqrt(sum(nx.^2+ny.^2));
        nx = nx.*quiver_scale;
        ny = ny.*quiver_scale;
        q1 = quiver(xPoints(1:end-1),yPoints(1:end-1),nx(1:end-1),ny(1:end-1),0.5,'linewidth',1);
        q2 = quiver(xPoints(end),yPoints(end),nx(end),ny(end),0.2,'linewidth',1,'color','magenta');
        
%         fprintf("5th [x y u v] = [%.2f,%.2f,%.2f,%.2f]\n",xPoints(5),yPoints(5),nx(5),ny(5));
%         fprintf("Last [x y u v] = [%.2f,%.2f,%.2f,%.2f]\n",xPoints(end),yPoints(end),nx(end),ny(end));       
        
        %% Extra poses and head & tail of each trial are marked
        % Draw the extra poses together with each trial; The start of each
        % trial is marked as green cross, and the end is marked as magenta
        % cross. The last extra poses from the previous trial is the start
        % of the next trial.
        if ~isequal(extra, {})
            extra_x = cellfun(@(m) double(m.Pose.Pose.Position.X),extra);
            extra_y = cellfun(@(m) double(m.Pose.Pose.Position.Y),extra);
            scatter(extra_x,extra_y,'o','white');
%             display('end of extra: ')
%             display([n extra_x(end) extra_y(end)]);
        end
        
        scatter(xPoints(1),yPoints(1),'x','green','LineWidth',1.5);
        scatter(xPoints(end),yPoints(end),'x','magenta','LineWidth',1.5);
%         Only comment out for debugging of the quiver size
%         fprintf('head: ')
%         display([n xPoints(1) yPoints(1)]);
        
        hold off;
        
        %% Magnify option
        % When 'Magnify' flag is turned on, each trial's map is automatically
        % zoomed in to the trial. The boundary of the zoom-in version is
        % defined by the global max and min values of the poses. The
        % boundary buffer (currently 50%) can be readjusted based on
        % different maps. The panel displayed in GUI is currently a square.
        if Args.Magnify 
            % to find the global max and min
            all_T = obj.data.poseStructs(1:end);
            all_x = cellfun(@(m) double(m.Pose.Pose.Position.X),all_T);
            all_y = cellfun(@(m) double(m.Pose.Pose.Position.Y),all_T);
            
            % global_miny and global_maxy are picked as the square panel's boundary
            % because y(s) has a larger range than the x(s)
            % global_minx = min(all_x,[],'all'); 
            % global_maxx = max(all_x,[],'all');
            global_miny = min(all_y,[],'all');
            global_maxy = max(all_y,[],'all');
            
            xlim([global_miny*1.5+7 global_maxy*1.5+7])
            ylim([global_miny*1.5+2 global_maxy*1.5+2])
            
            % To fix the axis of the panels of all trials
            set(gca,'XTickLabel',[-4 -2 0 2 4 6]);
            set(gca,'YTickLabel',[-5 -4 -3 -2 -1 0 1 2 3 4 5]);

        end
        hold off;
        %end
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
end

function [xunit yunit] = circle(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = transpose(r * cos(th) + x);
    yunit = transpose(r * sin(th) + y);
end