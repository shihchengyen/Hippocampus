function [obj, varargout] = plot(obj,varargin)
%@eyelink/plot Plot function for eyelink object.
%   eyelink = plot(eyelink_obj) creates a histogram of the saccades and
%   fixations per trial
%   OBJ = plot (eyelink_obj, 'Trial') Plots the x and y positions of the eye wrt time per trial
%   OBJ = plot (eyelink_obj, 'Trial', 'EyeMovements'); Plots all the x and y
%   positions of the eye per trial 
%

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Trial',0, 'XY', 0, 'Calibration', 0, 'Cmds', ''); %commands after Tiral 
Args.flags = {'LabelsOff','ArgsOnly','Trial','XY', 'Calibration'};
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

     %convert n into the column (session) and 
        %row (trial) that is being accessed.
    noOfSessions = obj.data.noOfSessions;
    i=1;
    if(noOfSessions > 1) 
       for i = 1:noOfSessions
           l = obj.data.noOfTrials(1,i);
       %find the actual length of the first column
           if ( (n-l) > 0)
               n = n-l;
           else 
               %disp(n);
               %fprintf('and');
               %disp(i);
               break;
           end 
       end 
     end %to identify column (i) and row (n) that the item is in  
    
    if(Args.Trial)
            
    %Plot x vs t and y vs t positions per trial
        i = 1+(i-1)*3;  %Use this for finding which column (session) is being accessed 
        x = obj.data.trial_timestamps(:, i:i+2); 

        trial_start_time = double(obj.data.timestamps(x(n,1))); %convert to double for math manip 
        trial_cue_time = double (obj.data.timestamps(x(n,2))) - trial_start_time;
        trial_end_time = (obj.data.timestamps(x(n,3)));

        %timestamps is the x axis to be plotted 
        timestamps = obj.data.timestamps(x(n,1)-500 : x(n,3)); %start indexing from 500 ms before start
        timestamps = double (timestamps);
        timestamps = timestamps - trial_start_time;

        y = obj.data.eye_pos(x(n,1)-500:x(n,3), :);

        %start plotting all the data: 
        plot (timestamps, y(:,1), 'b-');
        %plot the y data and label all the axes 
        hold on
        title (strcat(['Eye Movements versus Time for Trial' ' ' num2str(n)]));
        xlabel ('Time (ms)');
        ylabel ('Position (screen pixels)');
        plot (timestamps, y(:,2), 'g-');
        legend ('X position', 'Y Position');

        %Plotting lines to mark the start, cue offset, and end/timeout
        %for the trial
        line([0 0], ylim, 'Color','g');
        line([trial_cue_time trial_cue_time], ylim, 'Color', 'm');
        timedOut = obj.data.timeouts == trial_end_time; 
        trial_end_time = trial_end_time - trial_start_time;
        timedOut = find(timedOut,1);
        if (isempty(timedOut)) %this means the trial didn't time out
            line([trial_end_time trial_end_time], ylim, 'Color', 'b');
        else %trial did timeout
            line([trial_end_time trial_end_time], ylim, 'Color', 'r');
        end 

        hold off  
        
    elseif (Args.XY)
        %Plots the x and y movement of the eye per trial
        %extract all the trials from one session 
        i = 1+(i-1)*3;
        x = obj.data.trial_timestamps(:, i:i+2); 
        y = obj.data.eye_pos(x(n,1):x(n,3), :);
            
        %plot all the data points 
        plot (y(:,1), y(:,2), 'bo');
        hold on
        title (['Eye Movements on the screen for trial' ' ' num2str(n)]);
        xlabel ('Gaze position X (screen pixels)');
        ylabel ('Gaze position Y (screen pixels)');
        hold off
    
    elseif (Args.Calibration)
        %index into the eyeposition 
        eyePos = obj.data.eyePos;
        indices = uint32(obj.data.indices);
        y = eyePos(indices(n,1):indices(n,3), :);
        
        plot(y(:,1), y(:,2), 'bo');
        hold on
        title ('Calibration Eye movements from session');
        xlabel ('Gaze Position X (screen pixels)');
        ylabel ('Gaze Position Y (screen pixels)');
        hold off 
    
    else %some other argument - code to plot yet another kind of plot
        
        %Histograms of fixations and saccades per session
        sacc_durations = obj.data.sacc_event(:,n);
        fix_durations = obj.data.fix_event(:,n);
        
        %empty data points that are zero 
        sacc_durations(sacc_durations == 0) = [];
        fix_durations (fix_durations == 0) = [];
        
        lower = min(fix_durations);
        upper = max(fix_durations);
        edges = [lower:25:upper];
        edges(edges > 1000) = [];
        
        histogram (sacc_durations, edges);
        hold on;
        title ('Distribution of Saccades and Fixations for the Session');
        histogram (fix_durations, edges);
        legend ('Saccades', 'Fixations');
        ylabel ('# of events');
        xlabel ('Duration (ms)');
        hold off;

    end
    
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
