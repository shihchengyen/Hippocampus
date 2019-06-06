function [obj, varargout] = plot(obj,varargin)
%@eyelink/plot Plot function for eyelink object.
%   eyelink = plot(eyelink_obj) creates a histogram of the saccades and
%   fixations per trial
%   OBJ = plot (eyelink_obj, 'Trial') Plots the x and y positions of the eye wrt time per trial
%   OBJ = plot (eyelink_obj, 'Trial', 'EyeMovements'); Plots all the x and y
%   positions of the eye per trial 
%

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Trial',0, 'EyeMovements', 0, 'Cmds', ''); %commands after Tiral 
Args.flags = {'LabelsOff','ArgsOnly','Trial','EyeMovements'};
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
    
    if(Args.Trial)
       
        %the following loop helps break down n into the column(session) and
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
                    disp(n);
                    fprintf('and');
                    disp(i);
                    break;
                end 
            end 
        end %to identify column (i) and row (n) that the item is in  
            
        if (Args.EyeMovements) %Plots the x and y movement of the eye per trial
            %extract all the trials from one session 
            y = obj.data.eye_pos(x(n,1):x(n,3), :);
            
            %plot all the data points 
            plot (y(:,1), y(:,2), 'bo');
            hold on
            title ('Eye Movements on the screen for Trial');
            xlabel ('x pos');
            ylabel ('y pos');
            hold off
            
        else %Plot x and y positions wrt time per trial
            i = 1+(i-1)*3;  %Use this for finding which column (session) is being accessed 
            x = obj.data.trial_timestamps(:, i:i+2); 
            %x = obj.data.trial_idx(:, i:i+2);
            
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
            title ('Eye Movements versus Time for Trial');
            xlabel ('Time');
            ylabel ('Position');
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
        end 
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
        title ('Distribution of saccades and fixations for per Session');
        histogram (fix_durations, edges);
        legend ('Saccades', 'Fixations');
        ylabel ('# of events');
        xlabel ('Duration (ms)');
        hold off;

        
	end	

else
    
    
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
