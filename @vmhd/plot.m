function [obj, varargout] = plot(obj,varargin)
%@dirfiles/plot Plot function for dirfiles object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','','plotmap',1);
Args.flags = {'LabelsOff','ArgsOnly','plotmap'};
[Args,varargin2] = getOptArgs(varargin,Args);

if Args.plotmap
    
    ax = gca;
    mapL = obj.data.maps_raw;
        % Convert to x,y %
    bin_ang = (2*pi) / length(mapL);
%     angles = (bin_ang/2 : bin_ang : (2*pi)-(bin_ang/2)); % CCW from x axis
    angles = (-bin_ang/2 : -bin_ang : (-2*pi)+(bin_ang/2)); % CW from x axis
    angles=angles+0.5*pi; % north is 0 deg
    angles(angles>(2*pi))=angles(angles>(2*pi)) - 2*pi;
    angles(angles<=0)=angles(angles<=0) + 2*pi;
    % Normalize map, then multiply by factor if need to be scaled to place/view maps
    mapLnorm = (mapL./(nanmax(abs(mapL)))) .* 20;
    [x, y] = pol2cart(angles, mapLnorm);
    x(end+1) = x(1);    % Make line meet itself at end
    y(end+1) = y(1);    %
    % Color the plot line
    c = mapL;
    c(end+1) = c(1);
    maxrate = max(mapL,[],'omitnan');
    if isnan(maxrate)
        maxrate = 1;
    end

    if isempty(varargin) % Plot HD map with interp-colored line
        % Plot %
        hMap = patch(x,y,c,'FaceColor','interp','EdgeColor','interp','LineWidth',3,'FaceAlpha',0.5,'EdgeAlpha',0.5);
    elseif size(varargin,2) == 2 && strcmp(varargin{1},'headdirection') % Patch HD field within HD map
        % Plot %
%         ax = gca;
%         hMap = plot(ax,x, y, 'Color',varargin{2}); % Single arc
        x(isnan(x)) = 0;
        y(isnan(y)) = 0;
        hMap = patch(ax,x, y, varargin{2},'FaceAlpha',0.3,'EdgeColor','none'); % patch
    
        % Move center of plot
        zplus = repelem(raise,length(x));
        xplus = x+20;
        yplus = y+20;
        % Plot
%         ax = gca;
        if size(varargin,2) == 1 % Plot HD map with interp-colored line within place or view map
            hMap = patch(xplus,yplus,zplus,c,'FaceColor','interp','EdgeColor','interp','LineWidth',3,'FaceAlpha',0.5,'EdgeAlpha',0.5);
        elseif size(varargin,2) == 2 % Patching HD field within a place or view map
%             hMap = plot3(ax,xplus, yplus,zplus, varargin{2});
            xplus(isnan(xplus)) = 20;
            yplus(isnan(yplus)) = 20;
            hMap = patch(ax,xplus, yplus, zplus,varargin{2},'FaceAlpha',0.3,'EdgeColor','none'); % patch
        end 
    end
       % This is all formatting of axes %
    if isempty(varargin) 
        axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim')])); % Make axes square around biggest value
        if axlim == 0
            axlim = 1;
        end
        axis(ax,[-1.1*axlim 1.1*axlim -1.1*axlim 1.1*axlim]);                         %   ..
        line('xdata',0.95*[-axlim axlim],'ydata',[0 0],'parent',ax);   
        line('xdata',[0 0],'ydata',0.95*[-axlim axlim],'parent',ax); % centre-crossing axes
        % axis(ax, 'square', 'off');
%         axis(ax, 'square', 'off', 'tight');
        ax.FontSize = 14;
        set(ax,'CLim',[0 maxrate],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
    end
    
    view(-35,20);
    colormap jet;
    % colorbar;

    mapG = [x' y'];
    mapGdummy = [[1:size(x,2)-1 1]' [1:size(y,2)-1 1]'];
    
    axis(ax,'tight');
    w = ax.OuterPosition(3)*0.02;
    x = ax.OuterPosition(1) + ax.OuterPosition(3);
    y = ax.OuterPosition(2) + ax.OuterPosition(4)/4;
    h = ax.OuterPosition(4)/2;

    c = colorbar(ax,'Position',[x y w h]);
end

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	if(Args.Type1)
		% code to plot 1 kind of plot
		
		% label the axis
		xlabel('X Axis')
		ylabel('Y Axis')
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

	% add an appropriate title
	sdstr = get(obj,'SessionDirs');
	title(getDataOrder('ShortName','DirString',sdstr{1}))
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
