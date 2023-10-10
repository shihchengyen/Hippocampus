function [mapG,mapGdummy,maxrate]= plotmap(mapL,objtype,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot rate map in either place, spatial view or head direction frames
% 
% Inputs:
% ax: plotting axis
% mapL: linear rate map e.g. vms.data.maps_adsm
% objtype: 'place' or 'spatialview' or 'headdirection' to specify frame 
% 
% Outputs:
% mapG: rate map in grid form
% mapGdummy: map of bin numbers in grid form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for linearity of input map
if min(size(mapL)) ~= 1
    error('Input map must be linear, not grid')
end
ax = gca;

if strcmp(objtype,'place') || strcmp(objtype,'view')
    axis(ax,'tight');

    samesizeaxis = false;
    if strcmp(objtype,'place') && ~samesizeaxis % Plot place without z axis
        mapG = flipud(reshape(mapL, 40, 40)');
        mapGdummy = flipud(reshape(1:1600, 40, 40)');

        % Set up surf frame for plotting
        floor_x = repmat(0:40, 41, 1);
        floor_y = flipud(repmat([0:40]', 1, 41));
        floor_z = zeros(41,41);

        floor = flipud(reshape(mapL, 40, 40)');

        % Plot floor
        surf(floor_x, floor_y, floor_z, floor);
        alpha 1; shading flat;
        hold on;

        axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim')])); % Make axes square around biggest value
        axis(ax,[0 axlim 0 axlim]);  
    else 
        if strcmp(objtype,'place')
            % Insert floor place map into larger 3D view setting
            mapLtemp = mapL;
            mapL = nan(1,5122);
            mapL(3:3+1600-1) = mapLtemp;
            mapG = flipud(reshape(mapLtemp, 40, 40)');
            mapGdummy = flipud(reshape(1:1600, 40, 40)');
        end
    
        mapLdummy = 1:length(mapL);
        % Set up surf frame for plotting
        floor_x = repmat(0:40, 41, 1);
        floor_y = flipud(repmat([0:40]', 1, 41));
        floor_z = zeros(41,41);
    
        ceiling_x = floor_x;
        ceiling_y = floor_y;
        ceiling_z = 40.*ones(41,41);
    
        walls_x = repmat([0.*ones(1,40) 0:39 40.*ones(1,40) 40:-1:0], 9, 1);
        walls_y = repmat([0:39 40.*ones(1,40) 40:-1:1 0.*ones(1,41)], 9, 1);
        walls_z = repmat([24:-1:16]', 1, 40*4 + 1);
    
        P1_x = repmat([24.*ones(1,8) 24:31 32.*ones(1,8) 32:-1:24], 6, 1);
        P1_y = repmat([8:15 16.*ones(1,8) 16:-1:9 8.*ones(1,9)], 6, 1);
        PX_z = repmat([21:-1:16]', 1, 8*4 + 1);
    
        P2_x = repmat([8.*ones(1,8) 8:15 16.*ones(1,8) 16:-1:8], 6, 1);
        P2_y = P1_y;
    
        P3_x = P1_x;
        P3_y = repmat([24:31 32.*ones(1,8) 32:-1:25 24.*ones(1,9)], 6, 1);
    
        P4_x = P2_x;
        P4_y = P3_y;
    
        floor = flipud(reshape(mapL(3:3+1600-1), 40, 40)');
        floordum = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');
    
        % ceiling follows floor mapping, top down view
        ceiling = flipud(reshape(mapL(1603:1603+1600-1), 40, 40)');
        ceilingdum = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');
    
        % from top down, slit walls at bottom left corner, open outwards.
        % start from row closest to ground, rightwards, then climb rows
        walls = flipud(reshape(mapL(3203:3203+1280-1), 40*4, 8)');
        wallsdum = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');
    
        % BL - bottom left, and so on, from top view, same slicing as walls
        % pillar width 8, height 5
        P1_BR = flipud(reshape(mapL(4483:4483+160-1), 8*4, 5)');
        P1_BRdum = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
        P2_BL = flipud(reshape(mapL(4643:4643+160-1), 8*4, 5)');
        P2_BLdum = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
        P3_TR = flipud(reshape(mapL(4803:4803+160-1), 8*4, 5)');
        P3_TRdum = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
        P4_TL = flipud(reshape(mapL(4963:4963+160-1), 8*4, 5)');
        P4_TLdum = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
        if strcmp(objtype,'view')
            mapG = { NaN; NaN; floor; ceiling; walls; P1_BR; P2_BL; P3_TR; P4_TL };
            mapGdummy = { NaN; NaN; floordum; ceilingdum; wallsdum; P1_BRdum; P2_BLdum; P3_TRdum; P4_TLdum };
        end
    
        % Pad with NaNs for surf plots
        P1_BR = [P1_BR; nan(1,size(P1_BR,2))];
        P1_BR = [P1_BR nan(size(P1_BR,1),1)];
    
        P2_BL = [P2_BL; nan(1,size(P2_BL,2))];
        P2_BL = [P2_BL nan(size(P2_BL,1),1)];        
    
        P3_TR = [P3_TR; nan(1,size(P3_TR,2))];
        P3_TR = [P3_TR nan(size(P3_TR,1),1)];                
    
        P4_TL = [P4_TL; nan(1,size(P4_TL,2))];
        P4_TL = [P4_TL nan(size(P4_TL,1),1)];
    
        % Plot floor
        surf(floor_x, floor_y, floor_z, floor);
        alpha 1; shading flat;
        hold on;
    
        % Plot ceiling and walls
        surf(ceiling_x, ceiling_y, ceiling_z, ceiling);
        alpha 1; shading flat;
        surf(walls_x, walls_y, walls_z, walls);      
        alpha 1; shading flat;
    
        % Plot pillars
        surf(P1_x, P1_y, PX_z, P1_BR);
        alpha 1; shading flat;
        surf(P2_x, P2_y, PX_z, P2_BL);
        alpha 1; shading flat;
        surf(P3_x, P3_y, PX_z, P3_TR);
        alpha 1; shading flat;
        surf(P4_x, P4_y, PX_z, P4_TL);

        axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim'),get(ax, 'zlim')])); % Make axes square around biggest value
        axis(ax,[0 axlim 0 axlim 0 axlim]);  

    end

    % Display parameters
    if strcmp(objtype,'view')
        if ~isnan(nanmax(mapL(3:end))) && nanmax(mapL(3:end)) ~= 0
            maxrate = nanmax(mapL(3:end));
        else
            maxrate = 1;
        end
    else
        if ~isnan(nanmax(mapL)) && nanmax(mapL) ~= 0
            maxrate = nanmax(mapL);
        else 
            maxrate = 1;
        end
    end
    set(ax,'CLim',[0 maxrate],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
    
    colormap jet;
    
    alpha 1; shading flat; 
    view(-35,20);
    
elseif strcmp(objtype,'headdirection')
    
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
    else
        % If plotting both primary and secondary maps in 3D
        if strcmp(varargin{1},'place')
            raise = 60; % 20;
        elseif strcmp(varargin{1},'view')
            raise = 60;
        end
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
    else
        if ~strcmp(varargin{1},'headdirection')
            if strcmp(varargin{1},'place')
                raise = 60; % 20;
            elseif strcmp(varargin{1},'view')
                raise = 60;
            end
            % axis(ax, 'tight');
            axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim')])); % Make axes square around biggest value
            if axlim == 0
                axlim = 1;
            end
            axsize = max([abs(x),abs(y)]); 
            if isnan(axsize)
                axsize = 20;
            end
            line('xdata',0.95*[20-axsize 20+axsize],'ydata',[20 20],'zdata',[raise raise],'parent',ax);   
            line('xdata',[20 20],'ydata',0.95*[20-axsize 20+axsize],'zdata',[raise raise],'parent',ax); % centre-crossing axes
            if strcmp(varargin{1},'headdirection')
                axis(ax,[1.1*(20-axsize) 1.1*axlim 1.1*(20-axsize) 1.1*axlim]); 
%                 axis(ax, 'square', 'off');
%                 axis(ax, 'square', 'off', 'tight');
            else
%                 axis(ax,'square','off');
            end
            ax.FontSize = 14;
            view(-35,20);
        else
            % axis(ax, 'tight');
            axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim')])); % Make axes square around biggest value
            axis(ax,[-1.1*axlim 1.1*axlim -1.1*axlim 1.1*axlim]);                         %   ..
            line('xdata',0.95*[-axlim axlim],'ydata',[0 0],'parent',ax);   
            line('xdata',[0 0],'ydata',0.95*[-axlim axlim],'parent',ax); % centre-crossing axes
    %         axis(ax, 'square', 'off');
%             axis(ax, 'square', 'off', 'tight');
        end
    end
    view(-35,20);
    colormap jet;
    % colorbar;

    mapG = [x' y'];
    mapGdummy = [[1:size(x,2)-1 1]' [1:size(y,2)-1 1]'];
    
    
end
axis(ax,'tight');
w = ax.OuterPosition(3)*0.02;
x = ax.OuterPosition(1) + ax.OuterPosition(3);
y = ax.OuterPosition(2) + ax.OuterPosition(4)/4;
h = ax.OuterPosition(4)/2;

c = colorbar(ax,'Position',[x y w h]);