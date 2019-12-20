function [out,surfx,surfy,surfz,surfmap] = plotspatialview2(bb,map,lin_grid_reference,ax,mode,shading_mode)

% Plot maps
switch mode
    
    case 0    
        switch bb
            case 1 % Cue
                % Plot 
                hold(ax,'on');
                surfx = ones(9,41);
                surfy = repmat((0:40),9,1);
                surfz = repmat((0:8)',1,41);
                surfmap = nan(9,40);
                surfmap(2,20) = map;

                surf(surfx,surfy,surfz,surfmap);
                if ~shading_mode; shading flat; else; shading faceted; end;
                hold(ax,'off');
            case 2
                % Plot 
                hold(ax,'on');
                surfx = ones(9,41);
                surfy = repmat((0:40),9,1);
                surfz = repmat((0:8)',1,41);
                surfmap = nan(9,40);
                surfmap(6,20) = map;

                surf(surfx,surfy,surfz,surfmap);
                if ~shading_mode; shading flat; else; shading faceted; end;
                hold(ax,'off');
            case 3 % Ground
                hold(ax,'on');
                surfx = repmat((0:40)',1,41);
                surfy = repmat(0:40,41,1);
                surfz = -1.*ones(41);
                surfmap = map;
                lin_grid_reference1 = lin_grid_reference;

                g = surf(surfx,surfy,surfz,surfmap,'Tag','g','UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;
                box off;
                hold(ax,'off');
                out = {g};
            case 4 % Ceiling
                hold(ax,'on');
                surfx = repmat((0:40)',1,41);
                surfy = repmat(0:40,41,1);
                surfz = 2.*ones(41);
                surfmap = map;
                lin_grid_reference1 = lin_grid_reference;

                c = surf(surfx,surfy,surfz,surfmap,'Tag','c','UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;
                box off;
                hold(ax,'off');
                out = {c};
            case 5 % Walls
                hold(ax,'on');
        %         if half == 1
                    % Left 
                    surfx = repmat(fliplr(0:40),9,1);
                    surfy = zeros(9,41);
                    surfz = repmat(((0:8)./8)',1,41);
                    surfmap = flipud(map(:,1:40));
                    lin_grid_reference1 = flipud(lin_grid_reference(:,1:40));

                    w_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                    if ~shading_mode; shading flat; else; shading faceted; end;

                    % Top
                    surfx = zeros(9,41);
                    surfy = repmat((0:40),9,1);
                    surfz = repmat(((0:8)./8)',1,41);
                    surfmap = flipud(map(:,41:80));
                    lin_grid_reference1 = flipud(lin_grid_reference(:,41:80));

                    w_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                    if ~shading_mode; shading flat; else; shading faceted; end;

        %             hold(ax,'off');
        %         elseif half == 2
                    % Right
                    surfx = repmat((0:40),9,1);
                    surfy = repmat(40,9,41);
                    surfz = repmat(((0:8)./8)',1,41);
                    surfmap = flipud(map(:,81:120));
                    lin_grid_reference1 = flipud(lin_grid_reference(:,81:120));

                    w_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                    if ~shading_mode; shading flat; else; shading faceted; end;

                    % Bottom
                    surfx = repmat(40,9,41);
                    surfy = repmat(fliplr(0:40),9,1);
                    surfz = repmat(((0:8)./8)',1,41);
                    surfmap = flipud(map(:,121:160));
                    lin_grid_reference1 = flipud(lin_grid_reference(:,121:160));

                    w_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                    if ~shading_mode; shading flat; else; shading faceted; end;

        %         end
                box off;
                hold(ax,'off');
                out = {w_l, w_t, w_r, w_b};

            case 6 % Pillar1 (Bottom Right)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(24:32),6,1);
                surfy = repmat(24,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,1:8));
                lin_grid_reference1 = flipud(lin_grid_reference(:,1:8));

                p1_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Top
                surfx = repmat(24,6,9);
                surfy = repmat((24:32),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,9:16));
                lin_grid_reference1 = flipud(lin_grid_reference(:,9:16));

                p1_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Right
                surfx = repmat((24:32),6,1);
                surfy = repmat(32,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,17:24));
                lin_grid_reference1 = flipud(lin_grid_reference(:,17:24));

                p1_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Bottom
                surfx = repmat(32,6,9);
                surfy = repmat(fliplr(24:32),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,25:32));
                lin_grid_reference1 = flipud(lin_grid_reference(:,25:32));

                p1_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;
                box off;
                hold(ax,'off');
                out = {p1_l, p1_t, p1_r, p1_b};

            case 7 % Pillar 2 (Bottom Left)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(24:32),6,1);
                surfy = repmat(8,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,1:8));
                lin_grid_reference1 = flipud(lin_grid_reference(:,1:8));

                p2_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Top
                surfx = repmat(24,6,9);
                surfy = repmat((8:16),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,9:16));
                lin_grid_reference1 = flipud(lin_grid_reference(:,9:16));

                p2_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Right
                surfx = repmat((24:32),6,1);
                surfy = repmat(16,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,17:24));
                lin_grid_reference1 = flipud(lin_grid_reference(:,17:24));

                p2_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Bottom
                surfx = repmat(32,6,9);
                surfy = repmat(fliplr(8:16),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,25:32));
                lin_grid_reference1 = flipud(lin_grid_reference(:,25:32));

                p2_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;
                box off;
                hold(ax,'off');
                out = {p2_l, p2_t, p2_r, p2_b};

            case 8 % Pillar3 (Top Right)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(8:16),6,1);
                surfy = repmat(24,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,1:8));
                lin_grid_reference1 = flipud(lin_grid_reference(:,1:8));

                p3_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Top
                surfx = repmat(8,6,9);
                surfy = repmat((24:32),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,9:16));
                lin_grid_reference1 = flipud(lin_grid_reference(:,9:16));

                p3_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Right
                surfx = repmat((8:16),6,1);
                surfy = repmat(32,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,17:24));
                lin_grid_reference1 = flipud(lin_grid_reference(:,17:24));

                p3_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Bottom
                surfx = repmat(16,6,9);
                surfy = repmat(fliplr(24:32),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,24:32));
                lin_grid_reference1 = flipud(lin_grid_reference(:,24:32));

                p3_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;
                box off;
                hold(ax,'off');
                out = {p3_l, p3_t, p3_r, p3_b};

            case 9 % Pillar4 (Top Left)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(8:16),6,1);
                surfy = repmat(8,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,1:8));
                lin_grid_reference1 = flipud(lin_grid_reference(:,1:8));

                p4_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Top
                surfx = repmat(8,6,9);
                surfy = repmat((8:16),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,9:16));
                lin_grid_reference1 = flipud(lin_grid_reference(:,9:16));

                p4_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Right
                surfx = repmat((8:16),6,1);
                surfy = repmat(16,6,9);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,17:24));
                lin_grid_reference1 = flipud(lin_grid_reference(:,17:24));

                p4_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;

                % Bottom
                surfx = repmat(16,6,9);
                surfy = repmat(fliplr(8:16),6,1);
                surfz = repmat(((0:5)./5)',1,9);
                surfmap = flipud(map(:,25:32));
                lin_grid_reference1 = flipud(lin_grid_reference(:,25:32));

                p4_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                if ~shading_mode; shading flat; else; shading faceted; end;
                box off;
                hold(ax,'off');
                out = {p4_l, p4_t, p4_r, p4_b};
        end
        
    case 1
        
        switch bb
            case 1 
                disp('check correct mode, last orgument');
            case 2
                disp('check correct mode, last orgument');
            case 3
                disp('check correct mode, last orgument');
            case 4
                disp('check correct mode, last orgument');
            case 5
                disp('check correct mode, last orgument');

            case 0
                hold(ax,'on');
                surfx = repmat(19:21,2,1);
                surfy = repmat(40,2,3);
                surfz = repmat(([-1 2])',1,3);
                surfmap = zeros(size(surfy)-1);
                lin_grid_reference1 = surfmap;
                compass = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                out = {compass};
                
            case 6 % Pillar1 (Bottom Right)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(24:32),2,1);
                surfy = repmat(24,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p1_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Top
                surfx = repmat(24,2,9);
                surfy = repmat((24:32),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p1_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Right
                surfx = repmat((24:32),2,1);
                surfy = repmat(32,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p1_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Bottom
                surfx = repmat(32,2,9);
                surfy = repmat(fliplr(24:32),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p1_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                box off;
                hold(ax,'off');
                out = {p1_l, p1_t, p1_r, p1_b};

            case 7 % Pillar 2 (Bottom Left)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(24:32),2,1);
                surfy = repmat(8,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p2_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Top
                surfx = repmat(24,2,9);
                surfy = repmat((8:16),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p2_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Right
                surfx = repmat((24:32),2,1);
                surfy = repmat(16,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p2_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Bottom
                surfx = repmat(32,2,9);
                surfy = repmat(fliplr(8:16),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p2_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                box off;
                hold(ax,'off');
                out = {p2_l, p2_t, p2_r, p2_b};

            case 8 % Pillar3 (Top Right)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(8:16),2,1);
                surfy = repmat(24,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p3_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Top
                surfx = repmat(8,2,9);
                surfy = repmat((24:32),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p3_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Right
                surfx = repmat((8:16),2,1);
                surfy = repmat(32,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p3_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Bottom
                surfx = repmat(16,2,9);
                surfy = repmat(fliplr(24:32),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p3_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                box off;
                hold(ax,'off');
                out = {p3_l, p3_t, p3_r, p3_b};

            case 9 % Pillar4 (Top Left)
                hold(ax,'on');
                % Left
                surfx = repmat(fliplr(8:16),2,1);
                surfy = repmat(8,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p4_l = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Top
                surfx = repmat(8,2,9);
                surfy = repmat((8:16),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p4_t = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Right
                surfx = repmat((8:16),2,1);
                surfy = repmat(16,2,9);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p4_r = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);

                % Bottom
                surfx = repmat(16,2,9);
                surfy = repmat(fliplr(8:16),2,1);
                surfz = repmat(((0:1)./5)'-1,1,9);
                surfmap = nan(size(surfy)-1);
                lin_grid_reference1 = surfmap;

                p4_b = surf(surfx,surfy,surfz,surfmap,'UserData',lin_grid_reference1);
                box off;
                hold(ax,'off');
                out = {p4_l, p4_t, p4_r, p4_b};
        end        
        
        
        
        
end
        
end
