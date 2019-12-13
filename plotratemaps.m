function [] = plotratemaps(objtype,sortsic,plottype,save,maptype,varargin)

% Plot type
% 1-3 = Place (batch). Maps are grouped by session. Whole session map.
% 2 = Place (batch). Maps are grouped by session. 1st half + 2nd half maps.
% 4 = Spatial view (batch). Maps are grouped by session
% 5 For single cells only: Place by spatial view. objtype should be 'place'
% 
% Varargin:
% 1st: cell directory

% Specify saved figure location
cwd = pwd;
if save
%     figdir = [ '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/MountainSort/SpatialviewOnPlace/' 'ratemaps' 'plottype' num2str(plottype) 'maptype' num2str(maptype)];
    figdir = [ '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/MountainSort/alpha1e2/SpatialView/postsfn/' 'ratemaps' 'plottype' num2str(plottype) 'maptype' num2str(maptype)];
    if exist(figdir,'dir') ~= 7
%         rmdir(figdir,'s');
%     end
        mkdir(figdir);
    end
end
video = 0;

if nargin > 6
    % Load cell list
    cellList = varargin(1);
    mapGrid = varargin{2};
    ax = varargin{3};
    cd(cellList{1});
    % Load object
    switch objtype
        case 'Place'
            objMain = load('vmplacecell.mat');
            objMain = objMain.vp;
            stepsMain = [objMain.data.gridSteps objMain.data.gridSteps];
            mapGrid = mapGrid{1};
        case 'Spatialview'
            objMain = load('spatialview.mat');
            objMain = objMain.sv;
            stepsMain = objMain.data.binDepths;
        case 'Direction'
            objMain = load('vmhd.mat');
            objMain = objMain.vmd;
            stepsMain = objMain.data.DirSteps;
    end
elseif nargin > 5
    cellList = varargin(1);
    cd(cellList{1});
    % Load object
    switch objtype
        case 'Place'
            objMain = load('vmplacecell.mat');
            objMain = objMain.vp;
            stepsMain = [objMain.data.gridSteps objMain.data.gridSteps];
%             mapGrid = mapGrid{1};
        case 'Spatialview'
            objMain = load('spatialview.mat');
            objMain = objMain.sv;
            stepsMain = objMain.data.binDepths;
        case 'Direction'
            objMain = load('vmhd.mat');
            objMain = objMain.vmd;
            stepsMain = objMain.data.DirSteps;
    end
else
    % Load cell list
    cwd = '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/MountainSort/alpha1e2/SpatialView/postsfn';
    cd(cwd);
    % load celllist 
    % celllist = textread('/Volumes/User/huimin/Desktop/condor_shuffle/cell_list copy.txt','%c');
    fid = fopen([cwd '/cell_list.txt'],'rt');
    cellList = textscan(fid,'%s','Delimiter','\n');
    cellList = cellList{1};
    % Make sure no empty cells
    notempty = ~cellfun(@isempty,cellList);
    cellList = cellList(notempty,:);
    % Load object
    objMain = load('sv.mat');
    objMain = objMain.sv;
end

% Generate unique identifier for each cell
s = regexp(cellList{1},'session');
identifiers = zeros(size(cellList,1),5);
cellid = cell(size(cellList,1),1);
missing = [];
for ii = 1:size(cellList,1)
    if exist(cellList{ii},'dir') == 7
        % Collect date, session, array, channel, cell
        identifiers(ii,:) = [str2double(cellList{ii}(s-9:s-2)) str2double(cellList{ii}(s+7:s+8)) ...
            str2double(cellList{ii}(s+15:s+16)) str2double(cellList{ii}(s+25:s+27)) str2double(cellList{ii}(s+33:s+34))];
        % Cell identifier
        cellid{ii} = horzcat(num2str(identifiers(ii,4)),'-',num2str(identifiers(ii,5)));
    else
        missing = [missing ii];
    end
end
% Remove missing cells
identifiers(missing,:) = [];
cellid(missing) = [];
cellList(missing) = [];
setsessions = unique(identifiers(:,1));
setcells = unique(cellid);

% Load firing rate maps depending on smoothing type
switch maptype
    case 'raw'
        maps = objMain.data.maps_raw;
    case 'boxcar'
        maps = objMain.data.maps_boxsmooth;
    case 'adaptive'
        maps = objMain.data.maps_adsmooth';
end

% Load nptdata objMainect 
nCells = objMain.data.numSets;
% Plot params
plotgridh = 6;
plotgridv = 4;
switch objtype
    case 'Place'
        viewangle = [75,-20];
    case 'Spatialview'
        if video
            viewangle = [75,-20];
        else
            viewangle = [75,-20];
        end
end

if plottype == 1 || plottype == 2  % Place maps
    % For each session, plot rate maps for each cell
    fig = 1;
    subpnum = 1;
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));
%         % Get shuffled SIC threshold for this session
%         sicsh = vpc.data.SICsh(2:end,cells_ind);
%         sicsh = reshape(sicsh,numel(sicsh),1);
%         shprc = prctile(sicsh,95);
        % Sort cells by descending SIC if required
        sic = objMain.data.SIC(cells_ind);
        [~,sici] = sort(sic,'descend');
        for jj = 1:length(cells_ind) % For each cell
            % Get cell index
            if sortsic
                cell_ind = cells_ind(sici(jj));
            else 
                cell_ind = cells_ind(jj);
            end
            
            if plottype == 2 || plottype == 3 
                
                %% % Compute first and second half correlation
                
                % Find figure number
                if jj*2 > plotgridh * plotgridv && mod((jj*2), (plotgridh * plotgridv)) == 2
                    % Save figure
                    if save
                        % Save previous figure
                        figtitle = [num2str(setsessions(ii)) '-' num2str(floor((jj*2)/(plotgridh * plotgridv)))];
                        saveas(h,figtitle,'png');
                        print('-painters',figtitle,'-dsvg');
                    end
                    fig = fig + 1;
                    subpnum = 1;
                end
                
                %%%%%%% CORRELATION %%%%%%%%%%%
                
                % Intra-session correlation %%%%%% NOTE: Should use boxcar
                % smoothed map
                map1 = objMain.data.maps_adsmooth1sthalf(:,cell_ind);
                vis1 = ~isnan(map1);
                SI1 = objMain.data.SIC1sthalf(cell_ind);
                map2 = objMain.data.maps_adsmooth2ndhalf(:,cell_ind);
                vis2 = ~isnan(map2);
                SI2 = objMain.data.SIC2ndhalf(cell_ind);
                vis = vis1 & vis2; % Correlate only visited bins;
                intracorr = corr2(map1(vis), map2(vis));
                % Get shuffled SI cutoff for this cell - 95th percentile
                SIthr = prctile(objMain.data.SICsh(2:end,cell_ind),95);
                
                % Plot
                
                for kk = 1:2
                    
                    % Get map
                    if kk == 1
                        mapLin = map1;
                    else
                        mapLin = map2;
                    end
                    % Plot
                    h = figure(fig);
                    hold on;
                    colormap(jet);
                    set(h,'Name',num2str(setsessions(ii)),'Units','Normalized','Position',[0 1 1 1]);
                    ax = subplot(plotgridv,plotgridh,subpnum);
                    % Restructure bins from linear to square grid
                    mapGrid = nan(objMain.data.Args.GridSteps);
                    for mm = 1:objMain.data.Args.GridSteps
                        mapGrid(mm,:) = mapLin( (mm-1)*objMain.data.Args.GridSteps+1:mm*objMain.data.Args.GridSteps );
                    end
                    pcolor([mapGrid nan(size(mapGrid,1),1); nan(1, size(mapGrid,2)+1)]);
                    shading flat;
                    if ~isnan(max(max(mapGrid))) && max(max(mapGrid)) ~= 0
                        set(ax,'CLim',[0 max(max(mapGrid))]);
                    else
                        set(ax,'CLim',[0 1]);
                    end
                    subpax = gca;
%                     set(subpax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{});
                    subpax.YLabel.String = horzcat(num2str(max(mapLin)),'Hz');
                    if kk == 1
                        subpax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' corr=',num2str(intracorr,2), ', SI=',num2str(SI1,2), '/',num2str(SIthr,2));
                        if SI1 > SIthr
                            subpax.Title.Color = 'r';
                        else
                            subpax.Title.Color = 'k';
                        end
                    else
                        subpax.Title.String = horzcat('SI=',num2str(SI2,2),'/',num2str(SIthr,2));
                        if SI2 > SIthr
                            subpax.Title.Color = 'r';
                        else
                            subpax.Title.Color = 'k';
                        end
                    end
                    subpnum = subpnum + 1;
                    hold off;
                end
                
            else 
              
                %% Plot 1 map for 1 cell
                
                % Find figure number
                if jj > plotgridh * plotgridv && mod(jj, (plotgridh * plotgridv)) == 1
                    % Save figure
                    if save
                        % Save previous figure
                        figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                        saveas(h,figtitle,'png');
                        print('-painters',figtitle,'-dsvg');
                    end
                    fig = fig + 1;
                    subpnum = 1;
                end
                
                % Get shuffled SI cutoff for this cell - 95th percentile
                SI = objMain.data.SICsh(1,cell_ind);
                SIthr = prctile(objMain.data.SICsh(2:end,cell_ind),95);
                % Get map
                if nargin <= 6 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
                    mapLin = maps(:,cell_ind);
                    % Restructure bins from linear to square grid
                    mapGrid = nan(objMain.data.Args.GridSteps);
                    for mm = 1:objMain.data.Args.GridSteps
                        mapGrid(mm,:) = mapLin( (mm-1)*objMain.data.Args.GridSteps+1:mm*objMain.data.Args.GridSteps );
                    end
                    h = figure(fig);
                    ax = subplot(plotgridv,plotgridh,subpnum);
                end

                % Setup main object
                h = gcf;
                hold on;
                colormap(jet);
%                 figname = horzcat('Place',': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', field',num2str(fieldnum));
                figname = horzcat('Place',': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)));
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
                ax.DataAspectRatioMode = 'manual';
                ax.DataAspectRatio = [1 1 1];
                set(ax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XColor','none','YColor','none','ZColor','none','GridLineStyle','none');
%                 ax.YLabel.String = horzcat(num2str(max(mapLin)),'Hz');
                ax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' SI=',num2str(SI,2),'/',num2str(SIthr,2),',',horzcat(num2str(max(max(mapGrid)),3),'Hz'));
                if objMain.data.SIC(cell_ind,1) > SIthr
                    ax.Title.Color = 'r';
                else
                    ax.Title.Color = 'k';
                end
                
                % Plot main object
                surfx = repmat((0:40)',1,41);
                surfy = repmat(0:40,41,1);
                surfz = zeros(41);
                surf(surfx,surfy,surfz,mapGrid);
                view(ax,viewangle);
                shading flat;
                if ~isnan(nanmax(nanmax(mapGrid))) && nanmax(nanmax(mapGrid)) ~= 0
                    set(ax,'CLim',[0 nanmax(nanmax(mapGrid))]);
                else
                    set(ax,'CLim',[0 1]);
                end
                
                subpnum = subpnum + 1;
                hold off;
                
            end
                
        end
        if save
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
            print('-painters',figtitle,'-dsvg');
        end
        fig = fig + 1;
        subpnum = 1;
    end
    
elseif strcmp(objtype,'Direction')
    
    % For each session, plot rate maps for each cell
    fig = 1;
    subpnum = 1;
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));
%         % Get shuffled SIC threshold for this session
%         sicsh = vpc.data.SICsh(2:end,cells_ind);
%         sicsh = reshape(sicsh,numel(sicsh),1);
%         shprc = prctile(sicsh,95);
        % Sort cells by descending SIC if required
        sic = objMain.data.SIC(cells_ind);
        [~,sici] = sort(sic,'descend');
        for jj = 1:length(cells_ind) % For each cell
            % Get cell index
            if sortsic
                cell_ind = cells_ind(sici(jj));
            else 
                cell_ind = cells_ind(jj);
            end
            
            if plottype == 2 || plottype == 3 
                
                %% % Compute first and second half correlation
                
                % Find figure number
                if jj*2 > plotgridh * plotgridv && mod((jj*2), (plotgridh * plotgridv)) == 2
                    % Save figure
                    if save
                        % Save previous figure
                        figtitle = [num2str(setsessions(ii)) '-' num2str(floor((jj*2)/(plotgridh * plotgridv)))];
                        saveas(h,figtitle,'png');
                        print('-painters',figtitle,'-dsvg');
                    end
                    fig = fig + 1;
                    subpnum = 1;
                end
                
                %%%%%%% CORRELATION %%%%%%%%%%%
                
                % Intra-session correlation %%%%%% NOTE: Should use boxcar
                % smoothed map
                map1 = objMain.data.maps_adsmooth1sthalf(:,cell_ind);
                vis1 = ~isnan(map1);
                SI1 = objMain.data.SIC1sthalf(cell_ind);
                map2 = objMain.data.maps_adsmooth2ndhalf(:,cell_ind);
                vis2 = ~isnan(map2);
                SI2 = objMain.data.SIC2ndhalf(cell_ind);
                vis = vis1 & vis2; % Correlate only visited bins;
                intracorr = corr2(map1(vis), map2(vis));
                % Get shuffled SI cutoff for this cell - 95th percentile
                SIthr = prctile(objMain.data.SICsh(2:end,cell_ind),95);
                
                % Plot
                
                for kk = 1:2
                    
                    % Get map
                    if kk == 1
                        mapLin = map1;
                    else
                        mapLin = map2;
                    end
                    % Plot
                    h = figure(fig);
                    hold on;
                    colormap(jet);
                    set(h,'Name',num2str(setsessions(ii)),'Units','Normalized','Position',[0 1 1 1]);
                    ax = subplot(plotgridv,plotgridh,subpnum);
                    % Restructure bins from linear to square grid
                    mapGrid = nan(objMain.data.Args.GridSteps);
                    for mm = 1:objMain.data.Args.GridSteps
                        mapGrid(mm,:) = mapLin( (mm-1)*objMain.data.Args.GridSteps+1:mm*objMain.data.Args.GridSteps );
                    end
                    pcolor([mapGrid nan(size(mapGrid,1),1); nan(1, size(mapGrid,2)+1)]);
                    shading flat;
                    if ~isnan(max(max(mapGrid))) && max(max(mapGrid)) ~= 0
                        set(ax,'CLim',[0 max(max(mapGrid))]);
                    else
                        set(ax,'CLim',[0 1]);
                    end
                    subpax = gca;
%                     set(subpax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{});
                    subpax.YLabel.String = horzcat(num2str(max(mapLin)),'Hz');
                    if kk == 1
                        subpax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' corr=',num2str(intracorr,2), ', SI=',num2str(SI1,2), '/',num2str(SIthr,2));
                        if SI1 > SIthr
                            subpax.Title.Color = 'r';
                        else
                            subpax.Title.Color = 'k';
                        end
                    else
                        subpax.Title.String = horzcat('SI=',num2str(SI2,2),'/',num2str(SIthr,2));
                        if SI2 > SIthr
                            subpax.Title.Color = 'r';
                        else
                            subpax.Title.Color = 'k';
                        end
                    end
                    subpnum = subpnum + 1;
                    hold off;
                end
                
            else 
              
                %% Plot 1 map for 1 cell
                
                % Find figure number
                if jj > plotgridh * plotgridv && mod(jj, (plotgridh * plotgridv)) == 1
                    % Save figure
                    if save
                        % Save previous figure
                        figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                        saveas(h,figtitle,'png');
                        print('-painters',figtitle,'-dsvg');
                    end
                    fig = fig + 1;
                    subpnum = 1;
                end
                
                % Get shuffled SI cutoff for this cell - 95th percentile
                SI = objMain.data.SICsh(1,cell_ind);
                SIthr = prctile(objMain.data.SICsh(2:end,cell_ind),95);
                % Get shuffled RV cutoff for this cell - 95th percentile
                RV = objMain.data.RV(1,cell_ind);
                RVthr = prctile(objMain.data.RVsh(2:end,cell_ind),95);
                % Get map
                if nargin <= 6 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
                    mapLin = maps(:,cell_ind);
                    h = figure(fig);
                    ax = subplot(plotgridv,plotgridh,subpnum);
                else
                    h = gcf;
                end

                % Setup main object
                hold on;
                colormap(jet);
                figname = horzcat('Place',': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)));
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
                ax.DataAspectRatioMode = 'manual';
                ax.DataAspectRatio = [1 1 1];
                set(ax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XColor','none','YColor','none','ZColor','none','GridLineStyle','none');
                ax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', RVthr = ',num2str(RVthr,2));
%                 if SI > SIthr
%                     ax.Title.Color = 'r';
%                 else
%                     ax.Title.Color = 'k';
%                 end
                
                % Plot main object
                % Convert to x,y %
                bin_ang = (2*pi) / length(mapLin);
                angles = (bin_ang/2 : bin_ang : (2*pi)-(bin_ang/2));
%                 angles=angles+opt.rotateDir;
%                 angles(angles>(2*pi))=angles(angles>(2*pi)) - 2*pi;
%                 angles(angles<=0)=angles(angles<=0) + 2*pi;
                [x, y] = pol2cart(angles', mapLin);
                x(end+1) = x(1);    % Make line meet itself at end
                y(end+1) = y(1);    %
                % Plot %
                plot(ax, x, y, 'k-');
                % Plot circular mean, if requested %
                if 0 % opt.circ_mean
                    cm=circ_mean(angles',mapLin);
                    [cmX cmY]=pol2cart(cm,max(mapLin));
                    line('xdata',[0 cmX],'ydata',[0 cmY],'parent',ax,'linewidth',2);
                end
                % Format axes %
                set(ax, 'color', 'w');  % For SCAn preview pane.
                axis(ax, 'tight');
                axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim')])); % Make axes square around biggest value
                axis(ax,[-axlim axlim -axlim axlim]);                         %   ..
                line('xdata',0.95*[-axlim axlim],'ydata',[0 0],'parent',ax);   line('xdata',[0 0],'ydata',0.95*[-axlim axlim],'parent',ax); % centre-crossing axes
                axis(ax, 'square', 'off', 'tight');
                % Write max rate and SI
                set(ax, 'tag', '');
                rate_string = num2str( max(max(mapLin)) , '%4.1f');
                SI_string = num2str(RV,'%1.3f');
                fs_val = 0.08;
                fu_val = 'normalized';
                pos_val = [0 1];
                hoz_val = 'left';
                ver_val = 'cap';
                text_color = [0 0 0];
                text('position', pos_val, 'units', 'normalized', 'HorizontalAlignment', hoz_val, 'string', rate_string, ...
                    'FontUnits', fu_val, 'VerticalAlignment', ver_val, 'fontsize', fs_val, 'color', text_color,'parent',ax);
                hoz_val = 'right';
                pos_val = [1 1];
                if RV > RVthr
                    text_color = [1 0 0];
                end
                text('position', pos_val, 'units', 'normalized', 'HorizontalAlignment', hoz_val, 'string', SI_string, ...
                    'FontUnits', fu_val, 'VerticalAlignment', ver_val, 'fontsize', fs_val, 'color', text_color,'parent',ax);
                
                subpnum = subpnum + 1;
                hold off;
                
            end
                
        end
        if save
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
            print('-painters',figtitle,'-dsvg');
        end
        fig = fig + 1;
        subpnum = 1;
    end
    
elseif plottype == 4 % Spatial view
    
    fig = 1;
    % Set up cell counts
    crossallthresh = 0;
    crosscellthresh = 0;
    crosspopthresh = 0;
    % Get full population SIC threshold
    sicfull = objMain.data.SICsh(2:end,:);
    sicfull = reshape(sicfull,numel(sicfull),1);
    shprcfull = prctile(sicfull,95);
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));
%         % Get shuffled SIC threshold for this session
%         sicsh = vpc.data.SICsh(2:end,cells_ind);
%         sicsh = reshape(sicsh,numel(sicsh),1);
%         shprc = prctile(sicsh,95);
        % Sort cells by descending SIC if required
        sic = objMain.data.SIC(cells_ind);
        [~,sici] = sort(sic,'descend');
        for jj = 1:length(cells_ind) % For each cell
            
            % Get cell index
            if sortsic
                cell_ind = cells_ind(sici(jj));
            else 
                cell_ind = cells_ind(jj);
            end
            % Get shuffled SI cutoff for this cell - 95th percentile
            SI = objMain.data.SICsh(1,cell_ind);
            SIthr = prctile(objMain.data.SICsh(2:end,cell_ind),95);
            
            if nargin <= 6 % If mapGrid is not already specified (i.e. if plotting for a batch of cells)
                % Get map
                mapLin = maps(:,cell_ind);
                % Restructure bins from linear to square grid
                stepsMain = objMain.data.binDepths;
                % Initialise empty cell array for grid maps
                mapGrid = cell(size(stepsMain,1),1);
                for bb = 1:size(stepsMain,1) % for each grid
                    % Initialise empty matrices
                    map = nan(stepsMain(bb,1),stepsMain(bb,2));
                    % Assign linear bin to grid bin
                    for mm = 1:stepsMain(bb,1)*stepsMain(bb,2) % For every point in linear map
                        if mod(mm,stepsMain(bb,2)) == 0
                            y = stepsMain(bb,2);
                        else
                            y = mod(mm,stepsMain(bb,2));
                        end
                        x = ceil(mm/stepsMain(bb,2));
                        indbins_lin = mm + sum(stepsMain(1:bb-1,1).*stepsMain(1:bb-1,2));
                        % Assign
                        map(x,y) = mapLin(indbins_lin);
                    end
                    % Collect output 
                    mapGrid{bb} = map;
                end
                % Set up figure
                h = figure(fig);
                
            end
            maxc = nan(1,size(stepsMain,1));
            for bb = 1:size(stepsMain)
                maxc(bb) = nanmax(nanmax(mapGrid{bb}));
            end
            maxC = nanmax(maxc);
            
            h = gcf;
            hold on;
            colormap(jet);
            figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', SI: ',num2str(SI,2),' of ',num2str(SIthr,2));
            set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
            ax = gca;
            set(ax,'XColor','none','YColor','none','ZColor','none');
            ax.Title.String = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' SI=',num2str(SI,2),'/',num2str(SIthr,2),'/',num2str(shprcfull,2),',',horzcat(num2str(maxC,3)),'Hz');
            ax.Title.FontSize = 14;
            ax.DataAspectRatioMode = 'manual';
            ax.DataAspectRatio = [1 1 1];
            
            hold(ax,'on')
            for bb = 1:size(stepsMain,1)

                plotspatialview(bb,plotgridv,plotgridh,mapGrid{bb});
                view(ax,viewangle);
                if ~isnan(maxC) && maxC ~= 0
                    set(ax,'CLim',[0 maxC],'ZLim',[0 40],'XLim',[0 40],'YLim',[0 40]);
                else
                    set(ax,'CLim',[0 1],'ZLim',[0 40],'XLim',[0 40],'YLim',[0 40]);
                end

            end
            for bb = 1:size(stepsMain,1)
                patchboundaries(bb);
            end
            if SI >= SIthr && SI >= shprcfull
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
            elseif SI >= SIthr && SI < shprcfull
                ax.Title.Color = 'b';
                crosscellthresh = crosscellthresh + 1;
            elseif SI >= shprcfull && SI < SIthr
                ax.Title.Color = 'g';
                crosspopthresh = crosspopthresh + 1;
            else
                ax.Title.Color = 'k';
            end
            
            if video
%                 ax.Title.String = horzcat(num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' ',horzcat(num2str(maxC,3)),'Hz');
                h.Units = 'normalized';
                h.Position = [0 0 0.5 1];
                ax.Position = [0 0 1 1];
                ax.Title.Color = 'none';
                ax.Color = 'none';
                ax.CameraViewAngle = 10;
                axis vis3d;
                videoname = ['Video ' 'FigNum' num2str(h.Number) ' ' objtype,' ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)) '.avi'];
                v = VideoWriter(videoname);
                open(v);
                for kstep = 1:360
                    viewanglemod = [viewangle(1)+(kstep-1) viewangle(2)];
                    disp(['kstep ' num2str(kstep) ' viewangle ' num2str(viewanglemod(1))]);
                    view(ax,viewanglemod);
                    frame = getframe(gcf);
                    writeVideo(v,frame);
                end
                close(v);
            end
            
            if save
                cd(figdir)
                % Save figure
                figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))) ' FigNum ',num2str(h.Number)];
                saveas(h,figtitle,'png');
                print('-painters',figtitle,'-dsvg');
                close(figure(fig));
                cd(cwd);
            end
            fig = fig + 1;
        end
        
    end
    disp(['Cross all thresh = ', num2str(crossallthresh),' cells']);
    disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
    disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);

end

cd(cwd);

function [surfx,surfy,surfz,surfmap] = plotspatialview(bb,plotgridv,plotgridh,map)

% Plot maps
switch bb
    case 1 % Cue
%         % Plot 
%         hold(ax,'on');
%         surfx = ones(9,41);
%         surfy = repmat((0:40),9,1);
%         surfz = repmat((0:8)',1,41);
%         surfmap = nan(9,40);
%         surfmap(2,20) = map;
%         
%         surf(surfx,surfy,surfz,surfmap);
%         shading flat;
%         hold(ax,'off');
    case 2
%         % Plot 
%         hold(ax,'on');
%         surfx = ones(9,41);
%         surfy = repmat((0:40),9,1);
%         surfz = repmat((0:8)',1,41);
%         surfmap = nan(9,40);
%         surfmap(6,20) = map;
%         
%         surf(surfx,surfy,surfz,surfmap);
%         shading flat;
%         hold(ax,'off');
    case 3 % Ground
        surfx = repmat((0:40)',1,41);
        surfy = repmat(0:40,41,1);
        surfz = zeros(41);
        surfmap = map;
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
    case 4 % Ceiling
        surfx = repmat((0:40)',1,41);
        surfy = repmat(0:40,41,1);
        surfz = repmat(40,41,41);
        surfmap = map;
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
    case 5 % Walls
        surfz = repmat((16:24)',1,41);
        % Left 
        surfx = repmat(fliplr(0:40),9,1);
        surfy = zeros(9,41);
        surfmap = flipud(map(:,1:40));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Top
        surfx = zeros(9,41);
        surfy = repmat((0:40),9,1);
        surfmap = flipud(map(:,41:80));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Right
        surfx = repmat((0:40),9,1);
        surfy = repmat(40,9,41);
        surfmap = flipud(map(:,81:120));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Bottom
        surfx = repmat(40,9,41);
        surfy = repmat(fliplr(0:40),9,1);
        surfmap = flipud(map(:,121:160));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
    case 6 % Pillar1 (Bottom Right)
        surfz = repmat((16:21)',1,9);
        % Left
        surfx = repmat(fliplr(24:32),6,1);
        surfy = repmat(24,6,9);
        surfmap = flipud(map(:,1:8));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Top
        surfx = repmat(24,6,9);
        surfy = repmat((24:32),6,1);
        surfmap = flipud(map(:,9:16));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Right
        surfx = repmat((24:32),6,1);
        surfy = repmat(32,6,9);
        surfmap = flipud(map(:,17:24));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Bottom
        surfx = repmat(32,6,9);
        surfy = repmat(fliplr(24:32),6,1);
        surfmap = flipud(map(:,25:32));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
    case 7 % Pillar 2 (Bottom Left)
        surfz = repmat((16:21)',1,9);
        % Left
        surfx = repmat(fliplr(24:32),6,1);
        surfy = repmat(8,6,9);
        surfmap = flipud(map(:,1:8));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Top
        surfx = repmat(24,6,9);
        surfy = repmat((8:16),6,1);
        surfmap = flipud(map(:,9:16));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Right
        surfx = repmat((24:32),6,1);
        surfy = repmat(16,6,9);
        surfmap = flipud(map(:,17:24));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Bottom
        surfx = repmat(32,6,9);
        surfy = repmat(fliplr(8:16),6,1);
        surfmap = flipud(map(:,25:32));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
    case 8 % Pillar3 (Top Right)
        surfz = repmat((16:21)',1,9);
        % Left
        surfx = repmat(fliplr(8:16),6,1);
        surfy = repmat(24,6,9);
        surfmap = flipud(map(:,1:8));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Top
        surfx = repmat(8,6,9);
        surfy = repmat((24:32),6,1);
        surfmap = flipud(map(:,9:16));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Right
        surfx = repmat((8:16),6,1);
        surfy = repmat(32,6,9);
        surfmap = flipud(map(:,17:24));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Bottom
        surfx = repmat(16,6,9);
        surfy = repmat(fliplr(24:32),6,1);
        surfmap = flipud(map(:,24:32));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
    case 9 % Pillar4 (Top Left)
        surfz = repmat((16:21)',1,9);
        % Left
        surfx = repmat(fliplr(8:16),6,1);
        surfy = repmat(8,6,9);
        surfmap = flipud(map(:,1:8));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Top
        surfx = repmat(8,6,9);
        surfy = repmat((8:16),6,1);
        surfmap = flipud(map(:,9:16));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Right
        surfx = repmat((8:16),6,1);
        surfy = repmat(16,6,9);
        surfmap = flipud(map(:,17:24));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        % Bottom
        surfx = repmat(16,6,9);
        surfy = repmat(fliplr(8:16),6,1);
        surfmap = flipud(map(:,25:32));
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
end

function [] = patchboundaries(ax)

% for ii = 1:size(ax,1)
    
%     % Get axis of subplot
%     axes(ax{ii});
%     hold(ax{ii},'on');
    switch ax
        case 1 % Cue
            
        case 2  % Hint
            
        case 3  % Ground
            
            % Outer boundaries
            patch([0 0 40 40],[0 40 40 0], [0 0 0 0] ,[1 1 1 1],'FaceColor','none');
            % Pillars
            patch([24 24 32 32],[8 16 16 8], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Bottom left
            patch([8 8 16 16],[8 16 16 8], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Top left
            patch([8 8 16 16],[24 32 32 24], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Top right
            patch([24 24 32 32 ],[24 32 32 24], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Bottom right
            
        case 4  % Ceiling
            
            % Outer boundaries
            patch([0 0 40 40],[0 40 40 0], [40 40 40 40] ,[1 1 1 1],'FaceColor','none');
            
        case 5  % Walls
            
            % Outer boundaries
            patch([0 0 40 40],[0 0 0 0], [16 24 24 16] ,[1 1 1 1],'FaceColor','none');
            patch([0 0 0 0],[0 0 40 40], [16 24 24 16] ,[1 1 1 1],'FaceColor','none');
            patch([0 0 40 40], [40 40 40 40], [16 24 24 16] ,[1 1 1 1],'FaceColor','none');
            patch([40 40 40 40],[0 0 40 40], [16 24 24 16] ,[1 1 1 1],'FaceColor','none');
            
        case 6  % Pillar1
            
            % Outer boundaries
            patch([24 24 32 32],[24 24 24 24], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 24 24],[24 24 32 32], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 32 32],[32 32 32 32], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([32 32 32 32],[24 24 32 32], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');

            % Outline Rabit poster on m_wall_25
            patch([26.88 26.88 29.12 29.12],[32 32 32 32], [17.8 19.2 19.2 17.8] ,[1 1 1 1],'FaceColor','none');    
            
        case 7  % Pillar2
            
            % Outer boundaries
            patch([24 24 32 32],[8 8 8 8], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 24 24],[8 8 16 16], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 32 32],[16 16 16 16], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([32 32 32 32],[8 8 16 16], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
        
            % Outline Cat poster on m_wall_10
            patch([32 32 32 32],[10.88 10.88 13.12 13.12], [17.8 19.2 19.2 17.8] ,[1 1 1 1],'FaceColor','none');
%             patch([32 32 32 32],[10.88 10.88 13.12 13.12], [1.8 3.2 3.2 1.8] ,[1 1 1 1],'FaceColor','none');

            % Outline Pig poster on m_wall_29
            patch([24 24 24 24],[10.88 10.88 13.12 13.12], [17.8 19.2 19.2 17.8] ,[1 1 1 1],'FaceColor','none');
            
        case 8  % Pillar3
            
            % Outer boundaries
            patch([8 8 16 16],[24 24 24 24], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 8 8],[24 24 32 32], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 16 16],[32 32 32 32], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([16 16 16 16],[24 24 32 32], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');

            % Outline Croc poster on m_wall_4
            patch([16 16 16 16],[26.88 26.88 29.12 29.12], [17.8 19.2 19.2 17.8] ,[1 1 1 1],'FaceColor','none');

            % Outline Donkey poster on m_wall_15
            patch([8 8 8 8],[26.88 26.88 29.12 29.12], [17.8 19.2 19.2 17.8] ,[1 1 1 1],'FaceColor','none');
            
        case 9  % Pillar4
            
            % Outer boundaries
            patch([8 8 16 16],[8 8 8 8], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 8 8],[8 8 16 16], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 16 16],[16 16 16 16], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');
            patch([16 16 16 16],[8 8 16 16], [16 21 21 16] ,[1 1 1 1],'FaceColor','none');

            % Outline Camel poster on m_wall_20
            patch([10.88 10.88 13.12 13.12],[8 8 8 8], [17.8 19.2 19.2 17.8] ,[1 1 1 1],'FaceColor','none');
    end
    
% end







