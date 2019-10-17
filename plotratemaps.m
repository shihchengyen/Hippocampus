function [] = plotratemaps(objtype,sortsic,plottype,save,maptype,varargin)

% Plot type
% 1-3 = Place rate maps for a GROUP of cells. Maps are grouped by session
% 4 = Spatial view rate maps for a GROUP of cells. Maps are grouped by session
% 5 = For single cells only: Place by spatial view. objtype should be 'place'
cwd = pwd;

if nargin > 5
    % Load cell list
    if nargin > 6
        cellList = varargin(1);
        map2G = varargin{2};
        map2L = varargin{3};
        binDepths = varargin{4};
        fieldnum = varargin{5};
        peakcoords = varargin{6};
        cd(cellList{1});
    else
        cellList = varargin(1);
        cd(cellList{1});
    end
    % Load object
    switch objtype
        case 'place'
            obj = load('vmplacecell.mat');
            obj = obj.vp;
            gSteps = obj.data.gridSteps;
        case 'spatialview'
            obj = load('spatialview.mat');
            obj = obj.sv;
            binDepths = obj.data.binDepths;
    end
else
    % Load cell list
%     cwd = '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/MountainSort/alpha1e3/Spatialview/';
%     cd(cwd);
    % load celllist 
    % celllist = textread('/Volumes/User/huimin/Desktop/condor_shuffle/cell_list copy.txt','%c');

    fid = fopen([cwd '/cell_list.txt'],'rt');
    cellList = textscan(fid,'%s','Delimiter','\n');
    cellList = cellList{1};
    % Make sure no empty cells
    notempty = ~cellfun(@isempty,cellList);
    cellList = cellList(notempty,:);
    % Load object
    switch objtype
        case 'place'
            obj = load('vpc.mat');
            obj = obj.vp;
        case 'spatialview'
            obj = load('sv.mat');
            obj = obj.sv;
    end
end

% Generate unique identifier for each cell
s = regexp(cellList{1},'session');
identifiers = zeros(size(cellList,1),5);
cellid = cell(size(cellList,1),1);
missing = [];
for ii = 1:size(cellList,1)
    if exist(cellList{ii},'dir') == 7 
%         cd(cellList{ii});
%         if exist('vmplacecell.mat','file') == 2
            % Collect date, session, array, channel, cell
            identifiers(ii,:) = [str2double(cellList{ii}(s-9:s-2)) str2double(cellList{ii}(s+7:s+8)) ...
                str2double(cellList{ii}(s+15:s+16)) str2double(cellList{ii}(s+25:s+27)) str2double(cellList{ii}(s+33:s+34))];
            % Cell identifier
            cellid{ii} = horzcat(num2str(identifiers(ii,4)),'-',num2str(identifiers(ii,5)));
%         else
%             missing = [missing ii];
%         end
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
cd(cwd);

% Load firing rate maps depending on smoothing type
switch maptype
    case 'raw'
        maps = obj.data.maps_raw;
    case 'boxcar'
        maps = obj.data.maps_boxsmooth;
    case 'adaptive'
        maps = obj.data.maps_adsmooth;
end

% Load nptdata object 
nCells = obj.data.numSets;

% Specify saved figure location
cwd = pwd;
if save
    if nargin > 6
        figdir = ['/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/MountainSort/alpha1e2/SpatialView/sfn/' 'ratemaps' 'plottype' num2str(plottype) 'maptype' num2str(maptype) '/' num2str(identifiers(1)) 'ch' num2str(identifiers(4)), 'c' num2str(identifiers(5)) ];
    else
        figdir = [ '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/MountainSort/alpha1e2/SpatialView/sfn/' 'ratemaps' 'plottype' num2str(plottype) 'maptype' num2str(maptype)];
    end
    if exist(figdir,'dir') == 7
        rmdir(figdir,'s');
    end
    mkdir(figdir);
end

% Get overall SIC threshold across all sessions
allSICsh = obj.data.SICsh(2:end,:);
allSICsh = reshape(allSICsh,[size(allSICsh,1)*size(allSICsh,2),1]);
shprc_all = prctile(allSICsh,95);
    
if plottype == 1 || plottype == 2 
    % For each session, plot rate maps for each cell
    plotgridh = 6;
    plotgridv = 4;
    fig = 1;
    subpnum = 1;
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));
%         % Get shuffled SIC threshold for this session
%         sicsh = vpc.data.SICsh(2:end,cells_ind);
%         sicsh = reshape(sicsh,numel(sicsh),1);
%         shprc = prctile(sicsh,95);
        % Sort cells by descending SIC if required
        sic = obj.data.SIC(cells_ind);
        [~,sici] = sort(sic,'descend');
        for jj = 1:length(cells_ind) % For each cell
            % Get cell index
            if sortsic
                cell_ind = cells_ind(sici(jj));
            else 
                cell_ind = cells_ind(jj);
            end
            
            if plottype == 2 
                
                %% % Compute first and second half correlation
                
                % Find figure number
                if jj*2 > plotgridh * plotgridv && mod((jj*2), (plotgridh * plotgridv)) == 2
                    % Save figure
                    if save
                        cd(figdir);
                        % Save previous figure
                        figtitle = [num2str(setsessions(ii)) '-' num2str(floor((jj*2)/(plotgridh * plotgridv)))];
                        saveas(h,figtitle,'png');
                        print('-painters',figtitle,'-dsvg');
                        cd(cwd);
                        close(h);
                    end
                    fig = fig + 1;
                    subpnum = 1;
                end
                
                %%%%%%% CORRELATION %%%%%%%%%%%
                
                % Intra-session correlation %%%%%% NOTE: Should use boxcar
                % smoothed map
                map1 = obj.data.maps_adsmooth1sthalf(:,cell_ind);
                vis1 = ~isnan(map1);
                SI1 = obj.data.SIC1sthalf(cell_ind);
                map2G = obj.data.maps_adsmooth2ndhalf(:,cell_ind);
                vis2 = ~isnan(map2G);
                SI2 = obj.data.SIC2ndhalf(cell_ind);
                vis = vis1 & vis2; % Correlate only visited bins;
                intracorr = corr2(map1(vis), map2G(vis));
                % Get shuffled SI cutoff for this cell - 95th percentile
                shprc = prctile(obj.data.SICsh(2:end,cell_ind),95);
                
                % Plot
                plotgridh = 6;
                plotgridv = 4;
                
                for kk = 1:2
                    
                    % Get map
                    if kk == 1
                        mapLin = map1;
                    else
                        mapLin = map2G;
                    end
                    % Plot
                    h = figure(fig);
                    hold on;
                    colormap(jet);
                    set(h,'Name',num2str(setsessions(ii)),'Units','Normalized','Position',[0 1 1 1]);
                    ax = subplot(plotgridv,plotgridh,subpnum);
                    % Restructure bins from linear to square grid
                    mapGrid = nan(obj.data.Args.GridSteps);
                    for mm = 1:obj.data.Args.GridSteps
                        mapGrid(obj.data.Args.GridSteps-(mm-1),:) = mapLin( (mm-1)*obj.data.Args.GridSteps+1:mm*obj.data.Args.GridSteps );
                    end
%                     pcolor([mapGrid nan(size(mapGrid,1),1); nan(1, size(mapGrid,2)+1)]);
                    surfx = repmat((0:40)',1,41);
                    surfy = repmat(0:40,41,1);
                    surfz = zeros(41);
                    surf(surfx,surfy,surfz,mapGrid);
                    shading flat;
                    view(90,90);
                    shading flat;
                    if ~isnan(max(max(mapGrid))) && max(max(mapGrid)) ~= 0
                        set(ax,'CLim',[0 max(max(mapGrid))]);
                    else
                        set(ax,'CLim',[0 1]);
                    end
                    subpax = gca;
                    set(subpax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XLim',[0 40],'YLim',[0 40],'GridLineStyle','none','Color','none','Box','on');
                    subpax.YLabel.String = horzcat(num2str(max(mapLin)),'Hz');
                    if kk == 1
                        subpax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' corr=',num2str(intracorr,2), ', SI=',num2str(SI1,2), '/',num2str(shprc,2),'/',num2str(shprc_all,2));
                        if SI1 > shprc && SI1 > shprc_all
                            subpax.Title.Color = 'r';
                        else
                            subpax.Title.Color = 'k';
                        end
                    else
                        subpax.Title.String = horzcat('SI=',num2str(SI2,2),'/',num2str(shprc,2),'/',num2str(shprc_all,2));
                        if SI2 > shprc && SI2 > shprc_all
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
                        cd(figdir);
                        % Save previous figure
                        figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv)))];
                        saveas(h,figtitle,'png');
                        print('-painters',figtitle,'-dsvg');
                        cd(cwd);
                    end
                    fig = fig + 1;
                    subpnum = 1;
                end
                
                % Get shuffled SI cutoff for this cell - 95th percentile
                SI = obj.data.SICsh(1,cell_ind);
                shprc = prctile(obj.data.SICsh(2:end,cell_ind),95);
                % Get map
                mapLin = maps(:,cell_ind);
                % Restructure bins from linear to square grid
                mapGrid = nan(obj.data.Args.GridSteps);
                for mm = 1:obj.data.Args.GridSteps
                    mapGrid(obj.data.Args.GridSteps-(mm-1),:) = mapLin( (mm-1)*obj.data.Args.GridSteps+1:mm*obj.data.Args.GridSteps );
                end
                
                % Plot params
                plotgridh = 6;
                plotgridv = 4;
                h = figure(fig);
                hold on;
                colormap(jet);
                set(h,'Name',num2str(setsessions(ii)),'Units','Normalized','Position',[0 0.9 0.9 0.9]);
                ax = subplot(plotgridv,plotgridh,subpnum);
                % Plot
%                 pcolor([mapGrid nan(size(mapGrid,1),1); nan(1, size(mapGrid,2)+1)]);
                surfx = repmat((0:40)',1,41);
                surfy = repmat(0:40,41,1);
                surfz = zeros(41);
                surf(surfx,surfy,surfz,mapGrid);
                shading flat;
                view(90,90);
                if ~isnan(max(max(mapGrid))) && max(max(mapGrid)) ~= 0
                    set(ax,'CLim',[0 max(max(mapGrid))]);
                else
                    set(ax,'CLim',[0 1]);
                end
                
                %             colorbar;
                %             imagesc(mapGrid,[0,max(maps(:,cell_ind))]);
                %             imagesc(reshape(vpc.data.meanFRs(:,cell_ind),gSteps,gSteps));
                subpax = gca;
                set(subpax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XLim',[0 40],'YLim',[0 40],'GridLineStyle','none','Color','none','Box','on');
                subpax.YLabel.String = horzcat(num2str(max(mapLin)),'Hz');
                subpax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' SI=',num2str(SI,2),'/',num2str(shprc,2),'/',num2str(shprc_all,2));
                if obj.data.SIC(cell_ind,1) > shprc && obj.data.SIC(cell_ind,1) > shprc_all
                    subpax.Title.Color = 'r';
                else
                    subpax.Title.Color = 'k';
                end
                subpnum = subpnum + 1;
                hold off;
                
            end
                
        end
        if save
            cd(figdir);
            % Save figure
            if plottype == 1
                figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv)))];
            else
                figtitle = [num2str(setsessions(ii)) '-' num2str(ceil((2*length(cells_ind))/(plotgridh * plotgridv)))];
            end
            saveas(h,figtitle,'png');
            print('-painters',figtitle,'-dsvg');
            cd(cwd);
        end
        close(h);
        fig = fig + 1;
        subpnum = 1;
    end
    
elseif plottype == 3 || plottype == 4
    
    % Plot params
    plotgridh = 3;
    plotgridv = 3;
    % Plot params
    axnum = [   1 NaN 0.05 0.65 % Cue
                1 NaN 0.05 0.65 % Hint
                8 NaN 0.325 0.05    % Ground
                2 NaN 0.325 0.65    % Ceiling
%                 4 6   0.05 0.3
                5 NaN 0.325 0.325     % Walls
                5 NaN 0.325 0.325     % Pillar1
                5 NaN 0.325 0.325     % Pillar2
                5 NaN 0.325 0.325     % Pillar3
                5 NaN 0.325 0.325];   % Pillar4
    fig = 1;
    
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));
%         % Get shuffled SIC threshold for this session
%         sicsh = vpc.data.SICsh(2:end,cells_ind);
%         sicsh = reshape(sicsh,numel(sicsh),1);
%         shprc = prctile(sicsh,95);
        % Sort cells by descending SIC if required
        sic = obj.data.SIC(cells_ind);
        [~,sici] = sort(sic,'descend');
        for jj = 1:length(cells_ind) % For each cell
            
            % Get cell index
            if sortsic
                cell_ind = cells_ind(sici(jj));
            else 
                cell_ind = cells_ind(jj);
            end
                
            % Get shuffled SI cutoff for this cell - 95th percentile
            SI = obj.data.SICsh(1,cell_ind);
            shprc = prctile(obj.data.SICsh(2:end,cell_ind),95);
            
            h = figure(fig);
            hold on;
            colormap(jet);
            figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', SI=',num2str(SI,2),'/',num2str(shprc,2),'/',num2str(shprc_all,2));
            set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
            ax = gca;
            set(ax,'Visible','Off');
            
            if plottype == 4 % Intracorr
                
                % Intra-session correlation %%%%%% NOTE: Should use boxcar
                % smoothed map
                map1 = obj.data.maps_adsmooth1sthalf(:,cell_ind);
                vis1 = ~isnan(map1);
                SI1 = obj.data.SIC1sthalf(cell_ind);
                map2 = obj.data.maps_adsmooth2ndhalf(:,cell_ind);
                vis2 = ~isnan(map2);
                SI2 = obj.data.SIC2ndhalf(cell_ind);
                vis = vis1 & vis2; % Correlate only visited bins;
                intracorr = corr2(map1(vis), map2(vis));
                
                % Restructure bins from linear to square grid
                binDepths = obj.data.binDepths;
                maxC = nanmax(maps(:,cell_ind));
                for corrn = 1:2
                    if corrn == 1
                        mapLin = map1;
                        shift = -0.1;
                    else
                        mapLin = map2;
                        shift = 0.3;
                    end
                    % Initialise empty cell array for grid maps
                    grid_map = cell(size(binDepths,1),1);
                    axhandles2 = cell(size(binDepths,1),1);
                    for bb = 1:size(binDepths,1) % for each grid
                        % Initialise empty matrices
                        map = nan(binDepths(bb,1),binDepths(bb,2));
                        % Assign linear bin to grid bin
                        for mm = 1:binDepths(bb,1)*binDepths(bb,2) % For every point in linear map
                            if mod(mm,binDepths(bb,2)) == 0
                                y = binDepths(bb,2);
                            else
                                y = mod(mm,binDepths(bb,2));
                            end
                            x = ceil(mm/binDepths(bb,2));
                            indbins_lin = mm + sum(binDepths(1:bb-1,1).*binDepths(1:bb-1,2));
                            % Assign
                            map(x,y) = mapLin(indbins_lin);
                        end
                        % Collect output 
                        grid_map{bb} = map;

                        % PLOT 3D MAP
                        for half = 1:2
                            if ~isnan(axnum(bb,half))
        %                         ax = subplot(plotgridv,plotgridh,axnum(bb,half));
                                if half == 1
                                    ax = axes('Position',[axnum(bb,3)+shift axnum(bb,4) 0.35 0.4]);
                                else
                                    ax = axes('Position',[axnum(bb,3)+0.55+shift axnum(bb,4) 0.35 0.4]);
                                end
                                axhandles2{bb} = ax;
                                % Plot rate maps
                                plotspatialview(bb,plotgridv,plotgridh,map,ax,half);
                                if ~isnan(maxC) && maxC ~= 0
                                    set(ax,'CLim',[0 maxC],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
                                else
                                    set(ax,'CLim',[0 1],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
                                end
                                % Rotate view
                                view(ax,[120 70]);
                            end
                        end
                    end

                    % Plot outlines of shapes
                    patchboundaries(axhandles2);
                end
                
            elseif plottype == 3 % Single maps
                
                % Get map
                mapLin = maps(:,cell_ind);

%                 h = figure(fig);
%                 hold on;
%                 colormap(jet);
%                 figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', SI=',num2str(SI,2),'/',num2str(shprc,2),'/',num2str(shprc_all,2));
%                 set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
%                 ax = gca;
%                 set(ax,'Visible','Off');

                % Restructure bins from linear to square grid
                binDepths = obj.data.binDepths;
                maxC = nanmax(mapLin);
                % Initialise empty cell array for grid maps
                grid_map = cell(size(binDepths,1),1);
                axhandles2 = cell(size(binDepths,1),1);
                for bb = 1:size(binDepths,1) % for each grid
                    % Initialise empty matrices
                    map = nan(binDepths(bb,1),binDepths(bb,2));
                    % Assign linear bin to grid bin
                    for mm = 1:binDepths(bb,1)*binDepths(bb,2) % For every point in linear map
                        if mod(mm,binDepths(bb,2)) == 0
                            y = binDepths(bb,2);
                        else
                            y = mod(mm,binDepths(bb,2));
                        end
                        x = ceil(mm/binDepths(bb,2));
                        indbins_lin = mm + sum(binDepths(1:bb-1,1).*binDepths(1:bb-1,2));
                        % Assign
                        map(x,y) = mapLin(indbins_lin);
                    end
                    % Collect output 
                    grid_map{bb} = map;

                    % PLOT 3D MAP
                    for half = 1:2
                        if ~isnan(axnum(bb,half))
    %                         ax = subplot(plotgridv,plotgridh,axnum(bb,half));
                            if half == 1
                                ax = axes('Position',[axnum(bb,3) axnum(bb,4) 0.35 0.4]);
                            else
                                ax = axes('Position',[axnum(bb,3)+0.55 axnum(bb,4) 0.35 0.4]);
                            end
                            axhandles2{bb} = ax;
                            % Plot rate maps
                            plotspatialview(bb,plotgridv,plotgridh,map,ax,half);
                            if ~isnan(maxC) && maxC ~= 0
                                set(ax,'CLim',[0 maxC],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
                            else
                                set(ax,'CLim',[0 1],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
                            end
                            % Rotate view
                            view(ax,[120 70]);
                        end
                    end
                end

                % Plot outlines of shapes
                patchboundaries(axhandles2);
            end
            
            % Plot SI info
            ax = axes('Position',[0.05 0.325 0.1 0.1]);
            set(ax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XLim',[0 40],'YLim',[0 40],'GridLineStyle','none','Color','none','Box','on');
            ax.YLabel.String = horzcat(num2str(max(mapLin)),'Hz');
            if plottype == 4
                ax.XLabel.String = horzcat('corr = ',num2str(intracorr,2));
            end
            ax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' SI=',num2str(SI,2),'/',num2str(shprc,2),'/',num2str(shprc_all,2));
            if obj.data.SIC(cell_ind,1) > shprc && obj.data.SIC(cell_ind,1) > shprc_all
                ax.Title.Color = 'r';
            else
                ax.Title.Color = 'k';
            end
            
            if save
                cd(figdir)
                % Save figure
                figtitle = horzcat(objtype,' ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)));
%                     num2str(ceil(length(cells_ind)/(plotgridh * plotgridv)))];
                saveas(h,figtitle,'png');
                print('-painters',figtitle,'-dsvg');
                close(figure(fig));
                cd(cwd);
            end
            fig = fig + 1;
        end
        
    end
    
elseif plottype == 5
    
    % MAIN OBJECT
    % Set figure parameters
    h1 = figure;
    fignum = h1.Number;
    cell_ind = 1;
    % Get shuffled SI cutoff for this cell - 95th percentile
    SI = obj.data.SIC;
    shprc = prctile(obj.data.SICsh(2:end),95);
    % Get map
    mapLin = obj.data.maps_adsmooth;
    % Restructure bins from linear to square grid
    mapGrid = nan(obj.data.Args.GridSteps);
    for mm = 1:obj.data.Args.GridSteps
        mapGrid(obj.data.Args.GridSteps-(mm-1),:) = mapLin( (mm-1)*obj.data.Args.GridSteps+1:mm*obj.data.Args.GridSteps );
    end
    
    % SECONDARY OBJECT
    % Load object
    switch objtype
        case 'place'
            obj2 = load('spatialview.mat');
            obj2 = obj2.sv;
        case 'spatialview'
            obj2 = load('vmplacecell.mat');
            obj2 = obj2.vp;
    end
    obj2SI = obj2.data.SIC;
    obj2shprc = prctile(obj2.data.SICsh(2:end),95);
    
    % Plot main object
    hold on;
    colormap(jet);
    figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', field',num2str(fieldnum));
    set(h1,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
    % Plot
    surfx = repmat((0:40)',1,41);
    surfy = repmat(0:40,41,1);
    surfz = zeros(41);
    surf(surfx,surfy,surfz,mapGrid);
    view(90,90);
    shading flat;
    ax = gca;
    if ~isnan(max(max(mapGrid))) && max(max(mapGrid)) ~= 0
        set(ax,'CLim',[0 max(max(mapGrid))]);
    else
        set(ax,'CLim',[0 1]);
    end
    % Patch the field concerned
    patch([peakcoords(1)-1 peakcoords(1)-1 peakcoords(1)+2 peakcoords(1)+2],[peakcoords(2)-2 peakcoords(2)+1 peakcoords(2)+1 peakcoords(2)-2], [0 0 0 0] ,[1 1 1 1],'FaceColor','none');
    % Label 
    set(ax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XLim',[0 40],'YLim',[0 40]);
    ax.YLabel.String = horzcat(num2str(max(mapLin)),'Hz');
    ax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' SI=',num2str(SI,2),'/',num2str(shprc,2),'/',num2str(shprc_all,2));
    if obj.data.SIC(cell_ind,1) > shprc && obj.data.SIC(cell_ind,1) > shprc_all
        ax.Title.Color = 'r';
    else
        ax.Title.Color = 'k';
    end
    hold off;

    if save
        % Save figure
        cd(figdir);
        figtitle = ['PLACE' num2str(setsessions(ii))];
        saveas(h1,figtitle,'png');
        print('-painters',figtitle,'-dsvg');
        cd(cwd);
    end

    % Plot secondary object split by pixels of main object
    fignum = fignum + 1;
    h2 = figure(fignum);
    figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', SI=',num2str(obj2SI,2),'/',num2str(obj2shprc,2),', field',num2str(fieldnum));
    set(h2,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
    ax = gca;
    set(ax,'Visible','Off');
    % Plot for each pixel of main object
    pxSize = 0.3;
    pxcount = 1;
    for ii = 1:size(map2G,1)
        for jj = 1:size(map2G,2)
            % Get spatial view maps (unsmoothed)
            map = map2G{ii,jj};
            maxes = nan(size(map,1),1);
            for kk = 1:size(map,1)
                maxes(kk) = nanmax(map{kk}(:));
            end
            maxC = nanmax(maxes);
            
            % Plot params
            plotgridh = 3;
            plotgridv = 3;
            % Plot params
            axnum = [   1 NaN 0.05 0.65 % Cue
                    1 NaN 0.05 0.65 % Hint
                    8 NaN 0.325 0.05    % Ground
                    2 NaN 0.325 0.65    % Ceiling
    %                 4 6   0.05 0.3
                    5 NaN 0.325 0.325     % Walls
                    5 NaN 0.325 0.325     % Pillar1
                    5 NaN 0.325 0.325     % Pillar2
                    5 NaN 0.325 0.325     % Pillar3
                    5 NaN 0.325 0.325];   % Pillar4
            hold on;
            colormap(jet);
            
            % Initialise axis handles
            axhandles2 = cell(size(binDepths,1),1);
            for bb = 1:size(binDepths,1)
                % PLOT 3D MAP
                for half = 1:2
                    if ~isnan(axnum(bb,half))
                        if half == 1
                            ax = axes('Position',[0.05+axnum(bb,3)*pxSize+(jj-1)*pxSize 0.05+(plotgridv-ceil(pxcount/plotgridh))*pxSize+axnum(bb,4)*pxSize 0.1 0.1]);
                        else
                            ax = axes('Position',[0.05+(axnum(bb,3)+0.6)*pxSize+(jj-1)*pxSize 0.05+(plotgridv-ceil(pxcount/plotgridh))*pxSize+axnum(bb,4)*pxSize 0.1 0.1]);
                        end
                        axhandles2{bb} = ax;
                        % Plot rate maps
                        if ~isempty(map)
                            plotspatialview(bb,plotgridv,plotgridh,map{bb},ax,half);
                            
                            if ~isnan(maxC) && maxC ~= 0
                                set(ax,'CLim',[0 maxC],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
                            else
                                set(ax,'CLim',[0 1],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
                            end
                        end
                        
                        % Rotate view
                        view(ax,[120 70]);
                        if bb > 2
                            patch([peakcoords(1)-1 peakcoords(1)-1 peakcoords(1)+2 peakcoords(1)+2],[peakcoords(2)-2 peakcoords(2)+1 peakcoords(2)+1 peakcoords(2)-2], [0 0 0 0] ,[1 1 1 1],'FaceColor','r');
                        end
                    end
                end
            end
            % Plot outlines of shapes
            patchboundaries(axhandles2);
            % Go on to next pixel of main object
            pxcount = pxcount + 1;
        end
    end
    
%     % Plot single map of secondary object representing all pixels of main object
%     fignum = fignum + 1;
%     h3 = figure(fignum);
%     figname = horzcat(objtype,': ',num2str(setsessions(cell_ind)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', SI=',num2str(obj2SI,2),'/',num2str(obj2shprc,2),', field',num2str(fieldnum));
%     set(h3,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
%     ax = gca;
%     set(ax,'Visible','Off');
%     % Plot
%     axhandles3 = cell(size(binDepths,1),1);
%     map2Lavg = sum(map2L,2)/size(map2L,2);
%     % Convert from linear to grid map
%     map2LavgG = nan(obj.data.Args.GridSteps);
%     %%%%%%%%%%%%%%%%%%???
%     for bb = 1:size(binDepths,1)        
%         % PLOT 3D MAP
%         for half = 1:2
%             if ~isnan(axnum(bb,half))
%                 if half == 1
%                     ax = axes('Position',[axnum(bb,3) axnum(bb,4) 0.35 0.4]);
%                 else
%                     ax = axes('Position',[axnum(bb,3)+0.55 axnum(bb,4) 0.35 0.4]);
%                 end
%                 axhandles3{bb} = ax;
%                 % Plot rate maps
%                 plotspatialview(bb,plotgridv,plotgridh,map2L{bb},ax,half);
%                 if ~isnan(maxC) && maxC ~= 0
%                     set(ax,'CLim',[0 maxC],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
%                 else
%                     set(ax,'CLim',[0 1],'Visible','Off','ZLim',[0 15],'XLim',[0 40],'YLim',[0 40]);
%                 end
%                 % Rotate view
%                 view(ax,[120 70]);
%             end
%         end
%     end
%     % Plot outlines of shapes
%     patchboundaries(axhandles2);
end

cd(cwd);





function [surfx,surfy,surfz,surfmap] = plotspatialview(bb,plotgridv,plotgridh,map,ax,half)

% Plot maps
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
        shading flat;
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
        shading flat;
        hold(ax,'off');
    case 3 % Ground
        hold(ax,'on');
        surfx = repmat((0:40)',1,41);
        surfy = repmat(0:40,41,1);
        surfz = zeros(41);
        surfmap = map;
        
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        hold(ax,'off');
    case 4 % Ceiling
        hold(ax,'on');
        surfx = repmat((0:40)',1,41);
        surfy = repmat(0:40,41,1);
        surfz = zeros(41);
        surfmap = map;
        
        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        hold(ax,'off');
    case 5 % Walls
        hold(ax,'on');
%         if half == 1
            % Left 
            surfx = repmat(fliplr(0:40),9,1);
            surfy = zeros(9,41);
            surfz = repmat((0:8)',1,41);
            surfmap = flipud(map(:,1:40));

            surf(surfx,surfy,surfz,surfmap);
            shading flat;
            
            % Top
            surfx = zeros(9,41);
            surfy = repmat((0:40),9,1);
            surfz = repmat((0:8)',1,41);
            surfmap = flipud(map(:,41:80));

            surf(surfx,surfy,surfz,surfmap);
            shading flat;
            
%             hold(ax,'off');
%         elseif half == 2
            % Right
            surfx = repmat((0:40),9,1);
            surfy = repmat(40,9,41);
            surfz = repmat((0:8)',1,41);
            surfmap = flipud(map(:,81:120));

            surf(surfx,surfy,surfz,surfmap);
            shading flat;
            
            % Bottom
            surfx = repmat(40,9,41);
            surfy = repmat(fliplr(0:40),9,1);
            surfz = repmat((0:8)',1,41);
            surfmap = flipud(map(:,121:160));

            surf(surfx,surfy,surfz,surfmap);
            shading flat;
            
%         end
        hold(ax,'off');
        
    case 6 % Pillar1 (Bottom Right)
        hold(ax,'on');
        % Left
        surfx = repmat(fliplr(24:32),6,1);
        surfy = repmat(24,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,1:8));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Top
        surfx = repmat(24,6,9);
        surfy = repmat((24:32),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,9:16));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Right
        surfx = repmat((24:32),6,1);
        surfy = repmat(32,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,17:24));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Bottom
        surfx = repmat(32,6,9);
        surfy = repmat(fliplr(24:32),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,25:32));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        hold(ax,'off');

    case 7 % Pillar 2 (Bottom Left)
        hold(ax,'on');
        % Left
        surfx = repmat(fliplr(24:32),6,1);
        surfy = repmat(8,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,1:8));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;

        % Top
        surfx = repmat(24,6,9);
        surfy = repmat((8:16),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,9:16));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Right
        surfx = repmat((24:32),6,1);
        surfy = repmat(16,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,17:24));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Bottom
        surfx = repmat(32,6,9);
        surfy = repmat(fliplr(8:16),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,25:32));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        hold(ax,'off');
        
    case 8 % Pillar3 (Top Right)
        hold(ax,'on');
        % Left
        surfx = repmat(fliplr(8:16),6,1);
        surfy = repmat(24,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,1:8));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Top
        surfx = repmat(8,6,9);
        surfy = repmat((24:32),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,9:16));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Right
        surfx = repmat((8:16),6,1);
        surfy = repmat(32,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,17:24));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Bottom
        surfx = repmat(16,6,9);
        surfy = repmat(fliplr(24:32),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,24:32));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        hold(ax,'off');
    case 9 % Pillar4 (Top Left)
        hold(ax,'on');
        % Left
        surfx = repmat(fliplr(8:16),6,1);
        surfy = repmat(8,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,1:8));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Top
        surfx = repmat(8,6,9);
        surfy = repmat((8:16),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,9:16));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Right
        surfx = repmat((8:16),6,1);
        surfy = repmat(16,6,9);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,17:24));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        % Bottom
        surfx = repmat(16,6,9);
        surfy = repmat(fliplr(8:16),6,1);
        surfz = repmat((0:5)',1,9);
        surfmap = flipud(map(:,25:32));

        surf(surfx,surfy,surfz,surfmap);
        shading flat;
        
        hold(ax,'off');
end


function [] = patchboundaries(axhandles)

for ii = 1:size(axhandles,1)
    
    % Get axis of subplot
    axes(axhandles{ii});
    hold(axhandles{ii},'on');
    switch ii
        case 1 % Cue
            
        case 2  % Hint
            
        case 3  % Ground
            
            % Outer boundaries
            patch([0 0 40 40],[0 40 40 0], [0 0 0 0] ,[1 1 1 1],'FaceColor','none','LineWidth',1);
            % Pillars
            patch([24 24 32 32],[8 16 16 8], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Bottom left
            patch([8 8 16 16],[8 16 16 8], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Top left
            patch([8 8 16 16],[24 32 32 24], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Top right
            patch([24 24 32 32 ],[24 32 32 24], [0 0 0 0] ,[1 1 1 1],'FaceColor','none'); % Bottom right
            
        case 4  % Ceiling
            
            % Outer boundaries
            patch([0 0 40 40],[0 40 40 0], [0 0 0 0] ,[1 1 1 1],'FaceColor','none');
            
        case 5  % Walls
            
            % Outer boundaries
            patch([0 0 40 40],[0 0 0 0], [0 8 8 0] ,[1 1 1 1],'FaceColor','none');
            patch([0 0 0 0],[0 0 40 40], [0 8 8 0] ,[1 1 1 1],'FaceColor','none');
            patch([0 0 40 40], [40 40 40 40], [0 8 8 0] ,[1 1 1 1],'FaceColor','none');
            patch([40 40 40 40],[0 0 40 40], [0 8 8 0] ,[1 1 1 1],'FaceColor','none');
            
        case 6  % Pillar1
            
            % Outer boundaries
            patch([24 24 32 32],[24 24 24 24], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 24 24],[24 24 32 32], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 32 32],[32 32 32 32], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([32 32 32 32],[24 24 32 32], [0    5 5 0] ,[1 1 1 1],'FaceColor','none');

            % Outline Rabit poster on m_wall_25
            patch([26.88 26.88 29.12 29.12],[32 32 32 32], [1.8 3.2 3.2 1.8] ,[1 1 1 1],'FaceColor','none');    
            
        case 7  % Pillar2
            
            % Outer boundaries
            patch([24 24 32 32],[8 8 8 8], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 24 24],[8 8 16 16], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([24 24 32 32],[16 16 16 16], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([32 32 32 32],[8 8 16 16], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
        
            % Outline Cat poster on m_wall_10
            patch([32 32 32 32],[10.88 10.88 13.12 13.12], [1.8 3.2 3.2 1.8] ,[1 1 1 1],'FaceColor','none');

            % Outline Pig poster on m_wall_29
            patch([24 24 24 24],[10.88 10.88 13.12 13.12], [1.8 3.2 3.2 1.8] ,[1 1 1 1],'FaceColor','none');
            
        case 8  % Pillar3
            
            % Outer boundaries
            patch([8 8 16 16],[24 24 24 24], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 8 8],[24 24 32 32], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 16 16],[32 32 32 32], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([16 16 16 16],[24 24 32 32], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');

            % Outline Croc poster on m_wall_4
            patch([16 16 16 16],[26.88 26.88 29.12 29.12], [1.8 3.2 3.2 1.8] ,[1 1 1 1],'FaceColor','none');

            % Outline Donkey poster on m_wall_15
            patch([8 8 8 8],[26.88 26.88 29.12 29.12], [1.8 3.2 3.2 1.8] ,[1 1 1 1],'FaceColor','none');
            
        case 9  % Pillar4
            
            % Outer boundaries
            patch([8 8 16 16],[8 8 8 8], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 8 8],[8 8 16 16], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([8 8 16 16],[16 16 16 16], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');
            patch([16 16 16 16],[8 8 16 16], [0 5 5 0] ,[1 1 1 1],'FaceColor','none');

            % Outline Camel poster on m_wall_20
            patch([10.88 10.88 13.12 13.12],[8 8 8 8], [1.8 3.2 3.2 1.8] ,[1 1 1 1],'FaceColor','none');
    end
    
    
end


