function [] = plotratemaps(objtype,criteria,save,maptype,varargin)

% Function to plot rate maps for place, spatialview and heading direction
% Without additional inputs, it will work from and save in current
% directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% objtype: 'place', 'spatialview', or 'direction'
% 
% criteria: 'sic' or 'ise'
% 
% save: 1 (save figures) or 0.
% NOTE: if using save, figure folder location must be hard-coded below.
% 
% maptype: Type of smoothing, 'raw' (unsmoothed), 'boxcar'
% (boxcar-smoothed) or 'adaptive' (adaptive-smooted) NOTE: Currently only
% raw and adaptive smoothing have been tested.
% 
% Varargin: If only need to plot 1 cell's maps
% 1st: cell directory
%
% Note: replace 'criteria' with 'sortsic' if need to sort cells by SIC/ISE
% sortsic: 1 (sort the plotting by descending SIC) or 0 (sort the plotting
% by date) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sortsic = 0;
video = 0;
filttype = 'FiltVel';
pix = 1;

cwd = '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis';

if nargin > 4 % If plotting a single cell and map is already given as input (only used in placebyspatialview.m) 
    % Load cell list
    cellList = varargin(1);
%     mapGrid = varargin{2};
%     ax = varargin{3};
%     mapLin = varargin{4};
%     binDepths = varargin{5};
%     fieldCount = varargin{6};
%     fieldCoords = varargin{7};
    cd([cellList{1} '/' filttype '/' num2str(pix) 'px/']);
    % Load object
    switch objtype
        case 'place'
            objMain = load('vmpc.mat');
            objMain = objMain.vmp;
        case 'spatialview'
            objMain = load('vmsv.mat');
            objMain = objMain.vms;
        case 'direction'
            objMain = load('vmhd.mat');
            objMain = objMain.vmd;
    end
    if save
        %%% FILL IN
    end
elseif nargin > 4 % If plotting a single cell and map needs to be loaded
    cellList = varargin(1);
    cd(cellList{1});
    % Load object
    switch objtype
        case 'place'
            objMain = load('vmpc.mat');
            objMain = objMain.vmp;
        case 'spatialview'
            objMain = load('vmsv.mat');
            objMain = objMain.vms;
        case 'direction'
            objMain = load('vmhd.mat');
            objMain = objMain.vmd;
    end
    if save
        %%% FILL IN
    end
else % If plotting a batch of cells
    if save
        figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/' objtype '_' criteria '_' num2str(maptype)];
        if exist(figdir,'dir') ~= 7
            mkdir(figdir);
        else
            rmdir(figdir,'s');
            mkdir(figdir);
        end
    end
    % Load cell list
    cd(cwd);
    fid = fopen([cwd '/cell_list.txt'],'rt');
    cellList = textscan(fid,'%s','Delimiter','\n');
    cellList = cellList{1};
    % Make sure no empty cells
    notempty = ~cellfun(@isempty,cellList);
    cellList = cellList(notempty,:);
    
    % Load object
    c_objdir = [cwd '/Combined Objects/' filttype '/' num2str(pix) 'px/'];
    if exist(c_objdir,'dir') ~= 7
        disp('Missing combined object');
    else
        cd(c_objdir);
    end
    
    switch objtype
        case 'place'
            objMain = load('c_vmpc.mat');
            objMain = objMain.vmp;
        case 'spatialview'
            objMain = load('c_vmsv.mat');
            objMain = objMain.vms;
        case 'direction'
    end
    cd(cwd);
    if size(cellList,1) ~= size(objMain.data.SIC,1)
        disp('Object has different number of cells than CellList');
    end
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
        if strcmp(objtype,'Direction') % For now, 1st/2nd half plots only verified for Direction
            maps1 = objMain.data.maps_raw1;
            maps2 = objMain.data.maps_raw2;
        end
    case 'boxcar'
        maps = objMain.data.maps_boxsmooth;
    case 'adaptive'
        switch objtype
            case 'place'
                maps = objMain.data.maps_adsm;
            case 'spatialview'
                maps = objMain.data.maps_adsm;
        end
        if strcmp(objtype,'Direction') % For now, 1st/2nd half plots only verified for Direction
            maps1 = objMain.data.maps_adsm1;
            maps2 = objMain.data.maps_adsm2;
        end
end

% Load nptdata objMainect 
nCells = objMain.data.numSets;
% Plot params
plotgridh = 3;
plotgridv = 3;

if strcmp(objtype, 'place')  % Place maps
    % For each session, plot rate maps for each cell
    fig = 1;
    subpnum = 1;
    % Set up cell counts
    crossallthresh = 0;
    crosscellthresh = 0;
    crosspopthresh = 0;
    crosseitherthresh = 0;
    
    % Get shuffled SIC threshold for all cells
    switch criteria
        case 'sic'
            thr_sh = [objMain.data.SIC; objMain.data.SICsh];
        case 'ise'
            thr_sh = [objMain.data.ISE; objMain.data.ISEsh];
    end
%     thr_sh = reshape(thr_sh,numel(thr_sh),1);
    thr_pop = prctile(thr_sh,95);
    z_pop = zscore(thr_sh);
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));

        % Sort cells by descending SIC if required
        switch criteria
            case 'sic'
                thr_batch = objMain.data.SIC(cells_ind);
            case 'ise'
                thr_batch = objMain.data.ISE(cells_ind);
        end
        [~,thri] = sort(thr_batch,'descend');
        for jj = 1:length(cells_ind) % For each cell
            % Get cell index
            if sortsic
                cell_ind = cells_ind(thri(jj));
            else 
                cell_ind = cells_ind(jj);
            end

            %% Plot 1 map for 1 cell

            % Find figure number
            if jj*3 > plotgridh * plotgridv && mod((jj*3), (plotgridh * plotgridv)) == 3
                % Save figure
                if save
                    cwd = pwd;
                    cd(figdir);
                    % Save previous figure
                    figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                    saveas(h,figtitle,'png');
                    print('-painters',figtitle,'-dsvg');
                    cd(cwd);
                    close(figure(fig));
                end
                
                fig = fig + 1;
                subpnum = 1;
            end

            % Get shuffled SI cutoff for this cell - 95th percentile
            switch criteria
                case 'sic'
                    crit = objMain.data.SIC(cell_ind,1);
                    thr_cell = prctile(objMain.data.SICsh( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles ,1 ) ,95);
                    z_cell = z_pop(cell_ind,1);
                case 'ise'
                    crit = objMain.data.ISEsh(cell_ind,1);
                    thr_cell = prctile(objMain.data.ISEsh( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles ,1 ),95);
                    z_cell = z_pop(cell_ind,1);
            end
            
            % Get map
            if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
                mapLin = maps(cell_ind,:);

                h = figure(fig);
                ax = subplot(plotgridv,plotgridh,subpnum);
            end

            % Setup main object
            h = gcf;
            hold on;
            figname = horzcat(objtype,': ',num2str(setsessions(ii)));
            set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);

            % Plot main object
            [mapGrid,~]= plotmap(mapLin,objtype);
            
            % Set up axes
            if ~isnan(nanmax(mapLin(3:end))) && nanmax(mapLin(3:end)) ~= 0
                maxC = nanmax(mapLin(3:end));
            else
                maxC = 1;
            end
            set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
            ax.Title.String = horzcat(num2str(setsessions(ii)), 'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' ', criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2),', ','z',num2str(z_cell),', ',horzcat(num2str(nanmax(mapLin),3),'Hz'));
            
            if crit >= thr_cell && crit >= thr_pop
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
            elseif crit >= thr_cell && crit < thr_pop
                ax.Title.Color = 'm';
                crosscellthresh = crosscellthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            elseif crit < thr_cell && crit >= thr_pop
                ax.Title.Color = 'b';
                crosspopthresh = crosspopthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            else
                ax.Title.Color = 'k';
            end

            % Patch standing point if placebyspatialview
            if nargin > 5 % If placebyspatialview
                    fieldCoords
                    patch([fieldCoords(1)-2 fieldCoords(1)-2 fieldCoords(1)+1 fieldCoords(1)+1],[fieldCoords(2)-2 fieldCoords(2)+1 fieldCoords(2)+1 fieldCoords(2)-2], [0 0 0 0] ,[1 1 1 1],'FaceColor','none');
            end
            
            % Intra-session correlation %%%%%% NOTE: Should use boxcar
            % smoothed map
            switch maptype
                case 'adaptive'
                    map1 = objMain.data.maps_adsm1(cell_ind,:);
                    map2 = objMain.data.maps_adsm2(cell_ind,:);
                case 'raw'
                    map1 = objMain.data.maps_raw1(cell_ind,:);
                    map2 = objMain.data.maps_raw2(cell_ind,:);
            end
            vis1 = ~isnan(map1);
            vis2 = ~isnan(map2);
            vis = vis1 & vis2; % Correlate only visited bins;
            intracorr = corr2(map1(vis), map2(vis));
            switch criteria
                case 'sic'
                    crit1 = objMain.data.SIC1(cell_ind);
                    crit2 = objMain.data.SIC2(cell_ind);
                case 'ise'
                    crit1 = objMain.data.ISE1(cell_ind);
                    crit2 = objMain.data.ISE2(cell_ind);
            end

            % Plot
            subpnum = subpnum + 1;
            for kk = 1:2
                
                % Get map
                if kk == 1
                    mapLin = map1;
                    crit = crit1;
                    half = '1st';
                else
                    mapLin = map2;
                    crit = crit2;
                    half = '2nd';
                end
                % Restructure bins from linear to square grid
                mapGrid = flipud(reshape(mapLin, 40, 40)');

                % Setup object
                h = gcf;
                ax = subplot(plotgridv,plotgridh,subpnum);
                hold on;
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                
                % Plot map
                [~,~]= plotmap(mapLin,objtype);
                
                % Set up axes
                if ~isnan(nanmax(mapLin(3:end))) && nanmax(mapLin(3:end)) ~= 0
                    maxC = nanmax(mapLin(3:end));
                else
                    maxC = 1;
                end
                set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                ax.Title.String = horzcat(half,' half: ','corr=',num2str(intracorr,2),' ', criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2),', ',horzcat(num2str(nanmax(mapLin),3),'Hz'));
                
                if crit >= thr_cell && crit >= thr_pop
                    ax.Title.Color = 'r';
                elseif crit >= thr_cell && crit < thr_pop
                    ax.Title.Color = 'm';
                elseif crit < thr_cell && crit >= thr_pop
                    ax.Title.Color = 'b';
                else
                    ax.Title.Color = 'k';
                end

                subpnum = subpnum + 1;
            end
            
            hold off;

                
        end
        if save
            cwd = pwd;
            cd(figdir);
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
            print('-painters',figtitle,'-dsvg');
            cd(cwd);
            close(figure(fig));
        end
        
        fig = fig + 1;
        subpnum = 1;
        
    end
    disp(['Cross all thresh = ', num2str(crossallthresh),' cells']);
    disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
    disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);
    disp(['Cross either thresh = ', num2str(crosseitherthresh),' cells']);
    
elseif strcmp(objtype, 'spatialview') 
    
    fig = 1;
    subpnum = 1;
    % Set up cell counts
    crossallthresh = 0;
    crosscellthresh = 0;
    crosspopthresh = 0;
    crosseitherthresh = 0;
    % Get shuffled SIC threshold for all cells
    switch criteria
        case 'sic'
            thr_sh = [objMain.data.SIC; objMain.data.SICsh];
        case 'ise'
            thr_sh = [objMain.data.ISE; objMain.data.ISEsh];
    end
%     thr_sh = reshape(thr_sh,numel(thr_sh),1);
    thr_pop = prctile(thr_sh,95);
    z_pop = zscore(thr_sh);
    for ii = 1:size(setsessions,1) % For each session
        
        cells_ind = find(identifiers(:,1) == setsessions(ii));
        % Sort cells by descending SIC if required
        switch criteria
            case 'sic'
                thr_batch = objMain.data.SIC(cells_ind);
            case 'ise'
                thr_batch = objMain.data.ISE(cells_ind);
        end
        [~,thri] = sort(thr_batch,'descend');
        for jj = 1:length(cells_ind) % For each cell
            
            % Get cell index
            if sortsic
                cell_ind = cells_ind(thri(jj));
            else 
                cell_ind = cells_ind(jj);
            end
            
            %% Plot 1 map for 1 cell

            % Find figure number
            if jj*3 > plotgridh * plotgridv && mod((jj*3), (plotgridh * plotgridv)) == 3
                % Save figure
                if save
                    cwd = pwd;
                    cd(figdir);
                    % Save previous figure
                    figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                    saveas(h,figtitle,'png');
                    print('-painters',figtitle,'-dsvg');
                    cd(cwd);
                    close(figure(fig));
                end
                
                fig = fig + 1;
                subpnum = 1;
            end

            % Get shuffled SI cutoff for this cell - 95th percentile
            switch criteria
                case 'sic'
                    crit = objMain.data.SIC(cell_ind,1);
                    thr_cell = prctile(objMain.data.SICsh( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles,1 ),95);
                    z_cell = z_pop(cell_ind,1);
                case 'ise'
                    crit = objMain.data.ISE(cell_ind,1);
                    thr_cell = prctile(objMain.data.ISEsh( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles,1 ),95);
                    z_cell = z_pop(cell_ind,1);
            end
            
            if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells)

                mapLin = maps(cell_ind,:);
                % Set up figure
                h = figure(fig);
                ax = subplot(plotgridv,plotgridh,subpnum);
                
            end
            
            mapLin = emptyinsidepillar(mapLin);
            
            % Set up figure
            h = gcf;
            hold on;
            ax = gca;
            
            % Plot map
            [mapGrid,~]= plotmap(mapLin,objtype);
            
            % Figure and axes properties
            figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', ', criteria, ': ',num2str(crit,2),' of ',num2str(thr_cell,2));
            set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
            
            % Set up axes
            if ~isnan(nanmax(mapLin)) && nanmax(mapLin) ~= 0
                maxC = nanmax(mapLin(3:end));
            else
                maxC = 1;
            end
            set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
            % Identify cells sensitive to cue or hint
            if nanmax(mapLin(1)) > maxC
                ax.Title.String = horzcat('Cue: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' ',criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2),', ','z',num2str(z_cell),', ',horzcat(num2str(nanmax(mapLin),3)),'Hz');
            elseif nanmax(mapLin(2)) > maxC
                ax.Title.String = horzcat('Hint: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' ',criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2),', ','z',num2str(z_cell),', ',horzcat(num2str(nanmax(mapLin),3)),'Hz');
            else
                ax.Title.String = horzcat(num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),' ',criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2),', ','z',num2str(z_cell),', ',horzcat(num2str(nanmax(mapLin),3)),'Hz');
            end
            % Patch environment boundaries
            patchenvbounds(objtype);
            
            % Denote if significant spatial information
            if crit >= thr_cell && crit >= thr_pop
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
            elseif crit >= thr_cell && crit < thr_pop
                ax.Title.Color = 'm';
                crosscellthresh = crosscellthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            elseif crit >= thr_pop && crit < thr_cell
                ax.Title.Color = 'b';
                crosspopthresh = crosspopthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
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
            
            % Intra-session correlation %%%%%% NOTE: Should use boxcar
            % smoothed map
            
            switch maptype
                case 'adaptive'
                    map1 = objMain.data.maps_adsm1(cell_ind,:);
                    map2 = objMain.data.maps_adsm2(cell_ind,:);
                case 'raw'
                    map1 = objMain.data.maps_raw1(cell_ind,:);
                    map2 = objMain.data.maps_raw2(cell_ind,:);
            end
            map1 = emptyinsidepillar(map1);
            vis1 = ~isnan(map1);
            map2 = emptyinsidepillar(map2);
            vis2 = ~isnan(map2);
            vis = vis1 & vis2; % Correlate only visited bins;
            intracorr = corr2(map1(vis), map2(vis));
            switch criteria
                case 'sic'
                    crit1 = objMain.data.SIC1(cell_ind);
                    crit2 = objMain.data.SIC2(cell_ind);
                case 'ise'
                    crit1 = objMain.data.ISE1(cell_ind);
                    crit2 = objMain.data.ISE2(cell_ind);
            end

            % Plot
            subpnum = subpnum + 1;
            for kk = 1:2
                
                % Get map
                if kk == 1
                    mapLin = map1;
                    crit = crit1;
                    half = '1st';
                else
                    mapLin = map2;
                    crit = crit2;
                    half = '2nd';
                end

                % Setup object
                h = gcf;
                ax = subplot(plotgridv,plotgridh,subpnum);
                hold on;
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                
                % Plot map
                [mapGrid,~]= plotmap(mapLin,objtype);
                patchenvbounds(objtype);
                
                % Set up axes
                if ~isnan(nanmax(mapLin)) && nanmax(mapLin) ~= 0
                    maxC = nanmax(mapLin(3:end));
                else
                    maxC = 1;
                end
                set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                ax.Title.String = horzcat(half,' half: ','corr=',num2str(intracorr,2),' ', criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2),', ',horzcat(num2str(nanmax(mapLin),3),'Hz'));
                
                if crit >= thr_cell && crit >= thr_pop
                    ax.Title.Color = 'r';
                elseif crit >= thr_cell && crit < thr_pop
                    ax.Title.Color = 'm';
                elseif crit < thr_cell && crit >= thr_pop
                    ax.Title.Color = 'b';
                else
                    ax.Title.Color = 'k';
                end

                subpnum = subpnum + 1;
            end
            
            hold off;
            
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
        
        subpnum = 1;
        
    end
    disp(['Cross all thresh = ', num2str(crossallthresh),' cells']);
    disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
    disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);
    disp(['Cross either thresh = ', num2str(crosseitherthresh),' cells']);
    
elseif strcmp(objtype,'direction')
    
    maps = maps';
    maps1 = maps1';
    maps2 = maps2';
   
    
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
        thr_batch = objMain.data.SIC(cells_ind);
        [~,thri] = sort(thr_batch,'descend');
        
        %% Plot full session map for each cell
        
        for jj = 1:length(cells_ind) % For each cell
            % Get cell index
            if sortsic
                cell_ind = cells_ind(thri(jj));
            else 
                cell_ind = cells_ind(jj);
            end
            
                % Plot 1 map for 1 cell
                
                % Find figure number
                if jj > plotgridh * plotgridv && mod(jj, (plotgridh * plotgridv)) == 1
                    % Save figure
                    if save
                        % Save previous figure
                        figtitle = ['Direction: ' 'Full ' num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                        saveas(h,figtitle,'png');
                        print('-painters',figtitle,'-dsvg');
                    end
                    fig = fig + 1;
                    subpnum = 1;
                end
                
                % Get shuffled SI cutoff for this cell - 95th percentile
                crit = objMain.data.SIC(cell_ind,1);
                thr_cell = prctile(objMain.data.SICsh( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles,1 ),95);
                % Get shuffled RV cutoff for this cell - 95th percentile
                RV = objMain.data.RV(cell_ind,1);
                RVthr = prctile(objMain.data.RVsh( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles,1 ),95);
                % Get map
                if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
                    mapLin = maps(:,cell_ind);
                    h = figure(fig);
                    ax = subplot(plotgridv,plotgridh,subpnum);
                else
                    h = gcf;
                end

                % Setup main object
                hold on;
                colormap(jet);
                figname = horzcat('Direction: ','Full ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)));
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
                ax.DataAspectRatioMode = 'manual';
                ax.DataAspectRatio = [1 1 1];
                set(ax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XColor','none','YColor','none','ZColor','none','GridLineStyle','none');
                ax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', RVthr = ',num2str(RVthr,2));
                
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
        if save
            % Save figure
            figtitle = ['Direction: ' 'Full ' num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
            print('-painters',figtitle,'-dsvg');
        end
        fig = fig + 1;
        subpnum = 1;

        %% Compute first and second half correlation
        
        for jj = 1:length(cells_ind) % For each cell
            % Get cell index
            if sortsic
                cell_ind = cells_ind(thri(jj));
            else 
                cell_ind = cells_ind(jj);
            end

            % Find figure number
            if jj*2 > plotgridh * plotgridv && mod((jj*2), (plotgridh * plotgridv)) == 2
                % Save figure
                if save
                    % Save previous figure
                    figtitle = ['Direction: ' 'Half ' num2str(setsessions(ii)) '-' num2str(floor((jj*2)/(plotgridh * plotgridv)))];
                    saveas(h,figtitle,'png');
                    print('-painters',figtitle,'-dsvg');
                end
                fig = fig + 1;
                subpnum = 1;
            end

            %%%%%%% CORRELATION %%%%%%%%%%%

            % Intra-session correlation %%%%%% NOTE: Should use boxcar
            % smoothed map
            
            map1 = maps1(:,cell_ind);
            vis1 = ~isnan(map1);
            crit1 = objMain.data.SIC1;
            RV1 = objMain.data.RV1;
            map2 = maps2(:,cell_ind);
            vis2 = ~isnan(map2);
            crit2 = objMain.data.SIC2;
            RV2 = objMain.data.RV2;
            vis = vis1 & vis2; % Correlate only visited bins;
            intracorr = corr2(map1(vis), map2(vis));
            % Get shuffled SI cutoff for this cell - 95th percentile
            thr_cell = prctile(objMain.data.SICsh(2:end,cell_ind),95);
            RVthr = prctile(objMain.data.RVsh(2:end,cell_ind),95);

            % Plot

            for kk = 1:2
                
                % Get maps
                if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
                    if kk == 1
                        mapLin = map1;
                        RV = RV1;
                    else 
                        mapLin = map2;
                        RV = RV2;
                    end
                    h = figure(fig);
                    ax = subplot(plotgridv,plotgridh,subpnum);
                else
                    h = gcf;
                end

                % Setup main object
                hold on;
                colormap(jet);
                figname = horzcat('Direction: ','Half ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)));
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
                ax.DataAspectRatioMode = 'manual';
                ax.DataAspectRatio = [1 1 1];
                set(ax,'XTickLabelMode','manual','XTickLabel',{},'YTickLabelMode','manual','YTickLabel',{},'XColor','none','YColor','none','ZColor','none','GridLineStyle','none');
                ax.Title.String = horzcat('ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', RVthr = ',num2str(RVthr,2));
                
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
                RV_string = num2str(RV,'%1.3f');
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
                text('position', pos_val, 'units', 'normalized', 'HorizontalAlignment', hoz_val, 'string', RV_string, ...
                    'FontUnits', fu_val, 'VerticalAlignment', ver_val, 'fontsize', fs_val, 'color', text_color,'parent',ax);
                subpnum = subpnum + 1;
                hold off;
            end
        end
        
        if save
            % Save figure
            figtitle = ['Direction: ' 'Half ' num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
            print('-painters',figtitle,'-dsvg');
        end
        fig = fig + 1;
        subpnum = 1;
        
        
    end
    


end

cd(cwd);

function [surfx,surfy,surfz,surfmap] = plotspatialview(bb,plotgridv,plotgridh,map)
%%% DEFUNCT


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

%%% DEFUNCT

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

% Plot rate map (origin: placebyspatialview.m)
function [mapG,mapGdummy]= plotmap(mapL,objtype)

% Insert floor place map into larger 3D view setting
if strcmp(objtype,'place')
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
if strcmp(objtype,'spatialview')
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

disp(sum(sum(find(ceiling==0))) + sum(sum(find(floor==0))) + sum(sum(find(P4_TL==0))));

% Plot pillars
surf(P1_x, P1_y, PX_z, P1_BR);
alpha 1; shading flat;
surf(P2_x, P2_y, PX_z, P2_BL);
alpha 1; shading flat;
surf(P3_x, P3_y, PX_z, P3_TR);
alpha 1; shading flat;
surf(P4_x, P4_y, PX_z, P4_TL);
alpha 1; shading flat; 
view(-35,20);
colormap jet;
colorbar;

% Patch environment boundaries (origin: placebyspatialview.m)
function patchenvbounds(objtype)

switch objtype
    case 'place'
        
        % Floor outer bounds
        patch([0 0 40 40],[0 40 40 0],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        % Floor Pillar edges
        patch([8 8 16 16],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        
    case 'spatialview'
        
        % Floor outer bounds
        patch([0 0 40 40],[0 40 40 0],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        % Floor Pillar edges
        patch([8 8 16 16],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none');
        % Wall outer bounds
        patch([0 0 40 40],[0 0 0 0],[16 24 24 16],[1 1 1 1],'FaceColor','none');
        patch([0 0 0 0],[0 0 40 40],[16 24 24 16],[1 1 1 1],'FaceColor','none');
        patch([0 0 40 40], [40 40 40 40],[16 24 24 16],[1 1 1 1],'FaceColor','none');
        patch([40 40 40 40],[0 0 40 40],[16 24 24 16],[1 1 1 1],'FaceColor','none');
        % Pillar 1 (KW 3, top right)
        patch([24 24 32 32],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 24 24],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([26.88 26.88 29.12 29.12],[32 32 32 32],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Rabbit Poster on m_wall_25
        % Pillar 2 (KW 1, bottom right)
        patch([24 24 32 32],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 24 24],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[10.88 10.88 13.12 13.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Cat poster on m_wall_10
        patch([24 24 24 24],[10.88 10.88 13.12 13.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Pig poster on m_wall_29
        % Pillar 3 (KW 4, top left)
        patch([8 8 16 16],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 8 8],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[26.88 26.88 29.12 29.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Croc poster on m_wall_4
        patch([8 8 8 8],[26.88 26.88 29.12 29.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Donkey poster on m_wall_15
        % Pillar 4 (KW 2, bottom left)
        patch([8 8 16 16],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 8 8],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([10.88 10.88 13.12 13.12],[8 8 8 8],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Camel poster on m_wall_20
        % Ceiling
        patch([0 0 40 40],[0 40 40 0],[40 40 40 40],[1 1 1 1],'FaceColor','none');
        
end


function smoothviewmap = emptyinsidepillar(smoothviewmap)

% Temporary relief only!!!
% Removes the data from the inside of pillars where there should be no data
% at all

pillarcoords = [331:338,... % Bottom left
                371:378,...
                411:418,...
                451:458,...
                491:498,...
                531:538,...
                571:578,...
                611:618,...
                347:354,... % Bottom right
                387:394,...
                427:434,...
                467:474,...
                507:514,...
                547:554,...
                587:594,...
                627:634,...
                971:978,... % Top left
                1011:1018,...
                1051:1058,...
                1091:1098,...
                1131:1138,...
                1171:1178,...
                1211:1218,...
                1251:1258,...
                987:994,... % Top right
                1027:1034,...
                1067:1074,...
                1107:1114,...
                1147:1154,...
                1187:1194,...
                1227:1234,...
                1267:1274
                ];

smoothviewmap(pillarcoords) = nan;

