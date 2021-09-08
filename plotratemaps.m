function [selcell_orig] = plotratemaps(objtype,criteria,save,maptype,varargin)

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
    % Load object corrected for independent influences of place/view
    objCorr = load('vmcorr.mat');
    objCorr = objCorr.vmcorr;
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
    % Load object corrected for independent influences of place/view
    objCorr = load('vmcorr.mat');
    objCorr = objCorr.vmcorr;
    if save
        %%% FILL IN
    end
else % If plotting a batch of cells
    % Load cell list
    cd(cwd);
    fid = fopen([cwd '/cell_list_singlecell.txt'],'rt');
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
        case 'mixsel'
            objMain = load('c_vmms.mat');
            objMain = objMain.vmms;
            objPlace = load('c_vmpc.mat');
            objPlace = objPlace.vmp;
            objView = load('c_vmsv.mat');
            objView = objView.vms;
        case 'allprop'
            objMain0 = load('c_vmms0.mat');
            objMain0 = objMain0.vmms;
            objMain = load('c_vmms1.mat'); % Treat this as the main obj
            objMain = objMain.vmms;
            objPlace = load('c_vmpc.mat');
            objPlace = objPlace.vmp;
            objView = load('c_vmsv.mat');
            objView = objView.vms;
        case 'direction'
    end
    % Load object corrected for independent influences of place/view
    objCorr = load('c_vmcorr.mat');
    objCorr = objCorr.vmcorr;
    cd(cwd);
    % Check that combined object is the same size as cell list
    if size(cellList,1) ~= size(objMain.data.origin,1)
        disp('Object has different number of cells than CellList');
    end
    % Saving figure directory
    if save
        if strcmp(objtype,'mixsel')
            figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/' objtype '/UseCorr' num2str(objMain.data.Args.UseCorr)];
        else
            figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/' objtype '_' criteria '_' num2str(maptype)];
        end
        if exist(figdir,'dir') ~= 7
            mkdir(figdir);
        else
            rmdir(figdir,'s');
            mkdir(figdir);
        end
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

%%%% PATCH - c_vmcorr may have fewer cells than c_vmpc or c_vmsv. 
%%%% Match the cells in c_vmcorr to cells in other objects
if size(objCorr.data.origin,1) ~= size(objMain.data.origin,1)
    disp('Corrected object is not the same size as main object');
%     for ii = 1:size(objCorr.data.origin,1)
%         [~,matchcorr] = ismember(objCorr.data.origin{1},objMain.data.origin);
%     end
% else
%     matchcorr = (1:size(objCorr.data.origin,1))';
end

% Load firing rate maps depending on smoothing type
switch maptype
    case 'raw'
        mapname = 'maps_raw';
    case 'boxcar'
        mapname = 'maps_boxsmooth';
    case 'adaptive'
        mapname = 'maps_adsm';
end
if strcmp(objtype,'place') 
    maps = objMain.data.(mapname);
    mapscorr = objCorr.data.([mapname '_corrp']);
elseif strcmp(objtype,'spatialview')
    maps = objMain.data.(mapname);
    mapscorr = objCorr.data.([mapname '_corrsv']);
elseif strcmp(objtype,'mixsel')
    maps_p = objPlace.data.(mapname);
    maps_v = objView.data.(mapname);
elseif strcmp(objtype,'allprop')
    maps_p = objPlace.data.(mapname);
    maps_v = objView.data.(mapname);
end

% Load nptdata objMainect 
nCells = objMain.data.numSets;
% Plot params
plotgridh = 5;
plotgridv = 3;
% Set up cell counts
crossallthresh = 0;
crosscellthresh = 0;
crosspopthresh = 0;
crosseitherthresh = 0;
selcell_orig = [];
selcell_corr = [];

if strcmp(objtype, 'place')  % Place maps
    % For each session, plot rate maps for each cell
    fig = 1;
    subpnum = 1;
    
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
            if jj*5 > plotgridh * plotgridv && mod((jj*5), (plotgridh * plotgridv)) == 5
                % Save figure
                if save
                    cwd = pwd;
                    cd(figdir);
                    % Save previous figure
                    figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                    saveas(h,figtitle,'png');
%                     print('-painters',figtitle,'-dsvg');
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
                % if corrected map exists, get it
                if any(ismember(objMain.data.origin{cell_ind},objCorr.data.origin))
                    [~,corr_ind] = ismember(objMain.data.origin{cell_ind},objCorr.data.origin);
                    mapLincorr = mapscorr(corr_ind,:);
                else
                    mapLincorr = nan(size(mapLin));
                    corr_ind = [];
                end
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
            ax.Title.String = {horzcat(num2str(setsessions(ii)), 'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', ',num2str(nanmax(mapLin),3),'Hz'), horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z-',criteria,'= ',num2str(z_cell))};
            
            if crit >= thr_cell && crit >= thr_pop
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
%                 selcell_orig(end+1,1) = cell_ind;
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
            if crit >= thr_pop
                selcell_orig(end+1,1) = cell_ind;
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

                set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                ax.Title.String = {horzcat(half,' half: ','corr=',num2str(intracorr,2),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                
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
            
            % Plot corrected maps
            
            if isempty(corr_ind)
                critcorr = nan;
            else
                switch criteria
                    case 'sic'
                        switch objtype
                            case 'place'
                                critcorr = objCorr.data.SIC_corrp(corr_ind);
                            case 'spatialview'
                                critcorr = objCorr.data.SIC_corrsv(corr_ind);
                        end
                    case 'ise'
                        switch objtype
                            case 'place'
                                critcorr = objCorr.data.ISE_corrp(corr_ind);
                            case 'spatialview'
                                critcorr = objCorr.data.ISE_corrsv(corr_ind);
                        end
                end
            end
            
            % Setup object
            h = gcf;
            ax = subplot(plotgridv,plotgridh,subpnum);
            hold on;

            % Plot map
            [~,~]= plotmap(mapLincorr,objtype);
            maxCcorr = nanmax(mapLincorr);
            
            % Set up axes
            set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
            if ~isempty(corr_ind)
                ax.Title.String = {horzcat('Corrected',objCorr.data.llhpicklabel{corr_ind},'of',num2str(size(objCorr.data.llh{corr_ind},1)),': ',num2str(nanmax(mapLincorr),3),'Hz'),horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                set(ax,'CLim',[0 maxCcorr]);
            else
                ax.Title.String = 'No corrected map';
            end
            
            if critcorr >= thr_cell && critcorr >= thr_pop
                    ax.Title.Color = 'r';
%                     selcell_corr(end+1,1) = cell_ind;
            elseif critcorr >= thr_cell && critcorr < thr_pop
                ax.Title.Color = 'm';
            elseif critcorr < thr_cell && critcorr >= thr_pop
                ax.Title.Color = 'b';
            else
                ax.Title.Color = 'k';
            end
            hold off;
            if critcorr >= thr_pop
                selcell_corr(end+1,1) = cell_ind;
            end
            
            subpnum = subpnum + 1;
            
            % Plot covariance matrix
            
            % Setup object
            h = gcf;
            ax = subplot(plotgridv,plotgridh,subpnum);
            hold on;

            % Plot map
            im = imagesc(objCorr.data.covmat_norm{corr_ind});
            set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
            set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
            colormap jet;
            colorbar;
            
            % Replace NaNs with zeros in covariance matrix for norm calculations
            covmat = objCorr.data.covmat{corr_ind};
            covmat(isnan(covmat)) = 0;
            % Calculate norms
            norml1 = norm(covmat,1); % maximum of column sum
            norml2 = norm(covmat,2); % maximum single value
            
            % Set up axes
            set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
            if ~isempty(corr_ind)
                ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(norml1,2)), horzcat('l2=', num2str(norml2,2))};
            else
                ax.Title.String = 'No corrected map';
            end
            
            subpnum = subpnum + 1;
            
            hold off;

                
        end
        if save
            cwd = pwd;
            cd(figdir);
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
%             print('-painters',figtitle,'-dsvg');
            cd(cwd);
            close(figure(fig));
        end
        
        fig = fig + 1;
        subpnum = 1;
        
    end
    disp(['Cross all thresh, ', objtype, ' only = ', num2str(crossallthresh),' cells']);
    disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
    disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);
    disp(['Cross either thresh = ', num2str(crosseitherthresh),' cells']);
    disp(['Total number of cells = ',num2str(size(cellList,1)),' cells']);
    
elseif strcmp(objtype, 'spatialview') 
    
    fig = 1;
    subpnum = 1;

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
            if jj*5 > plotgridh * plotgridv && mod((jj*5), (plotgridh * plotgridv)) == 5
                % Save figure
                if save
                    cwd = pwd;
                    cd(figdir);
                    % Save previous figure
                    figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                    saveas(h,figtitle,'png');
%                     print('-painters',figtitle,'-dsvg');
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
                % if corrected map exists, get it
                if any(ismember(objMain.data.origin{cell_ind},objCorr.data.origin))
                    [~,corr_ind] = ismember(objMain.data.origin{cell_ind},objCorr.data.origin);
                    mapLincorr = mapscorr(corr_ind,:);
                else
                    mapLincorr = nan(size(mapLin));
                    corr_ind = [];
                end
                % Set up figure
                h = figure(fig);
                ax = subplot(plotgridv,plotgridh,subpnum);
                
            end
            
            mapLin = emptyinsidepillar(mapLin);
            mapLincorr = emptyinsidepillar(mapLincorr);
            
            % Set up figure
            h = gcf;
            hold on;
            ax = gca;
            
            % Plot map
            [mapGrid,~]= plotmap(mapLin,objtype);
            
            % Figure and axes properties
            figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)));
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
                ax.Title.String = {horzcat('Cue: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z',num2str(z_cell))};
            elseif nanmax(mapLin(2)) > maxC
                ax.Title.String = {horzcat('Hint: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z',num2str(z_cell))};
            else
                ax.Title.String = {horzcat(num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z',num2str(z_cell))};
            end
            % Patch environment boundaries
            patchenvbounds(objtype);
            
            % Denote if significant spatial information
            if crit >= thr_cell && crit >= thr_pop
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
%                 selcell_orig(end+1,1) = cell_ind;
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
            if crit >= thr_pop
                selcell_orig(end+1,1) = cell_ind;
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
%                 if ~isnan(nanmax(mapLin)) && nanmax(mapLin) ~= 0
%                     maxC = nanmax(mapLin(3:end));
%                 else
%                     maxC = 1;
%                 end
                set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                ax.Title.String = {horzcat(half,' half: ','corr=',num2str(intracorr,2),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                
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
            
            % Plot corrected maps
            
            if isempty(corr_ind)
                critcorr = nan;
            else
                switch criteria
                    case 'sic'
                        switch objtype
                            case 'place'
                                critcorr = objCorr.data.SIC_corrp(corr_ind);
                            case 'spatialview'
                                critcorr = objCorr.data.SIC_corrsv(corr_ind);
                        end
                    case 'ise'
                        switch objtype
                            case 'place'
                                critcorr = objCorr.data.ISE_corrp(corr_ind);
                            case 'spatialview'
                                critcorr = objCorr.data.ISE_corrsv(corr_ind);
                        end
                end
            end
            
            % Setup object
            h = gcf;
            ax = subplot(plotgridv,plotgridh,subpnum);
            hold on;

            % Plot map
            [~,~]= plotmap(mapLincorr,objtype);
            % Patch environment boundaries
            patchenvbounds(objtype);
            maxCcorr = nanmax(mapLincorr);
            
            % Set up axes
            set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
            if ~isempty(corr_ind)
                ax.Title.String = {horzcat('Corrected',objCorr.data.llhpicklabel{corr_ind},'of',num2str(size(objCorr.data.llh{corr_ind},1)),': ',num2str(nanmax(mapLincorr),3),'Hz'),horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                set(ax,'CLim',[0 maxCcorr]);
            else 
                ax.Title.String = 'No corrected map';
            end
            
            if critcorr >= thr_cell && critcorr >= thr_pop
                    ax.Title.Color = 'r';
%                     selcell_corr(end+1,1) = cell_ind;
            elseif critcorr >= thr_cell && critcorr < thr_pop
                ax.Title.Color = 'm';
            elseif critcorr < thr_cell && critcorr >= thr_pop
                ax.Title.Color = 'b';
            else
                ax.Title.Color = 'k';
            end
            if critcorr >= thr_pop
                selcell_corr(end+1,1) = cell_ind;
            end
            
            subpnum = subpnum + 1;
            
            % Plot covariance matrix
            
            % Setup object
            h = gcf;
            ax = subplot(plotgridv,plotgridh,subpnum);
            hold on;

            % Plot map
            im = imagesc(objCorr.data.covmat_norm{corr_ind});
            set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
%             set(ax,'CLim',[-1 1]);
            set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
            colormap jet;
            colorbar;
            
            % Replace NaNs with zeros in covariance matrix for norm calculations
            covmat = objCorr.data.covmat{corr_ind};
            covmat(isnan(covmat)) = 0;
            % Calculate norms
            norml1 = norm(covmat,1); % maximum of column sum
            norml2 = norm(covmat,2); % maximum single value
            
            % Set up axes
            set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
            if ~isempty(corr_ind)
                ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(norml1,2)), horzcat('l2=', num2str(norml2,2))};
            else
                ax.Title.String = 'No corrected map';
            end
            
            subpnum = subpnum + 1;
            
            hold off;
            
        end
        
        if save
            cwd = pwd;
            cd(figdir);
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_ind)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
%             print('-painters',figtitle,'-dsvg');
            cd(cwd);
            close(figure(fig));
        end
        
        fig = fig + 1;
        subpnum = 1;
        
    end
    disp(['Cross all thresh, ', objtype, ' only = ', num2str(crossallthresh),' cells']);
    disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
    disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);
    disp(['Cross either thresh = ', num2str(crosseitherthresh),' cells']);
    disp(['Total number of cells = ',num2str(size(cellList,1)),' cells']);
    
elseif strcmp(objtype,'mixsel')
    
%     plotgridh = 4;
    msobj = objMain.data.Args.msobj;
    
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_ind) % For each cell
            
            
            fig = 1;
            % Get cell index - unsorted
            cell_ind = cells_ind(jj);
            cell_id = [num2str(identifiers(cell_ind,1)) 'ch' num2str(identifiers(cell_ind,4)) 'c' num2str(identifiers(cell_ind,5))];
            cell_indP = strcmp(objPlace.data.origin,cellList{cell_ind});
            cell_indV = strcmp(objView.data.origin,cellList{cell_ind});
            cell_indCorr = strcmp(objCorr.data.origin,cellList{cell_ind});
            disp(['Plotting: ' cell_id]);
            
            %% Plot 1 map for 1 cell - mixed or single selective
            if objMain.data.placesel(cell_ind) || objMain.data.spatialviewsel(cell_ind)
                
                % Where to save plots
                if objMain.data.mixsel(cell_ind)
                    figdir2 = [figdir '/mixsel/' cell_id];
                    sel = 'mixed';
                elseif objMain.data.placesel(cell_ind)
                    figdir2 = [figdir '/placesel/' cell_id];
                    sel = 'place only';
                elseif objMain.data.spatialviewsel(cell_ind)
                    figdir2 = [figdir '/spatialviewsel/' cell_id];
                    sel = 'view only';
                end
                
                % Plot pv raw map vs pc/sv raw map vs corr
                h = figure(fig);
                hold on;
                figname = horzcat(objtype,': ',cell_id,'Raw maps');
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                for oo = 1:size(msobj,2)
                    % Plot pv maps
                    basedata = objMain.data.(msobj{oo})(cell_ind);
                    ax = subplot(2,3,1+(oo-1)*3);
                    plotmap(basedata.pvmap,msobj{oo});
                    % If firing rate or spike count = 0, set to black
                    settoblack(basedata.pvmap,msobj{oo});
                    % Scale color map
                    switch msobj{oo}
                        case 'place'
                            maxC = max(basedata.pvmap);
                            fieldcolor = 'r';
                        case 'spatialview'
                            maxC = max(basedata.pvmap(3:end)); % leave out cue and hint
                            fieldcolor = [50/256 205/256 50/256];
                    end
                    % Patch boundaries of base fields
                    for ff = 1:size(basedata.rate_components,1)
                        % Patch base fields
                        for pp = 1:size(basedata.fieldcoord{ff},1)
                            [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                        end
                    end
                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    colormap(ax,cool);
                    ax.Title.String = ['PV raw ',msobj{oo} ' map'];
                    % Plot maps from pc/sv objects
                    ax = subplot(2,3,2+(oo-1)*3);
                    if strcmp(msobj{oo},'place')
                        basemap = objPlace.data.maps_raw(cell_indP,:);
                        maxC = max(objPlace.data.maps_raw(cell_indP,:));
                        title = ['PC raw place map: ' sel];
                    elseif strcmp(msobj{oo},'spatialview')
                        basemap = objView.data.maps_raw(cell_indV,:);
                        maxC = max(objView.data.maps_raw(cell_indV,3:end));
                        title = ['SV raw view map: ' sel];
                    end
                    plotmap(basemap,msobj{oo});
                    % If firing rate or spike count = 0, set to black
                    settoblack(basemap,msobj{oo});
                    % Patch boundaries of base fields
                    for ff = 1:size(basedata.rate_components,1)
                        % Patch base fields
                        for pp = 1:size(basedata.fieldcoord{ff},1)
                            [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                        end
                    end
                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    colormap(ax,cool);
                    ax.Title.String = title;
                    % Plot maps from corr objects
                    ax = subplot(2,3,3+(oo-1)*3);
                    if strcmp(msobj{oo},'place')
                        basemap = objCorr.data.maps_raw_corrpset{cell_indCorr}(1,:);
                        maxC = max(basemap);
                        title = ['Corr raw place map: ' sel];
                    elseif strcmp(msobj{oo},'spatialview')
                        basemap = objCorr.data.maps_raw_corrsvset{cell_indCorr}(1,:);
                        maxC = max(basemap(3:end));
                        title = ['Corr raw view map: ' sel];
                    end
                    plotmap(basemap,msobj{oo});
                    % If firing rate or spike count = 0, set to black
                    settoblack(basemap,msobj{oo});
                    % Patch boundaries of base fields
                    for ff = 1:size(basedata.rate_components,1)
                        % Patch base fields
                        for pp = 1:size(basedata.fieldcoord{ff},1)
                            [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                        end
                    end
                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    colormap(ax,cool);
                    ax.Title.String = title;
                end
                % Save
                figtitle = horzcat(cell_id,': Raw maps: ', sel);
                if save
                    mkdir(figdir2);
                    savefigure(h,figtitle,figdir2);
                end
                close(h);
                fig = fig + 1;
                
                % Plot pc/sv vs corr smoothed maps - uncorrected
                h = figure(fig);
                hold on;
                figname = horzcat(objtype,': ',cell_id,'Smoothed maps');
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                for oo = 1:size(msobj,2)
                    basedata = objMain.data.(msobj{oo})(cell_ind);
                    % Plot maps from pc/sv objects
                    ax = subplot(2,2,1+(oo-1)*2);
                    if strcmp(msobj{oo},'place')
                        basemap = objPlace.data.maps_adsm(cell_indP,:);
                        maxC = max(basemap);
                        title = ['PC smoothed place map: ' sel];
                        fieldcolor = 'r';
                    elseif strcmp(msobj{oo},'spatialview')
                        basemap = objView.data.maps_adsm(cell_indV,:);
                        maxC = max(basemap(3:end));
                        title = ['SV smoothed view map: ' sel];
                        fieldcolor = [50/256 205/256 50/256];
                    end
                    plotmap(basemap,msobj{oo});
                    % If firing rate or spike count = 0, set to black
                    settoblack(basemap,msobj{oo});
                    % Patch boundaries of base fields
                    for ff = 1:size(basedata.rate_components,1)
                        % Patch base fields
                        for pp = 1:size(basedata.fieldcoord{ff},1)
                            [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                        end
                    end
                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    colormap(ax,cool);
                    ax.Title.String = title;
                    % Plot maps from corr objects
                    ax = subplot(2,2,2+(oo-1)*2);
                    if strcmp(msobj{oo},'place')
                        basemap = objCorr.data.maps_adsm_corrpset{cell_indCorr}(1,:);
                        maxC = max(basemap);
                        title = ['Corr smoothed place map: ' sel];
                    elseif strcmp(msobj{oo},'spatialview')
                        basemap = objCorr.data.maps_adsm_corrsvset{cell_indCorr}(1,:);
                        maxC = max(basemap(3:end));
                        title = ['Corr smoothed view map: ' sel];
                    end
                    plotmap(basemap,msobj{oo});
                    % If firing rate or spike count = 0, set to black
                    settoblack(basemap,msobj{oo});
                    % Patch boundaries of base fields
                    for ff = 1:size(basedata.rate_components,1)
                        % Patch base fields
                        for pp = 1:size(basedata.fieldcoord{ff},1)
                            [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                        end
                    end
                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    colormap(ax,cool);
                    ax.Title.String = title;
                end
                % Save
                figtitle = horzcat(cell_id,': Smoothed maps: ', sel);
                if save
                    savefigure(h,figtitle,figdir2);
                end
                close(h);
                fig = fig + 1;
                
                % Plot base maps for both objects and their fields
                h = figure(fig);
                hold on;
                figname = horzcat(objtype,': ','Base maps: ',sel);
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                for oo = 1:size(msobj,2)
                    basedata = objMain.data.(msobj{oo})(cell_ind);
                    ax = subplot(1,2,oo);
                    plotmap(basedata.basemapLsm,msobj{oo});
                    colormap(ax,cool);
                    % If firing rate or spike count = 0, set to black
                    settoblack(basedata.basemapLsm,msobj{oo});
                    % Scale color map
                    switch msobj{oo}
                        case 'place'
                            maxC = max(basedata.basemapLsm);
                            fieldcolor = 'r';
                        case 'spatialview'
                            maxC = max(basedata.basemapLsm(3:end));
                            fieldcolor = [50/256 205/256 50/256];
                    end
                    ax.Title.String = {['SI: ' num2str(basedata.SI,2)];['thr: ' num2str(basedata.SIthr,2)]};
                    ax.Title.FontSize = 16;
                    if basedata.SI>basedata.SIthr
                        ax.Title.Color = 'r';
                    end
                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    % Draw env boundaries
                    patchenvbounds(msobj{oo});
                    % Patch boundaries of base fields
                    for ff = 1:size(basedata.rate_components,1)
                        % Patch base fields
                        for pp = 1:size(basedata.fieldcoord{ff},1)
                            [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1);
                        end
                    end
                end
                % Save
                figtitle = horzcat(cell_id,': Base maps: ',sel);
                if save
                    savefigure(h,figtitle,figdir2);
                end
                close(h);
                fig = fig + 1;
                
                % Plot pixel maps for each base field
                for oo = 1:size(msobj,2)
                    
                    basedata = objMain.data.(msobj{oo})(cell_ind);
                    secdata = objMain.data.(msobj{2-oo+1})(cell_ind);
                    
                    for ff = 1:size(basedata.rate_components,1) % In mixselpv we've only extracted data for max of 3 fields per map
                            % Plot base maps with secondary pixel maps
                            h = figure(fig);
                            hold on;
                            ax = gca;
                            ax.Title.String = cell(size(basedata.secfieldrates{ff},1),1);
                            linked = [];
                            for ss = 1:size(basedata.secfieldrates{ff},1)
                                % Check if infield firing in sec field is higher than outfield firing
                                if basedata.secfieldrates{ff}(ss,1) > prctile(basedata.secfieldrates_sh{ff}{ss,1},95)
                                    linked = horzcat(linked,',',num2str(ss));
                                    ax.Title.String{ss,1} = horzcat('LINKED:',msobj{2-oo+1},' Field ',num2str(ss),' Infield: ', num2str(basedata.secfieldrates{ff}(ss,1),2),'Hz; Outfield: ', num2str(prctile(basedata.secfieldrates_sh{ff}{ss,1},95),2),'Hz');
                                else
                                    ax.Title.String{ss,1} = horzcat(msobj{2-oo+1},' Field ',num2str(ss),' Infield: ', num2str(basedata.secfieldrates{ff}(ss,1),2),'Hz; Outfield: ', num2str(prctile(basedata.secfieldrates_sh{ff}{ss,1},95),2),'Hz');
                                end
                            end
                            if isempty(linked)
                                figname = horzcat(msobj{oo},' field ',num2str(ff),' not linked');
                            else
                                figname = horzcat(msobj{oo},' field ',num2str(ff),' sig linked to ',msobj{2-oo+1},' field ', linked);
%                                 ax.Title.Color = 'r';
                            end
                            set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                            % Plot secondary full pixel set
                            tempsecmap = nan(size(basedata.secmapLsm));
                            tempsecmap(basedata.set_sec_linbin{ff,1}(:,1)) = basedata.set_sec_linbin{ff,1}(:,4);
                            plotmap(tempsecmap,msobj{2-oo+1});
                            colormap(ax,cool);
                            set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color','none');
                            patchenvbounds('spatialview');
                            % If firing rate or spike count = 0, set to black
                            settoblack(tempsecmap,msobj{2-oo+1});
                            % Patch basemap field
                            for pp = 1:size(basedata.fieldcoord{ff},1)
                                [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                if strcmp(msobj{oo},'place')
                                    patch(x,y,z,[1 1 1 1],'EdgeColor','r','FaceColor','none','LineWidth',2);
                                elseif strcmp(msobj{oo},'spatialview')
                                    patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',2); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                end
                            end
                            % Patch secondary fields
                            for mm = 1:size(secdata.rate_components,1) % For each secondary field
                                % Plot secondary field in consolidated pixel map
                                sec = secdata.fieldcoord{mm}; % Get bins for sec field
                                % Patch secondary fields in the consolidated pixel map
                                for pp = 1:size(sec,1)
                                    [x,y,z] = converttosurf(secdata.gridnum(mm),sec(pp,1),sec(pp,2));
                                    patch(x,y,z,'r','EdgeColor',[169/256 169/256 169/256],'FaceColor','none','LineWidth',2); % Dark Grey
                                    if strcmp(msobj{2-oo+1},'place')
                                        patch(x,y,z,'r','EdgeColor','r','FaceColor','none','LineWidth',2); % Red
                                    elseif strcmp(msobj{2-oo+1},'spatialview')
                                        patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',2); % Lime green
                                    end
                                end
                            end
                            % Save 
                            figtitle = figname;
                            if save
                                savefigure(h,figtitle,figdir2);
                            end
                            close(h);
                            fig = fig + 1;
                    end
                    
                end
                
            end
        end
    end
    
elseif strcmp(objtype,'allprop') % If plotting place, view, corr and mixsel all in one page
    
    msobj = objMain0.data.Args.msobj;
    
    for ii = 1:size(setsessions,1) % For each session
        cells_ind = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_ind) % For each cell
            
            cell_ind = cells_ind(jj); % ind within the given cell list
            cell_id = horzcat(num2str(identifiers(cell_ind,1)),'ch',num2str(identifiers(cell_ind,4)),...
                'c',num2str(identifiers(cell_ind,5)));
            disp(cell_id);
            % Get cell index within combined objects
            cell_indP = strcmp(objPlace.data.origin,cellList{cell_ind});
            cell_indV = strcmp(objView.data.origin,cellList{cell_ind});
            cell_indCorr = strcmp(objCorr.data.origin,cellList{cell_ind});
            cell_indMS0 = strcmp(objMain0.data.origin,cellList{cell_ind});
            cell_indMS = strcmp(objMain.data.origin,cellList{cell_ind});
            
            % What to save
            if objMain0.data.mixsel(cell_indMS0)
                presel = 'M';
            elseif objMain0.data.placesel(cell_indMS0)
                presel = 'P';
            elseif objMain0.data.spatialviewsel(cell_indMS0)
                presel = 'V';
            else
                presel = 'NS';
            end
            if objMain.data.mixsel(cell_indMS)
                postsel = 'M';
            elseif objMain.data.placesel(cell_indMS)
                postsel = 'P';
            elseif objMain.data.spatialviewsel(cell_indMS)
                postsel = 'V';
            else
                postsel = 'NS';
            end
            pcSIthr = prctile(objPlace.data.SICsh,95);
            svSIthr = prctile(objView.data.SICsh,95);
            
            % Set up plot page: Find out how many sig fields there are
            axnum_col = 8;
            axnum_row = 8;
            if strcmp(presel,'NS') && strcmp(postsel,'NS') % If not selective at all, just plot overall maps
                numfigs = 1;
            else % Plot pixel maps
                numfields = [objMain0.data.place(cell_indMS0).sigfields ...
                        objMain0.data.spatialview(cell_indMS0).sigfields ...
                        objMain.data.place(cell_indMS).sigfields ...
                        objMain.data.spatialview(cell_indMS).sigfields];
                if strcmp(presel,'M') || strcmp(postsel,'M')
                    numfigs = nanmax([numfields(1) numfields(3)]);
                else 
                    numfigs = 1;
                end
            end
            link_pre = {};
            link_post = {};
            
            % Plot 
            for kk = 1:numfigs
                % Set up figure
                h = figure(kk);
                hold on;
                figtitle = [cell_id ' ' presel '-' postsel ' ' num2str(kk)];
                set(h,'Name',figtitle,'Units','Normalized','Position',[0 0 1 1]);
                ax = gca;
                ax.Visible = 'off';
                axwidth = 0.8*(0.9/axnum_col);
                col_left = 0.05:0.9/axnum_col:0.95-0.9/axnum_col;
                axheight = 0.8*(0.9/axnum_row);
                row_bot = fliplr(0.05:0.9/axnum_row:0.95-0.9/axnum_row);

                % Plot place and view maps
                % Orig place
                ax = axes('Position',[col_left(1) row_bot(2) axwidth axheight*2]); % Full session
                map = objPlace.data.maps_adsm(cell_indP,:);
                plotmap(map,'place');
                patchenvbounds('place');
                colorbar off;
                maxC = max(map);
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objPlace.data.SIC(cell_indP),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objPlace.data.SIC(cell_indP) > pcSIthr % objMain0.data.placesel(cell_indMS0)
                    si.Color = 'r';
                end
                ax.Title.String = {'Orig Place'; [num2str(objMain0.data.place(cell_indMS0).sigfields) ' fields']};
                ax = axes('Position',[col_left(2) row_bot(1) axwidth axheight]); % 1st half
                map = objPlace.data.maps_adsm1(cell_indP,:);
                plotmap(map,'place');
                patchenvbounds('place');
                colorbar off;
                maxC = max(map);
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objPlace.data.SIC1(cell_indP),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objPlace.data.SIC1(cell_indP) > pcSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '1st half';
                ax = axes('Position',[col_left(2) row_bot(2) axwidth axheight]); % 2nd half
                map = objPlace.data.maps_adsm2(cell_indP,:);
                plotmap(map,'place');
                patchenvbounds('place');
                colorbar off;
                maxC = max(map);
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objPlace.data.SIC2(cell_indP),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objPlace.data.SIC2(cell_indP) > pcSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '2nd half';
                % Orig view
                ax = axes('Position',[col_left(3) row_bot(2) axwidth axheight*2]); % Full session
                map = objView.data.maps_adsm(cell_indV,:);
                plotmap(map,'spatialview');
                patchenvbounds('spatialview');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objView.data.SIC(cell_indV),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objView.data.SIC(cell_indV) > svSIthr
                    si.Color = 'r';
                end
                ax.Title.String = {'Orig View'; [num2str(objMain0.data.spatialview(cell_indMS0).sigfields) ' fields']};
                ax = axes('Position',[col_left(4) row_bot(1) axwidth axheight]); % 1st half
                map = objView.data.maps_adsm1(cell_indV,:);
                plotmap(map,'spatialview');
                patchenvbounds('spatialview');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objView.data.SIC1(cell_indV),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objView.data.SIC1(cell_indV) > svSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '1st half';
                ax = axes('Position',[col_left(4) row_bot(2) axwidth axheight]); % 2nd half
                map = objView.data.maps_adsm2(cell_indV,:);
                plotmap(map,'spatialview');
                patchenvbounds('spatialview');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objView.data.SIC2(cell_indV),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objView.data.SIC2(cell_indV) > svSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '2nd half';
                % Corr Place
                ax = axes('Position',[col_left(5) row_bot(2) axwidth axheight*2]); % Full session
                map = objCorr.data.maps_adsm_corrp(cell_indCorr,:);
                plotmap(map,'place');
                patchenvbounds('place');
                colorbar off;
                maxC = max(map);
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.SIC_corrp(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.SIC_corrp(cell_indCorr) > pcSIthr
                    si.Color = 'r';
                end
                ax.Title.String = {['Corr Place (' num2str(objCorr.data.llhpick(cell_indCorr)) ')']; ...
                    [num2str(objMain.data.place(cell_indMS).sigfields) ' fields']};
                % Corr view
                ax = axes('Position',[col_left(7) row_bot(2) axwidth axheight*2]); % Full session
                map = objCorr.data.maps_adsm_corrsv(cell_indCorr,:);
                plotmap(map,'spatialview');
                patchenvbounds('spatialview');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.SIC_corrsv(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.SIC_corrsv(cell_indCorr) > svSIthr
                    si.Color = 'r';
                end
                ax.Title.String = {['Corr View (' num2str(objCorr.data.llhpick(cell_indCorr)) ')'];...
                    [num2str(objMain.data.spatialview(cell_indMS).sigfields) ' fields']};
                
                % Plot orig pixel maps
                if ~isnan(objMain0.data.(msobj{1})(cell_indMS0).sigfields) && ...
                        ~isnan(objMain0.data.(msobj{2})(cell_indMS0).sigfields)
                    
                    for oo = 1:size(msobj,2)
                        basedata = objMain0.data.(msobj{oo})(cell_indMS0);
                        secdata = objMain0.data.(msobj{2-oo+1})(cell_indMS0);
                        if oo == 1
                            if basedata.sigfields == 0
                                continue;
                            elseif basedata.sigfields < kk % If num fields smaller than available axes
                                continue;
                            else
                                ffcount = kk;
                                fffcount = 1:secdata.sigfields;
                            end
                        elseif oo == 2
                            if basedata.sigfields == 0
                                continue;
%                             elseif basedata.sigfields < kk % If num fields smaller than available axes
%                                 continue;
                            else
                                ffcount = 1:basedata.sigfields;
                                if secdata.sigfields > 0
                                    fffcount = kk;
                                else
                                    fffcount = 1:secdata.sigfields; % 0;
                                end
                            end
                        end
                        for ff = ffcount % 1:basedata.sigfields
                            if ff > basedata.sigfields
                                continue;
                            end
                            % Plot base field
                            if oo == 1
                                ax = axes('Position',[col_left(2*(oo-1)+1) row_bot(4) axwidth axheight*2]);
                            elseif oo == 2
                                ax = axes('Position',[col_left(2*(oo-1)+1) row_bot(4+(ff-1)*2) axwidth axheight*2]);
                            end
                            map = nan(size(basedata.basemapLsm));
                            plotmap(map,msobj{oo});
                            patchenvbounds(msobj{oo});
                            axdisplay(ax,1);
                            colorbar off;
                            if strcmp(msobj{oo},'place')
                                fieldcolor = 'r';
                            elseif strcmp(msobj{oo},'spatialview')
                                fieldcolor = [50/256 205/256 50/256];
                            end
                            % Patch base fields
                            for pp = 1:size(basedata.fieldcoord{ff},1)
                                [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1);
                            end
                            ax.Title.String = [msobj{oo}(1) num2str(ff)];
                            % Plot sec pixels and fields
                            if isempty(fffcount)
                                % Plot sec pixels for single selectivecells
                                ax = axes('Position',[col_left(2*(oo-1)+2) row_bot(4+2*(ff-1)) axwidth axheight*2]);
                                map = nan(size(basedata.secmapLsm));
                                map(basedata.set_sec_linbin{ff,1}(:,1)) = basedata.set_sec_linbin{ff,1}(:,4);
                                maxC = nanmax(map(3:end));
                                plotmap(map,msobj{2-oo+1});
                                patchenvbounds(msobj{2-oo+1});
                                settoblack(map,msobj{2-oo+1});
                                axdisplay(ax,maxC);
                                colormap(ax,'cool');
                                c = colorbar;
                                c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4) 0.9/150 axheight*2];
                            else
                                for fff = fffcount
%                                     if fff == 0 
%                                         continue; % mostly for single selectivity to skip plotting pixels
                                    if fff > secdata.sigfields
                                        continue;
                                    end
                                    % Plot sec pixel maps
                                    if oo == 1
                                        ax = axes('Position',[col_left(2*(oo-1)+2) row_bot(4+2*(fff-1)) axwidth axheight*2]);
                                    elseif oo == 2 
                                        ax = axes('Position',[col_left(2*(oo-1)+2) row_bot(4+2*(ff-1)) axwidth axheight*2]); % ???? 
                                    end
                                    map = nan(size(basedata.secmapLsm));
                                    map(basedata.set_sec_linbin{ff,1}(:,1)) = basedata.set_sec_linbin{ff,1}(:,4);
                                    maxC = nanmax(map(3:end));
                                    plotmap(map,msobj{2-oo+1});
                                    patchenvbounds(msobj{2-oo+1});
                                    settoblack(map,msobj{2-oo+1});
                                    axdisplay(ax,maxC);
                                    colormap(ax,'cool');
                                    c = colorbar;
                                    if oo == 1
                                        c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(fff-1)) 0.9/150 axheight*2];
                                    elseif oo == 2
                                        c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
                                    end
                                    if strcmp(msobj{oo},'spatialview')
                                        fieldcolor = 'r';
                                    elseif strcmp(msobj{oo},'place')
                                        fieldcolor = [50/256 205/256 50/256]; % Lime green
                                    end
                                    % Patch secondary fields
                                    sec = secdata.fieldcoord{fff}; % Get bins for sec field
                                    for pp = 1:size(sec,1)
                                        [x,y,z] = converttosurf(secdata.gridnum(fff),sec(pp,1),sec(pp,2));
                                        patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                    end
                                    % Check if fields are reciprocally linked
                                    ax.Title.String = {horzcat(msobj{2-oo+1}(1),num2str(fff),' Infield: ', num2str(basedata.secfieldrates{ff}(fff,1),2),'Hz'); ...
                                        horzcat('Outfield: ', num2str(prctile(basedata.secfieldrates_sh{ff}{fff,1},95),2),'Hz')};
                                    if basedata.secfieldrates{ff}(fff,1) > prctile(basedata.secfieldrates_sh{ff}{fff,1},95)
                                        if secdata.secfieldrates{fff}(ff,1) > prctile(secdata.secfieldrates_sh{fff}{ff,1},95) % If reciprocal
                                            ax.Title.Color = 'r';
                                            if oo == 1
                                                link_pre{end+1,1} = [msobj{oo}(1) num2str(ff) 'RL' msobj{2-oo+1}(1) num2str(fff)];
                                            end
                                        else 
                                            ax.Title.Color = 'b';
                                            link_pre{end+1,1} = [msobj{oo}(1) num2str(ff) 'L' msobj{2-oo+1}(1) num2str(fff)];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                % Plot corr pixel maps
                if ~isnan(objMain.data.(msobj{1})(cell_indMS).sigfields) && ~isnan(objMain.data.(msobj{2})(cell_indMS).sigfields)
                    for oo =1:size(msobj,2)
                        basedata = objMain.data.(msobj{oo})(cell_indMS);
                        secdata = objMain.data.(msobj{2-oo+1})(cell_indMS);
                        if oo == 1
                            if basedata.sigfields == 0
                                continue;
                            elseif basedata.sigfields < kk % If num fields smaller than available axes
                                continue;
                            else
                                ffcount = kk;
                                fffcount = 1:secdata.sigfields;
                            end
                        elseif oo == 2
                            if basedata.sigfields == 0
                                continue;
%                             elseif basedata.sigfields < kk % If num fields smaller than available axes
%                                 continue;
                            else
                                ffcount = 1:basedata.sigfields;
                                if secdata.sigfields > 0
                                    fffcount = kk;
                                else
                                    fffcount = 1:secdata.sigfields; % 0;
                                end
                            end
                        end
                        for ff = ffcount % 1:basedata.sigfields
                            if ff > basedata.sigfields
                                continue;
                            end
                            % Plot base fields
                            if oo == 1
                                ax = axes('Position',[col_left(2*(oo-1)+5) row_bot(4) axwidth axheight*2]);
                            elseif oo == 2
                                ax = axes('Position',[col_left(2*(oo-1)+5) row_bot(4+(ff-1)*2) axwidth axheight*2]);
                            end
                            map = nan(size(basedata.basemapLsm));
                            plotmap(map,msobj{oo});
                            patchenvbounds(msobj{oo});
                            colorbar off;
                            axdisplay(ax,1);
                            if strcmp(msobj{oo},'place')
                                fieldcolor = 'r';
                            elseif strcmp(msobj{oo},'spatialview')
                                fieldcolor = [50/256 205/256 50/256];
                            end
                            % Patch base fields
                            for pp = 1:size(basedata.fieldcoord{ff},1)
                                [x,y,z] = converttosurf(basedata.gridnum(ff),basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1);
                            end
                            ax.Title.String = [msobj{oo}(1) num2str(ff)];
                            if isempty(fffcount)
                                % Plot sec pixels for single selectivecells
                                ax = axes('Position',[col_left(2*(oo-1)+6) row_bot(4+2*(ff-1)) axwidth axheight*2]);
                                map = nan(size(basedata.secmapLsm));
                                map(basedata.set_sec_linbin{ff,1}(:,1)) = basedata.set_sec_linbin{ff,1}(:,4);
                                maxC = nanmax(map(3:end));
                                plotmap(map,msobj{2-oo+1});
                                patchenvbounds(msobj{2-oo+1});
                                settoblack(map,msobj{2-oo+1});
                                axdisplay(ax,maxC);
                                colormap(ax,'cool');
                                c = colorbar;
                                c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4) 0.9/150 axheight*2];
                            else 
                                for fff = fffcount
%                                     if fff == 0 
%                                         continue; % mostly for single selectivity to skip plotting pixels
                                    if fff > secdata.sigfields
                                        continue;
                                    end
                                    % Sec pixel maps
                                    if oo == 1
                                        ax = axes('Position',[col_left(2*(oo-1)+6) row_bot(4+2*(fff-1)) axwidth axheight*2]);
                                    elseif oo == 2 
                                        ax = axes('Position',[col_left(2*(oo-1)+6) row_bot(4+2*(ff-1)) axwidth axheight*2]);
                                    end
                                    map = nan(size(basedata.secmapLsm));
                                    map(basedata.set_sec_linbin{ff,1}(:,1)) = basedata.set_sec_linbin{ff,1}(:,4);
                                    maxC = nanmax(map(3:end));
                                    plotmap(map,msobj{2-oo+1});
                                    patchenvbounds(msobj{2-oo+1});
                                    settoblack(map,msobj{2-oo+1});
                                    axdisplay(ax,maxC);
                                    colormap(ax,'cool');
                                    c = colorbar;
                                    if oo == 1
                                        c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4+2*(fff-1)) 0.9/150 axheight*2];
                                    elseif oo == 2
                                        c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
                                    end
                                    if strcmp(msobj{oo},'spatialview')
                                        fieldcolor = 'r';
                                    elseif strcmp(msobj{oo},'place')
                                        fieldcolor = [50/256 205/256 50/256]; % Lime green
                                    end
                                    % Patch secondary fields
                                    sec = secdata.fieldcoord{fff}; % Get bins for sec field
                                    for pp = 1:size(sec,1)
                                        [x,y,z] = converttosurf(secdata.gridnum(fff),sec(pp,1),sec(pp,2));
                                        patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                    end
                                    % Check if fields are reciprocally linked
                                    ax.Title.String = {horzcat(msobj{2-oo+1}(1),num2str(fff),' Infield: ', num2str(basedata.secfieldrates{ff}(fff,1),2),'Hz'); ...
                                        horzcat('Outfield: ', num2str(prctile(basedata.secfieldrates_sh{ff}{fff,1},95),2),'Hz')};
                                    if basedata.secfieldrates{ff}(fff,1) > prctile(basedata.secfieldrates_sh{ff}{fff,1},95)
                                        if secdata.secfieldrates{fff}(ff,1) > prctile(secdata.secfieldrates_sh{fff}{ff,1},95) % If reciprocal
                                            ax.Title.Color = 'r';
                                            if oo == 1
                                                link_post{end+1,1} = [msobj{oo}(1) num2str(ff) 'RL' msobj{2-oo+1}(1) num2str(fff)];
                                            end
                                        else 
                                            ax.Title.Color = 'b';
                                            link_post{end+1,1} = [msobj{oo}(1) num2str(ff) 'L' msobj{2-oo+1}(1) num2str(fff)];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                %%% Patch - Figure already named above. Adding here the code for whether fields are linked
                if ~isempty(link_pre) || ~isempty(link_post)
                    if ~isempty(strfind(link_pre,'R')) || ~isempty(strfind(link_post,'R'))
                        figtitle = [cell_id ' ' presel '-' postsel '-RL' ' ' num2str(kk)];
                    else
                        figtitle = [cell_id ' ' presel '-' postsel '-L' ' ' num2str(kk)];
                    end
                    set(h,'Name',figtitle);
                end
                if save
                    savefigure(h,figtitle,figdir);
                end
                close(h);
            end
            if ~isempty(strfind(link_pre,'R'))
                selcell_orig(end+1,1) = cell_ind;
            end 
            if ~isempty(strfind(link_post,'R'))
                selcell_corr(end+1,1) = cell_ind;
            end
        end
    end
 
end

selcell_orig = cellList(selcell_orig);
selcell_corr = cellList(selcell_corr);
% Save list of selective cells 
if save && ~isempty(selcell_orig) && ~isempty(selcell_corr)
    cd(figdir);
    fid = fopen('selectivecells_orig.txt','w');
%     fprintf(fid,'%s\n',[num2str(size(cellList,1)) ' cells']);
    for ii = 1:size(selcell_orig,1)
        fprintf(fid,'%s\n',selcell_orig{ii,1});
    end
    fclose(fid);
    fid = fopen('selectivecells_corr.txt','w');
%     fprintf(fid,'%s\n',[num2str(size(cellList,1)) ' cells']);
    for ii = 1:size(selcell_corr,1)
        fprintf(fid,'%s\n',selcell_corr{ii,1});
    end
end
cd(cwd);


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

% disp(sum(sum(find(ceiling==0))) + sum(sum(find(floor==0))) + sum(sum(find(P4_TL==0))));

% Plot pillars
surf(P1_x, P1_y, PX_z, P1_BR);
alpha 1; shading flat;
surf(P2_x, P2_y, PX_z, P2_BL);
alpha 1; shading flat;
surf(P3_x, P3_y, PX_z, P3_TR);
alpha 1; shading flat;
surf(P4_x, P4_y, PX_z, P4_TL);

% Display parameters
alpha 1; shading flat; 
view(-35,20);
colormap jet;
colorbar;

% Set up axes
function axdisplay(ax,maxC)
% Troubleshoot
if maxC == 0
    maxC = 1;
end
% Patch boundaries of base fields
set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');


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
        % Pillar 1 (bottom right)
        patch([24 24 32 32],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 24 24],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[10.88 10.88 13.12 13.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Rabbit Poster on m_wall_25
        % Pillar 2 (bottom left)
        patch([8 8 16 16],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 8 8],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([10.88 10.88 13.12 13.12],[8 8 8 8],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Cat poster on m_wall_10
        patch([10.88 10.88 13.12 13.12],[16 16 16 16],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Pig poster on m_wall_29
        % Pillar 3 (top right)
        patch([24 24 32 32],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 24 24],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([26.88 26.88 29.12 29.12],[24 24 24 24],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Croc poster on m_wall_4
        patch([26.88 26.88 29.12 29.12],[32 32 32 32],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Donkey poster on m_wall_15
        % Pillar 4 (top left)
        patch([8 8 16 16],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 8 8],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 8 8],[26.88 26.88 29.12 29.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Camel poster on m_wall_20
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

function [] = savefigure(h,figtitle,figdir)
cwd = pwd;
cd(figdir);
saveas(h,figtitle,'png');
% saveas(h,figtitle,'fig');
% print('-painters',figtitle,'-dsvg');
cd(cwd);

function [x,y,z] = converttosurf(gridnum,bx,by)

switch gridnum

    case 3 % Floor
        
        x = [bx-1 bx-1 bx bx];
        y = [by-1 by by by-1];
        z = [0 0 0 0];

    case 4 % Ceiling
        
        x = [bx-1 bx-1 bx bx];
        y = [by-1 by by by-1];
        z = [40 40 40 40];

    case 5 % Walls
        
        face = ceil(bx/40);
        startx = mod(bx,40);
        if startx == 0
            startx = 40;
        end

        switch face
            case 1 % Left
                x = [0 0 0 0];
                y = [0+startx-1 0+startx-1 0+startx 0+startx];
            case 2 % Top
                x = [0+startx-1 0+startx-1 0+startx 0+startx];
                y = [40 40 40 40];
            case 3 % Right
                x = [40 40 40 40];
                y = [0+(40-startx+1) 0+(40-startx+1) 0+(40-startx) 0+(40-startx)];
            case 4 % Bottom
                x = [0+(40-startx+1) 0+(40-startx+1) 0+(40-startx) 0+(40-startx)];
                y = [0 0 0 0];
        end
        z = [16+by-1 16+by 16+by 16+by-1];

    case 6 % Pillar 1 bottom right
        
        face = ceil(bx/8);
        startx = mod(bx,8);
        if startx == 0
            startx = 8;
        end

        switch face
            case 1 % Left
                x = [24 24 24 24];
                y = [8+startx-1 8+startx-1 8+startx 8+startx];
            case 2 % Top
                x = [24+startx-1 24+startx-1 24+startx 24+startx];
                y = [16 16 16 16];
            case 3 % Right
                x = [32 32 32 32];
                y = [8+(8-startx+1) 8+(8-startx+1) 8+(8-startx) 8+(8-startx)];
            case 4 % Bottom
                x = [24+(8-startx+1) 24+(8-startx+1) 24+(8-startx) 24+(8-startx)];
                y = [8 8 8 8];
        end
        z = [16+by-1 16+by 16+by 16+by-1];

    case 7 % Pillar 2 bottom left
        
        face = ceil(bx/8);
        startx = mod(bx,8);
        if startx == 0
            startx = 8;
        end
    
        switch face
            case 1 % Left
                x = [8 8 8 8];
                y = [8+startx-1 8+startx-1 8+startx 8+startx];
            case 2 % Top
                x = [8+startx-1 8+startx-1 8+startx 8+startx];
                y = [16 16 16 16];
            case 3 % Right
                x = [16 16 16 16];
                y = [8+(8-startx+1) 8+(8-startx+1) 8+(8-startx) 8+(8-startx)];
            case 4 % Bottom
                x = [8+(8-startx+1) 8+(8-startx+1) 8+(8-startx) 8+(8-startx)];
                y = [8 8 8 8];
        end
        z = [16+by-1 16+by 16+by 16+by-1];

    case 8 % Pillar 3 top right
        
        face = ceil(bx/8);
        startx = mod(bx,8);
        if startx == 0
            startx = 8;
        end
        
        switch face
            case 1 % Left
                x = [24 24 24 24];
                y = [24+startx-1 24+startx-1 24+startx 24+startx];
            case 2 % Top
                x = [24+startx-1 24+startx-1 24+startx 24+startx];
                y = [32 32 32 32];
            case 3 % Right
                x = [32 32 32 32];
                y = [24+(8-startx+1) 24+(8-startx+1) 24+(8-startx) 24+(8-startx)];
            case 4 % Bottom
                x = [24+(8-startx+1) 24+(8-startx+1) 24+(8-startx) 24+(8-startx)];
                y = [24 24 24 24];
        end
        z = [16+by-1 16+by 16+by 16+by-1];

    case 9 % Pillar 4 top left
        
        face = ceil(bx/8);
        startx = mod(bx,8);
        if startx == 0
            startx = 8;
        end
        
        switch face
            case 1 % Left
                x = [8 8 8 8];
                y = [24+startx-1 24+startx-1 24+startx 24+startx];
            case 2 % Top
                x = [8+startx-1 8+startx-1 8+startx 8+startx];
                y = [32 32 32 32];
            case 3 % Right
                x = [16 16 16 16];
                y = [24+(8-startx+1) 24+(8-startx+1) 24+(8-startx) 24+(8-startx)];
            case 4 % Bottom
                x = [8+(8-startx+1) 8+(8-startx+1) 8+(8-startx) 8+(8-startx)];
                y = [24 24 24 24];
        end
        z = [16+by-1 16+by 16+by 16+by-1];

end

function [gridnum,x,y] = findgrid(px,objtype)
% returns grid number and plot coords (x goes left to right, y goes bottom
% to top)

switch objtype
    case 'place'
        mapLdummy = 1:1600;
        gridnum = 3;
        temp = flipud(reshape(mapLdummy, 40, 40)');
    case 'spatialview'
        mapLdummy = 1:5122;
        if px == 1 % Cue
            gridnum = 1;
            x = 1;
            y = 1;
        elseif px == 2 % Hint
            gridnum = 2;
            x = 1; 
            y = 1;
        elseif px >= 3 && px <= 1602 % Floor
            gridnum = 3;
            temp = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');
        elseif px >= 1603 && px <= 3202 % Ceiling
            gridnum = 4;
            temp = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');
        elseif px >= 3203 && px <= 4482 % Walls
            gridnum = 5;
            temp = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');
        elseif px >= 4483 && px <= 4642 % Pillar 1
            gridnum = 6;
            temp = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
        elseif px >= 4643 && px <= 4802 % Pillar 2
            gridnum = 7;
            temp = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
        elseif px >= 4803 && px <= 4962 % Pillar 3
            gridnum = 8;
            temp = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
        elseif px >= 4963 && px <= 5122 % Pillar 4
            gridnum = 9;
            temp = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
        end
end
[y,x] = find(temp == px);
y = size(temp,1)-y+1;

% Set zero firing rate (occupied) pixels to black to differentiate from other low
% firing pixels
function [] = settoblack(basemap,baseobj)

blackpx = find(basemap == 0);
for pp = 1:size(blackpx,2)
    switch baseobj
        case 'place'
            gnum = 1;
            [gnum,xx,yy] = findgrid(blackpx(pp),baseobj);
        case 'spatialview'
            if blackpx(pp) > 3 % Make sure not cue or hint
               [gnum,xx,yy] = findgrid(blackpx(pp),baseobj);
            else
                continue;
            end

    end
    [x,y,z] = converttosurf(gnum,xx,yy);
    patch(x,y,z,[1 1 1 1],'FaceColor','k','FaceAlpha',0.7,'EdgeColor','none');
end
% 
% function [smoothedRate,smoothedSpk,smoothedDur]=adsmooth(dur,spk,alpha)
% % Adaptive smoothing of rate maps.
% %
% %       [smoothedRate,smoothedSpk,smoothedPos]=rates_adaptivesmooth(posMap,spkMap,alpha)
% %
% % Each bin is smoothed using a flat, circular kernal. The circle radius 
% % is set for each bin, indivdually, such that 
% %
% %   radius => alpha ./ ( sqrt(nSpike) .* nDwell )
% %
% % where nSpike and nDwell are the number of spikes, and the amount of dwell time (in s) within the kernel.
% %
% % smoothedRate, smoothedSpk, smoothedPos are the smoothed maps (spike and pos maps are smoothed 
% % with the same kernal group as for the rate map.
% 
% % Check for empty spk maps %
% if sum(sum(spk))==0
%     smoothedDur=dur;    smoothedDur(dur==0)=nan;
%     smoothedSpk=spk;    smoothedSpk(dur==0)=nan;
%     smoothedRate=spk;   smoothedRate(dur==0)=nan;
%     return
% end
% % Pre-assign output %
% smoothedDur=zeros(size(dur));
% smoothedSpk=zeros(size(dur));
% % Visited env template: use this to get numbers of visited bins in filter at edge of environemnt %
% vis=zeros(size(dur));
% vis(dur>0)=1;
% % Pre-assign map which records which bins have passed %
% smoothedCheck=false(size(dur));
% smoothedCheck(dur==0)=true; % Disregard unvisited - mark as already done.
% % Pre-assign list of radii used (this is for reporting purposes, not used for making maps) %
% radiiUsedList=nan(1,sum(sum(dur>0)));
% radiiUsedCount=1;
% 
% %%% Run increasing radius iterations %%%
% r=1; % Circle radius
% boundary=0; % IMFILTER boundary condition
% while any(any(~smoothedCheck))
%     % Check radius isn't getting too big (if >map/2, stop running) %
%     if r>max(size(dur))/2
% %     if r>20
%         smoothedSpk(~smoothedCheck)=nan;
%         smoothedDur(~smoothedCheck)=nan;
%         break
%     end
%     % Construct filter kernel ...
%     % Place: Flat disk, where r>=distance to bin centre %
%     f=fspecial('disk',r); 
%     f(f>=(max(max(f))/3))=1;
%     f(f~=1)=0;   
%     % Filter maps (get N spikes and pos sum within kernel) %
%     fSpk=imfilter(spk,f,boundary);
%     fDur=imfilter(dur,f,boundary);
%     fVis=imfilter(vis,f,boundary);
%     % Which bins pass criteria at this radius? %
%     warning('off', 'MATLAB:divideByZero');
%     binsPassed=alpha./(sqrt(fSpk).*fDur) <= r;
%     warning('on', 'MATLAB:divideByZero');
%     binsPassed=binsPassed & ~smoothedCheck; % Only get the bins that have passed in this iteration.
%     % Add these to list of radii used %
%     nBins=sum(binsPassed(:));
%     radiiUsedList(radiiUsedCount:radiiUsedCount+nBins-1)=r;
%     radiiUsedCount=radiiUsedCount+nBins;
%     % Assign values to smoothed maps %
%     smoothedSpk(binsPassed)=fSpk(binsPassed)./fVis(binsPassed);
%     smoothedDur(binsPassed)=fDur(binsPassed)./fVis(binsPassed);
%     % Record which bins were smoothed this iteration %
%     smoothedCheck(binsPassed)=true;
%     % Increase circle radius (half-bin steps) %
%     r=r+0.5; % Increase radius in 0.5 bin steps.
% end
% 
% % Assign Output %
% warning('off', 'MATLAB:divideByZero');
% smoothedRate=smoothedSpk./smoothedDur;
% warning('on', 'MATLAB:divideByZero');
% smoothedRate(dur==0)=nan;
% smoothedDur(dur==0)=nan;
% smoothedSpk(dur==0)=nan;
% 
% % Report radii sizes %
% if 0
%     hAllFigs = get(0, 'children');
%     hFig = findobj(hAllFigs, 'flat', 'tag', 'adaptiveSmoothPlotWindow');
%     if isempty(hFig);
%         hFig=figure;
%         set(hFig,'tag','adaptiveSmoothPlotWindow');
%     else
%         figure(hFig);
%     end
%     hist(radiiUsedList,1:10);
%     uiwait(hFig,1.5);
% end