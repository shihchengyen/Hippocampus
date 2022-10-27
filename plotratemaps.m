function [selcell_orig] = plotratemaps(objtype,criteria,save,maptype,varargin)

% Function to plot rate maps for place, view and heading direction
% Without additional inputs, it will work from and save in current
% directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% objtype: 'place', 'view', or 'direction'
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

sortcrit = 0;
video = 0;
filttype = 'FiltVel';
pix = 1;

cwd = '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis';

if nargin > 4 % If plotting a single cell and map is already given as input (only used in placebyview.m) 
    % Load cell list
    cellList = varargin(1);
    cd([cellList{1} '/' filttype '/' num2str(pix) 'px/']);
    % Load object
    switch objtype
        case 'place'
            objMain = load('vmpc.mat');
            objMain = objMain.vmp;
        case 'view'
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
        case 'view'
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
            objMS = load('c_vmms0_2var.mat');
            objMS = objMS.vmms;
        case 'view'
            objMain = load('c_vmsv.mat');
            objMain = objMain.vms;
            objMS = load('c_vmms0_2var.mat');
            objMS = objMS.vmms;
            objP = load('c_vmpc.mat');
            objP = objP.vmp;
        case 'mixsel0'
            objMain = load('c_vmms0.mat');
            objMain = objMain.vmms;
            objPlace = load('c_vmpc.mat');
            objPlace = objPlace.vmp;
            objView = load('c_vmsv.mat');
            objView = objView.vms;
            objHeaddirection = load('c_vmhd.mat');
            objHeaddirection = objHeaddirection.vmd;
        case 'mixsel1'
            objMain = load('c_vmms1.mat');
            objMain = objMain.vmms;
            objPlace = load('c_vmpc.mat');
            objPlace = objPlace.vmp;
            objView = load('c_vmsv.mat');
            objView = objView.vms;
            objHeaddirection = load('c_vmhd.mat');
            objHeaddirection = objHeaddirection.vmd;
        case 'allprop'
            objMain0 = load('c_vmms0.mat');
            objMain0 = objMain0.vmms;
            objMain = load('c_vmms1.mat'); % Treat this as the main obj
            objMain = objMain.vmms;
            objPlace = load('c_vmpc.mat');
            objPlace = objPlace.vmp;
            objView = load('c_vmsv.mat');
            objView = objView.vms;
            objHeaddirection = load('c_vmhd.mat');
            objHeaddirection = objHeaddirection.vmd;
        case 'headdirection'
            objMain = load('c_vmhd.mat');
            objMain = objMain.vmd;
            objMS = load('c_vmms0_2var.mat');
            objMS = objMS.vmms;
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
    if strcmp(objtype,'mixsel0') || strcmp(objtype,'mixsel1')
        figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/mixsel/UseCorr' objtype(end)];
    else
        figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/' objtype '_' criteria '_' num2str(maptype)];
    end
    if save
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

%%%% PATCH - c_vmcorr may have fewer cells than c_vmpc or c_vmv. 
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
        critname = [upper(criteria) '_adsm'];
        critshname = [upper(criteria) 'sh_adsm'];
    case 'smooth'
        mapname = 'maps_sm';
        critname = ['crit_sm'];
        critshname = ['critsh_sm'];
        secmapname = 'secsmoothmap_adsm';
    case 'dur'
        mapname = 'dur_raw';
        critname = [upper(criteria) '_adsm'];
        critshname = [upper(criteria) 'sh_adsm'];
    case 'boxcar'
        mapname = 'maps_bcsm';
        critname = [upper(criteria) '_bcsm'];
        critshname = [upper(criteria) 'sh_bcsm'];
        secmapname = 'secsmoothmap_bcsm';
    case 'disk'
        mapname = 'maps_dksm';
        critname = [upper(criteria) '_dksm'];
        critshname = [upper(criteria) 'sh_dksm'];
        secmapname = 'secsmoothmap_dksm';
    case 'adaptive'
        mapname = 'maps_adsm';
        critname = [upper(criteria) '_adsm'];
        critshname = [upper(criteria) 'sh_adsm'];
        secmapname = 'secsmoothmap_adsm';
end
if strcmp(objtype,'place') 
    maps = objMain.data.(mapname);
%     mapscorr = objCorr.data.pv.(['maps_sm' '_corrp']);
%     mapscorr = objCorr.data.([mapname '_corrp']);
    objtype_short = 'p';
elseif strcmp(objtype,'view')
    maps = objMain.data.(mapname);
%     mapscorr = objCorr.data.pv.(['maps_sm' '_corrp']);
%     mapscorr = objCorr.data.([mapname '_corrv']);
    objtype_short = 'v';
elseif strcmp(objtype,'headdirection')
    maps = objMain.data.(mapname);
%     mapscorr = objCorr.data.pv.(['maps_sm' '_corrp']);
%     mapscorr = objCorr.data.([mapname '_corrp']);
    objtype_short = 'h';
elseif strcmp(objtype,'mixsel0') || strcmp(objtype,'mixsel1')
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
    thr_sh = [objMain.data.(critname); objMain.data.(critshname)];
    thr_pop = prctile(thr_sh,95);
    z_pop = zscore(thr_sh);
    
    for ii = 1:size(setsessions,1) % For each session
        cells_indList = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_indList) % For each cell
            
            cell_indList = cells_indList(jj);
            cell_ind = find(strcmp(objMain.data.origin,cellList(cells_indList(jj))));
            okminspk = sum(objMain.data.spk_raw(cell_ind,:)) >= 100;
            if ~okminspk
                disp(cellList(cell_indList));
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
%                     print('-painters',figtitle,'-dvg');
                    cd(cwd);
                    close(figure(fig));
                end
                
                fig = fig + 1;
                subpnum = 1;
            end

            % Get shuffled SI cutoff for this cell - 95th percentile
            crit = objMain.data.(critname)(cell_ind,1);
            thr_cell = prctile(objMain.data.(critshname)( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles ,1 ) ,95);
            z_cell = z_pop(cell_ind,1);
            
            % Get map
            if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
                mapLin = maps(cell_ind,:);
                % if corrected map exists, get it
                if any(ismember(objMain.data.origin{cell_ind},objCorr.data.origin))
                    [~,corr_ind] = ismember(objMain.data.origin{cell_ind},objCorr.data.origin);
                    mapLincorr = objCorr.data.pv(corr_ind).(['maps_sm' '_corrp']);
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
            ax.Title.String = {horzcat(num2str(setsessions(ii)), 'ch',num2str(identifiers(cell_indList,4)),'c',num2str(identifiers(cell_indList,5)),', ',num2str(nanmax(mapLin),3),'Hz'), horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z-',criteria,'= ',num2str(z_cell))};
            
            if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
            elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'm';
                crosscellthresh = crosscellthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'b';
                crosspopthresh = crosspopthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            else
                ax.Title.Color = 'k';
            end
            if crit >= thr_pop && okminspk && maxC>=0.7
                selcell_orig(end+1,1) = cell_ind;
            end

            % Patch standing point if placebyview
            if nargin > 5 % If placebyview
                    fieldCoords
                    patch([fieldCoords(1)-2 fieldCoords(1)-2 fieldCoords(1)+1 fieldCoords(1)+1],[fieldCoords(2)-2 fieldCoords(2)+1 fieldCoords(2)+1 fieldCoords(2)-2], [0 0 0 0] ,[1 1 1 1],'FaceColor','none');
            end
            
            plothalves = false;
            if plothalves
                % Intra-session correlation %%%%%% NOTE: Should use boxcar
                % smoothed map
                map1 = objMain.data.([mapname '1'])(cell_ind,:);
                map2 = objMain.data.([mapname '2'])(cell_ind,:);
                vis1 = ~isnan(map1);
                vis2 = ~isnan(map2);
                vis = vis1 & vis2; % Correlate only visited bins;
                intracorr = corr2(map1(vis), map2(vis));

                crit1 = objMain.data.([critname '1'])(cell_ind);
                crit2 = objMain.data.([critname '2'])(cell_ind);

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
                    [~,~]= plotmap(mapLin,objtype);

                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    ax.Title.String = {horzcat(half,' half: ','corr=',num2str(intracorr,2),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};

                    if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                        ax.Title.Color = 'r';
                    elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
                        ax.Title.Color = 'm';
                    elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                        ax.Title.Color = 'b';
                    else
                        ax.Title.Color = 'k';
                    end

                    subpnum = subpnum + 1;
                end
            end
            
            % Plot corrected maps
            
            subpnum = subpnum + 1;
            critcorrset = [];
            for cc = 1:2 % pv and ph
                
                if cc == 1
                    temp = objCorr.data.pv;
                    msvar_short = 'pv';
                elseif cc == 2
                    temp = objCorr.data.ph;
                    msvar_short = 'ph';
                end
                if isempty(corr_ind)
                    critcorr = nan;
                else
                    critcorr = temp(corr_ind).(['SIC_sm' '_corr' objtype_short]);
    %                 switch objtype
    %                     case 'place'
    %                         critcorr = objCorr.data.([critname '_corrp'])(corr_ind);
    %                     case 'view'
    %                         critcorr = objCorr.data.([critname '_corrv'])(corr_ind);
    %                 end
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
                    ax.Title.String = {horzcat('Corrected',msvar_short,temp(corr_ind).llhpicklabel,'of',...
                        num2str(size(temp(corr_ind).llh,1)),': ',num2str(nanmax(mapLincorr),3),'Hz'), ...
                        horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                    set(ax,'CLim',[0 maxCcorr]);
                else
                    ax.Title.String = 'No corrected map';
                end

                if critcorr >= thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
                        ax.Title.Color = 'r';
                elseif critcorr >= thr_cell && critcorr < thr_pop && okminspk && maxCcorr>=0.7
                    ax.Title.Color = 'm';
                elseif critcorr < thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
                    ax.Title.Color = 'b';
                else
                    ax.Title.Color = 'k';
                end
                hold off;
                subpnum = subpnum + 1;
                critcorrset = [critcorrset critcorr];
            end
            if min(critcorrset) >= thr_pop && okminspk && maxCcorr>=0.7
                selcell_corr(end+1,1) = cell_ind;
            end
            
            
            % Plot predicted artefactual map
            
            for cc = 1:2 % pv and ph
                
                if cc == 1
                    temp = objCorr.data.pv;
                    msvar_short = 'pv';
                elseif cc == 2
                    temp = objCorr.data.ph;
                    msvar_short = 'ph';
                end
                % Setup object
                h = gcf;
                ax = subplot(plotgridv,plotgridh,subpnum);
                hold on;

                % Plot map
                distmap = temp(corr_ind).(['maps_dist_' objtype_short]);
                [~,~] = plotmap(distmap,objtype);
                maxCcorr = nanmax(distmap);

                % Set up axes
                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                if ~isempty(corr_ind)
                    if cc == 1
                        ax.Title.String = {['Distributed map ' msvar_short];...
                            ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)];...
                            ['dist ratio v= ' num2str(temp(corr_ind).distratio_v)]};
                    elseif cc == 2
                        ax.Title.String = {['Distributed map ' msvar_short];...
                            ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)];...
                            ['dist ratio h= ' num2str(temp(corr_ind).distratio_h)]};
                    end
                    set(ax,'CLim',[0 maxCcorr]);
                else
                    ax.Title.String = 'No distributed map';
                end

                hold off;

                subpnum = subpnum + 1;
            end
            
%             % Plot covariance matrix
%             
%             % Setup object
%             h = gcf;
%             ax = subplot(plotgridv,plotgridh,subpnum);
%             hold on;
% 
%             % Plot map
%             im = imagesc(objCorr.data.covmat_norm{corr_ind});
%             set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
%             set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
%             colormap jet;
%             colorbar;
%             
%             % Replace NaNs with zeros in covariance matrix for norm calculations
% %             covmat = objCorr.data.covmat{corr_ind};
%             l1norm = objCorr.data.l1norm(corr_ind);
%             l2norm = objCorr.data.l2norm(corr_ind);
% %             covmat(isnan(covmat)) = 0;
% %             % Calculate norms
% %             norml1 = norm(covmat,1); % maximum of column sum
% %             norml2 = norm(covmat,2); % maximum single value
%             
%             % Set up axes
%             set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%                 'XColor','none','YColor','none','ZColor','none',...
%                 'FontSize',14,'GridLineStyle','none','Color','none');
%             if ~isempty(corr_ind)
%                 ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(l1norm,2)), horzcat('l2=', num2str(l2norm,2))};
%             else
%                 ax.Title.String = 'No corrected map';
%             end
%             
%             subpnum = subpnum + 1;
%             
%             hold off;

                
        end
        if save
            cwd = pwd;
            cd(figdir);
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_indList)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
%             print('-painters',figtitle,'-dvg');
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
    
elseif strcmp(objtype, 'view') 
    
    fig = 1;
    subpnum = 1;

    thr_sh = [objMain.data.(critname); objMain.data.(critshname)];
    thr_pop = prctile(thr_sh,95);
    z_pop = zscore(thr_sh);
    
    for ii = 1:size(setsessions,1) % For each session
        
        cells_indList = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_indList) % For each cell
            
            cell_indList = cells_indList(jj);
            cell_ind = find(strcmp(objMain.data.origin,cellList(cells_indList(jj))));
            cell_indP = find(strcmp(objP.data.origin,cellList(cells_indList(jj))));
            okminspk = sum(objP.data.spk_raw(cell_indP,:)) >= 100;
            if ~okminspk
                disp(cellList(cell_indList));
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
%                     print('-painters',figtitle,'-dvg');
                    cd(cwd);
                    close(figure(fig));
                end
                
                fig = fig + 1;
                subpnum = 1;
            end

            % Get shuffled SI cutoff for this cell - 95th percentile
            crit = objMain.data.(critname)(cell_ind,1);
            thr_cell = prctile(objMain.data.(critshname)( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles,1 ),95);
            z_cell = z_pop(cell_ind,1);

            if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells)

                mapLin = maps(cell_ind,:);
                % if corrected map exists, get it
                if any(ismember(objMain.data.origin{cell_ind},objCorr.data.origin))
                    [~,corr_ind] = ismember(objMain.data.origin{cell_ind},objCorr.data.origin);
                    mapLincorr = objCorr.data.pv(corr_ind).(['maps_sm' '_corr' objtype_short]);
                else
                    mapLincorr = nan(size(mapLin));
                    corr_ind = [];
                end
                % Set up figure
                h = figure(fig);
                ax = subplot(plotgridv,plotgridh,subpnum);
                
            end
            
            % Set up figure
            h = gcf;
            hold on;
            ax = gca;
            
            % Plot map
            [mapGrid,~]= plotmap(mapLin,objtype);
            
            % Figure and axes properties
            figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),...
                'c',num2str(identifiers(cell_ind,5)));
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
                ax.Title.String = {horzcat('Cue: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),...
                    'c',num2str(identifiers(cell_ind,5)),', ',num2str(nanmax(mapLin),3),'Hz'),...
                    horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),...
                    horzcat('z',num2str(z_cell))};
            elseif nanmax(mapLin(2)) > maxC
                ax.Title.String = {horzcat('Hint: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z',num2str(z_cell))};
            else
                ax.Title.String = {horzcat(num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z',num2str(z_cell))};
            end
            % Patch environment boundaries
            patchenvbounds(objtype);
            
            % Denote if significant spatial information
            if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
            elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'm';
                crosscellthresh = crosscellthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            elseif crit >= thr_pop && crit < thr_cell && okminspk && maxC>=0.7
                ax.Title.Color = 'b';
                crosspopthresh = crosspopthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            else
                ax.Title.Color = 'k';
            end
            if crit >= thr_pop && okminspk && maxC>=0.7
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
            plothalves = false;
            if plothalves
                map1 = objMain.data.([mapname '1'])(cell_ind,:);
                map2 = objMain.data.([mapname '2'])(cell_ind,:);
    %             map1 = emptyinsidepillar(map1);
                vis1 = ~isnan(map1);
    %             map2 = emptyinsidepillar(map2);
                vis2 = ~isnan(map2);
                vis = vis1 & vis2; % Correlate only visited bins;
                intracorr = corr2(map1(vis), map2(vis));

                crit1 = objMain.data.([critname '1'])(cell_ind);
                crit2 = objMain.data.([critname '2'])(cell_ind);

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

                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    ax.Title.String = {horzcat(half,' half: ','corr=',num2str(intracorr,2),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};

                    if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                        ax.Title.Color = 'r';
                    elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
                        ax.Title.Color = 'm';
                    elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                        ax.Title.Color = 'b';
                    else
                        ax.Title.Color = 'k';
                    end

                    subpnum = subpnum + 1;
                end
            end
            
            % Plot corrected maps
            subpnum = subpnum + 1;
            critcorrset = [];
            for cc = 1 % pv & ph
                if cc == 1
                    temp = objCorr.data.pv;
                    msvar_short = 'pv';
                elseif cc == 2
                    temp = objCorr.data.ph;
                    msvar_short = 'ph';
                end
                if isempty(corr_ind)
                    critcorr = nan;
                else
                    critcorr = temp(corr_ind).(['crit_sm' '_corr' objtype_short]);
                end

                % Setup object
                h = gcf;
                ax = subplot(plotgridv,plotgridh,subpnum);
                hold on;

                % Plot map
                [~,~]= plotmap(mapLincorr,objtype);
                % Patch environment boundaries
                patchenvbounds(objtype);
                maxCcorr = nanmax(mapLincorr(3:end));

                % Set up axes
                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                if ~isempty(corr_ind)
                    ax.Title.String = {horzcat('Corrected',msvar_short,temp(corr_ind).llhpicklabel,'of',...
                        num2str(size(temp(corr_ind).llh,1)),': ',num2str(nanmax(mapLincorr),3),'Hz'), ...
                        horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                    set(ax,'CLim',[0 maxCcorr]);
                else 
                    ax.Title.String = 'No corrected map';
                end

                if critcorr >= thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
                        ax.Title.Color = 'r';
                elseif critcorr >= thr_cell && critcorr < thr_pop && okminspk && maxCcorr>=0.7
                    ax.Title.Color = 'm';
                elseif critcorr < thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
                    ax.Title.Color = 'b';
                else
                    ax.Title.Color = 'k';
                end
                critcorrset = [critcorrset critcorr];
                subpnum = subpnum + 1;
            end
            if min(critcorrset) >= thr_pop && okminspk && maxCcorr>=0.7
                selcell_corr(end+1,1) = cell_ind;
            end
            
            % Plot predicted artefactual map
            
            for cc = 1 % pv and ph
                if cc == 1
                    temp = objCorr.data.pv;
                    msvar_short = 'pv';
                elseif cc == 2
                    temp = objCorr.data.ph;
                    msvar_short = 'ph';
                end
                % Setup object
                h = gcf;
                ax = subplot(plotgridv,plotgridh,subpnum);
                hold on;

                % Plot map
                distmap = temp(corr_ind).(['maps_dist_' objtype_short]);
                [~,~] = plotmap(distmap,objtype);
                maxCcorr = nanmax(distmap);

                % Set up axes
                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                if ~isempty(corr_ind)
                    if cc == 1
                        ax.Title.String = {['Distributed map ' msvar_short];...
                            ['dist ratio v = ' num2str(temp(corr_ind).distratio_v)];...
                            ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)]};
                    elseif cc == 2
                        ax.Title.String = {['Distributed map ' msvar_short];...
                            ['dist ratio h = ' num2str(temp(corr_ind).distratio_h)];...
                            ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)]};

                    end
                    set(ax,'CLim',[0 maxCcorr]);
                else
                    ax.Title.String = 'No distributed map';
                end
                
            end
            hold off;
            
            
%             subpnum = subpnum + 1;
            % Patch
            if mod(subpnum,plotgridh) ~= 0
                subpnum = ceil(subpnum/plotgridh)*plotgridh+1;
            end
            
%             % Plot covariance matrix
%             
%             % Setup object
%             h = gcf;
%             ax = subplot(plotgridv,plotgridh,subpnum);
%             hold on;
% 
%             % Plot map
%             im = imagesc(objCorr.data.covmat_norm{corr_ind});
%             set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
% %             set(ax,'CLim',[-1 1]);
%             set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
%             colormap jet;
%             colorbar;
%             
%             % Replace NaNs with zeros in covariance matrix for norm calculations
%             l1norm = objCorr.data.l1norm(corr_ind);
%             l2norm = objCorr.data.l2norm(corr_ind);
% %             covmat = objCorr.data.covmat{corr_ind};
% %             covmat(isnan(covmat)) = 0;
% %             % Calculate norms
% %             l1norm = norm(covmat,1); % maximum of column sum
% %             l2norm = norm(covmat,2); % maximum single value
%             
%             % Set up axes
%             set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%                 'XColor','none','YColor','none','ZColor','none',...
%                 'FontSize',14,'GridLineStyle','none','Color','none');
%             if ~isempty(corr_ind)
%                 ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(l1norm,2)), horzcat('l2=', num2str(l2norm,2))};
%             else
%                 ax.Title.String = 'No corrected map';
%             end
%             
%             subpnum = subpnum + 1;
%             
%             hold off;
            
        end
        
        if save
            cwd = pwd;
            cd(figdir);
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_indList)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
%             print('-painters',figtitle,'-dvg');
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
    
elseif strcmp(objtype,'headdirection')
    
    % For each session, plot rate maps for each cell
    fig = 1;
    subpnum = 1;
    
    thr_sh = [objMain.data.(critname); objMain.data.(critshname)];
    thr_pop = prctile(thr_sh,95);
    z_pop = zscore(thr_sh);
    
    for ii = 1:size(setsessions,1) % For each session
        
        cells_indList = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_indList) % For each cell
            
            cell_indList = cells_indList(jj);
            cell_ind = find(strcmp(objMain.data.origin,cellList(cells_indList(jj))));
            okminspk = sum(objMain.data.spk_raw(cell_ind,:)) >= 100;
            if ~okminspk
                disp(cellList(cell_indList));
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
%                     print('-painters',figtitle,'-dvg');
                    cd(cwd);
                    close(figure(fig));
                end
                
                fig = fig + 1;
                subpnum = 1;
            end

            % Get shuffled SI cutoff for this cell - 95th percentile

            % Patch
            objMain.data.Args.NumShuffles = 10000;
            % End patch
            crit = objMain.data.(critname)(cell_ind,1);
            thr_cell = prctile(objMain.data.(critshname)( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles ,1 ) ,95);
            z_cell = z_pop(cell_ind,1);

            % Get map
            if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
                mapLin = maps(cell_ind,:);
                % if corrected map exists, get it
                if any(ismember(objMain.data.origin{cell_ind},objCorr.data.origin))
                    [~,corr_ind] = ismember(objMain.data.origin{cell_ind},objCorr.data.origin);
                    mapLincorr = objCorr.data.ph(corr_ind).(['maps_sm' '_corrh']);
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
            maxC = nanmax(mapLin);
 
            ax.Title.String = {horzcat(num2str(setsessions(ii)), 'ch',num2str(identifiers(cell_indList,4)),'c',num2str(identifiers(cell_indList,5)),', ',num2str(nanmax(mapLin),3),'Hz'), horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z-',criteria,'= ',num2str(z_cell))};
            
            if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'r';
                crossallthresh = crossallthresh + 1;
            elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'm';
                crosscellthresh = crosscellthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                ax.Title.Color = 'b';
                crosspopthresh = crosspopthresh + 1;
                crosseitherthresh = crosseitherthresh + 1;
            else
                ax.Title.Color = 'k';
            end
            if crit >= thr_pop && okminspk && maxC>=0.7
                selcell_orig(end+1,1) = cell_ind;
            end
            
            % Intra-session correlation %%%%%% NOTE: Should use boxcar
            % smoothed map
            map1 = objMain.data.([mapname '1'])(cell_ind,:);
            map2 = objMain.data.([mapname '2'])(cell_ind,:);

            vis1 = ~isnan(map1);
            vis2 = ~isnan(map2);
            vis = vis1 & vis2; % Correlate only visited bins;
            intracorr = corr2(map1(vis), map2(vis));
            
            crit1 = objMain.data.([critname '1'])(cell_ind);
            crit2 = objMain.data.([critname '2'])(cell_ind);

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
                [~,~]= plotmap(mapLin,objtype);
                maxC = nanmax(mapLin);

                ax.Title.String = {horzcat(half,' half: ','corr=',num2str(intracorr,2),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                
                if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                    ax.Title.Color = 'r';
                elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
                    ax.Title.Color = 'm';
                elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
                    ax.Title.Color = 'b';
                else
                    ax.Title.Color = 'k';
                end

                subpnum = subpnum + 1;
            end
            
            % Plot corrected map - PLACEHOLDER

%             % Setup object
%             h = gcf;
%             ax = subplot(plotgridv,plotgridh,subpnum);
%             hold on;

            
            critcorrset = [];
            for cc = 2 % pv and ph
                
                if cc == 1
                    temp = objCorr.data.pv;
                    msvar_short = 'pv';
                elseif cc == 2
                    temp = objCorr.data.ph;
                    msvar_short = 'ph';
                end
                if isempty(corr_ind)
                    critcorr = nan;
                else
                    critcorr = temp(corr_ind).(['crit_sm' '_corr' objtype_short]);
    %                 switch objtype
    %                     case 'place'
    %                         critcorr = objCorr.data.([critname '_corrp'])(corr_ind);
    %                     case 'view'
    %                         critcorr = objCorr.data.([critname '_corrv'])(corr_ind);
    %                 end
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
                    ax.Title.String = {horzcat('Corrected',msvar_short,temp(corr_ind).llhpicklabel,'of',...
                        num2str(size(temp(corr_ind).llh,1)),': ',num2str(nanmax(mapLincorr),3),'Hz'), ...
                        horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                    set(ax,'CLim',[0 maxCcorr]);
                else
                    ax.Title.String = 'No corrected map';
                end

                if critcorr >= thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
                        ax.Title.Color = 'r';
                elseif critcorr >= thr_cell && critcorr < thr_pop && okminspk && maxCcorr>=0.7
                    ax.Title.Color = 'm';
                elseif critcorr < thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
                    ax.Title.Color = 'b';
                else
                    ax.Title.Color = 'k';
                end
                hold off;
                subpnum = subpnum + 1;
                critcorrset = [critcorrset critcorr];
            end
            if min(critcorrset) >= thr_pop && okminspk && maxCcorr>=0.7
                selcell_corr(end+1,1) = cell_ind;
            end
            
            
            % Plot predicted artefactual map - PLACEHOLDER
            
%             % Setup object
%             h = gcf;
%             ax = subplot(plotgridv,plotgridh,subpnum);
%             hold on;
             
            for cc = 2 % pv and ph
                
                if cc == 1
                    temp = objCorr.data.pv;
                    msvar_short = 'pv';
                elseif cc == 2
                    temp = objCorr.data.ph;
                    msvar_short = 'ph';
                end
                % Setup object
                h = gcf;
                ax = subplot(plotgridv,plotgridh,subpnum);
                hold on;

                % Plot map
                distmap = temp(corr_ind).(['maps_dist_' objtype_short]);
                [~,~] = plotmap(distmap,objtype);
                maxCcorr = nanmax(distmap);

                % Set up axes
                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                if ~isempty(corr_ind)
                    if cc == 1
                        ax.Title.String = {['Distributed map ' msvar_short];...
                            ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)];...
                            ['dist ratio v= ' num2str(temp(corr_ind).distratio_v)]};
                    elseif cc == 2
                        ax.Title.String = {['Distributed map ' msvar_short];...
                            ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)];...
                            ['dist ratio h= ' num2str(temp(corr_ind).distratio_h)]};
                    end
                    set(ax,'CLim',[0 maxCcorr]);
                else
                    ax.Title.String = 'No distributed map';
                end

                hold off;

                subpnum = subpnum + 1;
            end
            
%             % Plot covariance matrix
%             
%             % Setup object
%             h = gcf;
%             ax = subplot(plotgridv,plotgridh,subpnum);
%             hold on;
% 
%             % Plot map
%             im = imagesc(objCorr.data.covmat_norm{corr_ind});
%             set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
%             set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
%             colormap jet;
%             colorbar;
%             
%             % Replace NaNs with zeros in covariance matrix for norm calculations
% %             covmat = objCorr.data.covmat{corr_ind};
%             l1norm = objCorr.data.l1norm(corr_ind);
%             l2norm = objCorr.data.l2norm(corr_ind);
% %             covmat(isnan(covmat)) = 0;
% %             % Calculate norms
% %             norml1 = norm(covmat,1); % maximum of column sum
% %             norml2 = norm(covmat,2); % maximum single value
%             
%             % Set up axes
%             set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%                 'XColor','none','YColor','none','ZColor','none',...
%                 'FontSize',14,'GridLineStyle','none','Color','none');
%             if ~isempty(corr_ind)
%                 ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(l1norm,2)), horzcat('l2=', num2str(l2norm,2))};
%             else
%                 ax.Title.String = 'No corrected map';
%             end
%             
%             subpnum = subpnum + 1;
%             
%             hold off;

                
        end
        if save
            cwd = pwd;
            cd(figdir);
            % Save figure
            figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_indList)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
            saveas(h,figtitle,'png');
%             print('-painters',figtitle,'-dvg');
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
    
elseif strcmp(objtype,'mixsel0') || strcmp(objtype,'mixsel1')

%     spatialvarpairs = objCorr.data.spatialvarpairs;
    spatialvarpairs = objMain.data.Args.spatialvarpairs;
    
    for ii = 1:size(setsessions,1) % For each session
        cells_indList = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_indList) % For each cell
            
            fig = 1;
            % Get cell index - unsorted
            cell_ind = cells_indList(jj); % within cell List
            cell_id = [num2str(identifiers(cell_ind,1)) 'ch' num2str(identifiers(cell_ind,4)) 'c' num2str(identifiers(cell_ind,5))];
            cell_indP = strcmp(objPlace.data.origin,cellList{cell_ind});
            cell_indV = strcmp(objView.data.origin,cellList{cell_ind});
            cell_indH = strcmp(objHeaddirection.data.origin,cellList{cell_ind});
            cell_indCorr = strcmp(objCorr.data.origin,cellList{cell_ind});
            cell_indMain = strcmp(objMain.data.origin,cellList{cell_ind});
            
%             if objMain.data.mixsel(cell_indMain) % Plotting only mixed selective cells
            if objMain.data.filtspkcount(cell_indMain) > 100
                disp(['Plotting: ' cell_id]);

                %% Where to save plots
                
                if objMain.data.placesel(cell_indMain) && objMain.data.viewsel(cell_indMain) && objMain.data.headdirectionsel(cell_indMain)
                    figdir2 = [figdir '/3sel/' cell_id];
                    sel = '3sel';
                elseif objMain.data.placesel(cell_indMain) && objMain.data.viewsel(cell_indMain) || ...
                        objMain.data.placesel(cell_indMain) && objMain.data.headdirectionsel(cell_indMain) || ...
                        objMain.data.headdirectionsel(cell_indMain) && objMain.data.viewsel(cell_indMain)
                    if objMain.data.placesel(cell_indMain) && objMain.data.viewsel(cell_indMain)
                        figdir2 = [figdir '/pvsel/' cell_id];
                        sel = 'pvsel';
                    elseif objMain.data.placesel(cell_indMain) && objMain.data.headdirectionsel(cell_indMain)
                        figdir2 = [figdir '/phsel/' cell_id];
                        sel = 'phsel';
                    elseif objMain.data.headdirectionsel(cell_indMain) && objMain.data.viewsel(cell_indMain)
                        figdir2 = [figdir '/hvsel/' cell_id];
                        sel = 'hvsel';
                    end
                elseif objMain.data.placesel(cell_indMain)
                    figdir2 = [figdir '/placesel/' cell_id];
                    sel = 'place only';
                elseif objMain.data.viewsel(cell_indMain)
                    figdir2 = [figdir '/viewsel/' cell_id];
                    sel = 'view only';
                elseif objMain.data.discard(cell_indMain)
                    figdir2 = [figdir '/discard/' cell_id];
                    sel = 'discard';
                else
                    figdir2 = [figdir '/ns/' cell_id];
                    sel = 'ns';
                end

                %% Plot vmms raw and smoothed maps for each variable pair
                
                for pair = 1:size(spatialvarpairs,2)
                    
                    msvar = spatialvarpairs{pair};
                    msvar_short = [msvar{1}(1) msvar{2}(1)];
                    % Set up
                    h = figure(fig);
                    hold on;
                    figname = horzcat(num2str(fig), '- Raw maps ',msvar_short,': ',cell_id,' ',sel);
                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                    % Plot 
                    for kk = 1:size(msvar,2)
                        basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{kk});
                        for mm = 1:2 % raw or smooth
                            if mm == 1
                                map = basedata.pvmap;
                                type = 'raw';
                            else
                                map = basedata.basemapLsm;
                                type = 'sm';
                            end
                            ax = subplot(2,2,2*(kk-1)+mm);
                            plotmap(map,msvar{kk});
                            % If firing rate or spike count = 0, set to black
                            settoblack(map,msvar{kk});
                            % Scale color map
                            switch msvar{kk}
                                case 'place'
                                    maxC = nanmax(map);
                                    fieldcolor = 'r';
                                case 'view'
                                    maxC = nanmax(map(3:end)); % leave out cue and hint
                                    fieldcolor = [50/256 205/256 50/256];
                                case 'headdirection'
                                    maxC = nanmax(map);
                            end
                            % Patch boundaries of base fields
                            if ~strcmp(msvar{kk},'headdirection')
                                for ff = 1:size(basedata.seccomponents_perbasepx,1)
                                    % Patch base fields
                                    for pp = 1:size(basedata.fieldcoord{ff},1)
                                        % PATCH
                                        if strcmp(msvar{kk},'place')
                                            plotgrid = 3;
                                        else
                                            plotgrid = basedata.gridnum(ff);
                                        end
                                        [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                        patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                                    end
                                end
                            else % For head direction, mark start and end of field
                                for ff = 1:size(basedata.seccomponents_perbasepx,1)
                                    hold on;
                                    section = nan(size(map));
                                    section(basedata.fieldlinbin{ff}) = 1.1*maxC;
                                    plotmap(section,msvar{kk},msvar{kk},'b');
                                end
                            end
                            set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color','none');
                            colormap(ax,cool);
                            ax.Title.String = ['vmms ' msvar_short '- ' type ' ' msvar{kk} ' map'];
                        end
                    end
                    % Save
                    figtitle = horzcat(num2str(fig),'- vmms maps ',msvar_short,': ',cell_id);
                    if save
                        mkdir(figdir2);
                        savefigure(h,figtitle,figdir2);
                    end
                    close(h);
                    fig = fig + 1;
                end
                
                %% Plot base maps (all 3 vars)
                
                % Set up
                h = figure(fig);
                hold on;
                figname = horzcat(num2str(fig),'- Base maps ',msvar_short,': ',cell_id,' ',sel);
                set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                var = {'place','view','headdirection'};
                varobj = {'pv','pv','hv'};
                for oo = 1:size(var,2)
                    basedata = objMain.data.(varobj{oo})(cell_indMain).(var{oo});
                    basemap = basedata.basemapLsm;
                    ax = subplot(1,3,oo);
                    plotmap(basedata.basemapLsm,var{oo});
                    % If firing rate or spike count = 0, set to black
                    settoblack(basedata.basemapLsm,var{oo});
                    maxC = nanmax(basemap);
                    % Rate and SIC 
                    rate = text(ax,1,1.1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',30,'HorizontalAlignment','right');
                    si = text(ax,0,1.1,1,num2str(basedata.SI,2),'Units','Normalized','FontSize',30,'HorizontalAlignment','left');
                    if basedata.SI > basedata.SIthr
                        si.Color = 'r';
                    end
                    set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                        'XColor','none','YColor','none','ZColor','none',...
                        'FontSize',14,'GridLineStyle','none','Color','none');
                    titlestring = text(ax,0.5,1.2,1,[' Base ' varobj{oo} ' sm map ' num2str(basedata.sigfields) ' sig fields'],...
                        'Units','Normalized','FontSize',30,'HorizontalAlignment','center');
                    % Draw env boundaries
                    if ~strcmp(var{oo},'headdirection')
                        patchenvbounds(var{oo});
                    end
                end
                % Save
                figtitle = horzcat(num2str(fig),'- Base maps sm ',': ',cell_id,' ',sel);
                if save
                    savefigure(h,figtitle,figdir2);
                end
                close(h);
                fig = fig + 1;
                
                
                
                %% Plot pixel maps for each base field
                
                for pair = 1:size(spatialvarpairs,2)
                    msvar = spatialvarpairs{pair};
                    msvar_short = [msvar{1}(1) msvar{2}(1)];

                    for oo = 1:size(msvar,2)
                        basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                        secdata = objMain.data.(msvar_short)(cell_indMain).(msvar{2-oo+1});
                        if strcmp(msvar{oo},'headdirection') || strcmp(msvar{2-oo+1},'headdirection')
                            disp('stop')
                        end
                        
                        for ff = 1:size(basedata.seccomponents_perbasepx,1) % In mixselpv we've only extracted data for max of 3 fields per map
                                % Plot raw base maps with secondary pixel maps
                                h = figure(fig);
                                figname = '';
                                set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
                                hold on;
                                ax = subplot(1,3,1);
                                % Plot secondary full pixel set
                                tempsecmap = nan(size(basedata.secmapLsm));
                                tempsecmap(basedata.secmaps_raw{ff,1}(:,1)) = basedata.secmaps_raw{ff,1}(:,4);
                                plotmap(tempsecmap,msvar{2-oo+1},msvar{oo});
    %                             colormap(ax,cool);
                                % Patch secondary fields
                                if ~strcmp(msvar{2-oo+1},'headdirection')
                                    for mm = 1:size(secdata.seccomponents_perbasepx,1) % For each secondary field
                                        % Plot secondary field in consolidated pixel map
                                        sec = secdata.fieldcoord{mm}; % Get bins for sec field
                                        % Patch secondary fields in the consolidated pixel map
                                        for pp = 1:size(sec,1)
                                            % PATCH
                                            if strcmp(msvar{2-oo+1},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = secdata.gridnum(mm);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
                                            patch(x,y,z,'r','EdgeColor',[169/256 169/256 169/256],'FaceColor','none','LineWidth',1); % Dark Grey
                                            if strcmp(msvar{2-oo+1},'place')
                                                patch(x,y,z,'r','EdgeColor','r','FaceColor','none','LineWidth',2); % Red
                                            elseif strcmp(msvar{2-oo+1},'view')
                                                patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green
                                            end
                                        end
                                    end
                                else
                                    for mm = 1:size(secdata.seccomponents_perbasepx,1) % For each secondary field
                                        hold on;
                                        section = nan(size(tempsecmap));
                                        section(secdata.fieldlinbin{mm}) = 1.1*nanmax(tempsecmap);
                                        plotmap(section,msvar{2-oo+1},msvar{oo},'b');
                                    end
                                end
                                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                    'XColor','none','YColor','none','ZColor','none',...
                                    'FontSize',14,'GridLineStyle','none','Color','none');
                                % Patch env bounds
                                if sum(strcmp(msvar,'headdirection')) == 0
                                    patchenvbounds('view');
                                elseif sum(strcmp(msvar,'place')) == 1
                                    patchenvbounds('place');
                                elseif sum(strcmp(msvar,'view'))
                                    patchenvbounds('view');
                                end
                                % If firing rate or spike count = 0, set to black
                                settoblack(tempsecmap,msvar{2-oo+1});
                                % Patch basemap field
                                if oo == 2 && strcmp(msvar{oo},'headdirection')
                                    disp('stop')
                                end
                                if ~strcmp(msvar{oo},'headdirection')
                                    for pp = 1:size(basedata.fieldcoord{ff},1)
                                        % PATCH
                                        if strcmp(msvar{oo},'place')
                                            plotgrid = 3;
                                        else
                                            plotgrid = basedata.gridnum(ff);
                                        end
                                        [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                        if strcmp(msvar{oo},'place')
                                            patch(x,y,z,[1 1 1 1],'EdgeColor','r','FaceColor','none','LineWidth',1);
                                        elseif strcmp(msvar{oo},'view')
                                            patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                        end
                                    end
                                else
                                    hold on;
                                    section = nan(size(basedata.basemapLrw));
                                    section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                    plotmap(section,msvar{oo},msvar{2-oo+1},'b');
                                end
                                % Strech axis back to normal if hv and h as base
                                if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                    axis equal;
                                else
                                    disp('stop')
                                end
                                ax.Title.String = cell(size(basedata.sec_infieldrates{ff},1),1);
                                linked = [];
                                for ss = 1:size(basedata.sec_infieldrates{ff},1)
                                    % Check if infield firing in sec field is higher than outfield firing
                                    if basedata.sec_infieldrates{ff}(ss,1) > prctile(basedata.sec_outfieldrates{ff}{ss,1},95)
                                        linked = horzcat(linked,',',num2str(ss));
                                        ax.Title.String{ss,1} = horzcat('LINKED:',msvar{2-oo+1},' Field ',num2str(ss),' Infield: ', num2str(basedata.sec_infieldrates{ff}(ss,1),2),'Hz; Outfield: ', num2str(prctile(basedata.sec_outfieldrates{ff}{ss,1},95),2),'Hz');
                                    else
                                        ax.Title.String{ss,1} = horzcat(msvar{2-oo+1},' Field ',num2str(ss),' Infield: ', num2str(basedata.sec_infieldrates{ff}(ss,1),2),'Hz; Outfield: ', num2str(prctile(basedata.sec_outfieldrates{ff}{ss,1},95),2),'Hz');
                                    end
                                end
                                ax.Title.String{end+1,1} = ['Rawpx = ' num2str(sum(~isnan(tempsecmap)))];
                                if isempty(linked)
                                    figname = horzcat(num2str(fig),'- Pixel maps ', msvar_short, '-', msvar{oo},' field ',num2str(ff),' not linked');
                                else
                                    figname = horzcat(num2str(fig),'- Pixel maps ', msvar_short, '-', msvar{oo},' field ',num2str(ff),' sig linked to ',msvar{2-oo+1},' field ', linked);
                                end

                                % Plot smoothed sec pixel map
                                ax = subplot(1,3,2);
%                                     temp1 = basedata.(['sec' mapname]){ff};
                                temp1 = basedata.(['secmaps_sm' ]){ff};
                                switch msvar{2-oo+1}
                                    case 'place'
                                        maxC = nanmax(temp1);
                                    case 'view'
                                        maxC = nanmax(temp1(3:end));
                                    case 'headdirection'
                                        maxC = nanmax(temp1);
                                end
                                plotmap(temp1,msvar{2-oo+1},msvar{oo});
                                % If firing rate or spike count = 0, set to black
                                settoblack(temp1,msvar{2-oo+1});
                                % Label
                                pseudoSIthr = prctile(basedata.(['pseudosecSIC_adsm' ]){ff},95); 
                                ax.Title.String = {horzcat('SI: ',num2str(basedata.(['secSIC_sm' ])(ff),2),'/',num2str(pseudoSIthr,2)); ...
                                    horzcat('Smoothpx = ',num2str(sum(~isnan(temp1)))); ...
                                    horzcat('Pseudosmoothpx = ', num2str(mean(sum(~isnan(basedata.pseudosecmaps_adsm{ff}),2)),2))};
                                ax.Title.FontSize = 16;
                                if basedata.(['secSIC_sm' ])(ff)>pseudoSIthr && maxC >= 0.7 
                                    ax.Title.Color = 'r';
                                end
                                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                    'XColor','none','YColor','none','ZColor','none',...
                                    'FontSize',14,'GridLineStyle','none','Color','none');
                                % Patch env bounds
                                if sum(strcmp(msvar,'headdirection')) == 0
                                    patchenvbounds('view');
                                elseif sum(strcmp(msvar,'place')) == 1
                                    patchenvbounds('place');
                                elseif sum(strcmp(msvar,'view'))
                                    patchenvbounds('view');
                                end
                                if strcmp(msvar{oo},'headdirection')
                                    disp('stop')
                                end
                                if ~strcmp(msvar{oo},'headdirection')
                                    % Patch basemap field
                                    for pp = 1:size(basedata.fieldcoord{ff},1)
                                        % PATCH
                                        if strcmp(msvar{oo},'place')
                                            plotgrid = 3;
                                        else
                                            plotgrid = basedata.gridnum(ff);
                                        end
                                        [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                        if strcmp(msvar{oo},'place')
                                            patch(x,y,z,[1 1 1 1],'EdgeColor','r','FaceColor','none','LineWidth',1);
                                        elseif strcmp(msvar{oo},'view')
                                            patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                        end
                                    end
                                else
                                    hold on;
                                    section = nan(size(basedata.basemapLrw));
                                    section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                    plotmap(section,msvar{oo},msvar{2-oo+1},'b');
                                end
                                % Patch fields found from smoothed conditioned maps
                                tempdata = basedata.seccond_fieldcoord{ff};
                                if ~isempty(tempdata) % If there are secondary fields in conditioned maps
                                    overlap = {};
                                    for gg = 1:size(tempdata,1) % For each of the secondary fields
                                        if sum(basedata.seccond_fieldoverlapind{ff}{gg})>0
                                            fieldcolor = 'y';
                                            xxx = find(basedata.seccond_fieldoverlapind{ff}{gg});
                                            for xx = 1:size(xxx,1)
                                                % Does the sec field overlap with any of the base fields from the other var?
                                                overlap{end+1,1} = ['Base field ' num2str(ff) ' X Sec field '...
                                                    num2str(xxx(xx))];
                                            end
                                        else
                                            fieldcolor = 'b';
                                        end
                                    end
                                    ax.Title.String{end+1,1} = [num2str(size(tempdata,1)) ' cond fields - mean coverage ' ...
                                        num2str(mean(basedata.seccond_fieldsizepercent{ff}),2) '% per field'];
                                    ax.Title.String{end+1,1} = cell2mat(overlap);

%                                     % Linear or nonlinear mixed selectivity?
%                                     outfieldrate_basefield = basedata.tertrate_orig(ff,:);
%                                     outfieldrate_pseudofield = basedata.tertrate_pseudo(ff,:);
%                                     [hval,pval,ci,stats] = ttest(outfieldrate_basefield,outfieldrate_pseudofield,'Tail','right');
%                                     if hval~=0
%                                        ax.Title.String{end+1,1} = ['Linear pv: Outfield firing greater from ' ...
%                                            msobj{oo} ' field than random field'];
%                                     else
%                                        ax.Title.String{end+1,1} = ['Nonlinear pv'];
%                                     end
%                                     ax.Title.String{end+1,1} = ['h=' num2str(hval) 't=' num2str(stats.tstat,2) 'p=' num2str(pval,2)];

                                else 
                                    ax.Title.String{end+1,1} = ['No cond fields'];
                                end
                                % Strech axis back to normal if hv and h as base
                                if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                    axis equal;
                                end

                                % Plot predicted sec pixel map
                                ax = subplot(1,3,3);
                                temp2 = basedata.secmaps_dist{ff};
                                switch msvar{2-oo+1}
                                    case 'place'
                                        maxC = nanmax(temp2);
                                    case 'view'
                                        maxC = nanmax(temp2(3:end));
                                    case 'headdirection'
                                        maxC = nanmax(temp2);
                                end
                                plotmap(temp2,msvar{2-oo+1},msvar{oo});
                                % If firing rate or spike count = 0, set to black
                                settoblack(temp2,msvar{2-oo+1});
                                % Label
                                pseudoSIthr = prctile(basedata.(['pseudosecSIC_adsm' ]){ff},95); 
                                ax.Title.String = {['Predicted ' msvar{2-oo+1} ' map']; ...
                                    ['dist ratio_' basedata.varname_short ' = ' num2str(basedata.base_distratio(ff))];...
                                    ['dist ratio_' secdata.varname_short ' = ' num2str(basedata.sec_distratio(ff))]};
                                ax.Title.FontSize = 16;
                                if basedata.(['secSIC_sm' ])(ff)>pseudoSIthr && maxC >= 0.7 
                                    ax.Title.Color = 'r';
                                end
                                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                    'XColor','none','YColor','none','ZColor','none',...
                                    'FontSize',14,'GridLineStyle','none','Color','none');
                                % Patch env bounds
                                if sum(strcmp(msvar,'headdirection')) == 0
                                    patchenvbounds('view');
                                elseif sum(strcmp(msvar,'place')) == 1
                                    patchenvbounds('place');
                                elseif sum(strcmp(msvar,'view'))
                                    patchenvbounds('view');
                                end
                                % Patch base field
                                if ~strcmp(msvar{oo},'headdirection')
                                    % Patch basemap field
                                    for pp = 1:size(basedata.fieldcoord{ff},1)
                                        % PATCH
                                        if strcmp(msvar{oo},'place')
                                            plotgrid = 3;
                                        else
                                            plotgrid = basedata.gridnum(ff);
                                        end
                                        [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                        if strcmp(msvar{oo},'place')
                                            patch(x,y,z,[1 1 1 1],'EdgeColor','r','FaceColor','none','LineWidth',1);
                                        elseif strcmp(msvar{oo},'view')
                                            patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                        end
                                    end
                                else
                                    hold on;
                                    section = nan(size(basedata.basemapLrw));
                                    section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                    plotmap(section,msvar{oo},msvar{2-oo+1},'b');
                                end
                                % Strech axis back to normal if hv and h as base
                                if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                    axis equal;
                                end

% %                                 if objMain.data.mixsel(cell_indMain)
%                                 % Linear or nonlinear mixed anova per field
%                                 base_infield = inrate_placeset((ff-1)*500+1:ff*500);
%                                 if strcmp(msobj{oo},'place')
%                                     inplace_field = inrate_placeset((ff-1)*500+1:ff*500);
%                                     outplace_field = outrate_placeset((ff-1)*500+1:ff*500);
%                                     ind_field = randsample(1:size(inrate_viewset,1),500);
%                                     inview_field = inrate_viewset(ind_field);
%                                     outview_field = outrate_viewset(ind_field);
% 
%                                 else
%                                     inview_field = inrate_viewset((ff-1)*500+1:ff*500);
%                                     outview_field = outrate_viewset((ff-1)*500+1:ff*500);
%                                     ind_field = randsample(1:size(inrate_placeset,1),500);
%                                     inplace_field = inrate_placeset(ind_field);
%                                     outplace_field = outrate_placeset(ind_field);
%                                 end
%                                 [pval,tbl,stats] = anova2([inplace_field outplace_field;inview_field outview_field],500);
% %                                 [pval,tbl,stats] = anova2([inplace_field inview_field; outplace_field outview_field],500);
%                                 figanova = gcf;
%                                 figtitle = ['ANOVA ' msobj{oo} ' field ' num2str(ff)];
%                                 if save
%                                     savefigure(figanova,figtitle,figdir2);
%                                 end
%                                 close(figanova);
%                                 if pval(3) < 0.05
%                                     ax.Title.String{end+1,1} = ['Nonlinear interaction: F=' num2str(tbl{4,5},2) ', p=' num2str(pval(3))];
%                                 else
%                                     ax.Title.String{end+1,1} = ['Linear: F=' num2str(tbl{4,5},2) ', p=' num2str(pval(3))];
%                                 end
% %                                 end

                                % Save 
                                figtitle = figname;
                                if save
                                    savefigure(h,figtitle,figdir2);
                                end
                                close(h);
                                fig = fig + 1;

    %                             % Plot a selection of the pseudopopulation of smoothed conditioned maps
    %                             h = figure(fig);
    %                             figname = [msobj{oo} ' field ' num2str(ff) ' - pseudopopulation'];
    %                             set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
    %                             hold on;
    %                             temp = basedata.pseudosecSIC_adsm{ff};
    %                             setpoints = quantile(temp,[0.25,0.5,0.75,0.975])';
    %                             for hh = 1:size(setpoints,1)
    %                                 inds_pseudoplot = find(abs(temp-setpoints(hh)) == min(abs(temp-setpoints(hh))),1);
    %                                 % Plot raw map
    %                                 ax = subplot(2,5,hh);
    %                                 plotmap(basedata.pseudosecdataperfield{ff}{4}(inds_pseudoplot,:),msobj{2-oo+1});
    %                                 colormap(ax,cool);
    %                                 % If firing rate or spike count = 0, set to black
    %                                 settoblack(basedata.pseudosecdataperfield{ff}{4}(inds_pseudoplot,:),msobj{2-oo+1});
    % %                                 % Patch tertiary field 
    % % %                                 for pp = 1:size(basedata.pseudosecdataperfield{ff}{1}{inds_pseudoplot},1)
    % %                                     if strcmp(msobj{oo},'place')
    % %                                         sectemp = nan(1,1600);
    % %                                         fieldcolor = 'r';
    % %                                     elseif strcmp(msobj{oo},'view')
    % %                                         sectemp = nan(1,5122);
    % %                                         fieldcolor = [50/256 205/256 50/256];
    % %                                     end
    % %                                     for pp = 1:size(basedata.pseudosecdataperfield{ff}{1}{inds_pseudoplot},1)
    % %                                         px = basedata.pseudosecdataperfield{ff}{1}{inds_pseudoplot}(pp);
    % %                                         if strcmp(msobj{oo},'place')
    % %                                             g = 1;
    % %                                         elseif strcmp(msobj{oo},'view')
    % %                                             g = findgrid(px,msobj{oo});
    % %                                         end
    % %                                         coords = find(basedata.dummygrid{g} == px);
    % %                                     end
    % %                                     sectemp(basedata.pseudosecdataperfield{ff}{1}{inds_pseudoplot}) = 1;
    % % %                                 end
    %                                 numspk = sum(basedata.pseudosecdataperfield{ff}{3}(inds_pseudoplot,:));
    %                                 ax.Title.String = {'Raw map';[num2str(numspk),' spk']};
    %                                 ax.Title.FontSize = 16;
    %                                 set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                                     'XColor','none','YColor','none','ZColor','none',...
    %                                     'FontSize',14,'GridLineStyle','none','Color','none');
    %                                 patchenvbounds('view');
    %                                 % Plot smoothed map
    %                                 ax = subplot(2,5,hh+5);
    %                                 plotmap(basedata.pseudosecmaps_adsm{ff}(inds_pseudoplot,:),msobj{2-oo+1});
    %                                 colormap(ax,cool);
    %                                 % If firing rate or spike count = 0, set to black
    %                                 settoblack(basedata.pseudosecmaps_adsm{ff}(inds_pseudoplot,:),msobj{2-oo+1});
    %                                 ax.Title.String = {'Smoothed map';...
    %                                     ['SI = ' num2str(basedata.pseudosecSIC_adsm{ff}(inds_pseudoplot,1),2) '/' num2str(pseudoSIthr,2)]};
    %                                 ax.Title.FontSize = 16;
    %                                 if basedata.pseudosecSIC_adsm{ff}(inds_pseudoplot,1) > pseudoSIthr
    %                                     ax.Title.Color = 'r';
    %                                 end
    %                                 set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                                     'XColor','none','YColor','none','ZColor','none',...
    %                                     'FontSize',14,'GridLineStyle','none','Color','none');
    %                                 patchenvbounds('view');
    %                             end
    %                             % Plot actual sec map for comparison
    %                             ax = subplot(2,5,5);
    %                             temp = nan(size(basedata.secmapLsm,2),1);
    %                             temp(basedata.secmaps_raw{ff}(:,1),1) = basedata.secmaps_raw{ff}(:,4);
    %                             plotmap(temp,msobj{2-oo+1});
    %                             colormap(ax,cool);
    %                             % If firing rate or spike count = 0, set to black
    %                             settoblack(temp,msobj{2-oo+1});
    %                             numspk = sum(basedata.secmaps_raw{ff}(:,3));
    %                             ax.Title.String = {'Actual raw map';[num2str(numspk),' spk']};
    %                             ax.Title.FontSize = 16;
    %                             set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                                 'XColor','none','YColor','none','ZColor','none',...
    %                                 'FontSize',14,'GridLineStyle','none','Color','none');
    %                             patchenvbounds('view');
    %                             ax = subplot(2,5,10);
    %                             plotmap(basedata.(['sec' mapname]){ff},msobj{2-oo+1});
    %                             colormap(ax,cool);
    %                             % If firing rate or spike count = 0, set to black
    %                             settoblack(basedata.(['sec' mapname]){ff},msobj{2-oo+1});
    %                             ax.Title.String = {'Actual smoothed map';...
    %                                 ['SI = ' num2str(basedata.(['sec' critname])(ff),2) '/' num2str(pseudoSIthr,2)]};
    %                             ax.Title.FontSize = 16;
    %                             if basedata.(['sec' critname])(ff) > pseudoSIthr
    %                                 ax.Title.Color = 'r';
    %                             end
    %                             set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                                 'XColor','none','YColor','none','ZColor','none',...
    %                                 'FontSize',14,'GridLineStyle','none','Color','none');
    %                             patchenvbounds('view');
    %                             % Save 
    %                             figtitle = figname;
    %                             if save
    %                                 savefigure(h,figtitle,figdir2);
    %                             end
    %                             close(h);
    %                             fig = fig + 1;

    %                             % Patch secondary fields
    %                             for mm = 1:size(secdata.seccomponents_perbasepx,1) % For each secondary field
    %                                 % Plot secondary field in consolidated pixel map
    %                                 sec = secdata.fieldcoord{mm}; % Get bins for sec field
    %                                 % Patch secondary fields in the consolidated pixel map
    %                                 for pp = 1:size(sec,1)
    %                                     % PATCH
    %                                     if strcmp(msobj{2-oo+1},'place')
    %                                         plotgrid = 3;
    %                                     else
    %                                         plotgrid = secdata.gridnum(mm);
    %                                     end
    %                                     [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
    %                                     patch(x,y,z,'r','EdgeColor',[169/256 169/256 169/256],'FaceColor','none','LineWidth',1); % Dark Grey
    %                                     if strcmp(msobj{2-oo+1},'place')
    %                                         patch(x,y,z,'r','EdgeColor','r','FaceColor','none','LineWidth',1); % Red
    %                                     elseif strcmp(msobj{2-oo+1},'view')
    %                                         patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green
    %                                     end
    %                                 end
    %                             end 



                        end


                    end
                    
                end
                
                %% Plot pc/sv vs corr smoothed maps - uncorrected
                % For checking filtering across objects - optional

                for pair = 1:size(spatialvarpairs,2)
                    msvar = spatialvarpairs{pair};
                    msvar_short = [msvar{1}(1) msvar{2}(1)];
                    temp = objCorr.data.(msvar_short);

                    h = figure(fig);
                    hold on;
                    figname = horzcat(num2str(fig),'- vmcorr Smoothed maps ', msvar_short, '-', cell_id,'-',num2str(fig));
                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                    for oo = 1:size(msvar,2)
                        basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                        % Plot maps from pc/v objects
                        ax = subplot(2,2,1+(oo-1)*2);
                        if strcmp(msvar{oo},'place')
                            basemap = objPlace.data.maps_sm(cell_indP,:);
                            maxC = max(basemap);
                            title = ['vmpc smoothed place map: ' sel];
                            fieldcolor = 'r';
                        elseif strcmp(msvar{oo},'view')
                            basemap = objView.data.maps_sm(cell_indV,:);
                            maxC = max(basemap(3:end));
                            title = ['vmsv smoothed view map: ' sel];
                            fieldcolor = [50/256 205/256 50/256];
                        elseif strcmp(msvar{oo},'headdirection')
                            basemap = objHeaddirection.data.maps_sm(cell_indH,:);
                            maxC = nanmax(basemap);
                            title = ['vmhd smoothed headdirection map: ' sel];
                        end
                        plotmap(basemap,msvar{oo});
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
                        % Patch boundaries of base fields
                        if ~strcmp(msvar{oo},'headdirection')
                            for ff = 1:size(basedata.seccomponents_perbasepx,1)
                                % Patch base fields
                                for pp = 1:size(basedata.fieldcoord{ff},1)
                                    % PATCH
                                    if strcmp(msvar{oo},'place')
                                        plotgrid = 3;
                                    else
                                        plotgrid = basedata.gridnum(ff);
                                    end
                                    [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                    patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                                end
                            end
                        else
                            hold on;
                            for ff = 1:size(basedata.seccomponents_perbasepx,1)
                                section = nan(size(basemap));
                                section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                plotmap(section,msvar{oo},msvar{2-oo+1},'b');
                            end
                        end
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = title;
                        % Plot maps from corr objects
                        ax = subplot(2,2,2+(oo-1)*2);
                        basemap = temp(cell_indCorr).(['maps_sm_corr' msvar{oo}(1) 'set'])(1,:);
                        title = ['vmcorr uncorrected smoothed ' msvar{oo} ' map: ' sel];
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                        else 
                            maxC = max(basemap(3:end));
                        end
                        plotmap(basemap,msvar{oo});
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
                        % Patch boundaries of base fields
                        if ~strcmp(msvar{oo},'headdirection')
                            for ff = 1:size(basedata.seccomponents_perbasepx,1)
                                % Patch base fields
                                for pp = 1:size(basedata.fieldcoord{ff},1)
                                    % PATCH
                                    if strcmp(msvar{oo},'place')
                                        plotgrid = 3;
                                    else
                                        plotgrid = basedata.gridnum(ff);
                                    end
                                    [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                    patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1,'EdgeAlpha',1);
                                end
                            end
                        else
                            hold on;
                            for ff = 1:size(basedata.seccomponents_perbasepx,1)
                                section = nan(size(basemap));
                                section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                plotmap(section,msvar{oo},msvar{2-oo+1},'b');
                            end
                        end
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = title;
                    end
                    % Save
                    figtitle = horzcat(num2str(fig),'- vmcorr Smoothed maps ', msvar_short, '-', cell_id,'-',num2str(fig));
                    if save
                        savefigure(h,figtitle,figdir2);
                    end
                    close(h);
                    fig = fig + 1;
                end
                
 
                %% Plot predicted raw maps (vmcorr)
                for pair = 1:size(spatialvarpairs,2)
                    msvar = spatialvarpairs{pair};
                    msvar_short = [msvar{1}(1) msvar{2}(1)];
                    temp = objCorr.data.(msvar_short);
                    if strcmp(msvar_short,'hv')
                        disp('stop')
                    end
                    
                    h = figure(fig);
                    hold on;
                    figname = horzcat(num2str(fig),'- Predicted raw maps ',msvar_short,'-',num2str(fig));
                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                    for oo = 1:size(msvar,2)
                        % Plot vmcorr smoothed maps
                        basemap = temp(cell_indCorr).(['maps_sm_corr' msvar{oo}(1) 'set'])(1,:);
                        title = ['vmcorr smoothed ' msvar{oo}(1) ' map: ' sel];
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                            if strcmp(msvar{oo},'place')
                                fieldcolor = 'r';
                            else
                                fieldcolor = 'b';
                            end
                        else
                            maxC = max(basemap(3:end));
                            fieldcolor = [50/256 205/256 50/256];
                        end
                        ax = subplot(2,3,1+(oo-1)*3);
                        plotmap(basemap,msvar{oo});
                        % If selective, patch field
%                         if pair == 1 || (pair==2 && oo == 1)
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
                        basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                        % Patch boundaries of base fields
                        for ff = 1:size(basedata.seccomponents_perbasepx,1)
                            if ~strcmp(msvar{oo},'headdirection')
                                % Patch base fields
                                for pp = 1:size(basedata.fieldcoord{ff},1)
                                    % PATCH
                                    if strcmp(msvar{oo},'place')
                                        plotgrid = 3;
                                    else
                                        plotgrid = basedata.gridnum(ff);
                                    end
                                    [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                    patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.5,'EdgeAlpha',1);
                                end
                            else
                                hold on;
                                section = nan(size(basemap));
                                section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                plotmap(section,msvar{oo},msvar{oo},fieldcolor);
                            end
                        end
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = ['vmcorr original smoothed ' msvar{oo} ' map: ' sel];

                        % Plot vmcorr raw pv maps
                        ax = subplot(2,3,2+(oo-1)*3);
                        basemap = temp(cell_indCorr).(['maps_raw_corr' msvar{oo}(1) 'set'])(1,:);
                        title = ['vmcorr raw ' msvar{oo} ' map: ' sel];
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                            if strcmp(msvar{oo},'place')
                                fieldcolor = 'r';
                            else
                                fieldcolor = 'b';
                            end
                        else
                            maxC = max(basemap(3:end));
                            fieldcolor = [50/256 205/256 50/256];
                        end
                        plotmap(basemap,msvar{oo});
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
%                         basedata = objMain.data.(msvar{oo})(cell_indMain);
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = ['vmcorr original raw ' msvar{oo} ' map: ' sel];

                        % Plot predicted maps from corr objects
                        ax = subplot(2,3,3+(oo-1)*3);
                        basemap = temp(cell_indCorr).(['maps_dist_' msvar{oo}(1)]);
                        title = {horzcat('vmcorr predicted raw ', msvar{oo}, ' map: ',sel); ...
                                horzcat('distratio_', msvar{oo}(1), ' = ',num2str(temp(cell_indCorr).(['distratio_' msvar{oo}(1)])));...
                                horzcat('distratio_', msvar{2-oo+1}(1), ' = ',num2str(temp(cell_indCorr).(['distratio_' msvar{2-oo+1}(1)])))};
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                        else
                            maxC = max(basemap(3:end));
                        end
                        plotmap(basemap,msvar{oo});
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = title;
                    end
                    % Save
                    figtitle = horzcat(num2str(fig),'- Predicted raw maps ',msvar_short,'-',num2str(fig));
                    if save
                        mkdir(figdir2);
                        savefigure(h,figtitle,figdir2);
                    end
                    close(h);
                    fig = fig + 1;
                end
                
                %% Plot predicted smoothed maps (vmcorr)
                for pair = 1:size(spatialvarpairs,2)
                    msvar = spatialvarpairs{pair};
                    msvar_short = [msvar{1}(1) msvar{2}(1)];
                    temp = objCorr.data.(msvar_short);
                    
                    h = figure(fig);
                    hold on;
                    figname = horzcat(num2str(fig),'- Predicted sm maps ',msvar_short,'-',num2str(fig));
                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                    
                    for oo = 1:size(msvar,2)
                        % Plot vmcorr smoothed original maps
                        ax = subplot(2,3,1+(oo-1)*3);
                        basemap = temp(cell_indCorr).(['maps_sm_corr' msvar{oo}(1) 'set'])(1,:);
                        title = ['vmcorr smoothed original ' msvar{oo}(1) ' map: ' sel];
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                            if strcmp(msvar{oo},'place')
                                fieldcolor = 'r';
                            else
                                fieldcolor = 'b';
                            end
                        else
                            maxC = max(basemap(3:end));
                            fieldcolor = [50/256 205/256 50/256];
                        end
                        plotmap(basemap,msvar{oo});
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
                        basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                        % Patch boundaries of base fields
                        for ff = 1:size(basedata.seccomponents_perbasepx,1)
                            if ~strcmp(msvar{oo},'headdirection')
                                % Patch base fields
                                for pp = 1:size(basedata.fieldcoord{ff},1)
                                    % PATCH
                                    if strcmp(msvar{oo},'place')
                                        plotgrid = 3;
                                    else
                                        plotgrid = basedata.gridnum(ff);
                                    end
                                    [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                    patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.5,'EdgeAlpha',1);
                                end
                            else
                                hold on;
                                section = nan(size(basemap));
                                section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                            end
                        end
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = ['vmcorr original smoothed ' msvar{oo} ' map: ' sel];

                        % Plot vmcorr smoothed original duration
                        ax = subplot(2,3,2+(oo-1)*3);
                        basemap = temp(cell_indCorr).(['dur_sm_corr' msvar{oo}(1)]);
                        basemap(basemap==0) = nan;
                        title = ['vmcorr smoothed original ' msvar{oo}(1) ' duration: ' sel];
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                        else
                            maxC = max(basemap(3:end));
                        end
                        plotmap(basemap,msvar{oo});
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = ['vmcorr original smoothed ' msvar{oo} ' duration: ' sel];

                        % Plot vmcorr smoothed predicted maps
                        ax = subplot(2,3,3+(oo-1)*3);
                        basemap = temp(cell_indCorr).(['maps_dist_' msvar{oo}(1) '_adsm']);
                        title = {['vmcorr predicted smoothed ' msvar{oo}(1) ' map: ' sel]; ...
                                horzcat('distratio_', msvar{oo}(1), '= ',num2str(temp(cell_indCorr).(['distratio_' msvar{oo}(1)])));...
                                horzcat('distratio_', msvar{2-oo+1}(1), ' = ',num2str(temp(cell_indCorr).(['distratio_' msvar{2-oo+1}(1)])))};
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                        else
                            maxC = max(basemap(3:end));
                        end
                            plotmap(basemap,msvar{oo});
                        % If firing rate or spike count = 0, set to black
                        settoblack(basemap,msvar{oo});
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = title;
                    end
                    % Save
                    figtitle = horzcat(num2str(fig),'- Predicted sm maps ',msvar_short,'-',num2str(fig));
                    if save
                        savefigure(h,figtitle,figdir2);
                    end
                    close(h);
                    fig = fig + 1;
                end
                
                
                % Plot linearised differences between actual and predicted
                % maps
                for pair = 1:size(spatialvarpairs,2)
                    msvar = spatialvarpairs{pair};
                    msvar_short = [msvar{1}(1) msvar{2}(1)];
                    temp = objCorr.data.(msvar_short);
                    
                    h = figure(fig);
                    hold on;
                    figname = horzcat(num2str(fig),'- Diff map ', msvar_short, ' - pred (red) vs raw (black)');
                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                    for oo = 1:size(msvar,2)
                        ax = subplot(1,2,oo);
                        % Plot vmcorr predicted raw map
                        basemap = temp(cell_indCorr).(['maps_dist_' msvar{oo}(1)]);
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                        else
                            maxC = max(basemap(3:end));
                        end
                        plot(1:size(basemap,2),log(basemap),'ro');
                        hold on;
                        % Plot vmcorr actual raw map
                        basemap = temp(cell_indCorr).(['maps_raw_corr' msvar{oo}(1) 'set'])(1,:);
                        title = {['vmcorr raw vs predicted ' msvar{oo}(1) ' map: ' sel]; ...
                                ['dist ratio_' msvar{oo}(1) ' = ' num2str(temp(cell_indCorr).(['distratio_' msvar{oo}(1)]))];...
                                ['dist ratio_' msvar{2-oo+1}(1) ' = ' num2str(temp(cell_indCorr).(['distratio_' msvar{2-oo+1}(1)]))]};
                        if ~strcmp(msvar{oo},'view')
                            maxC = max(basemap);
                        else
                            maxC = max(basemap(3:end));
                        end
                        plot(1:size(basemap,2),log(basemap),'ko');
                        hold off;
                        set(ax,'FontSize',14,'GridLineStyle','none','Color','none');
                        colormap(ax,cool);
                        ax.Title.String = title;
                    end
                    % Save
                    figtitle = figname;
                    if save
                        savefigure(h,figtitle,figdir2);
                    end
                    close(h);
                    fig = fig + 1;
                end
                
                % Linear or non linear mixed selectivity?
%                 if objMain.data.placesel(cell_indMain) && objMain.data.viewsel(cell_indMain)
%                 for pair = 1:size(spatialvarpairs,2)
%                     
%                     msvar = spatialvarpairs{pair};
%                     msvar_short = [msvar{1}(1) msvar{2}(1)];
%                     
%                     inrate_baseset = [];
%                     inrate_secset = [];
%                     outrate_baseset = [];
%                     outrate_secset = [];
%                     
%                     for oo = 1:size(msvar,2)
%                         basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
%                         secdata = objMain.data.(msvar_short)(cell_indMain).(msvar{2-oo+1});
% 
%                         for ff = 1:size(basedata.seccomponents_perbasepx,1)
%                             inrate_base = [];
%                             outrate_base = [];
%                             inrate_sec = [];
%                             outrate_sec = [];
%                             for xx = 1:objMain.data.Args.NumShuffles % For each shuffle
%                                 
%                                 basebin_size = length(basedata.fieldlinbin{ff});
%                                 inds = randsample(1:basebin_size,ceil(basebin_size/2));
%                                 indbins_base = basedata.fieldlinbin{ff}(inds);
%                                 inrate_base(xx) = nanmean(nanmean(basedata.basemapLsm(indbins_base)));
% 
%                                 pseudobasebins = basedata.pseudosecdataperfield{ff}{1}{xx};
%                                 nonoverlap = setdiff(pseudobasebins,basedata.fieldlinbin{ff});
%                                 pseudobasebin_size = size(nonoverlap,1);
%                                 inds = randsample(1:pseudobasebin_size,ceil(basebin_size/2));
%                                 indbins_pseudobase = nonoverlap(inds);
%                                 outrate_base(xx) = nanmean(basedata.basemapLsm(indbins_pseudobase));
% %                                 if strcmp(msvar{oo},'place')
% %                                     basebin_size = size(basedata.fieldlinbin{ff},1);
% %                                     inds = randsample(1:basebin_size,ceil(basebin_size/2));
% %                                     indbins_base = basedata.fieldlinbin{ff}(inds);
% %                                     inrate_base(xx) = nanmean(nanmean(basedata.basemapLsm(indbins_base)));
% % 
% %                                     pseudobasebins = basedata.pseudosecdataperfield{ff}{1}{xx};
% %                                     nonoverlap = setdiff(pseudobasebins,basedata.fieldlinbin{ff});
% %                                     pseudobasebin_size = size(nonoverlap,1);
% %                                     inds = randsample(1:pseudobasebin_size,ceil(basebin_size/2));
% %                                     indbins_pseudobase = nonoverlap(inds);
% %                                     outrate_base(xx) = nanmean(basedata.basemapLsm(indbins_pseudobase));
% %                                 elseif strcmp(msvar{oo},'view')
% %                                     basebin_size = ceil(size(basedata.fieldlinbin{ff},1));
% %                                     inds = randsample(1:basebin_size,ceil(basebin_size/2));
% %                                     indbins_base = basedata.fieldlinbin{ff}(inds);
% %                                     inrate_sec(xx) = nanmean(nanmean(basedata.basemapLsm(indbins_base)));
% % 
% %                                     pseudobasebins = basedata.pseudosecdataperfield{ff}{1}{xx};
% %                                     nonoverlap = setdiff(pseudobasebins,basedata.fieldlinbin{ff});
% %                                     pseudobasebin_size = size(nonoverlap,1);
% %                                     inds = randsample(1:pseudobasebin_size,ceil(basebin_size/2));
% %                                     indbins_pseudobase = nonoverlap(inds);
% %                                     outrate_sec(xx) = nanmean(basedata.basemapLsm(indbins_pseudobase));
% %                                 end
% 
%                             end
%                             
%                             inrate_baseset = [inrate_baseset;inrate_base'];
%                             outrate_baseset = [outrate_baseset;outrate_base'];
%                             inrate_secset = [inrate_secset;inrate_sec'];
%                             outrate_secset = [outrate_secset;outrate_sec'];
% 
% %                             if strcmp(msvar{oo},'place')
% %                                 inrate_baseset = [inrate_baseset;inrate_base'];
% %                                 outrate_baseset = [outrate_baseset;outrate_base'];
% %                             else
% %                                 inrate_secset = [inrate_secset;inrate_sec'];
% %                                 outrate_secset = [outrate_secset;outrate_sec'];
% %                             end 
%                         end
%                     end
%                     
%                     minobs = min([size(inrate_baseset,1),size(inrate_secset,1),size(outrate_baseset,1),size(outrate_secset,1)]);
%                     indsplace = randsample(1:size(inrate_baseset,1),minobs);
%                     indview = randsample(1:size(inrate_secset,1),minobs);
%                     inplace = inrate_baseset(indsplace);
%                     inview = inrate_secset(indview);
%                     outplace = outrate_baseset(indsplace);
%                     outview = outrate_secset(indview);
% 
%                     [pval,tbl,stats] = anova2([inplace outplace;inview outview],minobs);
% %                     [pval,tbl,stats] = anova2([inplace inview; outplace outview],minobs);                    
%                     figanova = gcf;
%                     figtitle = [numstr(fig) '- ANOVA ' msvar_short];
%                     if save
%                         savefigure(figanova,figtitle,figdir2);
%                     end
%                     close(figanova);
%                     h = figure(fig);
%                     ax = gca;
%                     ax.Title.String = {['Main p1 = ' num2str(pval(1)),2]; ...
%                         ['Main p2 = ' num2str(pval(2)),2]; ...
%                         ['Inter p = ' num2str(pval(3)),2]; ...
%                         ['Inter F = ' num2str(tbl{4,5},2)]};
%                     % Save 
%                     figtitle = [num2str(fig) '- Cell Stats' msvar_short];
%                     if save
%                         savefigure(h,figtitle,figdir2);
%                     end
%                     close(h);
%                     fig = fig + 1;
%                 end
                
%                 if objMain.data.placesel(cell_indMain) || objMain.data.viewsel(cell_indMain)

                    
      
%                 end
            end
        end
    end
    
    
    
elseif strcmp(objtype,'allprop') % If plotting place, view, corr and mixsel all in one page
    
    msvar = objMain0.data.Args.msobj;
    
    % Get SI/ISE threshold for original maps. ### TEmporary. should come from vmms object
    cells_indDiscardP = ismember(objPlace.data.origin,objMain0.data.origin(objMain0.data.discard));
    cells_indDiscardV = ismember(objView.data.origin,objMain0.data.origin(objMain0.data.discard));
    pcSIset = objPlace.data.(critshname);
    vSIset = objView.data.(critshname);
    cell_numDiscardP = find(cells_indDiscardP);
    cell_numDiscardV = find(cells_indDiscardV);
    for dd = 1:size(cell_numDiscardP,1)
        ind = cell_numDiscardP(dd);
        pcSIset((ind-1)*objPlace.data.Args.NumShuffles+1:ind*objPlace.data.Args.NumShuffles) = nan;
    end
    for dd = 1:size(cell_numDiscardV,1)
        ind = cell_numDiscardV(dd);
        vSIset((ind-1)*objView.data.Args.NumShuffles+1:ind*objView.data.Args.NumShuffles) = nan;
    end
    pcSIthr = prctile([objPlace.data.(critname)(~cells_indDiscardP); pcSIset],95);
    vSIthr = prctile([objView.data.(critname)(~cells_indDiscardV); vSIset],95);
    
    % Getting SI/ISE threshold for pseudopopulation of conditioned maps
    disp('Getting sel threshold for pseudopop of conditioned maps');
    pseudothr_orig = nan(size(msvar));
    pseudothr_corr = nan(size(msvar));
    for oo = 1:size(msvar,2)
        agg_orig = [];
        agg_corr = [];
        for aa = 1:size(objMain0.data.(msvar{oo}),1)
            if isempty(objMain0.data.(msvar{oo})(aa).pseudosecSIC_adsm)
                continue;
            end
            for bb = 1:size(objMain0.data.(msvar{oo})(aa).pseudosecdataperfield,1)
                agg_orig = [agg_orig;objMain0.data.(msvar{oo})(aa).pseudosecSIC_adsm{bb}];
            end
        end
        for aa = 1:size(objMain.data.(msvar{oo}),1)
            if isempty(objMain.data.(msvar{oo})(aa).pseudosecSIC_adsm)
                continue;
            end
            for bb = 1:size(objMain.data.(msvar{oo})(aa).pseudosecdataperfield,1)
                agg_corr = [agg_corr;objMain.data.(msvar{oo})(aa).pseudosecSIC_adsm{bb}];
            end
        end
        pseudothr_orig(oo) = prctile(agg_orig,95);
        pseudothr_corr(oo) = prctile(agg_corr,95);
    end
    
    for ii = 1:size(setsessions,1) % For each session
        cells_indList = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_indList) % For each cell
            
            cell_ind = cells_indList(jj); % ind within the given cell list
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
            if objMain0.data.discard(cell_indMS0)
                continue; % Skip cells with <100 spikes 
            elseif objMain0.data.mixsel(cell_indMS0)
                presel = 'M';
            elseif objMain0.data.placesel(cell_indMS0)
                presel = 'P';
            elseif objMain0.data.viewsel(cell_indMS0)
                presel = 'V';
            else
                presel = 'NS';
            end
            if objMain.data.mixsel(cell_indMS)
                postsel = 'M';
            elseif objMain.data.placesel(cell_indMS)
                postsel = 'P';
            elseif objMain.data.viewsel(cell_indMS)
                postsel = 'V';
            else
                postsel = 'NS';
            end
            
            % Set up plot page: Find out how many sig fields there are
            axnum_col = 8;
            axnum_row = 8;
            if strcmp(presel,'NS') && strcmp(postsel,'NS') % If not selective at all, just plot overall maps
                numfigs = 1;
            else % Plot pixel maps
                numfields = [objMain0.data.place(cell_indMS0).sigfields ...
                        objMain0.data.view(cell_indMS0).sigfields ...
                        objMain.data.place(cell_indMS).sigfields ...
                        objMain.data.view(cell_indMS).sigfields];
                if strcmp(presel,'V') || strcmp(postsel,'V')
                    numfigs = nanmax([numfields(2) numfields(4)]);
                else 
                    numfigs = nanmax([numfields(1) numfields(3)]);
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
                si = text(ax,0,1,1,num2str(objPlace.data.(critname)(cell_indP),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objPlace.data.(critname)(cell_indP) > pcSIthr % objMain0.data.placesel(cell_indMS0)
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
                si = text(ax,0,1,1,num2str(objPlace.data.([critname '1'])(cell_indP),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objPlace.data.([critname '1'])(cell_indP) > pcSIthr
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
                si = text(ax,0,1,1,num2str(objPlace.data.([critname '2'])(cell_indP),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objPlace.data.([critname '2'])(cell_indP) > pcSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '2nd half';
                % Orig view
                ax = axes('Position',[col_left(3) row_bot(2) axwidth axheight*2]); % Full session
                map = objView.data.maps_adsm(cell_indV,:);
                plotmap(map,'view');
                patchenvbounds('view');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objView.data.(critname)(cell_indV),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objView.data.(critname)(cell_indV) > vSIthr
                    si.Color = 'r';
                end
                ax.Title.String = {'Orig View'; [num2str(objMain0.data.view(cell_indMS0).sigfields) ' fields']};
                ax = axes('Position',[col_left(4) row_bot(1) axwidth axheight]); % 1st half
                map = objView.data.maps_adsm1(cell_indV,:);
                plotmap(map,'view');
                patchenvbounds('view');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objView.data.([critname '1'])(cell_indV),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objView.data.([critname '1'])(cell_indV) > vSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '1st half';
                ax = axes('Position',[col_left(4) row_bot(2) axwidth axheight]); % 2nd half
                map = objView.data.maps_adsm2(cell_indV,:);
                plotmap(map,'view');
                patchenvbounds('view');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objView.data.([critname '2'])(cell_indV),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objView.data.([critname '2'])(cell_indV) > vSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '2nd half';
                % Corr Place
                ax = axes('Position',[col_left(5) row_bot(2) axwidth axheight*2]); % Full session
                map = objCorr.data.([mapname '_corrp'])(cell_indCorr,:);
                plotmap(map,'place');
                patchenvbounds('place');
                colorbar off;
                maxC = max(map);
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.([critname '_corrp'])(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.([critname '_corrp'])(cell_indCorr) > pcSIthr
                    si.Color = 'r';
                end
                ax.Title.String = {['Corr Place (' num2str(objCorr.data.llhpick(cell_indCorr)) ')']; ...
                    [num2str(objMain.data.place(cell_indMS).sigfields) ' fields']};
                ax = axes('Position',[col_left(6) row_bot(1) axwidth axheight]); % 1st half
                map = objCorr.data.([mapname '_corrp1'])(cell_indCorr,:);
                plotmap(map,'place');
                patchenvbounds('place');
                colorbar off;
                maxC = max(map);
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.([critname '_corrp1'])(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.([critname '_corrp1'])(cell_indCorr) > pcSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '1st half';
                ax = axes('Position',[col_left(6) row_bot(2) axwidth axheight]); % 2nd half
                map = objCorr.data.([mapname '_corrp2'])(cell_indCorr,:);
                plotmap(map,'place');
                patchenvbounds('place');
                colorbar off;
                maxC = max(map);
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.([critname '_corrp2'])(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.([critname '_corrp2'])(cell_indCorr) > pcSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '2nd half';
                % Corr view
                ax = axes('Position',[col_left(7) row_bot(2) axwidth axheight*2]); % Full session
                map = objCorr.data.([mapname '_corrv'])(cell_indCorr,:);
                plotmap(map,'view');
                patchenvbounds('view');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.([critname '_corrv'])(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.([critname '_corrv'])(cell_indCorr) > vSIthr
                    si.Color = 'r';
                end
                ax.Title.String = {['Corr View (' num2str(objCorr.data.llhpick(cell_indCorr)) ')'];...
                    [num2str(objMain.data.view(cell_indMS).sigfields) ' fields']};
                ax = axes('Position',[col_left(8) row_bot(1) axwidth axheight]); % 1st half
                map = objCorr.data.([mapname '_corrv1'])(cell_indCorr,:);
                plotmap(map,'view');
                patchenvbounds('view');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.([critname '_corrv1'])(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.([critname '_corrv1'])(cell_indCorr) > vSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '1st half';
                ax = axes('Position',[col_left(8) row_bot(2) axwidth axheight]); % 2nd half
                map = objCorr.data.([mapname '_corrv2'])(cell_indCorr,:);
                plotmap(map,'view');
                patchenvbounds('view');
                colorbar off;
                maxC = max(map(3:end));
                axdisplay(ax,maxC);
                rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',14,'HorizontalAlignment','right');
                si = text(ax,0,1,1,num2str(objCorr.data.([critname '_corrv2'])(cell_indCorr),2),'Units','Normalized','FontSize',14,'HorizontalAlignment','left');
                if objCorr.data.([critname '_corrv2'])(cell_indCorr) > vSIthr
                    si.Color = 'r';
                end
                ax.Title.String = '1st half';
                
                % Plot orig pixel maps
                if ~isnan(objMain0.data.(msvar{1})(cell_indMS0).sigfields) && ...
                        ~isnan(objMain0.data.(msvar{2})(cell_indMS0).sigfields)
                    
                    for oo = 1:size(msvar,2)
                        basedata = objMain0.data.(msvar{oo})(cell_indMS0);
                        secdata = objMain0.data.(msvar{2-oo+1})(cell_indMS0);
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
%                                 ax = axes('Position',[col_left(2*(oo-1)+1) row_bot(4) axwidth axheight*2]);
                                ax = axes('Position',[col_left(2*(oo-1)+1) row_bot(4+(ff-1)*2) axwidth axheight*2]);
                            elseif oo == 2
                                ax = axes('Position',[col_left(2*(oo-1)+1) row_bot(4+(ff-1)*2) axwidth axheight*2]);
                            end
                            map = nan(size(basedata.basemapLsm));
                            plotmap(map,msvar{oo});
                            patchenvbounds(msvar{oo});
                            axdisplay(ax,1);
                            colorbar off;
                            if strcmp(msvar{oo},'place')
                                fieldcolor = 'r';
                            elseif strcmp(msvar{oo},'view')
                                fieldcolor = [50/256 205/256 50/256];
                            end
                            % Patch base fields
                            for pp = 1:size(basedata.fieldcoord{ff},1)
                                % PATCH
                                if strcmp(msvar{oo},'place')
                                    plotgrid = 3;
                                else
                                    plotgrid = basedata.gridnum(ff);
                                end
                                [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1);
                            end
                            ax.Title.String = [msvar{oo}(1) num2str(ff)];
                            % Plot sec pixels and fields
                            if isempty(fffcount)
                                % Plot sec pixels for single selectivecells
                                ax = axes('Position',[col_left(2*(oo-1)+2) row_bot(4+2*(ff-1)) axwidth axheight*2]);
%                                 ax = axes('Position',[col_left(2*(oo-1)+2) row_bot(4) axwidth axheight*2]);
                                if strcmp(maptype,'raw')
                                    map = nan(size(basedata.secmapLsm));
                                    map(basedata.secmaps_raw{ff,1}(:,1)) = basedata.secmaps_raw{ff,1}(:,4);
                                else
                                    map = basedata.(['sec' mapname]){ff,1};
                                end
                                switch objtype
                                    case 'place'
                                        maxC = nanmax(map);
                                    case 'view'
                                        maxC = nanmax(map(3:end));
                                end
                                plotmap(map,msvar{2-oo+1});
                                patchenvbounds(msvar{2-oo+1});
                                settoblack(map,msvar{2-oo+1});
                                axdisplay(ax,maxC);
                                colormap(ax,'cool');
                                c = colorbar;
                                c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
                                if ~strcmp(maptype,'raw')
                                    overlap = {};
                                    for dd = 1:size(basedata.cond_fieldoverlap{ff},1)
                                        if basedata.cond_fieldoverlap{ff}(dd)
                                            fieldcolor = 'y';
                                            xxx = find(basedata.seccond_fieldoverlapind{ff}{dd});
                                            xxx(xxx>3) = [];
                                            if isempty(xxx)
                                                continue;
                                            end
                                            for zz = 1:size(xxx,1)
                                                xx = xxx(zz);
                                                overlap{end+1,1} = ['c' msvar{2-oo+1}(1) num2str(dd) ' X ' msvar{2-oo+1}(1) ...
                                                    num2str(xx)];
                                                if ~isempty(secdata.cond_fieldoverlap{xx})
                                                    for ss = 1:size(secdata.cond_fieldoverlap{xx},1)
                                                        if any(secdata.seccond_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.seccond_fieldoverlapind{xx}{ss});
                                                            if yy == ff
                                                                link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'RL' msvar{2-oo+1}(1) num2str(xx)];
                                                            else
                                                                link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                            end
                                                        end
                                                    end
                                                else 
                                                    link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                end
                                            end
                                        else
                                            fieldcolor = 'b';
                                        end
                                        sec = basedata.seccond_fieldcoord{ff}{dd};
                                        for pp = 1:size(sec,1)
                                            % PATCH
                                            if strcmp(msvar{2-oo+1},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = basedata.cond_gridnum{ff}(dd);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
                                            patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                        end
                                    end
                                    % Does conditioned sec field overlap with base field of other var?
                                    ax.Title.FontSize = 16;
                                    ax.Title.String = {['SI=' num2str(basedata.secSIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_orig(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.secSIC_adsm(ff) > pseudothr_orig(oo) && maxC >= 0.7
                                        ax.Title.Color = 'r';
                                    end
                                end
                            else
                                if strcmp(maptype,'raw')
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
                                        map(basedata.secmaps_raw{ff,1}(:,1)) = basedata.secmaps_raw{ff,1}(:,4);
                                        maxC = nanmax(map(3:end));
                                        plotmap(map,msvar{2-oo+1});
                                        patchenvbounds(msvar{2-oo+1});
                                        settoblack(map,msvar{2-oo+1});
                                        axdisplay(ax,maxC);
                                        colormap(ax,'cool');
                                        c = colorbar;
                                        if oo == 1
                                            c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(fff-1)) 0.9/150 axheight*2];
                                        elseif oo == 2
                                            c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
                                        end
%                                         if strcmp(maptype,'raw')
                                        if strcmp(msvar{oo},'view')
                                            fieldcolor = 'r';
                                        elseif strcmp(msvar{oo},'place')
                                            fieldcolor = [50/256 205/256 50/256]; % Lime green
                                        end
                                        % Patch secondary fields
                                        sec = secdata.fieldcoord{fff}; % Get bins for sec field
                                        for pp = 1:size(sec,1)
                                            % PATCH
                                            if strcmp(msvar{2-oo+1},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = secdata.gridnum(fff);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
                                            patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                        end
                                        % Check if fields are reciprocally linked
                                        ax.Title.String = {horzcat(msvar{2-oo+1}(1),num2str(fff),' Infield: ', num2str(basedata.sec_infieldrates{ff}(fff,1),2),'Hz'); ...
                                            horzcat('Outfield: ', num2str(prctile(basedata.sec_outfieldrates{ff}{fff,1},95),2),'Hz')};
                                        if basedata.sec_infieldrates{ff}(fff,1) > prctile(basedata.sec_outfieldrates{ff}{fff,1},95)
                                            if secdata.secfieldrates{fff}(ff,1) > prctile(secdata.secfieldrates_sh{fff}{ff,1},95) % If reciprocal
                                                ax.Title.Color = 'r';
                                                if oo == 1
                                                    link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'RL' msvar{2-oo+1}(1) num2str(fff)];
    %                                                 link_pre{ff,fff} = [msobj{oo}(1) num2str(ff) 'RL' msobj{2-oo+1}(1) num2str(fff)];
                                                end
                                            else 
                                                ax.Title.Color = 'b';
                                                link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(fff)];
    %                                             link_pre{ff,fff} = [msobj{oo}(1) num2str(ff) 'L' msobj{2-oo+1}(1) num2str(fff)];
                                            end
                                        else 
                                            link_pre{end+1,1} = '';
    %                                         link_pre{ff,fff} = '';
                                        end

                                    end
                                else
                                    % Plot sec pixel maps
%                                     if oo == 1
                                        ax = axes('Position',[col_left(2*(oo-1)+2) row_bot(4+2*(ff-1)) axwidth axheight*2]);
%                                     elseif oo == 2 
%                                         ax = axes('Position',[col_left(2*(oo-1)+2) row_bot(4+2*(ff-1)) axwidth axheight*2]); % ???? 
%                                     end
                                    map = basedata.(['sec' mapname]){ff,1};
                                    switch msvar{oo}
                                        case 'place'
                                            maxC = nanmax(map);
                                        case 'view'
                                            maxC = nanmax(map(3:end));
                                    end
                                    plotmap(map,msvar{2-oo+1});
                                    patchenvbounds(msvar{2-oo+1});
                                    settoblack(map,msvar{2-oo+1});
                                    axdisplay(ax,maxC);
                                    colormap(ax,'cool');
                                    c = colorbar;
%                                     if oo == 1
%                                         c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(fff-1)) 0.9/150 axheight*2];
%                                     elseif oo == 2
                                        c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
%                                     end
                                    overlap = {};
                                    % e.g. if place object, ff = base place field num, dd = cond sec view  field num, xx = matching view base field num
                                    for dd = 1:size(basedata.cond_fieldoverlap{ff},1) 
                                        if basedata.cond_fieldoverlap{ff}(dd)
                                            fieldcolor = 'y';
                                            xxx = find(basedata.seccond_fieldoverlapind{ff}{dd});
                                            xxx(xxx>3) = [];
                                            if isempty(xxx)
                                                continue;
                                            end
                                            for zz = 1:size(xxx,1)
                                                xx = xxx(zz);
                                                overlap{end+1,1} = ['c' msvar{2-oo+1}(1) num2str(dd) ' X ' msvar{2-oo+1}(1) ...
                                                    num2str(xx)];
                                                if ~isempty(secdata.cond_fieldoverlap{xx})
                                                    for ss = 1:size(secdata.cond_fieldoverlap{xx},1)
                                                        if any(secdata.seccond_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.seccond_fieldoverlapind{xx}{ss});
                                                            if yy == ff
                                                                link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'RL' msvar{2-oo+1}(1) num2str(xx)];
                                                            else
                                                                link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                            end
                                                        end
                                                    end
                                                else 
                                                    link_pre{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                end
                                            end
                                        else
                                            fieldcolor = 'b';
                                        end
                                        sec = basedata.seccond_fieldcoord{ff}{dd};
                                        for pp = 1:size(sec,1)
                                            % PATCH
                                            if strcmp(msvar{2-oo+1},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = basedata.cond_gridnum{ff}(dd);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
                                            patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                        end
                                    end
                                    % Does conditioned sec field overlap with base field of other var?
                                    ax.Title.FontSize = 16;
                                    ax.Title.String = {['SI=' num2str(basedata.secSIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_orig(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.secSIC_adsm(ff) > pseudothr_orig(oo) && maxC >= 0.7
                                        ax.Title.Color = 'r';
                                    end

                                end
                            end
                        end
                    end
                end

                % Plot corr pixel maps
                if ~isnan(objMain.data.(msvar{1})(cell_indMS).sigfields) && ~isnan(objMain.data.(msvar{2})(cell_indMS).sigfields)
                    for oo =1:size(msvar,2)
                        basedata = objMain.data.(msvar{oo})(cell_indMS);
                        secdata = objMain.data.(msvar{2-oo+1})(cell_indMS);
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
%                                 ax = axes('Position',[col_left(2*(oo-1)+5) row_bot(4) axwidth axheight*2]);
                                ax = axes('Position',[col_left(2*(oo-1)+5) row_bot(4+(ff-1)*2) axwidth axheight*2]);
                            elseif oo == 2
                                ax = axes('Position',[col_left(2*(oo-1)+5) row_bot(4+(ff-1)*2) axwidth axheight*2]);
                            end
                            map = nan(size(basedata.basemapLsm));
                            plotmap(map,msvar{oo});
                            patchenvbounds(msvar{oo});
                            colorbar off;
                            axdisplay(ax,1);
                            if strcmp(msvar{oo},'place')
                                fieldcolor = 'r';
                            elseif strcmp(msvar{oo},'view')
                                fieldcolor = [50/256 205/256 50/256];
                            end
                            % Patch base fields
                            for pp = 1:size(basedata.fieldcoord{ff},1)
                                % PATCH
                                if strcmp(msvar{oo},'place')
                                    plotgrid = 3;
                                else
                                    plotgrid = basedata.gridnum(ff);
                                end
                                [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1);
                            end
                            ax.Title.String = [msvar{oo}(1) num2str(ff)];
                            if isempty(fffcount)
                                % Plot sec pixels for single selectivecells
                                ax = axes('Position',[col_left(2*(oo-1)+6) row_bot(4+2*(ff-1)) axwidth axheight*2]);
                                if strcmp(maptype,'raw')
                                    map = nan(size(basedata.secmapLsm));
                                    map(basedata.secmaps_raw{ff,1}(:,1)) = basedata.secmaps_raw{ff,1}(:,4);
                                else
                                    map = basedata.(['sec' mapname]){ff,1};
                                end
                                switch objtype
                                    case 'place'
                                        maxC = nanmax(map);
                                    case 'view'
                                        maxC = nanmax(map(3:end));
                                end
                                plotmap(map,msvar{2-oo+1});
                                patchenvbounds(msvar{2-oo+1});
                                settoblack(map,msvar{2-oo+1});
                                axdisplay(ax,maxC);
                                colormap(ax,'cool');
                                c = colorbar;
                                c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4+(ff-1)*2) 0.9/150 axheight*2];
                                if ~strcmp(maptype,'raw')
                                    overlap = {};
                                    for dd = 1:size(basedata.cond_fieldoverlap{ff},1)
                                        if basedata.cond_fieldoverlap{ff}(dd)
                                            fieldcolor = 'y';
                                            xxx = find(basedata.seccond_fieldoverlapind{ff}{dd});
                                            xxx(xxx>3) = [];
                                            if isempty(xxx)
                                                continue;
                                            end
                                            for zz = 1:size(xxx,1)
                                                xx = xxx(zz);
                                                overlap{end+1,1} = ['c' msvar{2-oo+1}(1) num2str(dd) ' X ' msvar{2-oo+1}(1) ...
                                                    num2str(xx)];
                                                if ~isempty(secdata.cond_fieldoverlap{xx})
                                                    for ss = 1:size(secdata.cond_fieldoverlap{xx},1)
                                                        if any(secdata.seccond_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.seccond_fieldoverlapind{xx}{ss});
                                                            if yy == ff
                                                                link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'RL' msvar{2-oo+1}(1) num2str(xx)];
                                                            else
                                                                link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                            end
                                                        end
                                                    end
                                                else 
                                                    link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                end
                                            end
                                        else
                                            fieldcolor = 'b';
                                        end
                                        sec = basedata.seccond_fieldcoord{ff}{dd};
                                        for pp = 1:size(sec,1)
                                            % PATCH
                                            if strcmp(msvar{2-oo+1},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = basedata.cond_gridnum{ff}(dd);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
                                            patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                        end
                                    end
                                    % Does conditioned sec field overlap with base field of other var?
                                    ax.Title.FontSize = 16;
                                    ax.Title.String = {['SI=' num2str(basedata.secSIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_corr(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.secSIC_adsm(ff) > pseudothr_corr(oo) && maxC >= 0.7
                                        ax.Title.Color = 'r';
                                    end
                                end
                            else 
                                if strcmp(maptype,'raw')
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
    %                                         if kk == ff
    %                                             ax = axes('Position',[col_left(2*(oo-1)+6) row_bot(4+2*(kk-1)) axwidth axheight*2]);
    %                                         else
    %                                             continue;
    %                                         end
                                        end
%                                         if strcmp(maptype,'raw')
                                            map = nan(size(basedata.secmapLsm));
                                            map(basedata.secmaps_raw{ff,1}(:,1)) = basedata.secmaps_raw{ff,1}(:,4);
%                                         else
%                                             map = basedata.(['sec' mapname]){ff,1};
%                                         end
                                        maxC = nanmax(map(3:end));
                                        plotmap(map,msvar{2-oo+1});
                                        patchenvbounds(msvar{2-oo+1});
                                        settoblack(map,msvar{2-oo+1});
                                        axdisplay(ax,maxC);
                                        colormap(ax,'cool');
                                        c = colorbar;
                                        if oo == 1
                                            c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4+2*(fff-1)) 0.9/150 axheight*2];
                                        elseif oo == 2
                                            c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
                                        end
                                        if strcmp(msvar{oo},'view')
                                            fieldcolor = 'r';
                                        elseif strcmp(msvar{oo},'place')
                                            fieldcolor = [50/256 205/256 50/256]; % Lime green
                                        end
                                        % Patch secondary fields
                                        sec = secdata.fieldcoord{fff}; % Get bins for sec field
                                        for pp = 1:size(sec,1)
                                            % PATCH
                                            if strcmp(msvar{2-oo+1},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = secdata.gridnum(fff);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
                                            patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                        end
                                        % Check if fields are reciprocally linked
                                        ax.Title.String = {horzcat(msvar{2-oo+1}(1),num2str(fff),' Infield: ', num2str(basedata.sec_infieldrates{ff}(fff,1),2),'Hz'); ...
                                            horzcat('Outfield: ', num2str(prctile(basedata.sec_outfieldrates{ff}{fff,1},95),2),'Hz')};
                                        if basedata.sec_infieldrates{ff}(fff,1) > prctile(basedata.sec_outfieldrates{ff}{fff,1},95)
                                            if secdata.sec_infieldrates{fff}(ff,1) > prctile(secdata.sec_outfieldrates{fff}{ff,1},95) % If reciprocal
                                                ax.Title.Color = 'r';
                                                if oo == 1
                                                    link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'RL' msvar{2-oo+1}(1) num2str(fff)];
    %                                                 link_post{ff,fff} = [msobj{oo}(1) num2str(ff) 'RL' msobj{2-oo+1}(1) num2str(fff)];
                                                end
                                            else 
                                                ax.Title.Color = 'b';
                                                link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(fff)];
    %                                             link_post{ff,fff} = [msobj{oo}(1) num2str(ff) 'L' msobj{2-oo+1}(1) num2str(fff)];
                                            end
                                        else 
                                            link_post{end+1,1} = '';
    %                                         link_post{ff,fff} = '';
                                        end
                                    end
                                else
                                    
                                    % Plot sec pixel maps
                                    ax = axes('Position',[col_left(2*(oo-1)+6) row_bot(4+2*(ff-1)) axwidth axheight*2]);
                                    map = basedata.(['sec' mapname]){ff,1};
                                    switch msvar{oo}
                                        case 'place'
                                            maxC = nanmax(map);
                                        case 'view'
                                            maxC = nanmax(map(3:end));
                                    end
                                    plotmap(map,msvar{2-oo+1});
                                    patchenvbounds(msvar{2-oo+1});
                                    settoblack(map,msvar{2-oo+1});
                                    axdisplay(ax,maxC);
                                    colormap(ax,'cool');
                                    c = colorbar;
%                                     if oo == 1
%                                         c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(fff-1)) 0.9/150 axheight*2];
%                                     elseif oo == 2
                                        c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
%                                     end
                                    overlap = {};
                                    for dd = 1:size(basedata.cond_fieldoverlap{ff},1)
                                        if basedata.cond_fieldoverlap{ff}(dd)
                                            fieldcolor = 'y';
                                            xxx = find(basedata.seccond_fieldoverlapind{ff}{dd});
                                            xxx(xxx>3) = [];
                                            if isempty(xxx)
                                                continue;
                                            end
                                            for zz = 1:size(xxx,1)
                                                xx = xxx(zz);
                                                overlap{end+1,1} = ['c' msvar{2-oo+1}(1) num2str(dd) ' X ' msvar{2-oo+1}(1) ...
                                                    num2str(xx)];
                                                if ~isempty(secdata.cond_fieldoverlap{xx})
                                                    for ss = 1:size(secdata.cond_fieldoverlap{xx},1)
                                                        if any(secdata.seccond_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.seccond_fieldoverlapind{xx}{ss});
                                                            if yy == ff
                                                                link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'RL' msvar{2-oo+1}(1) num2str(xx)];
                                                            else
                                                                link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                            end
                                                        end
                                                    end
                                                else 
                                                    link_post{end+1,1} = [msvar{oo}(1) num2str(ff) 'L' msvar{2-oo+1}(1) num2str(xx)];
                                                end
                                            end
                                        else
                                            fieldcolor = 'b';
                                        end
                                        sec = basedata.seccond_fieldcoord{ff}{dd};
                                        for pp = 1:size(sec,1)
                                            % PATCH
                                            if strcmp(msvar{2-oo+1},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = basedata.cond_gridnum{ff}(dd);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
                                            patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',0.1); 
                                        end
                                    end
                                    % Does conditioned sec field overlap with base field of other var?
                                    ax.Title.FontSize = 16;
                                    ax.Title.String = {['SI=' num2str(basedata.secSIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_corr(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.secSIC_adsm(ff) > pseudothr_corr(oo) && maxC >= 0.7
                                        ax.Title.Color = 'r';
                                    end
                                end
                            end
                        end
                    end
                end
                %%% Patch - Figure already named above. Adding here the code for whether fields are linked
                if ~isempty(link_pre) || ~isempty(link_post)
%                     link_all = {link_pre{:};link_post{:}};
                    if sum(contains(link_pre,'R')) > 0 || sum(contains(link_post,'R'))
                        if sum(contains(link_pre,['p' num2str(kk)])) > 0 || sum(contains(link_post,['p' num2str(kk)]))
                            figtitle = [cell_id ' ' presel '-' postsel ' ' num2str(kk) '-RL'];
                        end
                    elseif sum(contains(link_pre,'L')) > 0 || sum(contains(link_post,'L'))
                        if sum(contains(link_pre,['p' num2str(kk)])) > 0 || sum(contains(link_post,['p' num2str(kk)]))
                            figtitle = [cell_id ' ' presel '-' postsel ' ' num2str(kk) '-L'];
                        end
                    end
                    set(h,'Name',figtitle);
                end
                if save
                    savefigure(h,figtitle,figdir);
                end
                close(h);
            end
            if ~isempty(cell2mat(strfind(link_pre,'R')))
                selcell_orig(end+1,1) = cell_ind;
            end 
            if ~isempty(cell2mat(strfind(link_post,'R')))
                selcell_corr(end+1,1) = cell_ind;
            end
        end
    end
 
end

selcell_orig = cellList(selcell_orig);
selcell_corr = cellList(selcell_corr);
% Save list of selective cells 
if save && ~isempty(selcell_orig) 
    cd(figdir);
    fid = fopen('selectivecells_orig.txt','w');
%     fprintf(fid,'%s\n',[num2str(size(cellList,1)) ' cells']);
    for ii = 1:size(selcell_orig,1)
        fprintf(fid,'%s\n',selcell_orig{ii,1});
    end
    fclose(fid);
end
if save && ~isempty(selcell_corr)
    fid = fopen('selectivecells_corr.txt','w');
%     fprintf(fid,'%s\n',[num2str(size(cellList,1)) ' cells']);
    for ii = 1:size(selcell_corr,1)
        fprintf(fid,'%s\n',selcell_corr{ii,1});
    end
    fclose(fid);
end
cd(cwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up axes
function axdisplay(ax,maxC)
% Troubleshoot
if maxC == 0 || isnan(maxC)
    maxC = 1;
end
% Patch boundaries of base fields
set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');


% Patch environment boundaries (origin: placebyview.m)
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
        
    case 'view'
        
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
% print('-painters',figtitle,'-dvg');
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
    case 'view'
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

if size(basemap,1) > size(basemap,2) 
    basemap = basemap';
end
blackpx = find(basemap == 0);
if ~isempty(blackpx) && ~strcmp(baseobj,'headdirection')
    for pp = 1:size(blackpx,2)
        switch baseobj
            case 'place'
                gnum = 1;
                [gnum,xx,yy] = findgrid(blackpx(pp),baseobj);
            case 'view'
                if blackpx(pp) > 3 % Make sure not cue or hint
                   [gnum,xx,yy] = findgrid(blackpx(pp),baseobj);
                else
                    continue;
                end
        end
        [x,y,z] = converttosurf(gnum,xx,yy);
        patch(x,y,z,[1 1 1 1],'FaceColor','k','FaceAlpha',0.7,'EdgeColor','none');
    end
end
