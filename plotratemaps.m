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
pix = 100;

cwd = '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis';
% Saving figure directory
if strcmp(objtype,'mixsel0') || strcmp(objtype,'mixsel1')
    figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/mixsel/UseCorr' objtype(end)];
else
    figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/' objtype '_' criteria '_' num2str(maptype) '_' 'test'];
    % figdir = [cwd '/Figures/' filttype '/' num2str(pix) 'px' '/RateMaps/' objtype '_' criteria '_' num2str(maptype) '_' 'yy1px'];
end

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
else % If plotting a batch of cells
    if save
        if exist(figdir,'dir') ~= 7
            mkdir(figdir);
        else
            rmdir(figdir,'s');
            mkdir(figdir);
        end
    end
    % Load cell list
    cd(cwd);
    fid = fopen([cwd '/cell_list_11.txt'],'rt');
    cellList = textscan(fid,'%s','Delimiter','\n');
    cellList = cellList{1};
    % Make sure no empty cells
    notempty = ~cellfun(@isempty,cellList);
    cellList = cellList(notempty,:);
    % Make sure no duplicates
    cellList = unique(cellList);
    
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
        case 'view'
            objMain = load('c_vmsv.mat');
            objMain = objMain.vms;
        case 'headdirection'
            objMain = load('c_vmhd.mat');
            objMain = objMain.vmd;
        case 'mixsel0'
            tic;
            objMain = load('c_vmms0.mat');
            objMain = objMain.vmms;
            disp(['Loading c_vmms0 obj took ' num2str(toc) 's']); tic;
            objPlace = load('c_vmpc.mat');
            objPlace = objPlace.vmp;
            disp(['Loading c_vmpc obj took ' num2str(toc) 's']); tic;
            objView = load('c_vmsv.mat');
            objView = objView.vms;
            disp(['Loading c_vmsv obj took ' num2str(toc) 's']); tic;
            objHeaddirection = load('c_vmhd.mat');
            objHeaddirection = objHeaddirection.vmd;
            disp(['Loading c_vmhd obj took ' num2str(toc) 's']);
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
        case 'trajectory'
            objMain = load('c_vmpc.mat');
            objMain = objMain.vmp;
    end
    % Load object corrected for independent influences of place/view
    objCorr = load('c_vmcorr.mat');
    objCorr = objCorr.vmcorr;
    cd(cwd);
    % Check that combined object is the same size as cell list
    if size(cellList,1) ~= size(objMain.data.origin,1)
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
    objtype_short = 'p';
elseif strcmp(objtype,'view')
    maps = objMain.data.(mapname);
    objtype_short = 'v';
elseif strcmp(objtype,'headdirection')
    maps = objMain.data.(mapname);
    objtype_short = 'h';
elseif strcmp(objtype,'mixsel0') || strcmp(objtype,'mixsel1')
    maps_p = objPlace.data.(mapname);
    maps_v = objView.data.(mapname);
elseif strcmp(objtype,'allprop')
    maps_p = objPlace.data.(mapname);
    maps_v = objView.data.(mapname);
end

% Load nptdata objMain
nCells = objMain.data.numSets;
% Plot params
plotgridh = 5; % cols
plotgridv = 3; % rows
fig = 1;
subpnum = 1;
% Set up cell counts
crossallthresh = 0;
crosscellthresh = 0;
crosspopthresh = 0;
crosseitherthresh = 0;
selcell_orig = [];
selcell_corr = [];

if ismember(objtype,{'place','view','headdirection'}) % One row per cell

    %%% Pick (up to 5) things to plot
    panels = {'full','half1','half2','corr_pv','pred_hv'}; % 'pred','corr_hv',

    % Get shuffled SIC threshold for all cells
    % thr_sh = [objMain.data.(critname); objMain.data.(critshname)]; % all shuffles of population
    % thr_pop = prctile(thr_sh,95); % thr_place = 0.057; thr_view = ?
    % z_pop = zscore(thr_sh);

    % Get population threshold for selectivity
    cell_indsP = ismember(objMain.data.origin,cellList);
    leaveout_P = find(objMain.data.discard | ~cell_indsP);
    leaveout_P = repmat(1:objMain.data.Args.NumShuffles,size(leaveout_P)) + (leaveout_P*objMain.data.Args.NumShuffles-objMain.data.Args.NumShuffles);
    thr_sh = objMain.data.critsh_sm(setdiff(1:size(objMain.data.critsh_sm,1),leaveout_P));
    thr_pop = prctile(thr_sh,95);
    z_pop = zscore(thr_sh);
    
    for ii = 1:size(setsessions,1) % For each session
        cells_indList = find(identifiers(:,1) == setsessions(ii));

        for jj = 1:length(cells_indList) % For each cell
            
            cell_indList = cells_indList(jj);
            cell_ind = find(strcmp(objMain.data.origin,cellList(cells_indList(jj))));
            okminspk = objMain.data.filtspknum(cell_ind) >= 100;
            if ~okminspk
                disp(['Fewer than 100 spikes: ' cellList(cell_ind)]);
            end

            % Save previous figure
            if jj*5 > plotgridh * plotgridv && mod((jj*5), (plotgridh * plotgridv)) == 5
                % Save previous figure 
                figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
                savefigure(save,h,figtitle,figdir);
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
                    if strcmp(objtype,'place') || strcmp(objtype,'view')
                        mapLincorr = objCorr.data.pv(corr_ind).(['maps_sm' '_corr' objtype_short]);
                    elseif strcmp(objtype,'headdirection')
                        mapLincorr = objCorr.data.ph(corr_ind).(['maps_sm' '_corr' objtype_short]);
                    end
                else
                    mapLincorr = nan(size(mapLin));
                    corr_ind = [];
                end
            end

            % Setup figure
            h = figure(fig);
            hold on;
            figname = horzcat(objtype,': ',num2str(setsessions(ii)));
            set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);

            % Plot full session rate map
            if ismember('full',panels)
                % Find subplot number
                ax = subplot(plotgridv,plotgridh,(mod(jj-1,plotgridv))*plotgridh+find(strcmp('full',panels)));
                % Plot map
                [mapGrid,~,maxrate]= plotmap(mapLin,objtype);
                % Set up axes
                if strcmp(objtype,'view') && maxrate < nanmax(mapLin)
                    ax.Title.String = {horzcat('Cue: ',num2str(mapLin(1),2),'Hz'),...
                        horzcat('Hint: ',num2str(mapLin(2),2),'Hz')};
                    ax.Title.Color = 'r';
                else
                    ax.Title.String = {horzcat(num2str(setsessions(ii)), 'ch',num2str(identifiers(cell_indList,4)),...
                        'c',num2str(identifiers(cell_indList,5)),', ',num2str(maxrate,3),'Hz'), ...
                        horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),...
                        horzcat('z-',criteria,'= ',num2str(z_cell))};
                    if crit >= thr_cell && crit >= thr_pop && okminspk && maxrate>=0.7
                        ax.Title.Color = 'r';
                        crossallthresh = crossallthresh + 1;
                    elseif crit >= thr_cell && crit < thr_pop && okminspk && maxrate>=0.7
                        ax.Title.Color = 'm';
                        crosscellthresh = crosscellthresh + 1;
                        crosseitherthresh = crosseitherthresh + 1;
                    elseif crit < thr_cell && crit >= thr_pop && okminspk && maxrate>=0.7
                        ax.Title.Color = 'b';
                        crosspopthresh = crosspopthresh + 1;
                        crosseitherthresh = crosseitherthresh + 1;
                    end
                    if crit >= thr_pop && okminspk && maxrate>=0.7
                        selcell_orig(end+1,1) = cell_ind;
                    end
                end
                
                % Patch standing point if placebyview
                if nargin > 5 % If placebyview
                        fieldCoords
                        patch([fieldCoords(1)-2 fieldCoords(1)-2 fieldCoords(1)+1 fieldCoords(1)+1],[fieldCoords(2)-2 fieldCoords(2)+1 fieldCoords(2)+1 fieldCoords(2)-2], [0 0 0 0] ,[1 1 1 1],'FaceColor','none');
                end
            end
            
            % Plot session halves
            if ismember('half1',panels) || ismember('half2',panels)
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
                for kk = 1:2
                    if ismember(['half' num2str(kk)],panels)
                        ax = subplot(plotgridv,plotgridh,(mod(jj-1,plotgridv))*plotgridh++find(strcmp(['half' num2str(kk)],panels)));
                        mapLin = eval(['map' num2str(kk)]);
                        crit = eval(['crit' num2str(kk)]);
                        half = num2str(kk);
                        % Plot map
                        [~,~,maxrate_half]= plotmap(mapLin,objtype);
                        % Set up axes
                        ax.Title.String = {horzcat('half',half,': ','corr=',num2str(intracorr,2),', ',...
                            num2str(maxrate_half,3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                        if crit >= thr_cell && crit >= thr_pop && okminspk && maxrate>=0.7
                            ax.Title.Color = 'r';
                        elseif crit >= thr_cell && crit < thr_pop && okminspk && maxrate>=0.7
                            ax.Title.Color = 'm';
                        elseif crit < thr_cell && crit >= thr_pop && okminspk && maxrate>=0.7
                            ax.Title.Color = 'b';
                        end
                    end
                end
            end
            
            % Plot corrected maps 
            if cell2mat(regexp(panels,'corr'))
                howmanycorrpanels = find(~cellfun(@isempty,regexp(panels,'corr')));
                for kk = 1:size(howmanycorrpanels,2)
                    ax = subplot(plotgridv,plotgridh,(mod(jj-1,plotgridv))*plotgridh++howmanycorrpanels(kk));
                    corrobj = strsplit(panels{howmanycorrpanels(kk)},'_');
                    if ~isempty(regexp(corrobj{2},objtype_short))
                        temp = objCorr.data.(corrobj{2});
                        if isempty(corr_ind)
                            critcorr = nan;
                        else
                            critcorr = temp(corr_ind).(['crit_sm' '_corr' objtype_short]);
                        end
                        % Plot map
                        [~,~,maxratecorr]= plotmap(mapLincorr,objtype);
                        % Set up axes
                        if ~isempty(corr_ind)
                            ax.Title.String = {horzcat('Corrected',corrobj{2},temp(corr_ind).llhpicklabel,'of',...
                                num2str(size(temp(corr_ind).llh,1)),': ',num2str(maxratecorr,3),'Hz'), ...
                                horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
                        else
                            ax.Title.String = 'No corrected map';
                        end
                        if critcorr >= thr_cell && critcorr >= thr_pop && okminspk && maxratecorr>=0.7
                                ax.Title.Color = 'r';
                        elseif critcorr >= thr_cell && critcorr < thr_pop && okminspk && maxratecorr>=0.7
                            ax.Title.Color = 'm';
                        elseif critcorr < thr_cell && critcorr >= thr_pop && okminspk && maxratecorr>=0.7
                            ax.Title.Color = 'b';
                        end
                        if critcorr >= thr_pop && okminspk && maxratecorr>=0.7
                            selcell_corr(end+1,1) = cell_ind;
                        end
                    else
                        ax.Title.String = ['Requested corr obj does not have "' objtype '"'];
                    end
                end
            end
            
            
            % Plot predicted artefactual map
            if cell2mat(regexp(panels,'pred'))
                howmanypredpanels = find(~cellfun(@isempty,regexp(panels,'pred')));
                for kk = 1:size(howmanypredpanels,2)
                    ax = subplot(plotgridv,plotgridh,(mod(jj-1,plotgridv))*plotgridh++howmanypredpanels(kk));
                    predobj = strsplit(panels{howmanypredpanels(kk)},'_');
                    if ~isempty(regexp(predobj{2},objtype_short))
                        temp = objCorr.data.(predobj{2});
                        % Plot map
                        distmap = temp(corr_ind).(['maps_dist_' objtype_short]);
                        [~,~,~] = plotmap(distmap,objtype);
                        % Set up axes
                        if ~isempty(corr_ind)
                            ax.Title.String = {['Distributed map ' predobj{2}];...
                                ['dist ratio ' predobj{2}(1) ' = ' num2str(temp(corr_ind).(['distratio_' predobj{2}(1)]))];...
                                ['dist ratio ' predobj{2}(2) ' = ' num2str(temp(corr_ind).(['distratio_' predobj{2}(2)]))]};
                        else
                            ax.Title.String = 'No distributed map';
                        end
                    else 
                        ax.Title.String = ['Requested pred obj does not have "' objtype '"'];
                    end
                end
            end
            
            % Plot covariance matrix
            if ismember('covar',panels)
                % Find subplot number
                ax = subplot(plotgridv,plotgridh,(mod(jj-1,plotgridv))*plotgridh++find(strcmp('covar',panels)));
                % Plot map
                im = imagesc(objCorr.data.covmat_norm{corr_ind});
                set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
                set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
                colormap jet;
                colorbar;
                % Replace NaNs with zeros in covariance matrix for norm calculations
                l1norm = objCorr.data.l1norm(corr_ind);
                l2norm = objCorr.data.l2norm(corr_ind);
                % Set up axes
                if ~isempty(corr_ind)
                    ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(l1norm,2)), horzcat('l2=', num2str(l2norm,2))};
                else
                    ax.Title.String = 'No corrected map';
                end
            end
 
        end

        % Save figure
        figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_indList)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
        savefigure(save,h,figtitle,figdir)
            
        fig = fig + 1;
    end
    disp(['Cross all thresh, ', objtype, ' only = ', num2str(crossallthresh),' cells']);
    disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
    disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);
    disp(['Cross either thresh = ', num2str(crosseitherthresh),' cells']);
    disp(['Total number of cells = ',num2str(size(cellList,1)),' cells']);
    
% elseif strcmp(objtype, 'view') 
% 
%     fig = 1;
%     subpnum = 1;
% 
%     thr_sh = [objMain.data.(critname); objMain.data.(critshname)];
%     thr_pop = prctile(thr_sh,95);
%     z_pop = zscore(thr_sh);
% 
%     for ii = 1:size(setsessions,1) % For each session
% 
%         cells_indList = find(identifiers(:,1) == setsessions(ii));
% 
%         for jj = 1:length(cells_indList) % For each cell
% 
%             cell_indList = cells_indList(jj);
%             cell_ind = find(strcmp(objMain.data.origin,cellList(cells_indList(jj))));
%             cell_indP = find(strcmp(objP.data.origin,cellList(cells_indList(jj))));
%             okminspk = sum(objP.data.spk_raw(cell_indP,:)) >= 100;
%             if ~okminspk
%                 disp(cellList(cell_indList));
%             end
% 
%             %% Plot 1 map for 1 cell
% 
%             % Find figure number
%             if jj*5 > plotgridh * plotgridv && mod((jj*5), (plotgridh * plotgridv)) == 5
%                 % Save previous figure
%                 figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
%                 savefigure(save,h,figtitle,figdir);
%                 fig = fig + 1;
%                 subpnum = 1;
%             end
% 
%             % Get shuffled SI cutoff for this cell - 95th percentile
%             crit = objMain.data.(critname)(cell_ind,1);
%             thr_cell = prctile(objMain.data.(critshname)( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles,1 ),95);
%             z_cell = z_pop(cell_ind,1);
% 
%             if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells)
% 
%                 mapLin = maps(cell_ind,:);
%                 % if corrected map exists, get it
%                 if any(ismember(objMain.data.origin{cell_ind},objCorr.data.origin))
%                     [~,corr_ind] = ismember(objMain.data.origin{cell_ind},objCorr.data.origin);
%                     mapLincorr = objCorr.data.pv(corr_ind).(['maps_sm' '_corr' objtype_short]);
%                 else
%                     mapLincorr = nan(size(mapLin));
%                     corr_ind = [];
%                 end
%                 % Set up figure
%                 h = figure(fig);
%                 ax = subplot(plotgridv,plotgridh,subpnum);
% 
%             end
% 
%             % Set up figure
%             h = gcf;
%             hold on;
%             ax = gca;
% 
%             % Plot map
%             [mapGrid,~]= plotmap(mapLin,objtype);
% 
%             % Figure and axes properties
%             figname = horzcat(objtype,': ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_indList,4)),...
%                 'c',num2str(identifiers(cell_indList,5)));
%             set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
% 
%             % Identify cells sensitive to cue or hint
%             if nanmax(mapLin(1)) > maxC
%                 ax.Title.String = {horzcat('Cue: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_indList,4)),...
%                     'c',num2str(identifiers(cell_indList,5)),', ',num2str(nanmax(mapLin),3),'Hz'),...
%                     horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),...
%                     horzcat('z',num2str(z_cell))};
%             elseif nanmax(mapLin(2)) > maxC
%                 ax.Title.String = {horzcat('Hint: ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_indList,4)),'c',num2str(identifiers(cell_indList,5)),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z',num2str(z_cell))};
%             else
%                 ax.Title.String = {horzcat(num2str(setsessions(ii)),'ch',num2str(identifiers(cell_indList,4)),'c',num2str(identifiers(cell_indList,5)),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z',num2str(z_cell))};
%             end
%             % Patch environment boundaries
%             patchenvbounds(objtype);
% 
%             % Denote if significant spatial information
%             if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
%                 ax.Title.Color = 'r';
%                 crossallthresh = crossallthresh + 1;
%             elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
%                 ax.Title.Color = 'm';
%                 crosscellthresh = crosscellthresh + 1;
%                 crosseitherthresh = crosseitherthresh + 1;
%             elseif crit >= thr_pop && crit < thr_cell && okminspk && maxC>=0.7
%                 ax.Title.Color = 'b';
%                 crosspopthresh = crosspopthresh + 1;
%                 crosseitherthresh = crosseitherthresh + 1;
%             end
%             if crit >= thr_pop && okminspk && maxC>=0.7
%                 selcell_orig(end+1,1) = cell_ind;
%             end
% 
%             if video
%                 h.Units = 'normalized';
%                 h.Position = [0 0 0.5 1];
%                 ax.Position = [0 0 1 1];
%                 ax.Title.Color = 'none';
%                 ax.Color = 'none';
%                 ax.CameraViewAngle = 10;
%                 axis vis3d;
%                 videoname = ['Video ' 'FigNum' num2str(h.Number) ' ' objtype,' ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5)) '.avi'];
%                 v = VideoWriter(videoname);
%                 open(v);
%                 for kstep = 1:360
%                     viewanglemod = [viewangle(1)+(kstep-1) viewangle(2)];
%                     disp(['kstep ' num2str(kstep) ' viewangle ' num2str(viewanglemod(1))]);
%                     view(ax,viewanglemod);
%                     frame = getframe(gcf);
%                     writeVideo(v,frame);
%                 end
%                 close(v);
%             end
% 
%             % Intra-session correlation %%%%%% NOTE: Should use boxcar
%             % smoothed map
%             plothalves = false;
%             if plothalves
%                 map1 = objMain.data.([mapname '1'])(cell_ind,:);
%                 map2 = objMain.data.([mapname '2'])(cell_ind,:);
%     %             map1 = emptyinsidepillar(map1);
%                 vis1 = ~isnan(map1);
%     %             map2 = emptyinsidepillar(map2);
%                 vis2 = ~isnan(map2);
%                 vis = vis1 & vis2; % Correlate only visited bins;
%                 intracorr = corr2(map1(vis), map2(vis));
% 
%                 crit1 = objMain.data.([critname '1'])(cell_ind);
%                 crit2 = objMain.data.([critname '2'])(cell_ind);
% 
%                 % Plot
%                 subpnum = subpnum + 1;
%                 for kk = 1:2
% 
%                     % Get map
%                     if kk == 1
%                         mapLin = map1;
%                         crit = crit1;
%                         half = '1st';
%                     else
%                         mapLin = map2;
%                         crit = crit2;
%                         half = '2nd';
%                     end
% 
%                     % Setup object
%                     h = gcf;
%                     ax = subplot(plotgridv,plotgridh,subpnum);
%                     hold on;
%                     set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
% 
%                     % Plot map
%                     [mapGrid,~]= plotmap(mapLin,objtype);
%                     patchenvbounds(objtype);
% 
%                     ax.Title.String = {horzcat(half,' half: ','corr=',num2str(intracorr,2),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
%                     if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
%                         ax.Title.Color = 'r';
%                     elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
%                         ax.Title.Color = 'm';
%                     elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
%                         ax.Title.Color = 'b';
%                     end
% 
%                     subpnum = subpnum + 1;
%                 end
%             end
% 
%             % Plot corrected maps
%             subpnum = subpnum + 1;
%             critcorrset = [];
%             for cc = 1 % pv & ph
%                 if cc == 1
%                     temp = objCorr.data.pv;
%                     msvar_short = 'pv';
%                 elseif cc == 2
%                     temp = objCorr.data.ph;
%                     msvar_short = 'ph';
%                 end
%                 if isempty(corr_ind)
%                     critcorr = nan;
%                 else
%                     critcorr = temp(corr_ind).(['crit_sm' '_corr' objtype_short]);
%                 end
% 
%                 % Setup object
%                 h = gcf;
%                 ax = subplot(plotgridv,plotgridh,subpnum);
%                 hold on;
% 
%                 % Plot map
%                 [~,~]= plotmap(mapLincorr,objtype);
%                 % Patch environment boundaries
%                 patchenvbounds(objtype);
% 
%                 % Set up axes
%                 if ~isempty(corr_ind)
%                     ax.Title.String = {horzcat('Corrected',msvar_short,temp(corr_ind).llhpicklabel,'of',...
%                         num2str(size(temp(corr_ind).llh,1)),': ',num2str(nanmax(mapLincorr),3),'Hz'), ...
%                         horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
%                 else 
%                     ax.Title.String = 'No corrected map';
%                 end
% 
%                 if critcorr >= thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
%                         ax.Title.Color = 'r';
%                 elseif critcorr >= thr_cell && critcorr < thr_pop && okminspk && maxCcorr>=0.7
%                     ax.Title.Color = 'm';
%                 elseif critcorr < thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
%                     ax.Title.Color = 'b';
%                 end
%                 critcorrset = [critcorrset critcorr];
%                 subpnum = subpnum + 1;
%             end
%             if min(critcorrset) >= thr_pop && okminspk && maxCcorr>=0.7
%                 selcell_corr(end+1,1) = cell_ind;
%             end
% 
%             % Plot predicted artefactual map
% 
%             for cc = 1 % pv and ph
%                 if cc == 1
%                     temp = objCorr.data.pv;
%                     msvar_short = 'pv';
%                 elseif cc == 2
%                     temp = objCorr.data.ph;
%                     msvar_short = 'ph';
%                 end
%                 % Setup object
%                 h = gcf;
%                 ax = subplot(plotgridv,plotgridh,subpnum);
%                 hold on;
% 
%                 % Plot map
%                 distmap = temp(corr_ind).(['maps_dist_' objtype_short]);
%                 [~,~] = plotmap(distmap,objtype);
% 
%                 if ~isempty(corr_ind)
%                     if cc == 1
%                         ax.Title.String = {['Distributed map ' msvar_short];...
%                             ['dist ratio v = ' num2str(temp(corr_ind).distratio_v)];...
%                             ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)]};
%                     elseif cc == 2
%                         ax.Title.String = {['Distributed map ' msvar_short];...
%                             ['dist ratio h = ' num2str(temp(corr_ind).distratio_h)];...
%                             ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)]};
% 
%                     end
%                 else
%                     ax.Title.String = 'No distributed map';
%                 end
% 
%             end
%             hold off;
% 
%             % Patch
%             if mod(subpnum,plotgridh) ~= 0
%                 subpnum = ceil(subpnum/plotgridh)*plotgridh+1;
%             end
% 
% %             % Plot covariance matrix
% %             
% %             % Setup object
% %             h = gcf;
% %             ax = subplot(plotgridv,plotgridh,subpnum);
% %             hold on;
% % 
% %             % Plot map
% %             im = imagesc(objCorr.data.covmat_norm{corr_ind});
% %             set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
% % %             set(ax,'CLim',[-1 1]);
% %             set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
% %             colormap jet;
% %             colorbar;
% %             
% %             % Replace NaNs with zeros in covariance matrix for norm calculations
% %             l1norm = objCorr.data.l1norm(corr_ind);
% %             l2norm = objCorr.data.l2norm(corr_ind);
% % %             covmat = objCorr.data.covmat{corr_ind};
% % %             covmat(isnan(covmat)) = 0;
% % %             % Calculate norms
% % %             l1norm = norm(covmat,1); % maximum of column sum
% % %             l2norm = norm(covmat,2); % maximum single value
% %             
% %             % Set up axes
% %             if ~isempty(corr_ind)
% %                 ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(l1norm,2)), horzcat('l2=', num2str(l2norm,2))};
% %             else
% %                 ax.Title.String = 'No corrected map';
% %             end
% %             
% %             subpnum = subpnum + 1;
% %             
% %             hold off;
% 
%         end
%         % Save figure
%         figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_indList)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
%         savefigure(save,h,figtitle,figdir);
% 
%         fig = fig + 1;
%         subpnum = 1;
% 
%     end
%     disp(['Cross all thresh, ', objtype, ' only = ', num2str(crossallthresh),' cells']);
%     disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
%     disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);
%     disp(['Cross either thresh = ', num2str(crosseitherthresh),' cells']);
%     disp(['Total number of cells = ',num2str(size(cellList,1)),' cells']);
% 
% elseif strcmp(objtype,'headdirection')
% 
%     % For each session, plot rate maps for each cell
%     fig = 1;
%     subpnum = 1;
% 
%     thr_sh = [objMain.data.(critname); objMain.data.(critshname)];
%     thr_pop = prctile(thr_sh,95);
%     z_pop = zscore(thr_sh);
% 
%     for ii = 1:size(setsessions,1) % For each session
% 
%         cells_indList = find(identifiers(:,1) == setsessions(ii));
% 
%         for jj = 1:length(cells_indList) % For each cell
% 
%             cell_indList = cells_indList(jj);
%             cell_ind = find(strcmp(objMain.data.origin,cellList(cells_indList(jj))));
%             okminspk = sum(objMain.data.spk_raw(cell_ind,:)) >= 100;
%             if ~okminspk
%                 disp(cellList(cell_indList));
%             end
% 
%             %% Plot 1 map for 1 cell
% 
%             % Find figure number
%             if jj*5 > plotgridh * plotgridv && mod((jj*5), (plotgridh * plotgridv)) == 5
%                 % Save figure
%                 if save
%                     cwd = pwd;
%                     cd(figdir);
%                     % Save previous figure
%                     figtitle = [num2str(setsessions(ii)) '-' num2str(floor(jj/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
%                     saveas(h,figtitle,'png');
% %                     print('-painters',figtitle,'-dvg');
%                     cd(cwd);
%                     close(figure(fig));
%                 end
% 
%                 fig = fig + 1;
%                 subpnum = 1;
%             end
% 
%             % Get shuffled SI cutoff for this cell - 95th percentile
% 
%             % Patch
%             objMain.data.Args.NumShuffles = 10000;
%             % End patch
%             crit = objMain.data.(critname)(cell_ind,1);
%             thr_cell = prctile(objMain.data.(critshname)( (cell_ind-1)*objMain.data.Args.NumShuffles+1:cell_ind*objMain.data.Args.NumShuffles ,1 ) ,95);
%             z_cell = z_pop(cell_ind,1);
% 
%             % Get map
%             if nargin <= 5 % If mapGrid is not already specified (i.e. if plotting for a batch of cells
%                 mapLin = maps(cell_ind,:);
%                 % if corrected map exists, get it
%                 if any(ismember(objMain.data.origin{cell_ind},objCorr.data.origin))
%                     [~,corr_ind] = ismember(objMain.data.origin{cell_ind},objCorr.data.origin);
%                     mapLincorr = objCorr.data.ph(corr_ind).(['maps_sm' '_corrh']);
%                 else
%                     mapLincorr = nan(size(mapLin));
%                     corr_ind = [];
%                 end
%                 h = figure(fig);
%                 ax = subplot(plotgridv,plotgridh,subpnum);
%             end
% 
%             % Setup main object
%             h = gcf;
%             hold on;
%             figname = horzcat(objtype,': ',num2str(setsessions(ii)));
%             set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
% 
%             % Plot main object
%             [mapGrid,~]= plotmap(mapLin,objtype);
%             maxC = nanmax(mapLin);
% 
%             ax.Title.String = {horzcat(num2str(setsessions(ii)), 'ch',num2str(identifiers(cell_indList,4)),'c',num2str(identifiers(cell_indList,5)),', ',num2str(nanmax(mapLin),3),'Hz'), horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2)),horzcat('z-',criteria,'= ',num2str(z_cell))};
% 
%             if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
%                 ax.Title.Color = 'r';
%                 crossallthresh = crossallthresh + 1;
%             elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
%                 ax.Title.Color = 'm';
%                 crosscellthresh = crosscellthresh + 1;
%                 crosseitherthresh = crosseitherthresh + 1;
%             elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
%                 ax.Title.Color = 'b';
%                 crosspopthresh = crosspopthresh + 1;
%                 crosseitherthresh = crosseitherthresh + 1;
%             else
%                 ax.Title.Color = 'k';
%             end
%             if crit >= thr_pop && okminspk && maxC>=0.7
%                 selcell_orig(end+1,1) = cell_ind;
%             end
% 
%             % Intra-session correlation %%%%%% NOTE: Should use boxcar
%             % smoothed map
%             map1 = objMain.data.([mapname '1'])(cell_ind,:);
%             map2 = objMain.data.([mapname '2'])(cell_ind,:);
% 
%             vis1 = ~isnan(map1);
%             vis2 = ~isnan(map2);
%             vis = vis1 & vis2; % Correlate only visited bins;
%             intracorr = corr2(map1(vis), map2(vis));
% 
%             crit1 = objMain.data.([critname '1'])(cell_ind);
%             crit2 = objMain.data.([critname '2'])(cell_ind);
% 
%             % Plot
%             subpnum = subpnum + 1;
%             for kk = 1:2
% 
%                 % Get map
%                 if kk == 1
%                     mapLin = map1;
%                     crit = crit1;
%                     half = '1st';
%                 else
%                     mapLin = map2;
%                     crit = crit2;
%                     half = '2nd';
%                 end
% 
%                 % Setup object
%                 h = gcf;
%                 ax = subplot(plotgridv,plotgridh,subpnum);
%                 hold on;
%                 set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
% 
%                 % Plot map
%                 [~,~]= plotmap(mapLin,objtype);
%                 maxC = nanmax(mapLin);
% 
%                 ax.Title.String = {horzcat(half,' half: ','corr=',num2str(intracorr,2),', ',num2str(nanmax(mapLin),3),'Hz'),horzcat(criteria, '=',num2str(crit,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
% 
%                 if crit >= thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
%                     ax.Title.Color = 'r';
%                 elseif crit >= thr_cell && crit < thr_pop && okminspk && maxC>=0.7
%                     ax.Title.Color = 'm';
%                 elseif crit < thr_cell && crit >= thr_pop && okminspk && maxC>=0.7
%                     ax.Title.Color = 'b';
%                 else
%                     ax.Title.Color = 'k';
%                 end
% 
%                 subpnum = subpnum + 1;
%             end
% 
%             % Plot corrected map - PLACEHOLDER
% 
% %             % Setup object
% %             h = gcf;
% %             ax = subplot(plotgridv,plotgridh,subpnum);
% %             hold on;
% 
% 
%             critcorrset = [];
%             for cc = 2 % pv and ph
% 
%                 if cc == 1
%                     temp = objCorr.data.pv;
%                     msvar_short = 'pv';
%                 elseif cc == 2
%                     temp = objCorr.data.ph;
%                     msvar_short = 'ph';
%                 end
%                 if isempty(corr_ind)
%                     critcorr = nan;
%                 else
%                     critcorr = temp(corr_ind).(['crit_sm' '_corr' objtype_short]);
%     %                 switch objtype
%     %                     case 'place'
%     %                         critcorr = objCorr.data.([critname '_corrp'])(corr_ind);
%     %                     case 'view'
%     %                         critcorr = objCorr.data.([critname '_corrv'])(corr_ind);
%     %                 end
%                 end
% 
%                 % Setup object
%                 h = gcf;
%                 ax = subplot(plotgridv,plotgridh,subpnum);
%                 hold on;
% 
%                 % Plot map
%                 [~,~]= plotmap(mapLincorr,objtype);
%                 maxCcorr = nanmax(mapLincorr);
% 
%                 % Set up axes
%                 set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%                     'XColor','none','YColor','none','ZColor','none',...
%                     'FontSize',14,'GridLineStyle','none','Color','none');
%                 if ~isempty(corr_ind)
%                     ax.Title.String = {horzcat('Corrected',msvar_short,temp(corr_ind).llhpicklabel,'of',...
%                         num2str(size(temp(corr_ind).llh,1)),': ',num2str(nanmax(mapLincorr),3),'Hz'), ...
%                         horzcat(criteria, '=',num2str(critcorr,2),'/',num2str(thr_cell,2),'/',num2str(thr_pop,2))};
%                     set(ax,'CLim',[0 maxCcorr]);
%                 else
%                     ax.Title.String = 'No corrected map';
%                 end
% 
%                 if critcorr >= thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
%                         ax.Title.Color = 'r';
%                 elseif critcorr >= thr_cell && critcorr < thr_pop && okminspk && maxCcorr>=0.7
%                     ax.Title.Color = 'm';
%                 elseif critcorr < thr_cell && critcorr >= thr_pop && okminspk && maxCcorr>=0.7
%                     ax.Title.Color = 'b';
%                 else
%                     ax.Title.Color = 'k';
%                 end
%                 hold off;
%                 subpnum = subpnum + 1;
%                 critcorrset = [critcorrset critcorr];
%             end
%             if min(critcorrset) >= thr_pop && okminspk && maxCcorr>=0.7
%                 selcell_corr(end+1,1) = cell_ind;
%             end
% 
% 
%             % Plot predicted artefactual map - PLACEHOLDER
% 
% %             % Setup object
% %             h = gcf;
% %             ax = subplot(plotgridv,plotgridh,subpnum);
% %             hold on;
% 
%             for cc = 2 % pv and ph
% 
%                 if cc == 1
%                     temp = objCorr.data.pv;
%                     msvar_short = 'pv';
%                 elseif cc == 2
%                     temp = objCorr.data.ph;
%                     msvar_short = 'ph';
%                 end
%                 % Setup object
%                 h = gcf;
%                 ax = subplot(plotgridv,plotgridh,subpnum);
%                 hold on;
% 
%                 % Plot map
%                 distmap = temp(corr_ind).(['maps_dist_' objtype_short]);
%                 [~,~] = plotmap(distmap,objtype);
%                 maxCcorr = nanmax(distmap);
% 
%                 % Set up axes
%                 set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%                     'XColor','none','YColor','none','ZColor','none',...
%                     'FontSize',14,'GridLineStyle','none','Color','none');
%                 if ~isempty(corr_ind)
%                     if cc == 1
%                         ax.Title.String = {['Distributed map ' msvar_short];...
%                             ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)];...
%                             ['dist ratio v= ' num2str(temp(corr_ind).distratio_v)]};
%                     elseif cc == 2
%                         ax.Title.String = {['Distributed map ' msvar_short];...
%                             ['dist ratio p = ' num2str(temp(corr_ind).distratio_p)];...
%                             ['dist ratio h= ' num2str(temp(corr_ind).distratio_h)]};
%                     end
%                     set(ax,'CLim',[0 maxCcorr]);
%                 else
%                     ax.Title.String = 'No distributed map';
%                 end
% 
%                 hold off;
% 
%                 subpnum = subpnum + 1;
%             end
% 
% %             % Plot covariance matrix
% %             
% %             % Setup object
% %             h = gcf;
% %             ax = subplot(plotgridv,plotgridh,subpnum);
% %             hold on;
% % 
% %             % Plot map
% %             im = imagesc(objCorr.data.covmat_norm{corr_ind});
% %             set(im,'AlphaData',~isnan(objCorr.data.covmat_norm{corr_ind}));
% %             set(ax,'CLim',[-nanstd(nanstd(objCorr.data.covmat_norm{corr_ind})) nanstd(nanstd(objCorr.data.covmat_norm{corr_ind}))]);
% %             colormap jet;
% %             colorbar;
% %             
% %             % Replace NaNs with zeros in covariance matrix for norm calculations
% % %             covmat = objCorr.data.covmat{corr_ind};
% %             l1norm = objCorr.data.l1norm(corr_ind);
% %             l2norm = objCorr.data.l2norm(corr_ind);
% % %             covmat(isnan(covmat)) = 0;
% % %             % Calculate norms
% % %             norml1 = norm(covmat,1); % maximum of column sum
% % %             norml2 = norm(covmat,2); % maximum single value
% %             
% %             % Set up axes
% %             set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
% %                 'XColor','none','YColor','none','ZColor','none',...
% %                 'FontSize',14,'GridLineStyle','none','Color','none');
% %             if ~isempty(corr_ind)
% %                 ax.Title.String = {'Covariance place-view:',horzcat('l1=', num2str(l1norm,2)), horzcat('l2=', num2str(l2norm,2))};
% %             else
% %                 ax.Title.String = 'No corrected map';
% %             end
% %             
% %             subpnum = subpnum + 1;
% %             
% %             hold off;
% 
% 
%         end
%         if save
%             cwd = pwd;
%             cd(figdir);
%             % Save figure
%             figtitle = [num2str(setsessions(ii)) '-' num2str(ceil(length(cells_indList)/(plotgridh * plotgridv))),' FigNum ',num2str(h.Number)];
%             saveas(h,figtitle,'png');
% %             print('-painters',figtitle,'-dvg');
%             cd(cwd);
%             close(figure(fig));
%         end
% 
%         fig = fig + 1;
%         subpnum = 1;
% 
%     end
%     disp(['Cross all thresh, ', objtype, ' only = ', num2str(crossallthresh),' cells']);
%     disp(['Cross cell thresh only = ', num2str(crosscellthresh),' cells']);
%     disp(['Cross population thresh only = ', num2str(crosspopthresh),' cells']);
%     disp(['Cross either thresh = ', num2str(crosseitherthresh),' cells']);
%     disp(['Total number of cells = ',num2str(size(cellList,1)),' cells']);
    
elseif strcmp(objtype,'mixsel0') || strcmp(objtype,'mixsel1')

    %%% Select which plots you want to see (no limit to number selected)
    %%% summary, vmms_rawvsmooth, basemap, pixelmap, checkfilter, predmap,
    %%% linearizeddiff, lnlmixsel, fields
    whattoplot = {'fields'}; 
    
    % Basic settings
    colorset = setcolor('coolwarm64');
    fieldcolor = [255 215 0]./255; % gold
    axcolor = 'none';
    figdirpng = [figdir '/png']; mkdir(figdirpng);
    figdireps = [figdir '/eps']; mkdir(figdireps);
    spatialvarpairs = objMain.data.Args.spatialvarpairs;
    cellcount = 0;

    % Get population threshold for selectivity
    cell_indsP = ismember(objPlace.data.origin,cellList);
    cell_indsV = ismember(objView.data.origin,cellList);
    cell_indsH = ismember(objHeaddirection.data.origin,cellList);
    leaveout_P = find(objPlace.data.discard | ~cell_indsP);
    leaveout_P = repmat(1:objPlace.data.Args.NumShuffles,size(leaveout_P)) + (leaveout_P*objPlace.data.Args.NumShuffles-objPlace.data.Args.NumShuffles);
    leaveout_V = find(objView.data.discard | ~cell_indsV);
    leaveout_V = repmat(1:objView.data.Args.NumShuffles,size(leaveout_V)) + (leaveout_V*objView.data.Args.NumShuffles-objView.data.Args.NumShuffles);
    leaveout_H = find(objHeaddirection.data.discard | ~cell_indsH);
    leaveout_H = repmat(1:objHeaddirection.data.Args.NumShuffles,size(leaveout_H)) + (leaveout_H*objHeaddirection.data.Args.NumShuffles-objHeaddirection.data.Args.NumShuffles);
    pSIthr = prctile(objPlace.data.critsh_sm(setdiff(1:size(objPlace.data.critsh_sm,1),leaveout_P)),95);
    vSIthr = prctile(objView.data.critsh_sm(setdiff(1:size(objView.data.critsh_sm,1),leaveout_V)),95);
    hSIthr = prctile(objHeaddirection.data.critsh_sm(setdiff(1:size(objHeaddirection.data.critsh_sm,1),leaveout_H)),95);

    for ii = 1:size(setsessions,1) % For each session
        cells_indList = find(identifiers(:,1) == setsessions(ii));
        for jj = 1:length(cells_indList) % For each cell
            
            % Get cell index - unsorted
            cellcount = cellcount+1;
            cell_ind = cells_indList(jj); % within cell List
            cell_id = [num2str(identifiers(cell_ind,1)) 'ch' num2str(identifiers(cell_ind,4)) 'c' num2str(identifiers(cell_ind,5))];
            cell_indP = ismember(objPlace.data.origin,cellList{cell_ind}); % within combined vmpc object
            cell_indV = ismember(objView.data.origin,cellList{cell_ind}); % within combined vmsv object
            cell_indH = ismember(objHeaddirection.data.origin,cellList{cell_ind}); % within combined vmhd object
            cell_indCorr = ismember(objCorr.data.origin,cellList{cell_ind}); % within combined vmcorr object
            cell_indMain = ismember(objMain.data.origin,cellList{cell_ind}); % within combined vmms0/1 object
            
            if ~objMain.data.discard(cell_indMain) 

                % Report plot progress
                disp(['Plotting ' num2str(cellcount) ' of ' num2str(size(cellList,1)) ' cells: ' cell_id]);
                fig = 1;
                
                %% Where to save plots
                
                % Cell Selectivity
                pSI = objPlace.data.crit_sm(cell_indP);
                vSI = objView.data.crit_sm(cell_indV);
                hSI = objHeaddirection.data.crit_sm(cell_indH);
                rateokP = objPlace.data.rateok(cell_indP);
                rateokV = objView.data.rateok(cell_indV);
                rateokH = objHeaddirection.data.rateok(cell_indH);
                
                if pSI>pSIthr && vSI>vSIthr && hSI>hSIthr && rateokP && rateokV && rateokH
                    figdir2 = [figdir '/3sel/' cell_id];
                    sel = '3sel';
                elseif pSI>pSIthr && vSI>vSIthr && rateokP && rateokV || ...
                        pSI>pSIthr && hSI>hSIthr && rateokP && rateokH || ...
                        hSI>hSIthr && vSI>vSIthr && rateokH && rateokV
                    if pSI>pSIthr && vSI>vSIthr 
                        figdir2 = [figdir '/pvsel/' cell_id];
                        sel = 'pvsel';
                    elseif pSI>pSIthr && hSI>hSIthr 
                        figdir2 = [figdir '/phsel/' cell_id];
                        sel = 'phsel';
                    elseif hSI>hSIthr && vSI>vSIthr 
                        figdir2 = [figdir '/hvsel/' cell_id];
                        sel = 'hvsel';
                    end
                elseif pSI>pSIthr && rateokP
                    figdir2 = [figdir '/placesel/' cell_id];
                    sel = 'place only';
                elseif vSI>vSIthr && rateokV
                    figdir2 = [figdir '/viewsel/' cell_id];
                    sel = 'view only';
                elseif hSI>hSIthr && rateokH
                    figdir2 = [figdir '/discard/' cell_id];
                    sel = 'discard';
                else
                    figdir2 = [figdir '/ns/' cell_id];
                    sel = 'ns';
                end
                
                %% Plot SFN 2022 figure
                if ismember('summary',whattoplot)
                    disp('Plotting summary ...');
                    varnames = {'place','view','headdirection'};
                    varobj = {'pv','pv','hv'};
                    % Retrieve vmms datasets for each variable infield/outfield firing
                    inrate_set = cell(size(varnames,2),1);
                    outrate_set = cell(size(varnames,2),1);
                    for pair = 1:size(spatialvarpairs,2)
    
                        msvar = spatialvarpairs{pair};
                        msvar_short = [msvar{1}(1) msvar{2}(1)];
    
                        inrate_var = cell(2,1);
                        outrate_var = cell(2,1);
    
                        for oo = 1:size(msvar,2)
                            basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                            secdata = objMain.data.(msvar_short)(cell_indMain).(msvar{2-oo+1});
    
                            inrate_var{oo} = nan(objMain.data.Args.NumShuffles(cell_indMain)*size(basedata.condbase_componentsperpx,1),1);
                            outrate_var{oo} = nan(objMain.data.Args.NumShuffles(cell_indMain)*size(basedata.condbase_componentsperpx,1),1);
    
                            for ff = 1:size(basedata.condbase_componentsperpx,1)
                                for xx = 1:1000 % objMain.data.Args.NumShuffles(cell_indMain) % For each shuffle
    
                                    basebin_size = length(basedata.fieldlinbin{ff});
                                    pseudobasebins = basedata.pseudosecdataperfield{ff}{1}{xx};
                                    nonoverlap = setdiff(pseudobasebins,basedata.fieldlinbin{ff});
                                    pseudobasebin_size = size(nonoverlap,1);
                                    draw = round((min([basebin_size pseudobasebin_size]))/2); % or draw from ceil(basebin_size/2)
    
                                    inds = randsample(1:basebin_size,draw);
                                    indbins_base = basedata.fieldlinbin{ff}(inds);
                                    % draw rates from raw map
                                    inrate_var{oo}((ff-1)*objMain.data.Args.NumShuffles(cell_indMain)+xx) = nanmean(nanmean(basedata.basemapLrw(indbins_base)));
                                    inds = randsample(1:pseudobasebin_size,draw);
                                    indbins_pseudobase = nonoverlap(inds);
                                    % draw rates from raw map
                                    outrate_var{oo}((ff-1)*objMain.data.Args.NumShuffles(cell_indMain)+xx) = nanmean(basedata.basemapLrw(indbins_pseudobase));
    
                                end
                            end
                            ind_var = strcmp(varnames,msvar{oo});
                            inrate_set{ind_var} = [inrate_set{ind_var}; inrate_var{oo}];
                            outrate_set{ind_var} = [outrate_set{ind_var}; outrate_var{oo}];
    
                        end
                    end
                    
                    minobs = nan(size(varnames,2),1);
                    for vv = 1:size(varnames,2)
                        minobs(vv) = size(inrate_set{vv},1);
                    end
                    % Hedge's g for each var
                    hedgesg = nan(size(varnames,2),1);
                    for vv = 1:size(varnames,2)
                        if size(inrate_set{vv},1) > 0 % If var has a field
                            % Hedge's G
                            pooledstd = sqrt(( (minobs(vv)-1)*std(inrate_set{vv})*std(inrate_set{vv})+(minobs(vv)-1)*std(outrate_set{vv})*std(outrate_set{vv}) ) ...
                                /( (minobs(vv)-1)+(minobs(vv)-1) ));
                            hedgesg(vv) = ( mean(inrate_set{vv})-mean(outrate_set{vv}) )/(pooledstd);
                        end
                    end
                    
                    % ANOVA
                    if sum(minobs>0) > 1 % 2-way
                        ind = find(minobs>0);
                        anovavar = '';
                        minobs = min(minobs(ind));
                        anovaset = nan(minobs*2,size(ind,1));
                        for mm = 1:size(ind,1)
                            inds_var = randsample(1:size(inrate_set{ind(mm)},1),minobs);
                            anovavar = [anovavar '_' varnames{ind(mm)}(1)];
                            anovaset(1:minobs,mm) = inrate_set{ind(mm)}(inds_var); % in-field
                            anovaset(minobs+1:end,mm) = outrate_set{ind(mm)}(inds_var); % out-field
                        end
                        [pval,tbl,stats] = anova2(anovaset,minobs);
                        anovaf = tbl{4,5};
                        anovap = tbl{4,6};
                        close all;
                    elseif sum(minobs>0) == 1 % 1 way
                        ind = minobs>0;
                        anovavar = ['_' varnames{ind}(1)];
                        anovaset = [inrate_set{ind} outrate_set{ind}];
                        [p,tbl,stats] = anova1(anovaset);
                        anovaf = tbl{2,5};
                        anovap = p;
                        close all;
                    else
                        anovavar = '';
                        anovaf = nan;
                        anovap = nan;
                    end
                    
                    % Set up plot figures and axes
                    samesizeaxis = true;
                    if samesizeaxis
                        hset = cell(1,4);
                        fignameset = cell(1,4);
                        axlistset = cell(1,4);
                        for fig = 1:4
                            hset{fig} = figure(fig);
                            hold on;
                            fignameset{fig} = horzcat(cell_id,' ',sel, '- SFN 2022 summary figure ',num2str(fig));
                            colpx = 270;
                            % Set up axes
    %                         if fig ~=4 % Fig 4 is small plot for corrected map
                                numrow = 4;
                                maxnumcol = 4;
    %                         else
    %                             numrow = 4;
    %                             maxnumcol = 2;
    %                         end
                            set(hset{fig},'Name',fignameset{fig},'Units','pixels','InnerPosition',[0 0 maxnumcol*colpx numrow*colpx]);
                            numcol = 1;
                            % Place the variables in the correct column
                            condcol = nan(size(varnames));
                            for oo = 1:size(varnames,2)
                                basedata = objMain.data.(varobj{oo})(cell_indMain).(varnames{oo});
                                condcol(oo) = numcol + 1;
                                numcol = numcol + size(basedata.condbase_componentsperpx,1);
                            end
                            axmargin = 0.05;
                            axwidth = (1-2*axmargin)/maxnumcol;
                            axheight = axwidth;
                            plotwidth = axwidth-axmargin;
                            plotheight = axheight-axmargin;
                            % Set up grid of axes
                            axlistset{fig} = cell(numrow,maxnumcol);
                            for nn = 1:numrow
                                for mm = 1:maxnumcol
                                    axlistset{fig}{nn,mm} = subplot(numrow,maxnumcol,(nn-1)*maxnumcol+mm);
                                    axis(axlistset{fig}{nn,mm},'off');
                                end
                            end
                        end
                    else % square figure, adaptive axis size
                        h = figure(fig);
                        hold on;
                        figname = horzcat(cell_id,' ',sel, '- SFN 2022 figure');
                        set(h,'Name',figname,'Units','pixels','InnerPosition',[0 0 1080 1080]);
                        % Set up axes
                        numrow = 4;
                        numcol = 1;
                        varnames = {'place','view','headdirection'};
                        varobj = {'pv','pv','hv'};
                        condcol = nan(size(varnames));
                        for oo = 1:size(varnames,2)
                            basedata = objMain.data.(varobj{oo})(cell_indMain).(varnames{oo});
                            condcol(oo) = numcol + 1;
                            numcol = numcol + size(basedata.condbase_componentsperpx,1);
                        end
                        axmargin = 0.05;
                        axwidth = (1-2*axmargin)/numcol;
                        axheight = (1-2*axmargin)/numrow;
                        plotwidth = axwidth-axmargin;
                        plotheight = axheight-axmargin;
                        % Set up grid of axes
                        axlist = cell(numrow,numcol);
                        for nn = 1:numrow
                            for mm = 1:numcol
                                axlist{nn,mm} = subplot(numrow,numcol,(nn-1)*numcol+mm);
                                axis(axlist{nn,mm},'off');
                            end
                        end
                    end
                    
                    % Plot base maps - 1 set for pixel maps, 1 set for corrected maps
                    varplot = [varnames,'view'];
                    vartitle = {'Place','View','Head direction','Dist View'};
                    for bb = 1:4
                        if bb ~= 1 && bb ~= 4
                            continue;
                        end
                        h = figure(hset{bb});
                        for oo = 1:size(varplot,2)
                            if oo == 4
                                basedata = objCorr.data.pv(cell_indCorr); % artefactual view map for full data
                                basemap = basedata.maps_dist_v;
                            else
                                basedata = objMain.data.(varobj{oo})(cell_indMain).(varplot{oo}); % regular data
                                basemap = basedata.basemapLsm;
                            end
                            % Remove stray view pixels inside pillars
                            if strcmp(varplot{oo},'view') 
                                basemap = emptyinsidepillar(basemap);
                            end
                            ax = subplot(axlistset{bb}{oo,1});
                            plotmap(basemap,varplot{oo});
                            colormap(ax,colorset);
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,varplot{oo});
                            maxC = nanmax(basemap);
                            % Rate and SIC 
                            if oo < 4
                                rate = text(ax,1,1.05,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',16,'HorizontalAlignment','right');
                                si = text(ax,0,1.05,1,num2str(basedata.SI,2),'Units','Normalized','FontSize',16,'HorizontalAlignment','left');
                                if basedata.SI > basedata.SIthr
                                    si.Color = 'r';
                                end
                                if ~isnan(hedgesg(oo))
                                    gtext = text(ax,0.5,-0.05,1,['g=' num2str(hedgesg(oo),2)],'Units','Normalized','FontSize',16,'HorizontalAlignment','center');
                                end
                            else
                                dr = text(ax,0.5,-0.05,1,['dr = ' num2str(basedata.distratio_v,2)],'Units','Normalized','FontSize',16,'HorizontalAlignment','center');
                            end
                            % ANOVA
                            if oo == 1
                                if isnan(anovap)
                                    astring = '';
                                elseif anovap<0.001
                                    astring = ['F' anovavar '=' num2str(anovaf,2) ', p<0.001'];
                                else
                                    astring = ['F' anovavar '=' num2str(anovaf,2) ', p=' num2str(anovap,2)];
                                end
                                anovatext = text(ax,0.5,1.2,1,astring,'Units','Normalized','FontSize',18,'FontAngle','italic','HorizontalAlignment','center');
                            end
                            % Tweak color scale to discount outliers
    %                         set(ax,'CLim',[0 2*nanmean(basemap)],'Units','normalized',...
    %                             'Position',[axmargin axmargin+(numrow-oo)*axheight plotwidth plotheight],...
    %                             'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                             'XColor','none','YColor','none','ZColor','none',...
    %                             'FontSize',14,'GridLineStyle','none','Color',axcolor);
                            set(ax,'CLim',[0 3*nanstd(basemap)+nanmean(basemap)],'Units','normalized',...
                                'Position',[axmargin axmargin+(numrow-oo)*axheight plotwidth plotheight],...
                                'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
                            if oo == 1
                                colstring = text(ax,0.5,1.3,{'Full session'},...
                                    'Units','Normalized','FontSize',20,'FontWeight','bold',...
                                    'HorizontalAlignment','center','VerticalAlignment','bottom');
                            end
                            rowstring = text(ax,-0.2,0.5,{[vartitle{oo}]},...
                                'Units','Normalized','FontSize',20,'FontWeight','bold',...
                                'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
                            colorbar off;
                            % Draw env boundaries
                            if ~strcmp(varplot{oo},'headdirection')
                                patchenvbounds('view');
                            end
                            % Plot field outlines
                            if ~objMain.data.([varplot{oo} 'sel'])(cell_indMain) || oo == 4
                                continue;
                            end
    %                         for ff = 1:size(basedata.condbase_componentsperpx,1) 
    %                             % Patch basemap field % Actually nothing to
    %                             % plot but use this to make size of map
    %                             % equal across subplots
    %                             if ff > 2
    %                                 continue;
    %                             end
    %                             if ~strcmp(varplot{oo},'headdirection')
    %                                 if strcmp(varplot{oo},'place')
    %                                     color = [254 132 132]/255; % fieldcolor; % [254 132 132]/255;
    %                                 else
    %                                     color = [157 194 9]/255; % pistachio
    %                                 end
    %                                 for pp = 1:size(basedata.fieldcoord{ff},1)
    %                                     % PATCH
    %                                     if strcmp(varplot{oo},'place')
    %                                         plotgrid = 3;
    %                                     else
    %                                         plotgrid = basedata.gridnum(ff);
    %                                     end
    %                                     [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
    %                                     if strcmp(varplot{oo},'place')
    %                                         patch(x,y,z,[1 1 1 1],'EdgeColor',color,'FaceColor','none','LineWidth',1);
    %                                     elseif strcmp(varplot{oo},'view')
    %                                         patch(x,y,z,'r','EdgeColor',color,'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
    %                                     end
    %                                 end
    %                             else
    %                                 section = nan(size(basedata.basemapLrw));
    %                                 section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
    %                                 plotmap(section,varplot{oo},varplot{oo},fieldcolor);
    %                             end
    %                         end
                        end
                    end
    
                    % Plot pixel maps per base field
                    for pair = 1:size(spatialvarpairs,2)
                        
                        % Set up vars for objects
                        msvar = spatialvarpairs{pair};
                        msvar_short = [msvar{1}(1) msvar{2}(1)];
                        % Set up vars for mixed selectivity ANOVA
                        inrate_var = cell(2,1);
                        outrate_var = cell(2,1);
                        
                        for oo = 1:size(msvar,2)
                            
                            % Get objects
                            basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                            secdata = objMain.data.(msvar_short)(cell_indMain).(msvar{2-oo+1});
                            % Set up vars for mixed selectivity ANOVA
                            inrate_var{oo} = nan(objMain.data.Args.NumShuffles(cell_indMain)*size(basedata.condbase_componentsperpx,1),1);
                            outrate_var{oo} = nan(objMain.data.Args.NumShuffles(cell_indMain)*size(basedata.condbase_componentsperpx,1),1);
                            
                            for ff = 1:size(basedata.condbase_componentsperpx,1) 
                                
                                    % Get in-field/out-field data for mixsel ANOVA
                                    for xx = 1:1000 % objMain.data.Args.NumShuffles(cell_indMain) % For each shuffle
                                        basebin_size = length(basedata.fieldlinbin{ff});
                                        pseudobasebins = basedata.pseudosecdataperfield{ff}{1}{xx};
                                        nonoverlap = setdiff(pseudobasebins,basedata.fieldlinbin{ff});
                                        pseudobasebin_size = size(nonoverlap,1);
                                        draw = round((min([basebin_size pseudobasebin_size]))/2); % or draw from ceil(basebin_size/2)
    
                                        inds = randsample(1:basebin_size,draw);
                                        indbins_base = basedata.fieldlinbin{ff}(inds);
                                        % draw rates from raw map
                                        inrate_var{oo}((ff-1)*objMain.data.Args.NumShuffles(cell_indMain)+xx) = nanmean(nanmean(basedata.basemapLrw(indbins_base)));
                                        inds = randsample(1:pseudobasebin_size,draw);
                                        indbins_pseudobase = nonoverlap(inds);
                                        % draw rates from raw map
                                        outrate_var{oo}((ff-1)*objMain.data.Args.NumShuffles(cell_indMain)+xx) = nanmean(basedata.basemapLrw(indbins_pseudobase));
                                    end
                                
                                
                                    % Plot the field outline, no pixel data
                                    subrow = find(strcmp(varnames,msvar{oo}));
                                    subcol = condcol(strcmp(varnames,msvar{oo}))+ff-1;
                                    % Get the appropriate figure
                                    figh = ceil(subcol/maxnumcol);
                                    h = figure(hset{figh});
                                    if subcol>4
                                        subcol = rem(subcol,maxnumcol);
                                        if subcol == 0
                                            subcol = 4;
                                        end
                                    end
                                    ax = subplot(axlistset{figh}{subrow,subcol});
                                    hold on;
                                    tempmap = nan(size(basedata.basemapLsm));
                                    plotmap(tempmap,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    colormap(ax,colorset);
                                    % Patch env bounds
                                    patchenvbounds('view');
                                    % Patch basemap field % Actually nothing to
                                    % plot but use this to make size of map
                                    % equal across subplots
                                    if ~strcmp(msvar{oo},'headdirection')
                                        if strcmp(msvar{oo},'place')
                                            color = [254 132 132]/255;
                                        else
                                            color = [157 194 9]/255; % pistachio
                                        end
                                        for pp = 1:size(basedata.fieldcoord{ff},1)
                                            % PATCH
                                            if strcmp(msvar{oo},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = basedata.gridnum(ff);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                            if strcmp(msvar{oo},'place')
                                                patch(x,y,z,[1 1 1 1],'EdgeColor',color,'FaceColor','none','LineWidth',1);
                                            elseif strcmp(msvar{oo},'view')
                                                patch(x,y,z,'r','EdgeColor',color,'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                            end
                                        end
                                    else
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    end
                                    % Strech axis back to normal if hv and h as base
                                    if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                        axis equal;
                                    end
                                    set(ax,'CLim',[0 1],'Units','Normalized',...
                                        'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-subrow)*axheight plotwidth plotheight],...
                                        'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                        'XColor','none','YColor','none','ZColor','none',...
                                        'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                    % Outline in red the field-only plot
                                    dim = [ax.Position(1) ax.Position(2) plotwidth plotheight];
                                    annotation('rectangle',dim,'Color',fieldcolor,'LineWidth',1);
                                    % Title of axis
                                    titlestring = text(ax,0.5,1.05,{[msvar{oo} ' field ' num2str(ff)]},...
                                        'Units','Normalized','FontSize',16,...
                                        'HorizontalAlignment','center','VerticalAlignment','bottom');
                                    colorbar off;
                                    % Global column title
                                    if subrow == 1
                                        colstring = text(ax,0.5,1.3,{'Conditioned on';[msvar{oo} ' field ' num2str(ff)]},...
                                            'Units','Normalized','FontSize',20,'FontWeight','bold',...
                                            'HorizontalAlignment','center','VerticalAlignment','bottom');
                                    end
                                    
                                    % Plot raw base maps with secondary pixel maps
                                    subrow = find(strcmp(varnames,msvar{2-oo+1}));
                                    ax = subplot(axlistset{figh}{subrow,subcol});
                                    hold on;
                                    if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'view')
                                        disp('stop')
                                    end
                                    tempsecmap = nan(size(basedata.secmapLsm));
                                    tempsecmap(basedata.condbase_map_rw{ff,1}(:,1)) = basedata.condbase_map_rw{ff,1}(:,4);
                                    % Remove stray view pixels inside pillars
                                    if strcmp(msvar{2-oo+1},'view')
                                        tempsecmap = emptyinsidepillar(tempsecmap);
                                    end
                                    plotmap(tempsecmap,msvar{2-oo+1},msvar{oo});
                                    colormap(ax,colorset);
                                    % Tweak color scale to discount outliers
    %                                 set(ax,'CLim',[0 max([0.1 2*nanmean(tempsecmap)])],'Units','Normalized',...
    %                                     'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-subrow)*axheight plotwidth plotheight],...
    %                                     'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                                     'XColor','none','YColor','none','ZColor','none',...
    %                                     'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                    set(ax,'CLim',[0 max([0.1 3*nanstd(tempsecmap)+nanmean(tempsecmap)])],'Units','Normalized',...
                                        'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-subrow)*axheight plotwidth plotheight],...
                                        'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                        'XColor','none','YColor','none','ZColor','none',...
                                        'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                    % Patch env bounds
                                    patchenvbounds('view');
                                    % If firing rate or spike count = 0, set to black
                                    settoblack(tempsecmap,msvar{2-oo+1});
                                    % Patch base field
                                    if ~strcmp(msvar{oo},'headdirection')
                                        if strcmp(msvar{oo},'place')
                                            color = [254 132 132]/255;
                                        else
                                            color = [157 194 9]/255; % pistachio
                                        end
                                        for pp = 1:size(basedata.fieldcoord{ff},1)
                                            % PATCH
                                            if strcmp(msvar{oo},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = basedata.gridnum(ff);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                            if strcmp(msvar{oo},'place')
                                                patch(x,y,z,[1 1 1 1],'EdgeColor',color,'FaceColor','none','LineWidth',1);
                                            elseif strcmp(msvar{oo},'view')
                                                patch(x,y,z,'r','EdgeColor',color,'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                            end
                                        end
                                    else
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    end
                                    % Strech axis back to normal if hv and h as base
                                    if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                        axis equal;
                                    end
                                    % Resize colorbar
                                    c = colorbar;
                                    c.Position([1 3]) = [ax.Position(1)+ax.Position(3)+0.005 0.005];
                                    c.Ticks = [c.Ticks(1) c.Ticks(end)];
                                    % Global column title
                                    if subrow == 1
                                        colstring = text(ax,0.5,1.3,{'Conditioned on';[msvar{oo} ' field ' num2str(ff)]},...
                                            'Units','Normalized','FontSize',20,'FontWeight','bold',...
                                            'HorizontalAlignment','center','VerticalAlignment','bottom');
                                    end
    
                                    % Plot predicted sec pixel map
                                    if ~strcmp(msvar{2-oo+1},'view')
                                        continue;
                                    end
                                    subrow = 4; % find(strcmp(var,msvar{2-oo+1}));
                                    ax = subplot(axlistset{figh}{subrow,subcol});
                                    hold on;
                                    temp2 = basedata.secmaps_dist{ff};
                                    % Remove stray view pixels inside pillars
                                    if strcmp(msvar{oo},'view')
                                        temp2 = emptyinsidepillar(temp2);
                                    end
                                    switch msvar{2-oo+1}
                                        case 'place'
                                            maxC = nanmax(temp2);
                                        case 'view'
                                            maxC = nanmax(temp2(3:end));
                                        case 'headdirection'
                                            maxC = nanmax(temp2);
                                    end
    %                                 fieldcolor = 'y';
                                    plotmap(temp2,msvar{2-oo+1},msvar{oo});
                                    colormap(ax,colorset);
                                    % If firing rate or spike count = 0, set to black
                                    settoblack(temp2,msvar{2-oo+1});
                                    % Label
                                    pseudoSIthr = prctile(basedata.(['pseudosecSIC_adsm' ]){ff},95); 
                                    dr = text(ax,0.5,-0.05,1,['dr = ' num2str(basedata.sec_distratio(ff),2)],...
                                        'Units','Normalized','FontSize',16,'HorizontalAlignment','center');
                                    if basedata.(['condbase_SIC_sm' ])(ff)>pseudoSIthr && maxC >= 0.7 
                                        ax.Title.Color = 'r';
                                    end
    %                                 set(ax,'CLim',[0 2*nanmean(temp2)],'Units','Normalized',...
    %                                     'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-subrow)*axheight plotwidth plotheight],...
    %                                     'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                                     'XColor','none','YColor','none','ZColor','none',...
    %                                     'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                    set(ax,'CLim',[0 3*nanstd(temp2)+nanmean(temp2)],'Units','Normalized',...
                                        'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-subrow)*axheight plotwidth plotheight],...
                                        'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                        'XColor','none','YColor','none','ZColor','none',...
                                        'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                    % Patch env bounds
                                    patchenvbounds('view');
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
                                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1);
                                        end
                                    else
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    end
                                    % Strech axis back to normal if hv and h as base
                                    if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                        axis equal;
                                    end
                                    % Resize colorbar
                                    c = colorbar;
                                    c.Position([1 3]) = [ax.Position(1)+ax.Position(3)+0.005 0.005];
                                    c.Ticks = [c.Ticks(1) c.Ticks(end)];
                            end
                        end
                        
                        % Mixsel ANOVA
                        minobs = min([size(inrate_var{1},1),size(inrate_var{2},1),size(outrate_var{1},1),size(outrate_var{2},1)]);
                        
                    end
                    
                    
                    
                    % Plot corrected maps
                    figh = 4;
                    h = figure(hset{figh});
                    for pair = 1:size(spatialvarpairs,2) % 1
                        msvar = spatialvarpairs{pair};
                        msvar_short = [msvar{1}(1) msvar{2}(1)];
                        if strcmp(msvar_short,'hv')
                            disp('stop')
                        end
                        for oo = 1:size(msvar,2)
                            basedata = objCorr.data.(msvar_short)(cell_indCorr);
    %                         secdata = objCorr.data.(msvar_short)(cell_indCorr);
                            % Plot the corrected map
                            subrow = find(strcmp(varnames,msvar{oo}));
                            subcol = pair+1;
                            ax = subplot(axlistset{figh}{subrow,subcol}); 
                            hold on;
                            tempmap = basedata.(['maps_sm_corr' basedata.(['var' num2str(oo)])(1)]);
                            % Remove stray view pixels inside pillars
                            if strcmp(msvar{oo},'view') 
                                tempmap = emptyinsidepillar(tempmap);
                            end
                            plotmap(tempmap,msvar{oo});
                            colormap(ax,colorset);
                            % If firing rate or spike count = 0, set to black
                            settoblack(tempmap,msvar{oo});
                            maxC = nanmax(tempmap);
                            % Rate and SIC 
                            rate = text(ax,1,1.05,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',16,'HorizontalAlignment','right');
                            si = text(ax,0,1.05,1,num2str(basedata.(['crit_sm_corr' basedata.(['var' num2str(oo)])(1)]),2),'Units','Normalized','FontSize',16,'HorizontalAlignment','left');
                            if basedata.(['crit_sm_corr' basedata.(['var' num2str(oo)])(1)]) > objMain.data.(varobj{oo})(cell_indMain).(varplot{oo}).SIthr
                                si.Color = 'r';
                            end
    %                         maptype = text(ax,1,1.05,1,[msvar{oo}],'Units','Normalized','FontSize',16,'HorizontalAlignment','right');
                            % Tweak color scale to discount outliers
    %                         set(ax,'CLim',[0 2*nanmean(tempmap)],'Units','normalized',...
    %                             'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-subrow)*axheight plotwidth plotheight],...
    %                             'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                             'XColor','none','YColor','none','ZColor','none',...
    %                             'FontSize',14,'GridLineStyle','none','Color',axcolor);
                            set(ax,'CLim',[0 3*nanstd(tempmap)+nanmean(tempmap)],'Units','normalized',...
                                'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-subrow)*axheight plotwidth plotheight],...
                                'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
    %                         rowstring = text(ax,-0.2,0.5,{[vartitle{oo}]},...
    %                             'Units','Normalized','FontSize',20,'FontWeight','bold',...
    %                             'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
                            colorbar off;
                            % Draw env boundaries
                            if ~strcmp(msvar{oo},'headdirection')
                                patchenvbounds('view');
                            end
                            % Column title
                            if oo == 1
                                ax = subplot(axlistset{figh}{1,subcol}); 
                                colstring = text(ax,0.5,1.3,{'MLM-corrected';[spatialvarpairs{pair}{1} '-' spatialvarpairs{pair}{2}]},...
                                    'Units','Normalized','FontSize',20,'FontWeight','bold',...
                                    'HorizontalAlignment','center','VerticalAlignment','bottom');
                                set(ax,'Units','normalized',...
                                    'Position',[axmargin+(subcol-1)*axwidth axmargin+(numrow-1)*axheight plotwidth plotheight],...
                                    'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                    'XColor','none','YColor','none','ZColor','none',...
                                    'FontSize',14,'GridLineStyle','none','Color',axcolor);
                            end
                        end
                    end
                    
                    % Save
                    for fig = 1:4
                        figtitle = fignameset{fig};
                        h = figure(hset{fig});
                        if save && (fig<=ceil(numcol/maxnumcol) || fig ==4)
                            savefigure(h,figtitle,figdirpng);
                            print([figdireps '/' figtitle],'-depsc');
                        end
                        close(h);
                    end
                    fig = fig + 1;
                    
                    % Plot large base maps - 1 set for pixel maps, 1 set for corrected maps
    %                 varplot = [varnames,'view'];
    %                 vartitle = {'Place','View','Head direction','Dist View'};
    %                 for bb = 1:2
                        
                        for oo = 1:2 % size(varplot,2)
                            h = figure(fig);
                            figname = varplot{oo};
                            set(h,'Name',figname,'Units','pixels','InnerPosition',[0 0 1080 1080]);
                            if oo == 4
                                basedata = objCorr.data.pv(cell_indCorr); % artefactual view map for full data
                                basemap = basedata.maps_dist_v;
                            else
                                basedata = objMain.data.(varobj{oo})(cell_indMain).(varplot{oo}); % regular data
                                basemap = basedata.basemapLsm;
                            end
                            % Remove stray view pixels inside pillars
                            if strcmp(varplot{oo},'view') 
                                basemap = emptyinsidepillar(basemap);
                            end
    %                         ax = subplot(axlistset{bb}{oo,1});
                            ax = gca;
                            plotmap(basemap,varplot{oo});
                            colormap(ax,colorset);
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,varplot{oo});
                            maxC = nanmax(basemap);
                            % Rate and SIC 
                            if oo < 4
                                rate = text(ax,1,1.05,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',16,'HorizontalAlignment','right');
                                si = text(ax,0,1.05,1,num2str(basedata.SI,2),'Units','Normalized','FontSize',16,'HorizontalAlignment','left');
                                if basedata.SI > basedata.SIthr
                                    si.Color = 'r';
                                end
                                if ~isnan(hedgesg(oo))
                                    gtext = text(ax,0.5,-0.05,1,['g=' num2str(hedgesg(oo),2)],'Units','Normalized','FontSize',16,'HorizontalAlignment','center');
                                end
                            else
                                dr = text(ax,0.5,-0.05,1,['dr = ' num2str(basedata.distratio_v,2)],'Units','Normalized','FontSize',16,'HorizontalAlignment','center');
                            end
    %                         % ANOVA
    %                         if oo == 1
    %                             if isnan(anovap)
    %                                 astring = '';
    %                             elseif anovap<0.001
    %                                 astring = ['F' anovavar '=' num2str(anovaf,2) ', p<0.001'];
    %                             else
    %                                 astring = ['F' anovavar '=' num2str(anovaf,2) ', p=' num2str(anovap,2)];
    %                             end
    %                             anovatext = text(ax,0.5,1.2,1,astring,'Units','Normalized','FontSize',18,'FontAngle','italic','HorizontalAlignment','center');
    %                         end
                            % Tweak color scale to discount outliers
    %                         set(ax,'CLim',[0 2*nanmean(basemap)],'Units','normalized',...
    %                             'Position',[axmargin axmargin+(numrow-oo)*axheight plotwidth plotheight],...
    %                             'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    %                             'XColor','none','YColor','none','ZColor','none',...
    %                             'FontSize',14,'GridLineStyle','none','Color',axcolor);
                            set(ax,'CLim',[0 3*nanstd(basemap)+nanmean(basemap)],'Units','normalized',...
                                'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
                            if oo == 1
                                colstring = text(ax,0.5,1.3,{'Full session'},...
                                    'Units','Normalized','FontSize',20,'FontWeight','bold',...
                                    'HorizontalAlignment','center','VerticalAlignment','bottom');
                            end
                            rowstring = text(ax,-0.2,0.5,{[vartitle{oo}]},...
                                'Units','Normalized','FontSize',20,'FontWeight','bold',...
                                'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
                            colorbar off;
                            % Draw env boundaries
                            if ~strcmp(varplot{oo},'headdirection')
                                patchenvbounds('view');
                            end
                            % Plot field outlines
                            if ~objMain.data.([varplot{oo} 'sel'])(cell_indMain) || oo == 4
                                continue;
                            end
                            for ff = 1:size(basedata.condbase_componentsperpx,1) 
                                % Patch basemap field % Actually nothing to
                                % plot but use this to make size of map
                                % equal across subplots
                                if ff > 2
                                    continue;
                                end
                                if ~strcmp(varplot{oo},'headdirection')
                                    if strcmp(varplot{oo},'place')
                                        color = [254 132 132]/255; % fieldcolor; % [254 132 132]/255;
                                    else
                                        color = [157 194 9]/255; % pistachio
                                    end
                                    for pp = 1:size(basedata.fieldcoord{ff},1)
                                        % PATCH
                                        if strcmp(varplot{oo},'place')
                                            plotgrid = 3;
                                        else
                                            plotgrid = basedata.gridnum(ff);
                                        end
                                        [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                        if strcmp(varplot{oo},'place')
                                            patch(x,y,z,[1 1 1 1],'EdgeColor',color,'FaceColor','none','LineWidth',2);
                                        elseif strcmp(varplot{oo},'view')
                                            patch(x,y,z,'r','EdgeColor',color,'FaceColor','none','LineWidth',3); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                        end
                                    end
                                else
                                    section = nan(size(basedata.basemapLrw));
                                    section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                    plotmap(section,varplot{oo},varplot{oo},fieldcolor);
                                end
                            end
                            figtitle = horzcat(cell_id,' ',sel, '- SFN 2022 figure ',num2str(fig));
                            if save 
                                savefigure(h,figtitle,figdirpng);
                                print([figdireps '/' figtitle],'-depsc');
                            end
                            close(h);
                            fig = fig + 1;
                        end
    %                 end
                end
                
                %% Plot field outlines
                if ismember('fields',whattoplot)
                    disp('Plotting fields ...');
                    for pair = 1:size(spatialvarpairs,2)
                        msvar = spatialvarpairs{pair};
                        msvar_short = [msvar{1}(1) msvar{2}(1)];
                        % Set up
                        h = figure(fig);
                        hold on;
                        figname = horzcat('Fields ',num2str(fig),'-',msvar_short);
                        set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
                        % Plot
                        for kk = 1:size(msvar,2)
                            basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{kk});
                            map = basedata.basemapLsm;
                            maxC = max(map,[],'omitnan');
                            % Plot smoothed maps
                            for cycle = 1:2
                                ax = subplot(2,2,kk+(cycle-1)*2);
                                plotmap(map,msvar{kk});
                                % If firing rate or spike count = 0, set to black
                                settoblack(basedata.basemapLsm,msvar{kk});
                                if cycle == 1
                                    ax.Title.String = sel;
                                    ax.Title.FontSize = 30;
                                    % Rate and SIC 
                                    rate = text(ax,1,1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',30,'HorizontalAlignment','right');
                                    si = text(ax,0,1,1,num2str(eval([msvar_short(kk) 'SI']),2),'Units','Normalized','FontSize',30,'HorizontalAlignment','left');
                                    if eval([msvar_short(kk) 'SI']) > eval([msvar_short(kk) 'SIthr'])
                                        si.Color = 'r';
                                    end
                                elseif cycle == 2
                                    ax.Title.String = [num2str(basedata.sigfields) ' sig fields'];
                                    ax.Title.FontSize = 30;
                                    alpha(ax,0.2);
                                    % Patch fields
                                    color = {'r','m','b'}; % Order of field colors
                                    % Patch boundaries of base fields
                                    if ~strcmp(msvar{kk},'headdirection')
                                        for ff = 1:basedata.sigfields
                                            % Patch base fields
                                            for pp = 1:size(basedata.fieldcoord{ff},1)
                                                % PATCH
                                                if strcmp(msvar{kk},'place')
                                                    plotgrid = 3;
                                                else
                                                    plotgrid = basedata.gridnum(ff);
                                                end
                                                [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                                patch(x,y,z,[1 1 1 1],'EdgeColor',color{ff},'FaceColor','none','LineWidth',1,'EdgeAlpha',1);
                                            end
                                            fieldbintext = text(ax,1,1-(ff-1)*0.1,1,num2str(pp),'Color',color{ff}, 'Units','Normalized','FontSize',20,'HorizontalAlignment','right');
                                        end
                                    else % For head direction, mark start and end of field
                                        for ff = 1:basedata.sigfields
                                            hold on;
                                            section = nan(size(map));
                                            section(basedata.fieldlinbin{ff}) = 1.1*maxC;
                                            plotmap(section,msvar{kk},msvar{kk},color{ff});
                                            fieldbintext = text(ax,1,1-(ff-1)*0.1,1,num2str(pp),'Color',color{ff}, 'Units','Normalized','FontSize',20,'HorizontalAlignment','right');
                                        end
                                    end
                                end
                                % Draw env boundaries
                                if ~strcmp(msvar{kk},'headdirection')
                                    patchenvbounds(msvar{kk});
                                end
                            end
                        end
                        % Save
                        figtitle = figname;
                        if save
                            mkdir(figdir2);
                            savefigure(save,h,figtitle,figdir2);
                        end
                        fig = fig + 1;
                    end
                end
                
                %% Plot vmms raw and smoothed maps for each variable pair
                if ismember('vmms_rawvsmooth',whattoplot)
                    disp('Plotting vmms_rawsmooth ...');
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
                                % Remove stray view pixels inside pillars
                                if strcmp(msvar{kk},'view')
                                    map = emptyinsidepillar(map);
                                end
                                ax = subplot(2,2,2*(kk-1)+mm);
                                plotmap(map,msvar{kk});
                                colormap(ax,colorset);
                                % Scale color map
    %                             fieldcolor = 'y';
                                switch msvar{kk}
                                    case 'place'
                                        maxC = nanmax(map);
                                    case 'view'
                                        maxC = nanmax(map(3:end)); % leave out cue and hint
    %                                     fieldcolor = [50/256 205/256 50/256];
                                    case 'headdirection'
                                        maxC = nanmax(map);
                                end
                                % Patch boundaries of base fields
                                if ~strcmp(msvar{kk},'headdirection')
                                    for ff = 1:size(basedata.condbase_componentsperpx,1)
                                        if strcmp(msvar{kk},'place')
                                            color = [254 132 132]/255; % fieldcolor; % [254 132 132]/255;
                                        else
                                            color = [157 194 9]/255; % pistachio
                                        end
                                        % Patch base fields
                                        for pp = 1:size(basedata.fieldcoord{ff},1)
                                            % PATCH
                                            if strcmp(msvar{kk},'place')
                                                plotgrid = 3;
                                            else
                                                plotgrid = basedata.gridnum(ff);
                                            end
                                            [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                            patch(x,y,z,[1 1 1 1],'EdgeColor',color,'FaceColor','none','LineWidth',0.5,'EdgeAlpha',1);
                                        end
                                    end
                                else % For head direction, mark start and end of field
                                    for ff = 1:size(basedata.condbase_componentsperpx,1)
                                        hold on;
                                        section = nan(size(map));
                                        section(basedata.fieldlinbin{ff}) = 1.1*maxC;
                                        plotmap(section,msvar{kk},msvar{kk},fieldcolor);
                                    end
                                end
                                % Draw env boundaries
                                if ~strcmp(msvar{kk},'headdirection')
                                    patchenvbounds(msvar{kk});
                                end
                                if mm == 1
                                    % Tweak color scale to discount outliers
                                    set(ax,'CLim',[0 nanmean(map)+2*nanstd(map)]);
    %                                 set(ax,'CLim',[0 2*nanmean(map)]);
                                else 
                                    set(ax,'CLim',[0 maxC]);
                                end
                                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                    'XColor','none','YColor','none','ZColor','none',...
                                    'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                ax.Title.String = ['vmms ' msvar_short '- ' type ' ' msvar{kk} ' map'];
                            end
                        end
                        % Save
                        figtitle = horzcat(num2str(fig),'- vmms maps ',msvar_short,': ',cell_id);
                        if save
                            mkdir(figdir2);
                            savefigure(save, h,figtitle,figdir2);
    %                         saveas(h,figtitle,'png');
    %                         print('-painters',figtitle,'-dvg');
    %                         saveas(h,figtitle,'epsc');
                            print([figdir2 '/' figtitle],'-depsc');
                        end
                        close(h);
                        fig = fig + 1;
                    end
                end
                
                %% Plot base maps (all 3 vars)
                if ismember('basemap',whattoplot)
                    disp('Plotting basemap ...');
                    % Set up
                    h = figure(fig);
                    hold on;
                    figname = horzcat(num2str(fig),'- Base maps ',msvar_short,': ',cell_id,' ',sel);
                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                    varnames = {'place','view','headdirection'};
                    varobj = {'pv','pv','hv'};
                    for oo = 1:size(varnames,2)
                        basedata = objMain.data.(varobj{oo})(cell_indMain).(varnames{oo});
                        basemap = basedata.basemapLsm;
                        % Remove stray view pixels inside pillars
                        if strcmp(varnames{oo},'view')
                            basemap = emptyinsidepillar(basemap);
                        end
                        ax = subplot(1,3,oo);
                        plotmap(basedata.basemapLsm,varnames{oo});
                        colormap(ax,colorset);
                        % If firing rate or spike count = 0, set to black
                        settoblack(basedata.basemapLsm,varnames{oo});
                        maxC = nanmax(basemap);
                        % Rate and SIC 
                        rate = text(ax,1,1.1,1,[num2str(maxC,2) 'Hz'],'Units','Normalized','FontSize',30,'HorizontalAlignment','right');
                        si = text(ax,0,1.1,1,num2str(basedata.SI,2),'Units','Normalized','FontSize',30,'HorizontalAlignment','left');
                        if basedata.SI > basedata.SIthr
                            si.Color = 'r';
                        end
                        set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color',axcolor);
                        titlestring = text(ax,0.5,1.5,1,{[' Base ' varobj{oo} ' sm map ' ]; [num2str(basedata.sigfields) ' sig fields']},...
                            'Units','Normalized','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','bottom');
                        % Draw env boundaries
                        if ~strcmp(varnames{oo},'headdirection')
                            patchenvbounds(varnames{oo});
                        end
                    end
                    % Save
                    figtitle = horzcat(num2str(fig),'- Base maps sm ',': ',cell_id,' ',sel);
                    if save
                        savefigure(save, h,figtitle,figdir2);
                        print([figdir2 '/' figtitle],'-depsc');
                    end
                    close(h);
                    fig = fig + 1;
                end
                
                %% Plot pixel maps for each base field
                if ismember('pixelmap',whattoplot)
                    disp('Plotting pixelmap ...');
                    for pair = 1:size(spatialvarpairs,2)
                        msvar = spatialvarpairs{pair};
                        msvar_short = [msvar{1}(1) msvar{2}(1)];
    
                        for oo = 1:size(msvar,2)
                            basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                            secdata = objMain.data.(msvar_short)(cell_indMain).(msvar{2-oo+1});
                            
                            
                            for ff = 1:size(basedata.condbase_componentsperpx,1) % In mixselpv we've only extracted data for max of 3 fields per map
                                    if strcmp(msvar_short,'ph')
                                        disp('stop');
                                    end
                                    h = figure(fig);
                                    figname = '';
                                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 1 1]);
                                    hold on;
                                    
                                    % Plot the field outline
                                    ax = subplot(2,2,1);
                                    tempsecmap = nan(size(basedata.secmapLsm));
                                    % Remove stray view pixels inside pillars
                                    if strcmp(msvar{oo},'view')
                                        tempsecmap = emptyinsidepillar(tempsecmap);
                                    end
    %                                 fieldcolor = 'y';
                                    plotmap(tempsecmap,msvar{2-oo+1},msvar{oo});
                                    colormap(ax,colorset);
                                    % Patch env bounds
                                    if sum(strcmp(msvar,'headdirection')) == 0
                                        patchenvbounds('view');
                                    elseif sum(strcmp(msvar,'place')) == 1
                                        patchenvbounds('place');
                                    elseif sum(strcmp(msvar,'view'))
                                        patchenvbounds('view');
                                    end
                                    % Patch base fields
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
                                                patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1);
                                            elseif strcmp(msvar{oo},'view')
                                                patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1); 
                                            end
                                        end
                                    else
                                        hold on;
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    end
                                    % Patch basemap field % Actually nothing to
                                    % plot but use this to make size of map
                                    % equal across subplots
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
                                                patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1);
                                            elseif strcmp(msvar{oo},'view')
                                                patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                            end
                                        end
                                    else
                                        hold on;
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    end
                                    % Strech axis back to normal if hv and h as base
                                    if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                        axis equal;
    %                                 else
    %                                     disp('stop')
                                    end
                                    ax.Title.String = {'Field outline'};
                                    ax.Title.FontSize = 16;
                                    set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                        'XColor','none','YColor','none','ZColor','none',...
                                        'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                    
                                    
                                    % Plot raw base maps with secondary pixel maps
                                    ax = subplot(2,2,3);
                                    tempsecmap = nan(size(basedata.secmapLsm));
                                    tempsecmap(basedata.condbase_map_rw{ff,1}(:,1)) = basedata.condbase_map_rw{ff,1}(:,4);
                                    % Remove stray view pixels inside pillars
                                    if strcmp(msvar{oo},'view')
                                        tempsecmap = emptyinsidepillar(tempsecmap);
                                    end
    %                                 fieldcolor = 'y';
                                    plotmap(tempsecmap,msvar{2-oo+1},msvar{oo});
                                    colormap(ax,colorset);
                                    % Tweak color scale to discount outliers
                                    set(ax,'CLim',[0 max([1 2*nanstd(tempsecmap)+nanmean(tempsecmap)])]);
    %                                 % Patch secondary fields
    %                                 if ~strcmp(msvar{2-oo+1},'headdirection')
    %                                     for mm = 1:size(secdata.condbase_componentsperpx,1) % For each secondary field
    %                                         % Plot secondary field in consolidated pixel map
    %                                         sec = secdata.fieldcoord{mm}; % Get bins for sec field
    %                                         % Patch secondary fields in the consolidated pixel map
    %                                         for pp = 1:size(sec,1)
    %                                             % PATCH
    %                                             if strcmp(msvar{2-oo+1},'place')
    %                                                 plotgrid = 3;
    %                                             else
    %                                                 plotgrid = secdata.gridnum(mm);
    %                                             end
    %                                             [x,y,z] = converttosurf(plotgrid,sec(pp,1),sec(pp,2));
    %                                             patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1); % Dark Grey
    % %                                             if strcmp(msvar{2-oo+1},'place')
    % %                                                 patch(x,y,z,'r','EdgeColor','r','FaceColor','none','LineWidth',2); % Red
    % %                                             elseif strcmp(msvar{2-oo+1},'view')
    % %                                                 patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green
    % %                                             end
    %                                         end
    %                                     end
    %                                 else
    %                                     for mm = 1:size(secdata.condbase_componentsperpx,1) % For each secondary field
    %                                         hold on;
    %                                         section = nan(size(tempsecmap));
    %                                         section(secdata.fieldlinbin{mm}) = 1.1*nanmax(tempsecmap);
    %                                         plotmap(section,msvar{2-oo+1},msvar{oo},fieldcolor);
    %                                     end
    %                                 end
                                    set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                        'XColor','none','YColor','none','ZColor','none',...
                                        'FontSize',14,'GridLineStyle','none','Color',axcolor);
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
                                                patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1);
                                            elseif strcmp(msvar{oo},'view')
                                                patch(x,y,z,'r','EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
                                            end
                                        end
                                    else
                                        hold on;
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    end
                                    % Strech axis back to normal if hv and h as base
                                    if strcmp(msvar_short,'hv') && strcmp(msvar{oo},'headdirection')
                                        axis equal;
    %                                 else
    %                                     disp('stop')
                                    end
                                    ax.Title.String = cell(size(basedata.condbase_insecfieldrates{ff},1),1);
                                    linked = [];
                                    for ss = 1:size(basedata.condbase_insecfieldrates{ff},1)
                                        % Check if infield firing in sec field is higher than outfield firing
                                        if basedata.condbase_insecfieldrates{ff}(ss,1) > prctile(basedata.condbase_outsecfieldrates{ff}{ss,1},95)
                                            linked = horzcat(linked,',',num2str(ss));
                                            ax.Title.String{ss,1} = horzcat('LINKED:',msvar{2-oo+1},' Field ',num2str(ss),' Infield: ', num2str(basedata.condbase_insecfieldrates{ff}(ss,1),2),'Hz; Outfield: ', num2str(prctile(basedata.condbase_outsecfieldrates{ff}{ss,1},95),2),'Hz');
                                        else
                                            ax.Title.String{ss,1} = horzcat(msvar{2-oo+1},' Field ',num2str(ss),' Infield: ', num2str(basedata.condbase_insecfieldrates{ff}(ss,1),2),'Hz; Outfield: ', num2str(prctile(basedata.condbase_outsecfieldrates{ff}{ss,1},95),2),'Hz');
                                        end
                                    end
                                    ax.Title.String{end+1,1} = ['Rawpx = ' num2str(sum(~isnan(tempsecmap)))];
                                    if isempty(linked)
                                        figname = horzcat(num2str(fig),'- Pixel maps ', msvar_short, '-', msvar{oo},' field ',num2str(ff),' not linked');
                                    else
                                        figname = horzcat(num2str(fig),'- Pixel maps ', msvar_short, '-', msvar{oo},' field ',num2str(ff),' sig linked to ',msvar{2-oo+1},' field ', linked);
                                    end
    
                                    % Plot smoothed sec pixel map
                                    ax = subplot(2,2,2);
    %                                     temp1 = basedata.(['sec' mapname]){ff};
                                    temp1 = basedata.(['condbase_map_sm' ]){ff};
                                    % Remove stray view pixels inside pillars
                                    if strcmp(msvar{oo},'view')
                                        temp1 = emptyinsidepillar(temp1);
                                    end
                                    switch msvar{2-oo+1}
                                        case 'place'
                                            maxC = nanmax(temp1);
                                        case 'view'
                                            maxC = nanmax(temp1(3:end));
                                        case 'headdirection'
                                            maxC = nanmax(temp1);
                                    end
    %                                 fieldcolor = 'y';
                                    plotmap(temp1,msvar{2-oo+1},msvar{oo});
                                    colormap(ax,colorset);
                                    % If firing rate or spike count = 0, set to black
                                    settoblack(temp1,msvar{2-oo+1});
                                    % Label
                                    pseudoSIthr = prctile(basedata.(['pseudosecSIC_adsm' ]){ff},95); 
                                    ax.Title.String = {horzcat('SI: ',num2str(basedata.(['condbase_SIC_sm' ])(ff),2),'/',num2str(pseudoSIthr,2)); ...
                                        horzcat('Smoothpx = ',num2str(sum(~isnan(temp1)))); ...
                                        horzcat('Pseudosmoothpx = ', num2str(mean(sum(~isnan(basedata.pseudosecmaps_adsm{ff}),2)),2))};
                                    ax.Title.FontSize = 16;
                                    if basedata.(['condbase_SIC_sm' ])(ff)>pseudoSIthr && maxC >= 0.7 
                                        ax.Title.Color = 'r';
                                    end
                                    set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                        'XColor','none','YColor','none','ZColor','none',...
                                        'FontSize',14,'GridLineStyle','none','Color',axcolor);
                                    % Patch env bounds
                                    if sum(strcmp(msvar,'headdirection')) == 0
                                        patchenvbounds('view');
                                    elseif sum(strcmp(msvar,'place')) == 1
                                        patchenvbounds('place');
                                    elseif sum(strcmp(msvar,'view'))
                                        patchenvbounds('view');
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
                                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1);
    %                                         if strcmp(msvar{oo},'place')
    %                                             patch(x,y,z,[1 1 1 1],'EdgeColor','r','FaceColor','none','LineWidth',1);
    %                                         elseif strcmp(msvar{oo},'view')
    %                                             patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
    %                                         end
                                        end
                                    else
                                        hold on;
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
                                    end
                                    % Patch fields found from smoothed conditioned maps
                                    tempdata = basedata.condbase_fieldcoord{ff};
                                    if ~isempty(tempdata) % If there are secondary fields in conditioned maps
                                        overlap = {};
                                        for gg = 1:size(tempdata,1) % For each of the secondary fields
                                            if sum(basedata.condbase_fieldoverlapind{ff}{gg})>0
    %                                             fieldcolor = 'y';
                                                xxx = find(basedata.condbase_fieldoverlapind{ff}{gg});
                                                for xx = 1:size(xxx,1)
                                                    % Does the sec field overlap with any of the base fields from the other var?
                                                    %%%%%% PATCH
                                                    if ff<10
                                                        basefieldnum = ['0' num2str(ff)];
                                                    else basefieldnum = num2str(ff);
                                                    end
                                                    if xxx(xx)<10
                                                        secfieldnum = ['0' num2str(xxx(xx))];
                                                    else secfieldnum = num2str(xxx(xx));
                                                    end
                                                    overlap{end+1,1} = ['Base field ' basefieldnum ' X Sec field '...
                                                        secfieldnum];
    %                                                 if xxx(xx)>3
    %                                                     disp('stop')
    %                                                 end
                                                end
                                            else
    %                                             fieldcolor = 'b';
                                            end
                                        end
                                        ax.Title.String{end+1,1} = [num2str(size(tempdata,1)) ' cond fields - mean coverage ' ...
                                            num2str(mean(basedata.condbase_fieldsizepercent{ff}),2) '% per field'];
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
                                    ax = subplot(2,2,4);
                                    temp2 = basedata.secmaps_dist{ff};
                                    % Remove stray view pixels inside pillars
                                    if strcmp(msvar{oo},'view')
                                        temp2 = emptyinsidepillar(temp2);
                                    end
                                    switch msvar{2-oo+1}
                                        case 'place'
                                            maxC = nanmax(temp2);
                                        case 'view'
                                            maxC = nanmax(temp2(3:end));
                                        case 'headdirection'
                                            maxC = nanmax(temp2);
                                    end
    %                                 fieldcolor = 'y';
                                    plotmap(temp2,msvar{2-oo+1},msvar{oo});
                                    colormap(ax,colorset)
                                    % Tweak color scale to discount outliers
                                    set(ax,'CLim',[0 max([1 2*nanstd(temp2)+nanmean(temp2)])]);
                                    % If firing rate or spike count = 0, set to black
                                    settoblack(temp2,msvar{2-oo+1});
                                    % Label
                                    pseudoSIthr = prctile(basedata.(['pseudosecSIC_adsm' ]){ff},95); 
                                    ax.Title.String = {['Predicted ' msvar{2-oo+1} ' map']; ...
                                        ['Field ' num2str(ff) ' dist ratio_' basedata.varname_short ' = ' num2str(basedata.base_distratio(ff))];...
                                        ['Field ' num2str(ff) ' dist ratio_' secdata.varname_short ' = ' num2str(basedata.sec_distratio(ff))]};
                                    ax.Title.FontSize = 16;
                                    if basedata.(['condbase_SIC_sm' ])(ff)>pseudoSIthr && maxC >= 0.7 
                                        ax.Title.Color = 'r';
                                    end
                                    set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                        'XColor','none','YColor','none','ZColor','none',...
                                        'FontSize',14,'GridLineStyle','none','Color',axcolor);
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
                                            patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1);
    %                                         if strcmp(msvar{oo},'place')
    %                                             patch(x,y,z,[1 1 1 1],'EdgeColor','r','FaceColor','none','LineWidth',1);
    %                                         elseif strcmp(msvar{oo},'view')
    %                                             patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',1); % Lime green. Olive drab:[107/256 142/256 35/256], Medium sea green: [60/256 179/256 113/256]
    %                                         end
                                        end
                                    else
                                        hold on;
                                        section = nan(size(basedata.basemapLrw));
                                        section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                        plotmap(section,msvar{oo},msvar{2-oo+1},fieldcolor);
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
                                        print([figdir2 '/' figtitle],'-depsc');
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
        %                             temp(basedata.condbase_map_rw{ff}(:,1),1) = basedata.condbase_map_rw{ff}(:,4);
        %                             plotmap(temp,msobj{2-oo+1});
        %                             colormap(ax,cool);
        %                             % If firing rate or spike count = 0, set to black
        %                             settoblack(temp,msobj{2-oo+1});
        %                             numspk = sum(basedata.condbase_map_rw{ff}(:,3));
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
        %                             for mm = 1:size(secdata.condbase_componentsperpx,1) % For each secondary field
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
                end
                
                %% Plot pc/sv vs corr smoothed maps - uncorrected
                % For checking filtering across objects - optional
                if ismember('checkfilter',whattoplot)
                    disp('Plotting checkfilter ...');
                    for pair = 1:size(spatialvarpairs,2)
                        msvar = spatialvarpairs{pair};
                        msvar_short = [msvar{1}(1) msvar{2}(1)];
                        temp = objCorr.data.(msvar_short);
    
                        h = figure(fig);
                        hold on;
                        figname = horzcat(num2str(fig),'- Comparison of smoothed base maps ', msvar_short, '-', cell_id,'-',num2str(fig));
                        set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                        for oo = 1:size(msvar,2)
                            basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                            % Plot maps from pc/v objects
                            ax = subplot(2,2,1+(oo-1)*2);
                            if strcmp(msvar{oo},'place')
                                basemap = objPlace.data.maps_sm(cell_indP,:);
                                maxC = max(basemap);
                                title = ['vmpc smoothed place map: ' sel];
    %                             fieldcolor = 'r';
                            elseif strcmp(msvar{oo},'view')
                                basemap = objView.data.maps_sm(cell_indV,:);
                                maxC = max(basemap(3:end));
                                title = ['vmsv smoothed view map: ' sel];
    %                             fieldcolor = [50/256 205/256 50/256];
                            elseif strcmp(msvar{oo},'headdirection')
                                basemap = objHeaddirection.data.maps_sm(cell_indH,:);
                                maxC = nanmax(basemap);
                                title = ['vmhd smoothed headdirection map: ' sel];
                            end
                            % Remove stray view pixels inside pillars
                            if strcmp(msvar{oo},'view')
                                basemap = emptyinsidepillar(basemap);
                            end
    %                         fieldcolor = 'y';
                            plotmap(basemap,msvar{oo});
                            colormap(ax,colorset);
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,msvar{oo});
                            % Patch boundaries of base fields
                            if ~strcmp(msvar{oo},'headdirection')
                                for ff = 1:size(basedata.condbase_componentsperpx,1)
                                    % Patch base fields
                                    for pp = 1:size(basedata.fieldcoord{ff},1)
                                        % PATCH
                                        if strcmp(msvar{oo},'place')
                                            plotgrid = 3;
                                        else
                                            plotgrid = basedata.gridnum(ff);
                                        end
                                        [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                        patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1,'EdgeAlpha',1);
                                    end
                                end
                            else
                                hold on;
                                for ff = 1:size(basedata.condbase_componentsperpx,1)
                                    section = nan(size(basemap));
                                    section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                    plotmap(section,msvar{oo},msvar{oo},fieldcolor);
                                end
                            end
                            % Draw env boundaries
                            if ~strcmp(msvar{oo},'headdirection')
                                patchenvbounds(msvar{oo});
                            end
                            set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
    %                         colormap(ax,cool);
                            ax.Title.String = title;
                            
                            % Plot maps from corr objects
                            ax = subplot(2,2,2+(oo-1)*2);
                            basemap = temp(cell_indCorr).(['maps_sm_corr' msvar{oo}(1) 'set'])(1,:);
                            % Remove stray view pixels inside pillars
                            if strcmp(msvar{oo},'view')
                                basemap = emptyinsidepillar(basemap);
                            end
    %                         fieldcolor = 'y';
                            title = ['vmcorr uncorrected smoothed ' msvar{oo} ' map: ' sel];
                            if ~strcmp(msvar{oo},'view')
                                maxC = max(basemap);
                            else 
                                maxC = max(basemap(3:end));
                            end
                            plotmap(basemap,msvar{oo});
                            colormap(ax,colorset);
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,msvar{oo});
                            % Patch boundaries of base fields
                            if ~strcmp(msvar{oo},'headdirection')
                                for ff = 1:size(basedata.condbase_componentsperpx,1)
                                    % Patch base fields
                                    for pp = 1:size(basedata.fieldcoord{ff},1)
                                        % PATCH
                                        if strcmp(msvar{oo},'place')
                                            plotgrid = 3;
                                        else
                                            plotgrid = basedata.gridnum(ff);
                                        end
                                        [x,y,z] = converttosurf(plotgrid,basedata.fieldcoord{ff}(pp,1),basedata.fieldcoord{ff}(pp,2));
                                        patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1,'EdgeAlpha',1);
                                    end
                                end
                            else
                                hold on;
                                for ff = 1:size(basedata.condbase_componentsperpx,1)
                                    section = nan(size(basemap));
                                    section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                    plotmap(section,msvar{oo},msvar{oo},fieldcolor);
                                end
                            end
                            % Draw env boundaries
                            if ~strcmp(msvar{oo},'headdirection')
                                patchenvbounds(msvar{oo});
                            end
                            set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
    %                         colormap(ax,cool);
                            ax.Title.String = title;
                        end
                        % Save
                        figtitle = figname;
                        if save
                            savefigure(h,figtitle,figdir2);
                            print([figdir2 '/' figtitle],'-depsc');
                        end
                        close(h);
                        fig = fig + 1;
                    end
                end
                
 
                %% Plot predicted maps (vmcorr)
                if ismember('predmap',whattoplot)
                    disp('Plotting predmap ...');
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
                            % Remove stray view pixels inside pillars
                            if strcmp(msvar{oo},'view')
                                basemap = emptyinsidepillar(basemap);
                            end
    %                         title = {['vmcorr smoothed ' msvar{oo}(1) ' map: ' sel]};
                            if ~strcmp(msvar{oo},'view')
                                maxC = max(basemap);
    %                             if strcmp(msvar{oo},'place')
    % %                                 fieldcolor = 'r';
    %                             else
    % %                                 fieldcolor = 'b';
    %                             end
                            else
                                maxC = max(basemap(3:end));
    %                             fieldcolor = [50/256 205/256 50/256];
                            end
    %                         fieldcolor = 'y';
                            ax = subplot(2,4,1+(oo-1)*4);
                            plotmap(basemap,msvar{oo});
                            colormap(ax,colorset);
                            % If selective, patch field
    %                         if pair == 1 || (pair==2 && oo == 1)
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,msvar{oo});
                            basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                            % Patch boundaries of base fields
                            for ff = 1:size(basedata.condbase_componentsperpx,1)
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
                                        patch(x,y,z,[1 1 1 1],'EdgeColor',fieldcolor,'FaceColor','none','LineWidth',1,'EdgeAlpha',1);
                                    end
                                else
                                    hold on;
                                    section = nan(size(basemap));
                                    section(basedata.fieldlinbin{ff}) = 1; % *nanmax(tempsecmap);
                                    plotmap(section,msvar{oo},msvar{oo},fieldcolor);
                                end
                            end
                            % Draw env boundaries
                            if ~strcmp(msvar{oo},'headdirection')
                                patchenvbounds(msvar{oo});
                            end
                            set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
    %                         colormap(ax,cool);
                            ax.Title.String = {[msvar_short ' uncorrected smoothed '];[ msvar{oo} ' map: ' sel]};
    
                            % Plot vmcorr raw uncorrected pv maps
                            ax = subplot(2,4,2+(oo-1)*4);
                            basemap = temp(cell_indCorr).(['maps_raw_corr' msvar{oo}(1) 'set'])(1,:);
                            % Remove stray view pixels inside pillars
                            if strcmp(msvar{oo},'view')
                                basemap = emptyinsidepillar(basemap);
                            end
    %                         title = ['vmcorr raw ' msvar{oo} ' map: ' sel];
                            if ~strcmp(msvar{oo},'view')
                                maxC = max(basemap);
    %                             if strcmp(msvar{oo},'place')
    %                                 fieldcolor = 'r';
    %                             else
    %                                 fieldcolor = 'b';
    %                             end
                            else
                                maxC = max(basemap(3:end));
    %                             fieldcolor = [50/256 205/256 50/256];
                            end
    %                         fieldcolor = 'y';
                            plotmap(basemap,msvar{oo});
                            colormap(ax,colorset);
    %                         % Tweak color scale to discount outliers
    %                         set(ax,'CLim',[0 2*nanmean(basemap)]);
                            % Draw env boundaries
                            if ~strcmp(msvar{oo},'headdirection')
                                patchenvbounds(msvar{oo});
                            end
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,msvar{oo});
    %                         basedata = objMain.data.(msvar{oo})(cell_indMain);
                            set(ax,'CLim',[0 2*nanstd(basemap)+nanmean(basemap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
    %                         colormap(ax,cool);
                            ax.Title.String = {[msvar_short ' uncorrected raw '];[ msvar{oo} ' map: ' sel]};
    
                            % Plot predicted raw maps from corr objects
                            ax = subplot(2,4,3+(oo-1)*4);
                            basemap = temp(cell_indCorr).(['maps_dist_' msvar{oo}(1)]);
                            % Remove stray view pixels inside pillars
                            if strcmp(msvar{oo},'view')
                                basemap = emptyinsidepillar(basemap);
                            end
                            title = {horzcat(msvar_short,' predicted raw ', msvar{oo}, ' map: ',sel); ...
                                    horzcat('Cell distratio_', msvar{oo}(1), ' = ',num2str(temp(cell_indCorr).(['distratio_' msvar{oo}(1)])));...
                                    horzcat('Cell distratio_', msvar{2-oo+1}(1), ' = ',num2str(temp(cell_indCorr).(['distratio_' msvar{2-oo+1}(1)])))};
                            if ~strcmp(msvar{oo},'view')
                                maxC = max(basemap);
                            else
                                maxC = max(basemap(3:end));
                            end
                            plotmap(basemap,msvar{oo});
                            colormap(ax,colorset);
    %                         % Tweak color scale to discount outliers
    %                         set(ax,'CLim',[0 2*nanmean(basemap)]);
                            % Draw env boundaries
                            if ~strcmp(msvar{oo},'headdirection')
                                patchenvbounds(msvar{oo});
                            end
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,msvar{oo});
                            set(ax,'CLim',[0 2*nanstd(basemap)+nanmean(basemap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
    %                         colormap(ax,cool);
                            ax.Title.String = title;
                            
                            % Plot predicted smooothed maps from corr objects
                            ax = subplot(2,4,4+(oo-1)*4);
                            basemap = temp(cell_indCorr).(['maps_dist_' msvar{oo}(1) '_adsm']);
                            % Remove stray view pixels inside pillars
                            if strcmp(msvar{oo},'view')
                                basemap = emptyinsidepillar(basemap);
                            end
                            title = {[msvar_short ' predicted smoothed ' msvar{oo}(1) ' map: ' sel]; ...
                                    horzcat('Cell distratio_', msvar{oo}(1), '= ',num2str(temp(cell_indCorr).(['distratio_' msvar{oo}(1)])));...
                                    horzcat('Cell distratio_', msvar{2-oo+1}(1), ' = ',num2str(temp(cell_indCorr).(['distratio_' msvar{2-oo+1}(1)])))};
                            if ~strcmp(msvar{oo},'view')
                                maxC = max(basemap);
                            else
                                maxC = max(basemap(3:end));
                            end
                            plotmap(basemap,msvar{oo});
                            colormap(ax,colorset);
                            % Draw env boundaries
                            if ~strcmp(msvar{oo},'headdirection')
                                patchenvbounds(msvar{oo});
                            end
                            % If firing rate or spike count = 0, set to black
                            settoblack(basemap,msvar{oo});
                            set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                                'XColor','none','YColor','none','ZColor','none',...
                                'FontSize',14,'GridLineStyle','none','Color',axcolor);
    %                         colormap(ax,cool);
                            ax.Title.String = title;
                            
                        end
                        % Save
                        figtitle = horzcat(num2str(fig),'- Predicted raw maps ',msvar_short,'-',num2str(fig));
                        if save
    %                         mkdir(figdir2);
                            savefigure(h,figtitle,figdir2);
                            print([figdir2 '/' figtitle],'-depsc');
                        end
                        close(h);
                        fig = fig + 1;
                    end
                end
                
                %% Plot linearised differences between actual and predicted
                % maps
                if ismember('linearizeddiff',whattoplot)
                    disp('Plotting linearizeddiff ...');
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
                            % Plot vmcorr actual raw map
                            basemap = temp(cell_indCorr).(['maps_raw_corr' msvar{oo}(1) 'set'])(1,:);
                            title = {[msvar_short ' raw vs predicted ' msvar{oo}(1) ' map: ' sel]; ...
                                    ['Cell dist ratio_' msvar{oo}(1) ' = ' num2str(temp(cell_indCorr).(['distratio_' msvar{oo}(1)]))];...
                                    ['Cell dist ratio_' msvar{2-oo+1}(1) ' = ' num2str(temp(cell_indCorr).(['distratio_' msvar{2-oo+1}(1)]))]};
                            if ~strcmp(msvar{oo},'view')
                                maxC = max(basemap);
                            else
                                maxC = max(basemap(3:end));
                            end
                            plot(1:size(basemap,2),log(basemap),'ko-');
                            hold on;
                            % Plot vmcorr predicted raw map
                            basemap = temp(cell_indCorr).(['maps_dist_' msvar{oo}(1)]);
                            if ~strcmp(msvar{oo},'view')
                                maxC = max(basemap);
                            else
                                maxC = max(basemap(3:end));
                            end
                            plot(1:size(basemap,2),log(basemap),'ro-');
                            
                            hold off;
                            set(ax,'FontSize',14,'GridLineStyle','none','Color','none');
                            colormap(ax,cool);
                            ax.Title.String = title;
                        end
                        % Save
                        figtitle = figname;
                        if save
                            savefigure(h,figtitle,figdir2);
                            print([figdir2 '/' figtitle],'-depsc');
                        end
                        close(h);
                        fig = fig + 1;
                    end
                end
                
                % Linear or non linear mixed selectivity?
                if ismember('lnlmixsel',whattoplot)
                    disp('Plotting lnlmixsel ...');
                    for pair = 1:size(spatialvarpairs,2)
                        
                        msvar = spatialvarpairs{pair};
                        msvar_short = [msvar{1}(1) msvar{2}(1)];
                       
                        inrate_var = cell(2,1);
                        outrate_var = cell(2,1);
                        
                        for oo = 1:size(msvar,2)
                            basedata = objMain.data.(msvar_short)(cell_indMain).(msvar{oo});
                            secdata = objMain.data.(msvar_short)(cell_indMain).(msvar{2-oo+1});
                            
                            inrate_var{oo} = nan(objMain.data.Args.NumShuffles(cell_indMain)*size(basedata.condbase_componentsperpx,1),1);
                            outrate_var{oo} = nan(objMain.data.Args.NumShuffles(cell_indMain)*size(basedata.condbase_componentsperpx,1),1);
    
                            for ff = 1:size(basedata.condbase_componentsperpx,1)
    %                            
                                for xx = 1:1000 % objMain.data.Args.NumShuffles(cell_indMain) % For each shuffle
                                    
                                    basebin_size = length(basedata.fieldlinbin{ff});
                                    pseudobasebins = basedata.pseudosecdataperfield{ff}{1}{xx};
                                    nonoverlap = setdiff(pseudobasebins,basedata.fieldlinbin{ff});
                                    pseudobasebin_size = size(nonoverlap,1);
                                    draw = round((min([basebin_size pseudobasebin_size]))/2); % or draw from ceil(basebin_size/2)
                                    
                                    inds = randsample(1:basebin_size,draw);
                                    indbins_base = basedata.fieldlinbin{ff}(inds);
                                    % draw rates from raw map
                                    inrate_var{oo}((ff-1)*objMain.data.Args.NumShuffles(cell_indMain)+xx) = nanmean(nanmean(basedata.basemapLrw(indbins_base)));
                                    inds = randsample(1:pseudobasebin_size,draw);
                                    indbins_pseudobase = nonoverlap(inds);
                                    % draw rates from raw map
                                    outrate_var{oo}((ff-1)*objMain.data.Args.NumShuffles(cell_indMain)+xx) = nanmean(basedata.basemapLrw(indbins_pseudobase));
                                    
                                end
                                
                            end
                            
                            
                        end
                        
                        minobs = min([size(inrate_var{1},1),size(inrate_var{2},1),size(outrate_var{1},1),size(outrate_var{2},1)]);
                        
                        % 2-way ANOVA
                        if minobs > 0
                            
                            h = figure(fig);
                            hold on;
                            figname = horzcat(num2str(fig),'- Hedges g ', msvar_short);
                            set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                            if size(inrate_var{1},1) > 0 % If var 1 has a field
                                inds_var1 = randsample(1:size(inrate_var{1},1),minobs);
                                in_var1 = inrate_var{1}(inds_var1);
                                out_var1 = outrate_var{1}(inds_var1);
                                % Hedge's G
                                pooledstd1 = sqrt(( (minobs-1)*std(in_var1)*std(in_var1)+(minobs-1)*std(out_var1)*std(out_var1) ) ...
                                    /( (minobs-1)+(minobs-1) ));
                                g1 = ( mean(in_var1)-mean(out_var1) )/(pooledstd1);
                                ax = subplot(1,2,1); 
                                var1g = text(ax,0.5,0.5,1,{msvar{1}; ' in/out-field Hedges g = '; num2str(g1,2)},...
                                    'Units','Normalized','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle');
                            end
                            if size(inrate_var{2},1) > 0 % If var 2 has a field
                                inds_var2 = randsample(1:size(inrate_var{2},1),minobs);
                                in_var2 = inrate_var{2}(inds_var2);
                                out_var2 = outrate_var{2}(inds_var2);
                                % Hedge's G
                                pooledstd2 = sqrt(( (minobs-1)*std(in_var2)*std(in_var2)+(minobs-1)*std(out_var2)*std(out_var2) ) ...
                                    /( (minobs-1)+(minobs-1) ));
                                g2 = ( mean(in_var2)-mean(out_var2) )/(pooledstd2);
                                ax = subplot(1,2,2); 
                                var2g = text(ax,0.5,0.5,1,{msvar{2}; ' in/out-field Hedges g = '; num2str(g2,2)},...
                                    'Units','Normalized','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle');
                            end
                            % Save
                            figtitle = figname;
                            if save
                                savefigure(h,figtitle,figdir2);
                                print([figdir2 '/' figtitle],'-depsc');
                            end
                            close(h);
                            fig = fig + 1;
                        
                            % 2-way ANOVA test
                            [pval,tbl,stats] = anova2([in_var1 out_var1;in_var2 out_var2],minobs);
        %                     [pval,tbl,stats] = anova2([inplace inview; outplace outview],minobs); 
    %                         c1 = multcompare(stats);
                            figanova = gcf;
                            figtitle = [num2str(fig) '- 2way ANOVA Table ' msvar_short];
                            if save
                                savefigure(figanova,figtitle,figdir2);
                                print([figdir2 '/' figtitle],'-depsc');
                            end
                            close(figanova);
                            fig = fig + 1;
                            
                            h = figure(fig);
                            hold on;
                            set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                            ax = gca;
                            detail = {['Main p1 = ' num2str(pval(1)),2]; ...
                                ['Main p2 = ' num2str(pval(2)),2]; ...
                                ['Inter p = ' num2str(pval(3)),2]; ...
                                ['Inter F = ' num2str(tbl{4,5},2)];...
                                ['Sigma sq = ' num2str(stats.sigmasq)];...
                                ['Col means 1 = ' num2str(stats.colmeans(1)) ', Col means 2 = ' num2str(stats.colmeans(2))];...
                                ['Row means 1 = ' num2str(stats.rowmeans(1)) ', Row means 2 = ' num2str(stats.rowmeans(2))]};
                            text(ax,0.5,0.5,1,detail,...
                                'Units','Normalized','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle');
                            % Save 
                            figtitle = [num2str(fig) '- 2way ANOVA Cell Stats ' msvar_short];
                            if save
                                savefigure(h,figtitle,figdir2);
                                print([figdir2 '/' figtitle],'-depsc');
                            end
                            close(h);
                            fig = fig + 1;
                        else
                            h = figure(fig);
                            hold on;
                            figname = horzcat(num2str(fig),'- Hedges g ', msvar_short);
                            set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                            if size(inrate_var{1},1) > 0 % If var 1 has a field
                                in_var1 = inrate_var{1};
                                out_var1 = outrate_var{1};
                                numobs1 = size(inrate_var{1},1);
                                % Hedge's G
                                pooledstd1 = sqrt(( (numobs1-1)*std(in_var1)*std(in_var1)+(numobs1-1)*std(out_var1)*std(out_var1) ) ...
                                    /( (numobs1-1)+(numobs1-1) ));
                                g1 = ( mean(in_var1)-mean(out_var1) )/(pooledstd1);
                                ax = subplot(1,2,1); 
                                var1g = text(ax,0.5,0.5,1,{msvar{1}; ' in/out-field Hedges g = '; num2str(g1,2)},...
                                    'Units','Normalized','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle');
                            end
                            if size(inrate_var{2},1) > 0 % If var 2 has a field
                                in_var2 = inrate_var{2};
                                out_var2 = outrate_var{2};
                                numobs2 = size(inrate_var{2},1);
                                % Hedge's G
                                pooledstd2 = sqrt(( (numobs2-1)*std(in_var2)*std(in_var2)+(numobs2-1)*std(out_var2)*std(out_var2) ) ...
                                    /( (numobs2-1)+(numobs2-1) ));
                                g2 = ( mean(in_var2)-mean(out_var2) )/(pooledstd2);
                                ax = subplot(1,2,2); 
                                var2g = text(ax,0.5,0.5,1,{msvar{2}; ' in/out-field Hedges g = '; num2str(g2,2)},...
                                    'Units','Normalized','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle');
                            end
                            % Save
                            figtitle = figname;
                            if save
                                savefigure(h,figtitle,figdir2);
                                print([figdir2 '/' figtitle],'-depsc');
                            end
                            close(h);
                            fig = fig + 1;
                            
                            % 1-way ANOVA
                            for oo = 1:size(msvar,2)
                                if size(inrate_var{oo},1) > 0 % If var 1 has a field
                                    [pval,tbl,stats] = anova1([eval(['in_var' num2str(oo)]) eval(['out_var' num2str(oo)])]);
                                    figanova = figure(1);
                                    figtitle = [num2str(fig) '- 1way ANOVA ' msvar{oo} ' table'];
                                    if save
                                        savefigure(figanova,figtitle,figdir2);
                                        print([figdir2 '/' figtitle],'-depsc');
                                    end
                                    close(figanova);
                                    fig = fig + 1;
    
                                    figplot = figure(2);
                                    figtitle = [num2str(fig) '- 1way ANOVA ' msvar{oo} ' plot'];
                                    if save
                                        savefigure(figplot,figtitle,figdir2);
                                        print([figdir2 '/' figtitle],'-depsc');
                                    end
                                    close(figplot);
                                    fig = fig + 1;
    
                                    h = figure(fig);
                                    hold on;
                                    set(h,'Name',figname,'Units','Normalized','Position',[0 1 0.75 0.75]);
                                    ax = gca;
                                    detail = {['Main p1 = ' num2str(pval(1)),2]; ...
                                        ['Main F = ' num2str(tbl{2,5},2)];...
                                        ['Col means 1 = ' num2str(stats.means(1)) ', Col means 2 = ' num2str(stats.means(2))]};
                                    text(ax,0.5,0.5,1,detail,...
                                        'Units','Normalized','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle');
                                    % Save 
                                    figtitle = [num2str(fig) '- 1way ANOVA Stats' msvar{oo}];
                                    if save
                                        savefigure(h,figtitle,figdir2);
                                        print([figdir2 '/' figtitle],'-depsc');
                                    end
                                    close(h);
                                    fig = fig + 1;
                                end
                            end
                        end
                    end
                end
            else
                disp(['Cell ' num2str(cellcount) ' ' cell_id ' discarded']);
            end
        end
    end
    close all;
    
    
    
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
                                    map(basedata.condbase_map_rw{ff,1}(:,1)) = basedata.condbase_map_rw{ff,1}(:,4);
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
%                                 settoblack(map,msvar{2-oo+1});
                                axdisplay(ax,maxC);
                                colormap(ax,'cool');
                                c = colorbar;
                                c.Position = [col_left(2*(oo-1)+2)+1.05*axwidth row_bot(4+2*(ff-1)) 0.9/150 axheight*2];
                                if ~strcmp(maptype,'raw')
                                    overlap = {};
                                    for dd = 1:size(basedata.cond_fieldoverlap{ff},1)
                                        if basedata.cond_fieldoverlap{ff}(dd)
%                                             fieldcolor = 'y';
                                            xxx = find(basedata.condbase_fieldoverlapind{ff}{dd});
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
                                                        if any(secdata.condbase_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.condbase_fieldoverlapind{xx}{ss});
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
%                                             fieldcolor = 'b';
                                        end
                                        sec = basedata.condbase_fieldcoord{ff}{dd};
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
                                    ax.Title.String = {['SI=' num2str(basedata.condbase_SIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_orig(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.condbase_SIC_adsm(ff) > pseudothr_orig(oo) && maxC >= 0.7
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
                                        map(basedata.condbase_map_rw{ff,1}(:,1)) = basedata.condbase_map_rw{ff,1}(:,4);
                                        maxC = nanmax(map(3:end));
                                        plotmap(map,msvar{2-oo+1});
                                        patchenvbounds(msvar{2-oo+1});
%                                         settoblack(map,msvar{2-oo+1});
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
                                        ax.Title.String = {horzcat(msvar{2-oo+1}(1),num2str(fff),' Infield: ', num2str(basedata.condbase_insecfieldrates{ff}(fff,1),2),'Hz'); ...
                                            horzcat('Outfield: ', num2str(prctile(basedata.condbase_outsecfieldrates{ff}{fff,1},95),2),'Hz')};
                                        if basedata.condbase_insecfieldrates{ff}(fff,1) > prctile(basedata.condbase_outsecfieldrates{ff}{fff,1},95)
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
%                                     settoblack(map,msvar{2-oo+1});
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
%                                             fieldcolor = 'y';
                                            xxx = find(basedata.condbase_fieldoverlapind{ff}{dd});
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
                                                        if any(secdata.condbase_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.condbase_fieldoverlapind{xx}{ss});
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
%                                             fieldcolor = 'b';
                                        end
                                        sec = basedata.condbase_fieldcoord{ff}{dd};
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
                                    ax.Title.String = {['SI=' num2str(basedata.condbase_SIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_orig(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.condbase_SIC_adsm(ff) > pseudothr_orig(oo) && maxC >= 0.7
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
                                    map(basedata.condbase_map_rw{ff,1}(:,1)) = basedata.condbase_map_rw{ff,1}(:,4);
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
%                                 settoblack(map,msvar{2-oo+1});
                                axdisplay(ax,maxC);
                                colormap(ax,'cool');
                                c = colorbar;
                                c.Position = [col_left(2*(oo-1)+6)+1.05*axwidth row_bot(4+(ff-1)*2) 0.9/150 axheight*2];
                                if ~strcmp(maptype,'raw')
                                    overlap = {};
                                    for dd = 1:size(basedata.cond_fieldoverlap{ff},1)
                                        if basedata.cond_fieldoverlap{ff}(dd)
%                                             fieldcolor = 'y';
                                            xxx = find(basedata.condbase_fieldoverlapind{ff}{dd});
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
                                                        if any(secdata.condbase_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.condbase_fieldoverlapind{xx}{ss});
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
%                                             fieldcolor = 'b';
                                        end
                                        sec = basedata.condbase_fieldcoord{ff}{dd};
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
                                    ax.Title.String = {['SI=' num2str(basedata.condbase_SIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_corr(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.condbase_SIC_adsm(ff) > pseudothr_corr(oo) && maxC >= 0.7
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
                                            map(basedata.condbase_map_rw{ff,1}(:,1)) = basedata.condbase_map_rw{ff,1}(:,4);
%                                         else
%                                             map = basedata.(['sec' mapname]){ff,1};
%                                         end
                                        maxC = nanmax(map(3:end));
                                        plotmap(map,msvar{2-oo+1});
                                        patchenvbounds(msvar{2-oo+1});
%                                         settoblack(map,msvar{2-oo+1});
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
                                        ax.Title.String = {horzcat(msvar{2-oo+1}(1),num2str(fff),' Infield: ', num2str(basedata.condbase_insecfieldrates{ff}(fff,1),2),'Hz'); ...
                                            horzcat('Outfield: ', num2str(prctile(basedata.condbase_outsecfieldrates{ff}{fff,1},95),2),'Hz')};
                                        if basedata.condbase_insecfieldrates{ff}(fff,1) > prctile(basedata.condbase_outsecfieldrates{ff}{fff,1},95)
                                            if secdata.condbase_insecfieldrates{fff}(ff,1) > prctile(secdata.condbase_outsecfieldrates{fff}{ff,1},95) % If reciprocal
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
%                                     settoblack(map,msvar{2-oo+1});
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
%                                             fieldcolor = 'y';
                                            xxx = find(basedata.condbase_fieldoverlapind{ff}{dd});
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
                                                        if any(secdata.condbase_fieldoverlapind{xx}{ss})
                                                            yy = find(secdata.condbase_fieldoverlapind{xx}{ss});
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
%                                             fieldcolor = 'b';
                                        end
                                        sec = basedata.condbase_fieldcoord{ff}{dd};
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
                                    ax.Title.String = {['SI=' num2str(basedata.condbase_SIC_adsm(ff),2) '/'...
                                        num2str(pseudothr_corr(oo),2)]};
                                    ax.Title.String{end+1,1} = cell2mat(overlap);
                                    if basedata.condbase_SIC_adsm(ff) > pseudothr_corr(oo) && maxC >= 0.7
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
    
elseif strcmp(objtype,'trajectory')
    
    plotview = 1;
    video = 1;
    
    for ii = 1:size(setsessions,1) % For each session
        cells_indList = find(identifiers(:,1) == setsessions(ii));

        for kk = 1:size(cells_indList,1)
            cell_ind = cells_indList(kk);
            cd(cellList{cell_ind});
            % Load spike data
            spktime = load('spiketrain.mat');
            spktime = spktime.timestamps/1000;
            cd ..; cd ..; cd ..;
            % Load behavioral data
            pv = load([num2str(pix) 'vmpv.mat']);
            pv = pv.pv.data;
            pvdata = pv.sessionTimeC;
            % Load trial structure and markers
            uma = load('umaze.mat');
            uma = uma.uma.data;
            uData = uma.unityData;
            uInd = uma.unityTriggers;
            uTime = uma.unityTime;
            
            % Get place trajectory with either binned or continuous data
            usingbinnedplacedata = false;
            if usingbinnedplacedata
                %% Drawing trajectory using binned place data
                % Convert place bin numbers to surf coords 
                placecol = 2;
                hdcol = 3;
                viewcol = 4;
                placex = nan(size(pvdata,1),4);
                placey = nan(size(pvdata,1),4);
                placez = zeros(size(pvdata,1),4);
                tic;
                disp('extracting coords');
                for pp = 1:size(pvdata,1)
                    if pvdata(pp,2) <= 0 
                        continue;
                    end
                    [whichgrid,binx,biny] = findgrid(pvdata(pp,placecol),'place');
                    [placex(pp,:) placey(pp,:) placez(pp,:)] = converttosurf(whichgrid,binx,biny);
                end
                toc;
                % Convert place coords from patch to point
                placecoord = [mean(placex,2) mean(placey,2) mean(placez,2)];
            else
                %% Drawing trajectory based on continuous place x y data
                % Load view data
                if plotview
                    viewdata = h5read([num2str(pix) 'binData.hdf'],'/data');
                    viewdata(1,:) = viewdata(1,:)/1000;
                else
                    viewdata = [0;nan];
                end
                vst = viewdata';
                % Sort view data into unityData structure (copied from vmpv.m)
                % getting combined sessiontime
                pst = [uma.unityData(1,:); uma.unityData];
                pst(1,2) = 0; % Insert row for start of session for binning
                pst(:,8) = cumsum(pst(:,2))+uma.unityTime(1);
                % cst: cumtime place_x place_y hd viewbin cst1/2
                cst_1 = [pst(:,[8 3 4 5]) nan(size(pst,1),1) ones(size(pst,1),1)];
                cst_2 = [vst(:,1) nan(size(vst,1),3) vst(:,2) 2.*ones(size(vst,1),1)];
                cst = [cst_2; cst_1];
                placecol = [2 3];
                hdcol = 4;
                viewcol = 5;

                % hybrid insertion sort instead of default sortrows for speed
                % does insertion all at one go, at the end
                % can possibly be improved through some sort of linked-list?
                disp('sorting in');
                tic;
                index_in = nan(1,size(cst_1,1));
                idx1 = 1;
                history1 = 1;
                for r = size(cst_2,1)+1:size(cst,1) % for each place time
                    for i = history1:size(cst_2,1) % for each view time
                        if cst(r,1) == cst(i,1) || cst(r,1) < cst(i,1) % if place time less than view time
                            history1 = i;
                            index_in(idx1) = i + idx1 - 1; % where to insert the place row in between view rows
                            break;
                        end
                    end            
                    if isnan(index_in(idx1))
                        index_in(idx1) = idx1 + size(cst_2,1);
                        history1 = size(cst_2,1);
                    end
%                     if rem(idx1, 10000) == 0
%                         disp([num2str(idx1) '/' num2str(size(cst_1,1))]);
%                     end
                    idx1 = idx1 + 1;
                end
                to_insert = cst_1;
                cst(setdiff(1:size(cst,1),index_in),:) = cst(1:size(cst_2,1),:); % setdiff one length more than expected here
                cst(index_in,:) = to_insert;
                disp(['sorting in done, time elapsed: ' num2str(toc)]);

                % cst = sortrows(cst, [1 4]); % sort by time, then put place in front
                disp('filling (no progress marker)');
                cst(:,2) = fillmissing(cst(:,2),'previous'); % fills in view rows with place_x it current is at
                cst(:,3) = fillmissing(cst(:,3),'previous'); % fills in view rows with place_y it current is at
                cst(:,4) = fillmissing(cst(:,4),'previous'); % fills in view rows with hd it current is at
                place_rows = find(cst(:,6) == 1);
                disp(['filling done, time elapsed: ' num2str(toc)]);
                disp('place pruning start (no progress marker)');
                for row = length(place_rows):-1:1 
                    % in the event that place changed at the same time as view's new sample, 
                    % we want to remove the place row, as the place information has
                    % already been transferred down to the subsequent view rows.
                    % note that visual inspection might see it as the same when it
                    % isn't, due to rounding in GUI.
                    if place_rows(row) ~= size(cst,1) % last row will have 0 and 0 for place and hd bins, and nan for view
                        if cst(place_rows(row),1) ~= cst(place_rows(row)+1,1)
                            place_rows(row) = [];
                        end
                    end
                end
                cst(place_rows,:) = [];
                disp(['place pruning end, time elapsed: ' num2str(toc)]);

                % now we need to replace surviving place rows, with duplicates of
                % itself, corresponding to the number of views in the previous time
                % bin
                place_rows = find(cst(:,6) == 1);
                first_non_nan = find(isnan(cst(:,5))==0);
                first_non_nan = first_non_nan(1);
                place_rows(place_rows < first_non_nan) = [];
                place_x = cst(place_rows,2);
                place_y = cst(place_rows,3);
                hd = cst(place_rows,4);
                reference_rows = place_rows - 1; 
                reference_times = cst(reference_rows,1); % timestamps from which to pull view bins from
                actual_times = cst(place_rows,1);

                % iterating and inserting rows takes too much time, below method to
                % make it feasible. first we count the number of duplicated place rows
                % to be added, then preallocate the final-sized array. 
                disp('counting max array size');
                membership = ismember(cst(:,1),reference_times);
                pvdata = nan(size(cst,1)+sum(membership)-length(place_rows),6);
                searching_portion = cst(membership,:);
                % size of gaps needed is calculated here
                insertion_gaps = nan(length(reference_times),1);
                for idx = 1:length(insertion_gaps)
%                     if rem(idx, 10000) == 0
%                         disp([num2str(idx) '/' num2str(length(insertion_gaps))]);
%                     end
                    insertion_gaps(idx) = sum(searching_portion(:,1)==reference_times(idx));
                end
                disp(['counting max done, time elapsed: ' num2str(toc)]);
                % slotting in original cst into enlarged one
                full_start = 1;
                original_start = 1;
                disp('slotting');
                for idx = 1:length(reference_rows)
%                     if rem(idx, 10000) == 0
%                         disp([num2str(idx) '/' num2str(length(reference_rows))]);
%                     end
                    chunk_to_insert = cst(original_start:reference_rows(idx),:);
                    pvdata(full_start:full_start-1+size(chunk_to_insert,1),:) = chunk_to_insert;
                    original_start = reference_rows(idx)+2;
                    full_start = full_start-1+size(chunk_to_insert,1)+insertion_gaps(idx)+1;
                end
                pvdata(full_start:end,:) = cst(reference_rows(end)+2:end,:);
                disp(['slotting end, time elapsed: ' num2str(toc)]);
                % constructing the new portion
                disp('inserting (no progress marker)');
                inserting_portion = searching_portion;
                inserting_portion(:,1) = repelem(actual_times,insertion_gaps);
                inserting_portion(:,2) = repelem(place_x,insertion_gaps);
                inserting_portion(:,3) = repelem(place_y,insertion_gaps);
                inserting_portion(:,4) = repelem(hd,insertion_gaps);
                inserting_portion(:,6) = 6;
                pvdata(isnan(pvdata(:,1)),:) = inserting_portion;
                disp(['inserting end, time elapsed: ' num2str(toc)]);
                if ~isempty(find(diff(pvdata(:,1))<0))
                    error('combined sessiontime misaligned!');
                end
                pvdata = pvdata(:,1:5);
                disp('guaranteeing unique and ascending');
                dti = find(diff(pvdata(:,1)));
                dti(:,2) = [1; 1 + dti(1:end-1,1)];
                dti = dti(:,[2 1]);
                dti = [dti; [dti(end,2)+1 size(pvdata,1)]];
                to_remove = [];
                for chunk = 1:size(dti, 1)
                    pvdata(dti(chunk, 1): dti(chunk,2),:) = sortrows(pvdata(dti(chunk, 1): dti(chunk,2),:), [1 5]);
                    identify_dup = diff(pvdata(dti(chunk, 1): dti(chunk,2),:)); % duplicate consecutive place rows
                    if dti(chunk, 1) - dti(chunk, 2) ~= 0
                        if ~isempty(find(sum(identify_dup,2)==0))
                            to_remove = [to_remove; find(sum(identify_dup,2)==0)+dti(chunk,1)];
                        end
                    end
                end
                pvdata(to_remove,:) = [];
                disp(['duplicates found: ' num2str(length(to_remove))]);
                disp(['guaranteeing unique and ascending done, time elapsed: ' num2str(toc)]);
                % Downsample 
                placecoord = [pvdata(:,placecol(1)) pvdata(:,placecol(2)) zeros(size(pvdata,1),1)];
            end
            cumtime = pvdata(:,1);
            % Convert view data from bins to surf coords
            viewx = nan(size(pvdata,1),4);
            viewy = nan(size(pvdata,1),4);
            viewz = nan(size(pvdata,1),4);
            tic;
            disp('extracting coords');
            for pp = 1:size(pvdata,1)
                % Leave as NaNs the times when view is on cue, hint or nan
                if pvdata(pp,viewcol) < 3 || isnan(pvdata(pp,viewcol))
                    continue;
                end
                [whichgrid,binx,biny] = findgrid(pvdata(pp,viewcol),'view');
                [viewx(pp,:) viewy(pp,:) viewz(pp,:)] = converttosurf(whichgrid,binx,biny);
            end
            if ~usingbinnedplacedata
                viewx = (viewx/40)*25-12.5;
                viewy = (viewy/40)*25-12.5;
%                 viewz = (viewz/40)*25-12.5;
            end
            toc;
            % Convert head direction bins to rad
            if usingbinnedplacedata
                deg = mean([pvdata(:,hdcol)*6 (pvdata(:,hdcol)-1)*6],2);
                deg(deg<0) = 0;
            else
                deg = pvdata(:,hdcol);
            end
            deg = 360-deg+90; % Convert to increasing CCW and starting from x axis
            deg(deg>360) = deg(deg>360)-360;
            rad = deg2rad(deg); % Convert cartesian to polar coords
            [dx,dy] = pol2cart(rad,2);
            % Bin spikes
            spktime = spktime';
            spktime(spktime(:,1) > pvdata(end,1),:) = [];   
            spk = histcounts(spktime, pvdata(:,1))'; 
            spk = [spk; 0]; 
            % Fill in spikes for view radius > 1 px
            spk(:,2) = [0;diff(pvdata(:,1))]; 
            spk(:,3) = nan(size(spk,1),1); 
            spk(spk(:,2)~=0,3) = 0; 
            spk(spk(:,1)>0,3) = spk(spk(:,1)>0,1); 
            spk(:,3) = fillmissing(spk(:,3),'next'); 
            spk(:,1:2) = [];
            % Insert trial info into sessionTime structure
            trialnums = nan(size(pvdata,1),1);
            trialtarget = nan(size(pvdata,1),1);
            for pp = 1:size(uInd,1)
                starttime = uTime(uInd(pp,2)+1);
                endtime = uTime(uInd(pp,3));
                trialnums(pvdata(:,1)>=starttime & pvdata(:,1)<=endtime) = pp;
                trialtarget(pvdata(:,1)>=starttime & pvdata(:,1)<=endtime) = rem(uData(uInd(pp,2)),10);
            end
            combinedtarget = cell(2,size(uInd,1)); % trial target and trial direction
            trialdirlookup = {'O' 'N' 'E' 'N' 'N' 'N'; % Cat to others 
                              'S' 'O' 'E' 'E' 'E' 'S'; % Cam to others
                              'W' 'W' 'O' 'N' 'N' 'W'; % Rab to others
                              'S' 'W' 'S' 'O' 'S' 'S'; % Donk to others
                              'S' 'W' 'S' 'N' 'O' 'W'; % Croc to others
                              'S' 'N' 'E' 'N' 'E' 'O' % Pig to others
                              }; 
            for pp = 1:size(uInd,1)
                combinedtarget{1,pp} = rem(uData(uInd(pp,2)),10);
                if pp == 1
                    combinedtarget{2,pp} = trialdirlookup{1,rem(uData(uInd(pp,2)),10)};
                else
                    combinedtarget{2,pp} = trialdirlookup{rem(uData(uInd(pp-1,2)),10),rem(uData(uInd(pp,2)),10)};
                end
            end
            % Compute mean trial heading in deg
            meantrialheading = nan(size(uInd,1),1);
            for jj = 1:size(uInd,1)
                meantrialheading(jj) = nanmean(deg(trialnums == jj));
            end

            %% Set up figure 
            f = figure;
            set(f,'Units','Normalized','Position',[0 0 1 1]);
            ax = gca;
            h = animatedline(ax,'Tag','wholetrajectory');
            ptrailcolor = [111 111 111]/255;
            trailalpha = 0.3;
            vtrailcolor = [176 196 222]/255;
            h.Color = ptrailcolor;
            if plotview
                if usingbinnedplacedata
                    patchenvbounds('view');
                    ax.XLim = [0 40];
                    ax.YLim = [0 40];
                    ax.ZLim = [0 40];
                else
                    patchenvbounds('view_vu');
                    ax.XLim = [-12.5 12.5];
                    ax.YLim = [-12.5 12.5];
                    ax.ZLim = [0 40];
                end
                view(-35,20);
            else
                if usingbinnedplacedata
                    patchenvbounds('place');
                    ax.XLim = [0 40];
                    ax.YLim = [0 40];
                    ax.ZLim = [0 40];
                else
                    patchenvbounds('place_vu');
                    ax.XLim = [-12.5 12.5];
                    ax.YLim = [-12.5 12.5];
                    ax.ZLim = [0 40];
                end
            end
            axis square;
            axis off;
            hold on;
            
            % Set up marker for body location
            first = find(trialnums == 1,1)+1;
            boddir = plot([placecoord(first,1)-dx(first) placecoord(first,1)],...
                [placecoord(first,2)-dy(first) placecoord(first,2)],'b','LineWidth',1);
            bod = plot(placecoord(first,1),placecoord(first,2),'b.','MarkerSize',80);
%             bodtrajcoords = [placecoord(first,1) placecoord(first,2)];
            % Set up markers for trail of body location
            skip = 100; % downsampling - plot one every 'skip' points
            ptrail = plot(repmat(placecoord(first,1),10*skip,1),repmat(placecoord(first,2),10*skip,1),'k.','MarkerSize',10);
            % Set up marker for view fixation point
            vertices = [viewx(first,:)' viewy(first,:)' viewz(first,:)'];
            fixfaces = [1 2 3 4];
            % Set up markers for trail of view location
            vtrailvertices = vertices;
            vtrailfaces = fixfaces;
            % Patch first pixel of trajectory
            patch('Vertices',vtrailvertices,'Faces',vtrailfaces,...
                                'FaceColor',vtrailcolor,'EdgeColor',vtrailcolor,'Tag','totalviewtrail');
            % Patch first fixation point
            viewpt = patch('Vertices',vertices,'Faces',fixfaces,'FaceColor','b','Tag','Fixation');
%             uistack(findobj('Tag','Fixation'),'top');
            % Label trial number
            trialnum = annotation('textbox',[0.0 0.85 0.4 0.1],'String',['Trial 1'],...
                'FitBoxToText','on','HorizontalAlignment','center','FontSize',20);
            % Label time elapsed in trial
            trialtime = annotation('textbox',[0.0 0.8 0.4 0.1],'String',['Time = 0'],...
                'FitBoxToText','on','HorizontalAlignment','center','FontSize',20);
            % Label target identity
            target = annotation('textbox',[0.0 0.75 0.4 0.1],'String',[],...
                'FitBoxToText','on','HorizontalAlignment','center','FontSize',20);
            % Label number of spikes
            numspk = annotation('textbox',[0.0 0.70 0.4 0.1],'String',[],...
                'FitBoxToText','on','HorizontalAlignment','center','FontSize',20);
            targetnames = {'Cat','Camel','Rabbit','Donkey','Crocodile','Pig'};
            
            % Filter which trials to use
            trialfilter = 'subset';
            switch trialfilter
                case 'none'
                    includetrials = (1:size(uInd,1))';
                case 'subset'
                    includetrials = ([1:5])';
                case 'target_cat'
                    includetrials = unique(trialnums(trialtarget == 1));
                case 'target_camel'
                    includetrials = unique(trialnums(trialtarget == 2));
                case 'target_rabbit'
                    includetrials = unique(trialnums(trialtarget == 3));
                case 'target_donkey'
                    includetrials = unique(trialnums(trialtarget == 4));
                case 'target_crocodile'
                    includetrials = unique(trialnums(trialtarget == 5));
                case 'target_pig'
                    includetrials = unique(trialnums(trialtarget == 6));
                case 'north'
                    includetrials = find(strcmp(combinedtarget(2,:),'N'))';
                case 'west'
                    includetrials = find(strcmp(combinedtarget(2,:),'W'))';
                case 'south'
                    includetrials = find(strcmp(combinedtarget(2,:),'S'))';
                case 'east'
                    includetrials = find(strcmp(combinedtarget(2,:),'E'))';
            end
            % Restrict number of trials shown when recording movie
            if size(includetrials,1) > 30
                includetrials = randsample(includetrials,30);
            end
            
            % debug
            trialpoints = 0;
            % Record frames for movie
            if video
                videoname = ['Video ' objtype,' ',num2str(setsessions(ii)),'ch',num2str(identifiers(cell_ind,4)),'c',num2str(identifiers(cell_ind,5))];
                v = VideoWriter(videoname,'MPEG-4');
                open(v);
            end
            mov = struct('cdata',[],'colormap',[]);
            for jj = 1:size(includetrials,1) % For each trial
                t = includetrials(jj);
                
                % Dim spikes of last trial
                if exist('spkmarker_p','var')
                    ax = gca;
                    for pp = 1:spkcount
                        spkmarker_p{pp}.MarkerEdgeColor = [250 160 160]/255; % pink
                        spkmarker_p{pp}.MarkerSize = 20;
                        spkline{pp}.Color = [250 160 160]/255; % pink
                        spkmarker_v{pp}.Color = [255 255 194]/255; % sun
                        spkmarker_v{pp}.MarkerSize = 20;
                    end
                end
                % Change color of view trajectory from last trial to grey
                delete(findobj('Tag','currentviewtrail'));
                delete(findobj('Tag','totalviewtrail'));
                patch('Vertices',vtrailvertices,'Faces',vtrailfaces,'FaceColor',ptrailcolor,...
                    'FaceAlpha',trailalpha,'EdgeColor',ptrailcolor,'EdgeAlpha',trailalpha,'Tag','totalviewtrail');
                
                spkcount = 0;
                spkmarker_p = [];
                spkline = [];
                spkmarker_v = [];
                bodtrajcoords = [];
                % Get trial data
                inds = find(trialnums == t);
                % Update trial time
                trialtime.String = ['Time = 0'];
                % Update target identity
                targetnum = rem(uData(uInd(t,1),1),10);
                target.String = targetnames{targetnum};
                % Update number of spikes
                trialspk = sum(spk(inds));
                numspk.String = ['Trial spikes = ' num2str(trialspk)];
                % Highlight target pillar face
                delete(findobj('Tag','target'));
                patchenvbounds([lower(target.String) '_vu']);
                targface = findobj('Tag',target.String);
                targface.Tag = 'target';
                
                % Update trial number
                trialnum.String = ['Trial ' num2str(t)];
                if sum(cumtime(inds(end),1)-cumtime(inds(1)))>=25 % Mark incomplete trials
                    trialnum.Color = 'r';
                else
                    trialnum.Color = 'k';
                end
                % Mark current body location
                % Flash green to mark start of trial on trajectory
                bod.Color = [0 1 0];
                viewpt.FaceColor = [0 1 0];
                viewpt.EdgeColor = [0 1 0];
                % Change current body location
                bod.XData = placecoord(inds(1),1);
                bod.YData = placecoord(inds(1),2);
                if jj == 1
                    bodtrajcoords = [placecoord(inds(1),1) placecoord(inds(1),2)];
                end
                % Change current body direction
                boddir.XData = [placecoord(inds(1),1)-dx(inds(1)) placecoord(inds(1),1)];
                boddir.YData = [placecoord(inds(1),2)-dy(inds(1)) placecoord(inds(1),2)];
                % Reset animated line so it doesn't join old trajectory to new starting point
                if ~strcmp(trialfilter,'none') 
%                     % All this is to reduce number of children in axis so plotting is not slowed down
%                     % Remove previous trial animated line
%                     delete(findobj('Tag',['trial' num2str(includetrials(jj-1)) 'trajectory']));
%                     % Add previous trial animated line to whole trajectory
%                     addpoints(h,bodtrajcoords(:,1),bodtrajcoords(:,2));
%                     drawnow;
                    % Assign new animated line to current trial trajectory
                    h1 = animatedline(ax,'Tag',['trial' num2str(t) 'trajectory']);
                    h1.Color = ptrailcolor;
                else
                    h1 = findobj('Tag','wholetrajectory');
                end
                % Pad with points from start of this trial to pause drawing
                for kk = 1:20
                    addpoints(h1,placecoord(inds(1),1),placecoord(inds(1),2));
                    drawnow;
                end
                % Return body marker color to blue
                bod.Color = [0 0 1];
                % Return view marker color to blue
                viewpt.FaceColor = [0 0 1];
                viewpt.EdgeColor = [0 0 1];
                % Find how many points to plot before speeding up
                left = rem(size(inds,1),skip);
                if left == 0
                    if skip == 1
                        left = 1;
                    else
                        left = skip-1;
                    end
                end
%                 % Add first points of trial to body trajectory
%                 bodtrajcoords = [placecoord(inds(1:left-1),1) placecoord(inds(1:left-1),2)];

                % Downsample plotting
                count = 1;
                % Store frame for movie
                mov(1) = getframe(gcf);
                for pp = left:skip:size(inds,1)
                    % Update trial time
                    trialtime.String = ['Time = ' num2str(cumtime(inds(pp),1)-cumtime(inds(1),1),3)];
%                     trialtime.String = ['pp = ' num2str(pp)];
                    if pp == left
                        chunksize = left;
                    else
                        chunksize = skip;
                    end
                    inds_ds = inds(pp-chunksize+1:pp);
                    % Change current body location
                    bod.XData = placecoord(inds(pp),1);
                    bod.YData = placecoord(inds(pp),2);
                    % Change current body direction
                    boddir.XData = [placecoord(inds(pp),1)-dx(inds(pp)) placecoord(inds(pp),1)];
                    boddir.YData = [placecoord(inds(pp),2)-dy(inds(pp)) placecoord(inds(pp),2)];
                    % Add to total place trajectory
                    addpoints(h1,placecoord(inds_ds,1),placecoord(inds_ds,2));
                    bodtrajcoords = [bodtrajcoords; placecoord(inds_ds,1) placecoord(inds_ds,2)];
                    
                    % Find view trajectory for this trial
                    viewbins = pvdata(inds_ds,5);
                    viewinds = logical([1;diff(viewbins)~=0]);
%                     viewinds(find(viewinds,1,'last')) = false; % So that fixation point is not double-plotted
                    viewtrajx = viewx(inds_ds,:)'; 
                    viewtrajy = viewy(inds_ds,:)';
                    viewtrajz = viewz(inds_ds,:)';
                    % Downsample view trajectory to only unique fixation points
                    viewtrajx = viewtrajx(:,viewinds);
                    viewtrajy = viewtrajy(:,viewinds);
                    viewtrajz = viewtrajz(:,viewinds);
                    vtrailvertices_trial = [reshape(viewtrajx,size(viewtrajx,1)*size(viewtrajx,2),1) ...
                         reshape(viewtrajy,size(viewtrajy,1)*size(viewtrajy,2),1) ...
                         reshape(viewtrajz,size(viewtrajz,1)*size(viewtrajz,2),1)];
                    vtrailfaces_trial = reshape(1:size(vtrailvertices_trial,1),4,size(vtrailvertices_trial,1)/4)';
                    % Add to total view trajectory
                    vtrailvertices = [vtrailvertices; vtrailvertices_trial];
                    addfaces = reshape(vtrailfaces(end,end)+[1:size(vtrailvertices_trial,1)],4,size(vtrailvertices_trial,1)/4)';
                    vtrailfaces = [vtrailfaces; addfaces];
%                     viewpt = patch('Vertices',[viewx(inds(pp),:)' viewy(inds(pp),:)' viewz(inds(pp),:)'],...
%                         'Faces',fixfaces,'FaceColor','b','Tag','Fixation');
                    % Trail current place trajectory
                    ptrail.XData(1:end-size(inds_ds,1)) = ptrail.XData(size(inds_ds,1)+1:end);
                    ptrail.YData(1:end-size(inds_ds,1)) = ptrail.YData(size(inds_ds,1)+1:end);
                    ptrail.XData(end-size(inds_ds,1)+1:end) = placecoord(inds_ds,1);
                    ptrail.YData(end-size(inds_ds,1)+1:end) = placecoord(inds_ds,2);
                    % Trail current view trajectory
                    patch('Vertices',vtrailvertices_trial,'Faces',vtrailfaces_trial,'FaceColor',vtrailcolor,...
                        'EdgeColor',vtrailcolor,'Tag','currentviewtrail');
                    % Move current view fixation point 
                    viewpt.Vertices = [viewx(inds(pp),:)' viewy(inds(pp),:)' viewz(inds(pp),:)'];
                    uistack(viewpt,'top');
                    
                    % Overlay spikes of this trial
                    spkind = find(spk(inds_ds));
                    if size(spkind,1) > 0
                        % Plot spikes
                        for ss = 1:size(spkind,1) % For each time interval
                            ns = spk(inds_ds(spkind(ss)));
                            for sss = 1:ns % For each spike
                                spkcount = spkcount + 1;
                                % Spike overlaid on place
                                spkmarker_p{spkcount} = plot(placecoord(inds_ds(spkind(ss)),1),placecoord(inds_ds(spkind(ss)),2),'r.','MarkerSize',60);
                                spkline{spkcount} = plot([placecoord(inds_ds(spkind(ss)),1)-dx(inds_ds(spkind(ss))) placecoord(inds_ds(spkind(ss)),1)],...
                                    [placecoord(inds_ds(spkind(ss)),2)-dy(inds_ds(spkind(ss))) placecoord(inds_ds(spkind(ss)),2)],'r');
                                % Spike overlaid on view
                                spkx = viewx(inds_ds(spkind(ss)),:)';
                                spky = viewy(inds_ds(spkind(ss)),:)';
                                spkz = viewz(inds_ds(spkind(ss)),:)';
                                spkmarker_v{spkcount} = plot3(mean(spkx),mean(spky),mean(spkz),'.','MarkerSize',60,'Color',[255 232 124]/255);
                            end
                        end
                    end
                    count = count + 1;
                    drawnow;
                    % Store frame for movie
                    mov(end+1) = getframe(gcf);
                end
                trialpoints(t) = count;
                % Flash black to mark end of trial on trajectory
                bod.Color = [0 0 0];
                viewpt.FaceColor = [0 0 0];
                viewpt.EdgeColor = [0 0 0];
                % Pad with points from end of last trial to pause drawing
                for kk = 1:20
                    addpoints(h1,bod.XData,bod.YData);
                    drawnow;
                    % Store frame for movie
                    mov(end+1) = getframe(gcf);
                end
                if trialspk ~= spkcount
                    disp('error in number of spikes for trial');
                end
%                 % Remove prior fixation points (this is only necessary when
%                 % plotting the view trajectory with opaque points
%                 delete(findobj('Tag','Fixation'));
            end
            if video
                writeVideo(v,mov);
            end
        end
        
        
    end
 
end

% selcell_orig = cellList(selcell_orig);
selcell_orig = objMain.data.origin(selcell_orig);
% selcell_corr = cellList(selcell_corr);
selcell_corr = objMain.data.origin(selcell_corr);
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
        patch([0 0 40 40],[0 40 40 0],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Floor');
        % Floor Pillar edges
        patch([8 8 16 16],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([8 8 16 16],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 32 32],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 32 32],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        
    case 'place_vu' % center of maze at 0, size 25x25 virtual units
        % Floor outer bounds
        patch([-12.5 -12.5 12.5 12.5],[-12.5 12.5 12.5 -12.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Floor');
        % Floor Pillar edges
        patch([-7.5 -7.5 -2.5 -2.5],[-7.5 -2.5 -2.5 -7.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-7.5 -7.5 -2.5 -2.5],[2.5 7.5 7.5 2.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 7.5 7.5],[2.5 7.5 7.5 2.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 7.5 7.5],[-2.5 -7.5 -7.5 -2.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        
    case 'view'
        
        % Floor outer bounds
        patch([0 0 40 40],[0 40 40 0],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Floor');
        % Floor Pillar edges
        patch([8 8 16 16],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([8 8 16 16],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 32 32],[24 32 32 24],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 32 32],[8 16 16 8],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        % Wall outer bounds
        patch([0 0 40 40],[0 0 0 0],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        patch([0 0 0 0],[0 0 40 40],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        patch([0 0 40 40], [40 40 40 40],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        patch([40 40 40 40],[0 0 40 40],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        % Pillar 1 (bottom right)
        patch([24 24 32 32],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 24 24],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 32 32],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([32 32 32 32],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([32 32 32 32],[10.88 10.88 13.12 13.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none','Tag','Rabbit'); % Rabbit Poster on m_wall_25
        % Pillar 2 (bottom left)
        patch([8 8 16 16],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([8 8 8 8],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([8 8 16 16],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([16 16 16 16],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([10.88 10.88 13.12 13.12],[8 8 8 8],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none','Tag','Cat'); % Cat poster on m_wall_10
        patch([10.88 10.88 13.12 13.12],[16 16 16 16],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none','Tag','Pig'); % Pig poster on m_wall_29
        % Pillar 3 (top right)
        patch([24 24 32 32],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 24 24],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([24 24 32 32],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([32 32 32 32],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([26.88 26.88 29.12 29.12],[24 24 24 24],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none','Tag','Croc'); % Croc poster on m_wall_4
        patch([26.88 26.88 29.12 29.12],[32 32 32 32],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none','Tag','Donk'); % Donkey poster on m_wall_15
        % Pillar 4 (top left)
        patch([8 8 16 16],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([8 8 8 8],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([8 8 16 16],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([16 16 16 16],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([8 8 8 8],[26.88 26.88 29.12 29.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none','Tag','Camel'); % Camel poster on m_wall_20
        % Ceiling
        patch([0 0 40 40],[0 40 40 0],[40 40 40 40],[1 1 1 1],'FaceColor','none','Tag','Ceiling');
        
    case 'view_vu'
        % Floor outer bounds
        patch([-12.5 -12.5 12.5 12.5],[-12.5 12.5 12.5 -12.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Floor');
        % Floor Pillar edges
        patch([-7.5 -7.5 -2.5 -2.5],[-7.5 -2.5 -2.5 -7.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-7.5 -7.5 -2.5 -2.5],[2.5 7.5 7.5 2.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 7.5 7.5],[2.5 7.5 7.5 2.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 7.5 7.5],[-2.5 -7.5 -7.5 -2.5],[0 0 0 0],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        % Wall outer bounds
        patch([-12.5 -12.5 12.5 12.5],[-12.5 -12.5 -12.5 -12.5],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        patch([-12.5 -12.5 -12.5 -12.5],[-12.5 -12.5 12.5 12.5],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        patch([-12.5 -12.5 12.5 12.5], [12.5 12.5 12.5 12.5],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        patch([12.5 12.5 12.5 12.5],[-12.5 -12.5 12.5 12.5],[16 24 24 16],[1 1 1 1],'FaceColor','none','Tag','Wall');
        % Pillar 1 (bottom right)
        patch([2.5 2.5 7.5 7.5],[-7.5 -7.5 -7.5 -7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 2.5 2.5],[-7.5 -7.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 7.5 7.5],[-2.5 -2.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([7.5 7.5 7.5 7.5],[-7.5 -7.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
%         patch([32 32 32 32],[10.88 10.88 13.12 13.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Rabbit Poster on m_wall_25
        % Pillar 2 (bottom left)
        patch([-7.5 -7.5 -2.5 -2.5],[-7.5 -7.5 -7.5 -7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-7.5 -7.5 -7.5 -7.5],[-7.5 -7.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-7.5 -7.5 -2.5 -2.5],[-2.5 -2.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-2.5 -2.5 -2.5 -2.5],[-7.5 -7.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
%         patch([10.88 10.88 13.12 13.12],[8 8 8 8],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Cat poster on m_wall_10
%         patch([10.88 10.88 13.12 13.12],[16 16 16 16],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Pig poster on m_wall_29
        % Pillar 3 (top right)
        patch([2.5 2.5 7.5 7.5],[2.5 2.5 2.5 2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 2.5 2.5],[2.5 2.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([2.5 2.5 7.5 7.5],[7.5 7.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([7.5 7.5 7.5 7.5],[2.5 2.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
%         patch([26.88 26.88 29.12 29.12],[24 24 24 24],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Croc poster on m_wall_4
%         patch([26.88 26.88 29.12 29.12],[32 32 32 32],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Donkey poster on m_wall_15
        % Pillar 4 (top left)
        patch([-7.5 -7.5 -2.5 -2.5],[2.5 2.5 2.5 2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-7.5 -7.5 -7.5 -7.5],[2.5 2.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-7.5 -7.5 -2.5 -2.5],[7.5 7.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
        patch([-2.5 -2.5 -2.5 -2.5],[2.5 2.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','Tag','Pillar');
%         patch([8 8 8 8],[26.88 26.88 29.12 29.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Camel poster on m_wall_20
        % Ceiling
        patch([-12.5 -12.5 12.5 12.5],[-12.5 12.5 12.5 -12.5],[40 40 40 40],[1 1 1 1],'FaceColor','none','Tag','Ceiling');
    case 'cat_vu'
        % Pillar 2 (bottom left)
        patch([-7.5 -7.5 -2.5 -2.5],[-7.5 -7.5 -7.5 -7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','EdgeColor','g','Tag','Cat','LineWidth',2);
    case 'camel_vu'
        % Pillar 4 (top left)
        patch([-7.5 -7.5 -7.5 -7.5],[2.5 2.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','EdgeColor','g','Tag','Camel','LineWidth',2);
    case 'donkey_vu'
        % Pillar 3 (top right)
        patch([2.5 2.5 7.5 7.5],[7.5 7.5 7.5 7.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','EdgeColor','g','Tag','Donkey','LineWidth',2);
    case 'crocodile_vu'
        % Pillar 3 (top right)
        patch([2.5 2.5 7.5 7.5],[2.5 2.5 2.5 2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','EdgeColor','g','Tag','Crocodile','LineWidth',2);
    case 'pig_vu'
        % Pillar 2 (bottom left)
        patch([-7.5 -7.5 -2.5 -2.5],[-2.5 -2.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','EdgeColor','g','Tag','Pig','LineWidth',2);
    case 'rabbit_vu'
        % Pillar 1 (bottom right)
        patch([7.5 7.5 7.5 7.5],[-7.5 -7.5 -2.5 -2.5],[16 21 21 16],[1 1 1 1],'FaceColor','none','EdgeColor','g','Tag','Rabbit','LineWidth',2);
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

function [] = savefigure(saveoption,h,figtitle,figdir)
cwd = pwd;
cd(figdir);
if saveoption
    saveas(h,figtitle,'png');
    % saveas(h,figtitle,'fig');
    % print('-painters',figtitle,'-dvg');
    % print([figdir '/' figtitle],'-depsc');
end
cd(cwd);
close(figure(h));


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

% pick colormap
function [set] = setcolor(choice)

switch choice
    case 'twotonered'
        set = [ 0.4,    0.4,    0.4;...
                        0.425,  0.425,  0.425;...
                        0.45,   0.45,   0.45;...
                        0.475,  0.475,  0.475;...
                        0.5,    0.5,    0.5;...
                        0.525,  0.525,  0.525;...
                        0.55,   0.55,   0.55;...
                        0.575,  0.575,  0.575;...
                        0.6,    0.6,    0.6;...
                        0.625,  0.625,  0.625;...
                        0.65,   0.65,   0.65;...
                        0.675,  0.675,  0.675;...
                        0.7,    0.7,    0.7;...
                        0.725,  0.7,    0.7;...
                        0.75,   0.7,    0.7;...
                        0.775,  0.7,    0.7;...
                        0.8,    0.7,    0.7;...
                        0.825,  0.7,    0.7;...
                        0.85,   0.7,    0.7;...
                        0.875,  0.7,    0.7;...
                        0.9,    0.7,    0.7;...
                        0.925,  0.7,    0.7;...
                        0.95,   0.7,    0.7;...
                        0.975,  0.7,    0.7;...
                        1,      0.7,    0.7;...
                        1,      0.675,  0.675;
                        1,      0.65,   0.65;...
                        1,      0.625,  0.625;...
                        1,      0.6,    0.6;...
                        1,      0.575,  0.575;
                        1,      0.55,   0.55;...
                        1,      0.525,  0.525;...
                        1,      0.5,    0.5;...
                        1,      0.475,  0.475;...
                        1,      0.45,   0.45;...
                        1,      0.425,  0.425;...
                        1,      0.4,    0.4;...
                        1,      0.375,  0.375;...
                        1,      0.35,   0.35;...
        %                 1,      0.325,  0.325;...
                        1,      0.3,    0.3;...
        %                 1,      0.275,  0.275;...
                        1,      0.25,   0.25;...
        %                 1,      0.225,  0.225;...
                        1,      0.2,    0.2;...
        %                 1,      0.175,  0.175;...
                        1,      0.15,   0.15;...
        %                 1,      0.125,  0.125;...
                        1,      0.1,    0.1;...
        %                 1,      0.075,  0.075;...
                        1,      0.05,   0.05;...
        %                 1,      0.025,  0.025;...
                        1,      0,      0];
    case 'coolwarm32'

        set = [  59	76	192
                        68	91	205
                        78	105	216
                        88	118	226
                        99	132	235
                        110	144	242
                        121	156	248
                        133	168	252
                        144	178	254
                        155	188	255
                        167	196	254
                        178	204	251
                        188	210	247
                        198	215	241
                        208	218	234
                        217	220	226
                        225	219	215
                        233	214	203
                        239	207	190
                        243	200	178
                        246	191	165
                        247	180	152
                        247	169	139
                        245	157	126
                        242	144	113
                        237	130	101
                        230	115	89
                        223	99	78
                        214	82	67
                        204	63	57
                        192	41	47
                        180	4	38]./255;
    case 'coolwarm64'
            
        set = [  59	76	192
                        63	83	199
                        68	90	205
                        73	97	210
                        78	104	216
                        83	111	221
                        88	118	226
                        93	124	230
                        98	131	234
                        104	137	238
                        109	143	242
                        115	149	245
                        120	155	247
                        126	161	250
                        131	167	252
                        137	172	253
                        142	177	254
                        148	182	255
                        154	186	255
                        159	191	255
                        165	195	254
                        170	199	253
                        176	202	252
                        181	206	250
                        186	209	248
                        191	211	246
                        196	214	243
                        201	216	240
                        206	217	236
                        210	219	232
                        215	220	228
                        219	220	223
                        223	220	218
                        227	218	212
                        231	215	206
                        234	212	200
                        237	209	194
                        240	206	188
                        242	202	182
                        244	198	175
                        245	194	169
                        246	189	163
                        247	184	156
                        247	179	150
                        247	173	143
                        247	168	137
                        246	162	131
                        245	155	124
                        243	149	118
                        241	142	112
                        239	135	106
                        236	128	100
                        233	121	94
                        230	113	88
                        226	106	83
                        222	98	77
                        218	89	72
                        213	81	67
                        208	72	61
                        203	62	57
                        198	52	52
                        192	41	47
                        186	27	43
                        180	4	38]./255;
end
            
