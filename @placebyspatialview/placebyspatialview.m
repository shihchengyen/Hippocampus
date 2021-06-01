function [obj, varargout] = placebyspatialview(savefig,varargin)


peekinsidepillar = 0;
filttype = 'FiltVel';
pix = 1;

if nargin > 1 % Single cell
    cwd = varargin{1};
    cellList = {cwd};
else % Batch 
    % Load cell list
%     cwd = '/Users/yuhsuan/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis';
    cwd='/Users/yuhsuan/Desktop';
%     fid = fopen([cwd '/cell_list.txt'],'rt');
    fid = fopen([cwd '/cell_list copy.txt'],'rt');
%     fid = fopen([cwd '/cell_list_1pxFiltAll.txt'],'rt');
    cellList = textscan(fid,'%s','Delimiter','\n');
    cellList = cellList{1};
end

% Make sure no empty cells
notempty = ~cellfun(@isempty,cellList);
cellList = cellList(notempty,:);
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
% setcells = unique(cellid);

% Remove cell directories that already exist 
figdir = ['/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Figures/' filttype '/' num2str(pix) 'px/PlaceView/'];
if savefig    
    if exist(figdir,'dir') ~= 7
        mkdir(figdir);
    else
        rmdir(figdir,'s');
        mkdir(figdir);
    end
    % Create txt file to save list of selective cells
    cd(figdir);
    fid = fopen('sigcells.txt','a');
    cd(cwd);
end

sigcells = {};

% Find SI/ISE thresholds from population
cwd = pwd;
c_objdir = ['/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Combined Objects/' filttype '/' num2str(pix) 'px/'];
cd(c_objdir);
c_pc = load('c_vmpc.mat');
c_pc = c_pc.vmp.data;
c_sv = load('c_vmsv.mat');
c_sv = c_sv.vms.data;
pcSIthr = prctile([c_pc.SIC; c_pc.SICsh(:)],95);
svSIthr = prctile([c_sv.SIC; c_sv.SICsh(:)],95);
cd(cwd);

% Load vmpv object for each session
for ss = 1:size(setsessions,1)

    cd(['/Volumes/Hippocampus/Data/picasso-misc/' num2str(setsessions(ss)) '/session01']);
    pv = load([num2str(pix) 'vmpv.mat']);
    pv = pv.pv;
    setcells = find(identifiers(:,1) == setsessions(ss));
    
    % Plot places he's standing when viewing inside of pillar
    if peekinsidepillar
        cd(cellList{setcells(1)});
        placebyspatialview_insidepillar(pv,objtype,identifiers(setcells(1),:),savefig,figdir,filttype,pix);
        close all;
    end
    
    % Process placebyspatialview cell by cell

    for cc = 1:size(setcells,1)
        
        cd(cellList{setcells(cc)});
        disp(cellList{setcells(cc)});
        [figureout,data,celltype] = placebyspatialview_cell(pv,identifiers(setcells(cc),:),savefig,figdir,filttype,pix,pcSIthr,svSIthr);
        if figureout
            sigcells{end+1,1} = cellList{setcells(cc)};
            ID = [num2str(identifiers(setcells(cc),1)) 'ch' num2str(identifiers(setcells(cc),4)) 'c' num2str(identifiers(setcells(cc),5))];
            fprintf(fid,'%s\n',[ID celltype]);
        end
        close all;
        
    end
    
end
disp(sigcells);
    


function [figureout,data,celltype] = placebyspatialview_cell(pv,identifier,savefig,figdir,filttype,pix,pcSIthr,svSIthr)

ID = [num2str(identifier(1)) 'ch' num2str(identifier(4)) 'c' num2str(identifier(5))];
figureout = 0;
fieldsize = 5;
objtype = {'place','spatialview'};
data = struct;
celltype = '';

cwd = pwd;
objdir = [cwd '/' filttype '/' num2str(pix) 'px/'];
cd(objdir);

% Load vmpc object
pc = load('vmpc.mat');
pc = pc.vmp.data;

% Load vmsv object
sv = load('vmsv.mat');
sv = sv.vms.data;
cd(cwd);

% Find SI/ISE for this cell
pcSI = pc.SIC;
svSI = sv.SIC;

% Plot session place and view maps
if pcSI > pcSIthr || svSI > svSIthr

    % Record cell
    figureout = 1;
    
    if pcSI > pcSIthr && svSI > svSIthr
        celltype = 'pv';
    elseif pcSI > pcSIthr && svSI < svSIthr
        celltype = 'p';
    elseif pcSI < pcSIthr && svSI > svSIthr
        celltype = 'v';
    end
    figdir = [figdir ID celltype];
    if savefig
        if exist(figdir,'dir') ~= 7
            mkdir(figdir);
        else
            rmdir(figdir,'s');
            mkdir(figdir);
        end
    end

    % Load spike train
    spiketrain = load('spiketrain.mat');
    spiketrain = spiketrain.timestamps ./ 1000; % in seconds

    % Combine place and view info with spikes and make rate maps
    pvT = pv.data.sessionTimeC;
    pvT(:,4) = [diff(pvT(:,1)); 0];
    binned = histcounts(spiketrain, pvT(:,1))';
    pvT(:,5) = [binned; 0];
    % Filter out segments 
    switch filttype
        case 'FiltAll'
            pvT(~get(pv,'SpeedLimit',pc.Args.ThresVel),:) = []; % Velocity < threshold
            pvT(~ismember(pvT(:,2),pv.data.place_good_bins),:) = []; % Num place obs < MinObsPlace
            pvT(~ismember(pvT(:,3),pv.data.view_good_bins),:) = []; % Num view obs < MinObsView
        case 'FiltVel'
            pvT(~get(pv,'SpeedLimit',pc.Args.ThresVel),:) = []; % Velocity < threshold
        case 'FiltObs'
            pvT(~ismember(pvT(:,2),pv.data.place_good_bins),:) = []; % Num place obs < MinObsPlace
            pvT(~ismember(pvT(:,3),pv.data.view_good_bins),:) = []; % Num view obs < MinObsView
    end
    pvT(pvT(:,2)==0,:) = []; % ITI

    % Create base for backfilling 
    pvTfill = pvT;

    view_durations = NaN(5122,1600);
    view_spikes = NaN(5122,1600);
    place_durations = NaN(1,1600);
    place_spikes = NaN(1,1600);

    for i = 1:1600

        inds = pvT(:,2)==i;
        subsample = [pvT(inds, [3 4 5])];
        if ~isempty(subsample)
            disp(i);
        end

        % Get spikes and duration for place only
        place_durations(1,i) = sum(subsample(:,2));
        place_spikes(1,i) = sum(subsample(:,3));

        % back-filling spikes for view
        subsample(subsample(:,3)==0,3) = nan;
        subsample(:,4) = circshift(subsample(:,2)~=0 ,-1);
        subsample(isnan(subsample(:,3)) & subsample(:,4), 3) = 0;
        subsample(:,4) = [];
        subsample(:,3) = fillmissing(subsample(:,3), 'next');
        % back-filling time for view
        subsample(subsample(:,2)==0,2) = nan;
        subsample(:,2) = fillmissing(subsample(:,2), 'previous');
        % Put backfill into sessionTimeC array
        pvTfill(inds,[3 4 5]) = subsample;

        % padding with 5122 bin
        subsample = [subsample; [5122 0 0]];
        % remove bad view spots
        subsample(isnan(subsample(:,1)),:) = [];
        % sum durations
        view_durations(:,i) = accumarray(subsample(:,1), subsample(:,2),[],[],NaN);
        % sum spikes
        view_spikes(:,i) = accumarray(subsample(:,1), subsample(:,3),[],[],NaN); 
    end
%     gpd = accumarray(pvT(:,2),pvT(:,4),[1600 1])';
%     spk = accumarray(pvT(:,2),pvT(:,5),[1600 1])';
    % Plot full place maps
    h = figure(11);
    ax = gca;
    h.Name = [ID 'vmpvRawPlaceMap'];
    set(h,'Units','normalized','Position',[0 0 1 1]);
    % emptyplacegrids = all(isnan(full_rate),1);
    % rawplacemap1 = nansum(full_rate,1);
    % rawplacemap1(1,emptyplacegrids) = NaN;
    rawplacemap1 = place_spikes./place_durations;
    rawplacemap1(place_durations == 0) = NaN;
    plotmap(rawplacemap1,'place');
    ax.Title.String = ['vmpvRawPlaceMap: SI ' num2str(pcSI) '/' num2str(pcSIthr)];
    set(ax,'CLim',[0 max(rawplacemap1)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
        'XColor','none','YColor','none','ZColor','none',...
        'FontSize',14,'GridLineStyle','none','Color','none');
    patchenvbounds('place');
    if savefig
        savefigure(h,h.Name,figdir);
    end

    % Plot full place maps
    h = figure(12);
    ax = gca;
    h.Name = [ID 'vmpcRawPlaceMap'];
    set(h,'Units','normalized','Position',[0 0 1 1]);
    rawplacemap2 = pc.maps_raw;
    [rawplacemapG,~] = plotmap(rawplacemap2,'place');
    ax.Title.String = ['vmpcRawPlaceMap: SI ' num2str(pcSI) '/' num2str(pcSIthr)];
    set(ax,'CLim',[0 max(rawplacemap2)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
        'XColor','none','YColor','none','ZColor','none',...
        'FontSize',14,'GridLineStyle','none','Color','none');
    patchenvbounds('place');
    if savefig
        savefigure(h,h.Name,figdir);
    end

    h = figure(13);
    ax = gca;
    h.Name = [ID 'vmpcSmoothPlaceMap'];
    set(h,'Units','normalized','Position',[0 0 1 1]);
    smoothplacemap = pc.maps_adsm;
    [smoothplacemapG,placemapGdummy] = plotmap(smoothplacemap,'place');
    % [smoothplacemapG,placemapGdummy] = plotplacemap(smoothplacemap); % Do not use grid output as input for another plot as it will turn out rotated 90deg CCW
    ax.Title.String = ['vmpcSmoothPlaceMap: SI ' num2str(pcSI) '/' num2str(pcSIthr)];
    set(ax,'CLim',[0 max(smoothplacemap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
        'XColor','none','YColor','none','ZColor','none',...
        'FontSize',14,'GridLineStyle','none','Color','none');
    patchenvbounds('place');
    if savefig
        savefigure(h,h.Name,figdir);
    end
    data.place.smoothmapL = smoothplacemap;
    data.place.smoothmapG = smoothplacemapG;
    data.place.dummymapG = placemapGdummy;

    % Plot full view maps
    h = figure(21); 
    ax = gca;
    h.Name = [ID 'vmpvRawViewMap'];
    set(h,'Units','normalized','Position',[0 0 1 1]);
    spikes = nansum(view_spikes,2);
    durations = nansum(view_durations,2);
    rawviewmap1 = spikes./durations;
    rawviewmap1(durations==0) = nan;
    rawviewmap1 = emptyinsidepillar(rawviewmap1); % Temporary measure only! Remove data from inside of pillar where it should be empty
    plotmap(rawviewmap1,'spatialview');
    ax.Title.String = ['vmpvRawViewMap: SI ' num2str(svSI) '/' num2str(svSIthr)];
    set(ax,'CLim',[0 max(rawviewmap1(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
        'XColor','none','YColor','none','ZColor','none',...
        'FontSize',14,'GridLineStyle','none','Color','none');
    patchenvbounds('spatialview');
    if savefig
        savefigure(h,h.Name,figdir);
    end

    % Plot full view maps
    rawviewmap2 = sv.maps_raw;
    rawviewmap2 = emptyinsidepillar(rawviewmap2); % Temporary measure only! Remove data from inside of pillar where it should be empty
    h = figure(22);
    ax = gca;
    h.Name = [ID 'vmsvRawViewMap'];
    set(h,'Units','normalized','Position',[0 0 1 1]);
    [rawviewmapG,~] = plotmap(rawviewmap2,'spatialview');
    ax.Title.String = ['vmsvRawViewMap: SI ' num2str(svSI) '/' num2str(svSIthr)];
    set(ax,'CLim',[0 max(rawviewmap2(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
        'XColor','none','YColor','none','ZColor','none',...
        'FontSize',14,'GridLineStyle','none','Color','none');
    patchenvbounds('spatialview');
    if savefig
        savefigure(h,h.Name,figdir);
    end

    % Plot full view maps
    smoothviewmap = sv.maps_adsm;
    smoothviewmap = emptyinsidepillar(smoothviewmap); % Temporary measure only! Remove data from inside of pillar where it should be empty
    h = figure(23);
    ax = gca;
    h.Name = [ID 'vmsvSmoothViewMap'];
    set(h,'Units','normalized','Position',[0 0 1 1]);
    [smoothviewmapG,viewmapGdummy] = plotmap(smoothviewmap,'spatialview');
    ax.Title.String = ['vmsvSmoothViewMap: SI ' num2str(svSI) '/' num2str(svSIthr)];
    set(ax,'CLim',[0 max(smoothviewmap(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
        'XColor','none','YColor','none','ZColor','none',...
        'FontSize',14,'GridLineStyle','none','Color','none');
    patchenvbounds('spatialview');
    if savefig
        savefigure(h,h.Name,figdir);
    end
    data.spatialview.smoothmapL = smoothviewmap;
    data.spatialview.smoothmapG = smoothviewmapG;
    data.spatialview.dummymapG = viewmapGdummy;

    % Find fields in both place and view maps
    for oo = 1:size(objtype,2)

        cd(cwd);
        switch objtype{oo}
            case 'place'
                SI = pcSI;
                SIthr = pcSIthr;
                basename = 'Place';
                secname = 'View';
                basemapLsm = smoothplacemap;
                secmapLsm = smoothviewmap;
                basemapGsm = {smoothplacemapG};
                basemapGrw = {rawplacemapG};
                dummygrid = {placemapGdummy};
                peakrate_full = nanmax(smoothplacemap);
                peakrate_set = nanmax(smoothplacemap);
                prI = 1;
            case 'spatialview'
                SI = svSI;
                SIthr = svSIthr;
                basename = 'View';
                secname = 'Place';
                basemapLsm = smoothviewmap;
                secmapLsm = smoothplacemap;
                basemapGsm = smoothviewmapG;
                basemapGrw = rawviewmapG;
                dummygrid = viewmapGdummy;
                for ii = 1:size(smoothviewmapG,1)
                    maxset(ii) = nanmax(reshape(smoothviewmapG{ii},size(smoothviewmapG{ii},1)*size(smoothviewmapG{ii},2),1));
                end
                [peakrate_full,prI] = max(maxset);
                peakrate_set = max(maxset(3:end));
        end
        
        % Find 3 maxima of rate maps
        count = 0;
        for gg = 1:size(basemapGsm,1)
            % Skip cue and hint view maps 
            if size(basemapGsm{gg},1) == 1
                continue;
            end
            % Find fields with at least 1 pixel of > 70% peak rate
            ind_fields = basemapGsm{gg} > 0.7*peakrate_full;
            % Find separate fields
            [fieldlabel,fieldcount] = bwlabel(ind_fields,4);
            % For walls and pillars, if there are fields that wrap around
            % split, merge them.
            if gg > 4
               % Find possible split fields
               if any(fieldlabel(:,1)) && any(fieldlabel(:,end))
                   boundary = [fieldlabel(:,1) fieldlabel(:,end)];
                   [blabel,bcount] = bwlabel(boundary,8);
                   for bb = 1:bcount
                       label = unique(fieldlabel(blabel(:,1)==bcount(bb),1));
                       label = label(1);
                       labelset = unique(boundary);
                       labelset = labelset(labelset>0);
                       labelset = setdiff(labelset,label);
                       for ll = 1:length(labelset)
                          fieldlabel(fieldlabel == labelset(ll)) = label; 
                          fieldcount = fieldcount - 1;
                       end
                   end
               end
            end
            % Split fields that are too large i.e. > 1/4 of area
            for ii = 1:fieldcount
                inds = fieldlabel == ii;
                % Split fields that are too large
                if sum(inds(:)) >= (size(basemapLsm,2)/4) 
                    subinds = inds & basemapGsm{gg} > 0.8*peakrate_full;
                    [sublabel,subcount] = bwlabel(subinds,4);
                    sublabel(subinds) = sublabel(subinds) + fieldcount;
                    fieldcount = fieldcount + subcount;
                    fieldlabel(inds) = 0;
                    fieldlabel(inds) = sublabel(inds);
                end
            end
            % Find fields that are big enough (i.e. > 15 bins around peak, 5x5bins area)
            for ii = 1:fieldcount
                inds = fieldlabel == ii;
                % Make sure field size is > 15 bins 
                if sum(inds(:)) >= 15
                    count = count + 1;
                    fieldmaxrate_rw(count,1) = max(basemapGrw{gg}(inds)); 
                    fieldmaxrate_sm(count,1) = max(basemapGsm{gg}(inds));
                    switch objtype{oo}
                        case 'place'
                            gridnum(count,1) = 3;
                        case 'spatialview'
                            gridnum(count,1) = gg;
                    end
                    [fieldcoordx_mat, fieldcoordy_mat] = find(inds);
                    fieldcoord{count,1} = [fieldcoordy_mat size(basemapGrw{gg},1)-fieldcoordx_mat+1]; % In plot coords. x left to right, y bottom to top
                    linbin{count,1} = dummygrid{gg}(inds);
                end
            end
        end
        % If no significant fields, skip to next cell
        if count == 0
            return;
        end

        % Sort fields
        [fieldmaxrate_sm,I] = sort(fieldmaxrate_sm,'descend');
        fieldmaxrate_rw = fieldmaxrate_rw(I);
        gridnum = gridnum(I);
        fieldcoord = fieldcoord(I);
        linbin = linbin(I);

        for ii = 1:3 % Plot only for first 3 fields per base map

            if ii > size(fieldcoord,1)
                break;
            end

%             % Plot base map and field of interest
%             fignum = str2double([num2str(oo) '0' num2str(ii)]);
%             h = figure(fignum);
%             ax = gca;
%             h.Name = [ID ' Base' basename ' : Field ' num2str(ii)];
%             set(h,'Units','normalized','Position',[0 0 1 1]);
%             plotmap(data.(objtype{oo}).smoothmapL,objtype{oo});
%             % Patch base pixels
%             for pp = 1:size(fieldcoord{ii},1)
%                 [x,y,z] = converttosurf(gridnum(ii),fieldcoord{ii}(pp,1),fieldcoord{ii}(pp,2));
%                 patch(x,y,z,[1 1 1 1],'EdgeColor','k','FaceColor','none','LineWidth',1);
%             end
%             set(ax,'CLim',[0 max(data.(objtype{oo}).smoothmapL(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%                 'XColor','none','YColor','none','ZColor','none',...
%                 'FontSize',14,'GridLineStyle','none','Color','none');
%             ax.Title.String = ['Base' basename ': SI ' num2str(SI) '/' num2str(SIthr)];
%             % Patch environment boundaries
%             patchenvbounds(objtype{oo});
%             if savefig
%                 savefigure(h,h.Name,figdir);
%             end

            % Initialise variables
            rate_components_field = cell(size(fieldcoord{ii},1),1);
            maxrate_px = nan(size(fieldcoord{ii},1),1);
            secs = cell(size(fieldcoord{ii},1),1);
            usecpxs = [];
            % Find pixels of secondary map
            for pp = 1:size(fieldcoord{ii},1)
                %Get corresponding secondary pixels from pv object
                sec = [];
                switch objtype{oo}
                    case 'place'
                        ind_pv = pvT(:,2) == linbin{ii}(pp);
                        sec(:,1) = pvT(ind_pv,3); % view px
                        sec(:,2) = pvT(ind_pv,4); % dur
                        sec(:,3) = pvT(ind_pv,5); % spikes
                    case 'spatialview'
                        ind_pv = pvTfill(:,3) == linbin{ii}(pp);
                        sec(:,1) = pvTfill(ind_pv,2); % place px
                        sec(:,2) = pvTfill(ind_pv,4); % dur
                        sec(:,3) = pvTfill(ind_pv,5); % spikes
                end
                % Get firing rates 
                usecpx = unique(sec(:,1));
                usecpx(isnan(usecpx)) = [];
                rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
                rate_components_px(:,1) = usecpx;
                ind_timechange = find(sec(:,2) > 0);
                for cc = 1:size(ind_timechange,1) % For each instance of being in this base pixel
                    linbintemp = nan(size(usecpx,1),2); % Temp dur and spikes for this secondary pixel(s)
                    if cc == size(ind_timechange,1)
                        newind = ind_timechange(cc):size(sec,1);
                    else
                        newind = ind_timechange(cc):ind_timechange(cc+1)-1; % index into secondary pixel(s) for this instance
                    end
                    newview = sec(newind,1); 
                    newset = ismember(usecpx,newview);
                    linbintemp(newset,1) = sec(newind(1),2); % duration for this instance listed with first secondary pixel
                    linbintemp(newset,2) = sec(newind(end),3); % spikes for this instance listed with last secondary pixel
                    rate_components_px(:,2) = nansum( [rate_components_px(:,2) linbintemp(:,1)] ,2); % Sum duration for this sec pixel across instances
                    rate_components_px(:,3) = nansum( [rate_components_px(:,3) linbintemp(:,2)] ,2); % Sum spikes for this sec pixel across instances
                end
                rightfulnans = rate_components_px(:,2) == 0;
                rate_components_px(rightfulnans,:) = nan;
                rate_components_px(:,4) = rate_components_px(:,3)./rate_components_px(:,2); % Compute firing rates

                % Collect all sec rate information for the base pixels
                rate_components_field{pp} = ( rate_components_px );
                if ~isempty(rate_components_px)
                    maxrate_px(pp) = max(rate_components_px(:,4));
                else 
                    maxrate_px(pp) = NaN;
                end
                secs{pp} = sec;
                usecpxs = union(usecpxs,usecpx); % set of secondary pixels covered in whole session
                
            end
            disp(ii);

            session_seclinbin = usecpxs; 
            session_seclinbin(:,2:3) = NaN;
            % Consolidate occupancy and spikes across base field
            for pp = 1:size(fieldcoord{ii},1)
                tempbins = nan(size(session_seclinbin,1),2);
                tempbins(ismember(session_seclinbin(:,1),rate_components_field{pp}(:,1)),1) = rate_components_field{pp}(:,2); % Duration 
                tempbins(ismember(session_seclinbin(:,1),rate_components_field{pp}(:,1)),2) = rate_components_field{pp}(:,3); % Spikes 
                session_seclinbin(:,2) = nansum( [session_seclinbin(:,2) tempbins(:,1)] ,2); % Sum duration across session
                session_seclinbin(:,3) = nansum( [session_seclinbin(:,3) tempbins(:,2)] ,2); % Sum spikes across session
            end
            session_seclinbin(:,4) = session_seclinbin(:,3)./session_seclinbin(:,2);

            
            
            rate_components_full{ii,1} = rate_components_field;
            seclinbin_full{ii,1} = session_seclinbin;
%             secgridmap_full{ii,1} = session_secgridmap;

        end

        % Store data
        data.(objtype{oo}).SI = SI;
        data.(objtype{oo}).SIthr = SIthr;
        data.(objtype{oo}).basemapLsm = basemapLsm;
        data.(objtype{oo}).secmapLsm = secmapLsm;
        data.(objtype{oo}).basemapGsm = basemapGsm;
        data.(objtype{oo}).basemapGrw = basemapGrw;
        data.(objtype{oo}).dummygrid = dummygrid;
        data.(objtype{oo}).peakrate_full = peakrate_full;
        data.(objtype{oo}).peakrate_set = peakrate_set;
        data.(objtype{oo}).fieldmaxrate_sm = fieldmaxrate_sm;
        data.(objtype{oo}).fieldmaxrate_rw = fieldmaxrate_rw;
        data.(objtype{oo}).gridnum = gridnum;
        data.(objtype{oo}).fieldcoord = fieldcoord;
        data.(objtype{oo}).linbin = linbin;
        data.(objtype{oo}).set_sec_linbin = seclinbin_full;
%         data.(objtype{oo}).set_sec_gridmap = secgridmap_full;
        data.(objtype{oo}).rate_components = rate_components_full;

        clear fieldmaxrate_sm; clear fieldmaxrate_rw; clear gridnum; clear fieldcoord; clear linbin;
        clear seclinbin_full; % clear secgridmap_full; 
        clear rate_components_full; clear dummygrid; clear SI; clear SIthr;

        disp(ii);
    end

    % Figures and Stats
    for oo = 1:size(objtype,2)
        switch objtype{oo}
            case 'place'
                basename = 'Place';
                secname = 'View';
            case 'spatialview'
                basename = 'View';
                secname = 'Place';
        end
        
        % Plot base map
        basefignum = str2double([num2str(oo) '01']);
        h = figure(basefignum);
        ax = gca;
        h.Name = [ID ' Base' basename ' : ' num2str(size(data.(objtype{oo}).rate_components,1)) 'Fields'];
        set(h,'Units','normalized','Position',[0 0 1 1]);
        % Plot
        plotmap(data.(objtype{oo}).smoothmapL,objtype{oo});
        set(ax,'CLim',[0 max(data.(objtype{oo}).smoothmapL(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
            'XColor','none','YColor','none','ZColor','none',...
            'FontSize',14,'GridLineStyle','none','Color','none');
        ax.Title.String = ['Base' basename ': SI ' num2str(data.(objtype{oo}).SI) '/' num2str(data.(objtype{oo}).SIthr)];
        % Patch environment boundaries
        patchenvbounds(objtype{oo});
        % Patch base fields
        for ii = 1:size(data.(objtype{oo}).rate_components,1) % For each base field
            for pp = 1:size(data.(objtype{oo}).fieldcoord{ii},1)
                [x,y,z] = converttosurf(data.(objtype{oo}).gridnum(ii),data.(objtype{oo}).fieldcoord{ii}(pp,1),data.(objtype{oo}).fieldcoord{ii}(pp,2));
                patch(x,y,z,[1 1 1 1],'EdgeColor','k','FaceColor','none','LineWidth',1);
            end 

        end
        if savefig
            savefigure(h,h.Name,figdir);
        end
        
        % Plot consolidated pixel maps and secondary fields
        for ii = 1:size(data.(objtype{oo}).rate_components,1) % For each base field

            session_seclinbin = data.(objtype{oo}).set_sec_linbin{ii}; % Get set of secondary pixels covered in whole session
%             session_secgridmap = data.(objtype{oo}).set_sec_gridmap{ii};
            dummygridbase = data.(objtype{oo}).dummygrid;
            dummygridsec = data.(objtype{2-oo+1}).dummygrid;
            fieldcoord = data.(objtype{oo}).fieldcoord;
            gridnum = data.(objtype{oo}).gridnum;
            
            % Plot consolidated pixel maps
            fignum = str2double([num2str(oo) '0' num2str(ii) '1']);
            h = figure(fignum);
            h.Name = [ID ' Base' basename ' : Field ' num2str(ii) ' Sec Pixels Full'];
            set(h,'Units','normalized','Position',[0 0 1 1]);

            for ff = 1%:3
                ax = gca;
    %                 ax = subplot(1,3,ff);
                tempsecmap = nan(size(data.(objtype{oo}).secmapLsm));
                switch ff
                    case 1 % Rate
                        tempsecmap(1,session_seclinbin(:,1)) = session_seclinbin(:,4);
                        ax.Title.String = 'Rate';
                    case 2 % Spike
                        tempsecmap(1,session_seclinbin(:,1)) = session_seclinbin(:,3);
                        ax.Title.String = 'Spikes';
                    case 3 % Duration
                        tempsecmap(1,session_seclinbin(:,1)) = session_seclinbin(:,2);
                        caxis([min(session_seclinbin(session_seclinbin(:,2)>0,2)) max(session_seclinbin(:,2))]);
                        ax.Title.String = 'Duration';                        
                end
                % Plot 
                [session_secgridmap,~] = plotmap(tempsecmap,objtype{2-oo+1});
                if strcmp(objtype{2-oo+1},'place')
                    session_secgridmap = {session_secgridmap};
                end
                colormap(ax,cool);
                % If firing rate or spike count = 0, set to black
                blackpx = find(tempsecmap == 0);
                for pp = 1:size(blackpx,2)
                    switch objtype{2-oo+1}
                        case 'place'
                            gnum = 1;
                            [gnum,xx,yy] = findgrid(blackpx(pp),objtype{2-oo+1});
                        case 'spatialview'
                            if blackpx(pp) > 3 % Make sure not cue or hint
                               [gnum,xx,yy] = findgrid(blackpx(pp),objtype{2-oo+1});
                            else
                                continue;
                            end

                    end
                    [x,y,z] = converttosurf(gnum,xx,yy);
                    patch(x,y,z,[1 1 1 1],'FaceColor','k');
                end
                set(ax,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                ax.Title.FontSize = 16;
                % Patch environment boundaries
                patchenvbounds('spatialview');
                % Patch basemap field
                for pp = 1:size(fieldcoord{ii},1)
                    [x,y,z] = converttosurf(gridnum(ii),fieldcoord{ii}(pp,1),fieldcoord{ii}(pp,2));
                    if strcmp(objtype{oo},'place')
                        patch(x,y,z,[1 1 1 1],'EdgeColor','r','FaceColor','none','LineWidth',2);
                    elseif strcmp(objtype{oo},'spatialview')
    %                         patch(x,y,z,'r','EdgeColor',[107/256 142/256 35/256],'FaceColor','none','LineWidth',2); % Olive drab
    %                         patch(x,y,z,'r','EdgeColor',[60/256 179/256 113/256],'FaceColor','none','LineWidth',2); % Medium sea green
                        patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',2); % Lime green
    %                         patch(x,y,z,'r','EdgeColor','g','FaceColor','none','LineWidth',2);
                    end
                end
                % Patch secondary fields
                for jj = 1:size(data.(objtype{2-oo+1}).rate_components,1) % For each secondary field
                    % Plot secondary field in consolidated pixel map
                    sec = data.(objtype{2-oo+1}).fieldcoord{jj}; % Get bins for sec field
                    % Patch secondary fields in the consolidated pixel map
                    fignum = str2double([num2str(oo) '0' num2str(ii) '1']);
                    h = figure(fignum); 
                    for ff = 1%:3
                        ax = gca;
    %                     ax = subplot(1,3,ff);
                        for pp = 1:size(sec,1)
                            [x,y,z] = converttosurf(data.(objtype{2-oo+1}).gridnum(jj),sec(pp,1),sec(pp,2));
                            patch(x,y,z,'r','EdgeColor',[169/256 169/256 169/256],'FaceColor','none','LineWidth',2); % Dark Grey
                            if strcmp(objtype{2-oo+1},'place')
                                patch(x,y,z,'r','EdgeColor','r','FaceColor','none','LineWidth',2); % Red
                            elseif strcmp(objtype{2-oo+1},'spatialview')
                                patch(x,y,z,'r','EdgeColor',[50/256 205/256 50/256],'FaceColor','none','LineWidth',2); % Lime green
                            end
                        end
                    end
                    if savefig
                        savefigure(h,h.Name,figdir);
                    end
                end
                

            end
            
            % Stats
            for jj = 1:size(data.(objtype{2-oo+1}).rate_components,1) % For each secondary field
                
                linbin = data.(objtype{2-oo+1}).linbin{jj}; % pixels that make up the sec field. Not all of these will be sampled from this base field
                inds_infield = ismember(session_seclinbin(:,1),linbin);
                if isempty(inds_infield)
                    continue;
                end
                meanrate_infield = sum(session_seclinbin(inds_infield,3)) / sum(session_seclinbin(inds_infield,2));
                % Sample 1000 outfields
                meanrate_outfield = nan(10000,1);
                session_seclinbin_out = session_seclinbin(~inds_infield,:);
                
                % If all pixels are within sec field
                if isempty(session_seclinbin_out)
                    % Plot histogram
                    h = figure(str2double(['9' num2str(oo) num2str(ii) num2str(jj)]));
                    h.Name = [ID 'Base' basename ' field ' num2str(ii) ' insufficient px w Sec ' secname ' field ' num2str(jj)];
                else
                    % Generate a psuedopopulation of outfields same size as sec field 
                    for ff = 1:10000
                        % Start from a random pixel
                        startpx = randsample(1:size(session_seclinbin,1),1);
                        startpx = session_seclinbin(startpx,1);
                        % Find which gridmap this sec start px is on
                        switch objtype{2-oo+1}
                            case 'place'
                                gnum = 1;
                            case 'spatialview'
                                if startpx == 1 || startpx == 2 % Make sure not cue or hint
                                   ff = ff - 1;
                                   continue;
                                end
                                [gnum,~,~] = findgrid(startpx,objtype{2-oo+1});
                        end
                        % If random px overlaps with sec field, repeat
                        if ismember(startpx,linbin)
                            ff = ff - 1;
                            continue;
                        end
                        % Get the actual sampled sec grid map
                        tempmap = session_secgridmap{gnum};
                        % Get the grid coords of starting px
                        [startindx,startindy] = find(dummygridsec{gnum} == startpx);
                        % Expand radius around starting px until hit the requisite number of px
%                         while sum(sum(~isnan(tempmap(startindx,startindy)))) < sum(inds_infield)
                        while length(startindx)*length(startindy) < size(linbin,1)
                            if startindx(1) > 1 && startindx(end) < size(dummygridsec{gnum},1)
                                startindx = [startindx(1)-1 startindx startindx(end)+1];
                            elseif startindx(1) == 1 && startindx(end) < size(dummygridsec{gnum},1)
                                startindx = [startindx startindx(end)+1];
                            elseif startindx(1) > 1 && startindx(end) == size(dummygridsec{gnum},1)
                                startindx = [startindx(1)-1 startindx];
                            end
                            % If reach required num of px
                            if length(startindx)*length(startindy) > size(linbin,1)
                                break;
                            end
                            if startindy(1) > 1 && startindy(end) < size(dummygridsec{gnum},2)
                                startindy = [startindy(1)-1 startindy startindy(end)+1];
                            elseif startindy(1) == 1 && startindy(end) < size(dummygridsec{gnum},2)
                                startindy = [startindy startindy(end)+1];
                            elseif startindy(1) > 1 && startindy(end) == size(dummygridsec{gnum},2)
                                startindy = [startindy(1)-1 startindy];
                            end
                            % If exceed map bounds
                            if startindx(1) == 1 && startindx(end) == size(dummygridsec{gnum},1) && startindy(1) == 1 && startindy(end) == size(dummygridsec{gnum},2)
                                break;
                            end
                        end
                        % get another starting pixel if cannot sample enough px
                        if length(startindx)*length(startindy) < size(linbin,1)
                            ff = ff - 1;
                            continue;
                        end
                        % Remove excess px
                        pxsub = dummygridsec{gnum}(startindx,startindy);
                        inds_sampled = ~isnan(tempmap(startindx,startindy));
                        pxsub = reshape(pxsub,size(pxsub,1)*size(pxsub,2),1);
                        inds_sampled = reshape(inds_sampled,size(inds_sampled,1)*size(inds_sampled,2),1);
                        inds_keep = randsample(1:length(pxsub),size(linbin,1))';
                        pxsub = pxsub(inds_keep);
                        inds_sampled = inds_sampled(inds_keep);
                        % Start over if exactly the same px as field of interest
                        if isempty(setdiff(pxsub,linbin))
                            ff = ff - 1;
                            continue;
                        end
                        % Get sec bins actually sampled 
                        outfieldlinbin = pxsub(inds_sampled);

                        % draw out in the map my outfield px
%                         for pp = 1:size(outfieldlinbin,1)
%                             [gridnum,xx,yy] = findgrid(outfieldlinbin(pp));
%                             [x,y,z] = converttosurf(gridnum,xx,yy);
%                             patch(x,y,z,'r','EdgeColor','y','FaceColor','none','LineWidth',2);
%                         end

                        sess_inds_outfield = ismember(session_seclinbin(:,1),outfieldlinbin);
                        meanrate_outfield(ff,1) = sum(session_seclinbin(sess_inds_outfield,3)) / sum(session_seclinbin(sess_inds_outfield,2));
                        
                    end
                    meanrate_thr = prctile(meanrate_outfield(:),95);
                    % Plot histogram
                    h = figure(str2double(['9' num2str(oo) num2str(ii) num2str(jj)]));
                    if meanrate_infield > meanrate_thr % If significantly more infield firing than outfield
                        h.Name = [ID 'Base' basename ' field ' num2str(ii) ' sig linked w Sec ' secname ' field ' num2str(jj)];
                    else
                        h.Name = [ID 'Base' basename ' field ' num2str(ii) ' not linked w Sec ' secname ' field ' num2str(jj)];
                    end
                    hist(meanrate_outfield);
                    vline(meanrate_thr,'k:');
                    vline(meanrate_infield,'r');
                    ax = gca;
                    ax.Title.String = horzcat('Meanrate infield ',num2str(meanrate_infield),' / Meanrate outfield 95th prc ',num2str(meanrate_thr));
                    set(ax,'FontSize',14);
                end
                if savefig
                    savefigure(h,h.Name,figdir);
                end
                
            end

        end
    end
    % Save data
    if savefig
        cd(figdir);
        Args.classname = 'placebyspatialview';
        Args.matname = [Args.classname '.mat'];
        Args.matvarname = 'vmpxv';
        Args.SaveLevels = 1;
        n = nptdata(1,0,pwd);
        d.data = data;
        obj = class(d,Args.classname,n);
        saveObject(obj,'ArgsC',Args);
        cd(cwd);
    end
else
    disp('No significant place or view fields');
end

function [] = placebyspatialview_insidepillar(pv,objtype,identifier,savefig,figdir,filttype,pix)

ID = [num2str(identifier(1)) 'ch' num2str(identifier(4)) 'c' num2str(identifier(5))];

% Load vmpc object
pc = load('vmpc.mat');
pc = pc.vmp.data;

% Load vmsv object
sv = load('vmsv.mat');
sv = sv.vms.data;

% Load spike train
spiketrain = load('spiketrain.mat');
spiketrain = spiketrain.timestamps ./ 1000; % in seconds

% Combine place and view info with spikes and make rate maps
pvT = pv.data.sessionTimeC;
pvT(:,4) = [0; diff(pvT(:,1))];
binned = histcounts(spiketrain, pvT(:,1))';
pvT(:,5) = [binned; 0];
pvT(pvT(:,2)==-1,:) = [];
pvT(pvT(:,2)==0,:) = [];

% Create base for backfilling %%%%%%%%%% BUT CHECK good_rows!!!!
pvTfill = pvT;

view_durations = NaN(5122,1600);
view_spikes = NaN(5122,1600);
place_durations = NaN(1,1600);
place_spikes = NaN(1,1600);

for i = 1:1600
    
    inds = pvT(:,2)==i;
    subsample = [pvT(inds, [3 4 5])];
    if ~isempty(subsample)
        disp(i);
    end
    
    % Get spikes and duration for place only
    place_durations(1,i) = sum(subsample(:,2));
    place_spikes(1,i) = sum(subsample(:,3));

    % back-filling spikes for view
    subsample(subsample(:,3)==0,3) = nan;
    subsample(:,4) = circshift(subsample(:,2)~=0 ,-1);
    subsample(isnan(subsample(:,3)) & subsample(:,4), 3) = 0;
    subsample(:,4) = [];
    subsample(:,3) = fillmissing(subsample(:,3), 'next');
    % back-filling time for view
    subsample(subsample(:,2)==0,2) = nan;
    subsample(:,2) = fillmissing(subsample(:,2), 'previous');
    % Put backfill into sessionTimeC array
    pvTfill(inds,[3 4 5]) = subsample;

    % padding with 5122 bin
    subsample = [subsample; [5122 0 0]];
    % remove bad view spots
    subsample(isnan(subsample(:,1)),:) = [];
    % sum durations
    view_durations(:,i) = accumarray(subsample(:,1), subsample(:,2),[],[],NaN);
    % sum spikes
    view_spikes(:,i) = accumarray(subsample(:,1), subsample(:,3),[],[],NaN); 
end

% Find 3 maxima of rate maps
switch objtype
    case 'place'
        basemap = pc.maps_adsm;
        secmap = sv.maps_adsm;
        [smoothplacemapG,placemapGdummy] = plotmap(pc.maps_adsm,'place');
        close all;
        basemapG = {smoothplacemapG};
        SIthr = prctile(pc.SICsh(1,2:end),95);
        SI = pc.SIC;
        dummygrid = {placemapGdummy};
        peakrate = nanmax(reshape(basemapG{1},size(basemapG{1},1)*size(basemapG{1},2),1));
        prI = 1;
    case 'spatialview'
        basemap = sv.maps_adsm;
        secmap = pc.maps_adsm;
        [smoothviewmapG,viewmapGdummy] = plotmap(sv.maps_adsm,'spatialview');
        close all;
        basemapG = smoothviewmapG;
        SIthr = prctile(sv.SICsh(1,2:end),95);
        SI = sv.SIC;
        dummygrid = viewmapGdummy;
        for ii = 1:size(basemapG,1)
            maxset(ii) = nanmax(reshape(basemapG{ii},size(basemapG{ii},1)*size(basemapG{ii},2),1));
        end
        [peakrate,prI] = max(maxset);
end
   
threshrate = peakrate/2;
count = 0;
gg = 3;

xmapGset = [10 13 16 26 29 32];
ymapGset = [10 13 16 26 29 32];
% xmapGset = [10 13 16 ];
% ymapGset = [10 13 16 ];

for xx = 1:size(xmapGset,2)
    for yy = 1:size(xmapGset,2)
        
        xmapG = xmapGset(xx);
        ymapG = ymapGset(yy);
        
        xbounds_mapG = [xmapG+1 xmapG xmapG-1]; % Plotting coords (x is left to right, y is bottom to top)
        ybounds_mapG = [ymapG-1 ymapG ymapG+1];
        xbounds_mapG(xbounds_mapG > size(basemapG{gg},1) | xbounds_mapG < 1) = [];
        ybounds_mapG(ybounds_mapG > size(basemapG{gg},2) | ybounds_mapG < 1) = [];
        % Find the pixels surrounding peak of this field, starting from bottom left to top right, moving left to right
        xbounds_plot = ybounds_mapG;
        ybounds_plot = size(basemapG{gg},1)-xbounds_mapG+1;
        
        count = count + 1;
        gridnum(count) = gg;
        linbin = dummygrid{gg}(sort(xbounds_mapG),ybounds_mapG);
%         linbin(inds(sort(xbounds_mapG), ybounds_mapG)<1) = nan;
        linearobjbin(count) = {linbin}; % in plotting coords (x is left to right, y is bottom to top)
        fieldbounds_plot(count,1:2) = {xbounds_plot ybounds_plot};
        fieldbounds_mat(count,1:2) = {xbounds_mapG ybounds_mapG};

    end
end
            
for ii = 1:count
%     % Plot base map and field of interest
%     fignum = str2double(['10' num2str(ii)]);
%     h = figure(fignum);
%     ax = gca;
%     h.Name = [ID ' BaseMap: Pixel ' num2str(ii)];
%     switch objtype 
%         case 'place'
%             plotmap(basemap,'place');
%             % Patch base pixels
%             patch([fieldbounds_plot{ii,1}(1)-1 fieldbounds_plot{ii,1}(1)-1 fieldbounds_plot{ii,1}(end) fieldbounds_plot{ii,1}(end)],[fieldbounds_plot{ii,2}(1)-1 fieldbounds_plot{ii,2}(end) fieldbounds_plot{ii,2}(end) fieldbounds_plot{ii,2}(1)-1], [0 0 0 0] ,'EdgeColor','k','FaceColor','none');
%         case 'spatialview'
%             plotmap(basemap,'spatialview');
%             for pp = 1:size(fieldbounds_plot{ii,1},2) % For each pixel
%                 for mm = 1:size(fieldbounds_plot{ii,2},2)
%                     if ~isnan(linearobjbin{ii}(size(linearobjbin{ii},1)-mm+1,pp))
%                         [x,y,z] = converttosurf(gridnum(ii),fieldbounds_plot{ii,1}(pp),fieldbounds_plot{ii,2}(mm));
%                         patch(x,y,z,[1 1 1 1],'FaceColor','none');
%                     end
%                 end
%             end
%     end
%     set(ax,'CLim',[0 max(basemap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%         'XColor','none','YColor','none','ZColor','none',...
%         'FontSize',14,'GridLineStyle','none','Color','none');
%     ax.Title.String = 'BaseMap';
%     % Patch environment boundaries
%     patchenvbounds(objtype);

    % Initialise variables
    rate_components = cell(size(fieldbounds_plot{ii,2},2),size(fieldbounds_plot{ii,1},2));
    maxrate_px = nan(size(fieldbounds_plot{ii,2},2),size(fieldbounds_plot{ii,1},2));
    for yy = 1:size(fieldbounds_plot{ii,2},2) % Going left to right, bottom to top
        for xx = 1:size(fieldbounds_plot{ii,1},2)
            % Get base pixel
            if isnan(fieldbounds_plot{ii,1}(xx)) || isnan(fieldbounds_plot{ii,2}(yy)) || isnan(linearobjbin{ii}(size(linearobjbin{ii},1)-yy+1,xx))
                continue;
            end

            % Get corresponding secondary pixels from pv object
            sec = [];
            switch objtype
                case 'place'
                    ind_pv = pvT(:,2) == linearobjbin{ii}(size(linearobjbin{ii},1)-yy+1,xx);
                    sec(:,1) = pvT(ind_pv,3); % view px
                    sec(:,2) = pvT(ind_pv,4); % dur
                    sec(:,3) = pvT(ind_pv,5); % spikes
                case 'spatialview'
                    ind_pv = pvTfill(:,3) == linearobjbin{ii}(size(linearobjbin{ii},1)-yy+1,xx);
                    sec(:,1) = pvTfill(ind_pv,2); % place px
                    sec(:,2) = pvTfill(ind_pv,4); % dur
                    sec(:,3) = pvTfill(ind_pv,5); % spikes
            end

            % Get firing rates 
            usecpx = unique(sec(:,1));
            usecpx(isnan(usecpx)) = [];
            rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
            rate_components_px(:,1) = usecpx;
            ind_timechange = find(sec(:,2) > 0);
            for pp = 1:size(ind_timechange,1) % For each instance of being in this base pixel
                linbin = nan(size(usecpx,1),2); % Temp dur and spikes for this secondary pixel(s)
                if pp == size(ind_timechange,1)
                    newind = ind_timechange(pp):size(sec,1);
                else
                    newind = ind_timechange(pp):ind_timechange(pp+1)-1; % index into secondary pixel(s) for this instance
                end
                newview = sec(newind,1); 
                newset = ismember(usecpx,newview);
                linbin(newset,1) = sec(newind(1),2); % duration for this instance listed with first secondary pixel
                linbin(newset,2) = sec(newind(end),3); % spikes for this instance listed with last secondary pixel
                rate_components_px(:,2) = nansum( [rate_components_px(:,2) linbin(:,1)] ,2); % Sum duration for this sec pixel across instances
                rate_components_px(:,3) = nansum( [rate_components_px(:,3) linbin(:,2)] ,2); % Sum spikes for this sec pixel across instances
            end
            rightfulnans = rate_components_px(:,2) == 0;
            rate_components_px(rightfulnans,:) = nan;
            rate_components_px(:,4) = rate_components_px(:,3)./rate_components_px(:,2);

            % Collect all sec rate information for the 9 base pixels
            rate_components(3-yy+1,xx) = { rate_components_px };
    %                 if max(rate_components_px(:,4)) > maxrate_px 
            if ~isempty(rate_components_px)
                maxrate_px(3-yy+1,xx) = max(rate_components_px(:,4));
            else 
                maxrate_px(3-yy+1,xx) = NaN;
            end
    %                 end

        end
    end
    disp(ii);

    % Plot secondary maps by base pixel
%     fignum = str2double(['10' num2str(ii) '1']);
    fignum = ii;
    h = figure(fignum);
    h.Name = [ID ' SecondaryMaps: Pixel ' num2str(ii)];
    for yy = 1:size(fieldbounds_plot{ii,2},2) % Going left to right, bottom to top
        for xx = 1:size(fieldbounds_plot{ii,1},2)
            if isempty(rate_components{3-yy+1,xx}) || isnan(linearobjbin{ii}(size(linearobjbin{ii},1)-yy+1,xx))
                continue;
            end
            % Plot secondary map for each base pixel in field
            tempsecmap = nan(size(secmap));
            tempsecmap(1,rate_components{3-yy+1,xx}(:,1)) = rate_components{3-yy+1,xx}(:,4);
            ax = subplot('Position',[ 0.05+(xx-1)*(0.9/3) 0.05+(yy-1)*(0.9/3) 0.8/3 0.8/3 ]);
            switch objtype
                case 'place'
                    plotmap(tempsecmap,'spatialview');
                    % Patch environment boundaries
                    patchenvbounds('spatialview');
                    % Patch basemap pixel
                    patch([fieldbounds_plot{ii,1}(xx)-1 fieldbounds_plot{ii,1}(xx)-1 fieldbounds_plot{ii,1}(xx) fieldbounds_plot{ii,1}(xx)],[fieldbounds_plot{ii,2}(yy)-1 fieldbounds_plot{ii,2}(yy) fieldbounds_plot{ii,2}(yy) fieldbounds_plot{ii,2}(yy)-1], [0 0 0 0] ,'red','EdgeColor','k','FaceColor','r');
                    patch([fieldbounds_plot{ii,1}(xx)-1 fieldbounds_plot{ii,1}(xx)-1 fieldbounds_plot{ii,1}(xx) fieldbounds_plot{ii,1}(xx)],[fieldbounds_plot{ii,2}(yy)-1 fieldbounds_plot{ii,2}(yy) fieldbounds_plot{ii,2}(yy) fieldbounds_plot{ii,2}(yy)-1], [16 16 16 16] ,'red','EdgeColor','k','FaceColor','r');
                case 'spatialview'
                    plotmap(tempsecmap,'place');
                    % Patch environment boundaries
                    patchenvbounds('spatialview');
                    % Patch basemap pixel
                    switch objtype 
                        case 'place'
                            % Patch base pixels
                            patch([fieldbounds_plot{ii,1}(1)-1 fieldbounds_plot{ii,1}(1)-1 fieldbounds_plot{ii,1}(end) fieldbounds_plot{ii,1}(end)],[fieldbounds_plot{ii,2}(1)-1 fieldbounds_plot{ii,2}(end) fieldbounds_plot{ii,2}(end) fieldbounds_plot{ii,2}(1)-1], [0 0 0 0] ,'EdgeColor','k','FaceColor','none');
                        case 'spatialview'
                            % Patch base pixels
                            [x,y,z] = converttosurf(gridnum(ii),fieldbounds_plot{ii,1}(xx),fieldbounds_plot{ii,2}(yy));
                            patch(x,y,z,'r');
                    end    
            end
            ax.Title.String = horzcat('x',num2str(xx),'y',num2str(yy));
            if maxrate_px(3-yy+1,xx) == 0
                CLim = 0.001;
            else
                CLim = maxrate_px(3-yy+1,xx);
            end
            set(ax,'CLim',[0 CLim],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                'XColor','none','YColor','none','ZColor','none',...
                'FontSize',14,'GridLineStyle','none','Color','none');
            view(0,90);
            disp(ii);
        end
    end
    
end




% Find grid in spatial view frame for 1 pixel
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

% Plot rate map
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
colorbar('Location','east','FontSize',36);

% Patch environment boundaries
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
saveas(h,figtitle,'fig');
% print('-painters',figtitle,'-dsvg');
cd(cwd);
