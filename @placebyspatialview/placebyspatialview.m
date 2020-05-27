function [obj, varargout] = placebyspatialview(objtype)

% Load vmpc object
pc = vmpc('auto');
pc = pc.data;
% 
% Load vmsv object
% sv = vmsv('auto');
sv = vmsv('auto','MinOccDur',0,'MinOcc',0);
sv = sv.data;

% Load spike train
spiketrain = load('spiketrain.mat');
spiketrain = spiketrain.timestamps ./ 1000; % in seconds

% Combine place and view info with spikes and make rate maps
pv = vmpv('auto');
pv = pv.data.sessionTimeC;
pv(:,4) = [0; diff(pv(:,1))];
binned = histcounts(spiketrain, pv(:,1))';
pv(:,5) = [binned; 0];

pv(pv(:,2)==-1,:) = [];
pv(pv(:,2)==0,:) = [];

full_durations = NaN(5122,1600);
full_spikes = NaN(5122,1600);

for i = 1:1600

    disp(i);
    subsample = [pv(pv(:,2)==i, [3 4 5])];

    % filling spikes
    subsample(subsample(:,3)==0,3) = nan;
    subsample(:,4) = circshift(subsample(:,2)~=0 ,-1);
    subsample(isnan(subsample(:,3)) & subsample(:,4), 3) = 0;
    subsample(:,4) = [];
    subsample(:,3) = fillmissing(subsample(:,3), 'next');

    % filling time
    subsample(subsample(:,2)==0,2) = nan;
    subsample(:,2) = fillmissing(subsample(:,2), 'previous');

    % padding with 5122 bin
    subsample = [subsample; [5122 0 0]];

    % remove bad view spots
    subsample(isnan(subsample(:,1)),:) = [];

    % sum durations
    full_durations(:,i) = accumarray(subsample(:,1), subsample(:,2),[],[],NaN);

    % sum spikes
    full_spikes(:,i) = accumarray(subsample(:,1), subsample(:,3),[],[],NaN);

end
full_rate = full_spikes ./ full_durations;

% Plot full place maps
h = figure(11);
ax = gca;
h.Name = 'vmpvRawPlaceMap';
emptyplacegrids = all(isnan(full_rate),1);
rawplacemap1 = nansum(full_rate,1);
rawplacemap1(1,emptyplacegrids) = NaN;
plotplacemap(rawplacemap1);
ax.Title.String = 'vmpvRawPlaceMap';
set(ax,'CLim',[0 max(rawplacemap1)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
patchenvbounds('place');


h = figure(12);
ax = gca;
h.Name = 'vmpcRawPlaceMap';
rawplacemap2 = pc.maps_raw;
plotplacemap(rawplacemap2);
ax.Title.String = 'vmpcRawPlaceMap';
set(ax,'CLim',[0 max(rawplacemap2)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
patchenvbounds('place');


h = figure(13);
ax = gca;
h.Name = 'vmpcSmoothPlaceMap';
smoothplacemap = pc.maps_adsmooth;
smoothplacemapG = plotplacemap(smoothplacemap); % Do not use grid output as input for another plot as it will turn out rotated 90deg CCW
ax.Title.String = 'vmpcSmoothPlaceMap';
set(ax,'CLim',[0 max(smoothplacemap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
patchenvbounds('place');


% Plot full view maps
h = figure(21); 
ax = gca;
h.Name = 'vmpvRawViewMap';
emptyviewgrids = all(isnan(full_rate),2);
rawviewmap = nansum(full_rate,2);
rawviewmap(emptyviewgrids,1) = NaN;
plotviewmap(rawviewmap);
ax.Title.String = 'vmpvRawViewMap';
set(ax,'CLim',[0 max(rawviewmap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
patchenvbounds('spatialview');

smoothviewmap = sv.linear_map;
h = figure(23);
ax = gca;
h.Name = 'vmsvSmoothViewMap';
smoothviewmapG = plotviewmap(smoothviewmap); % Do not use grid output as input for another plot as it will turn out rotated 90deg CCW
ax.Title.String = 'vmsvSmoothViewMap';
set(ax,'CLim',[0 max(smoothviewmap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
patchenvbounds('spatialview');

% Find 3 maxima of rate maps
switch objtype
    case 'place'
        basemap = smoothplacemap;
        basemapG = smoothplacemapG;
        SIthr = prctile(pc.SICsh(1,2:end),95);
        SI = pc.SIC;
        pv_base = pv(:,2);
        pv_sec = pv(:,3);
    case 'spatialview'
        basemap = smoothviewmap;
        basemapG = smoothviewmapG;
        SIthr = prctile(sv.SICsh(1,2:end),95);
        SI = sv.SIC;
end

dummygrid = flipud(reshape((1:1600),40,40)'); % TEST VIEW GRIDS! 
if SI > SIthr
    % Find fields with at least 1 pixel of > half peak rate
    peakrate = max(max((basemapG)));
    threshrate = peakrate/2;
    ind_fields = basemapG > threshrate;
    % Find separate fields
    [fieldlabel,fieldcount] = bwlabel(ind_fields,8);
    % Find fields that are big enough (i.e. > 6 bins around peak)
    count = 0;
    for ii = 1:fieldcount
        inds = fieldlabel == ii;
        % Find peak of this field in grid matrix coords
        [xmapG,ymapG] = find( (basemapG == max(basemapG(inds))) ,1); 
        xbounds_mapG = [xmapG+1 xmapG xmapG-1]; % Plotting coords (x is left to right, y is bottom to top)
        ybounds_mapG = [ymapG-1 ymapG ymapG+1];
        xbounds_mapG(xbounds_mapG > 40 | xbounds_mapG < 1) = [];
        ybounds_mapG(ybounds_mapG > 40 | ybounds_mapG < 1) = [];
        % Find the pixels surrounding peak of this field, starting from bottom left to top right, moving left to right
        xbounds_plot = ybounds_mapG;
        ybounds_plot = size(basemapG,1)-xbounds_mapG+1;
        % Make sure at least 6 joined pixels
        if sum(sum(inds(sort(xbounds_mapG), ybounds_mapG))) >= 6
            count = count + 1;
            fieldmaxrate(count) = basemapG(xmapG,ymapG);
            linbin = dummygrid(sort(xbounds_mapG),ybounds_mapG);
            linbin(inds(sort(xbounds_mapG), ybounds_mapG)<1) = nan;
            linearobjbin(count) = {linbin};
            fieldbounds_plot(count,1:2) = {xbounds_plot ybounds_plot};
            fieldbounds_mat(count,1:2) = {xbounds_mapG ybounds_mapG};
        end
    end
    % Sort fields
    [fieldmaxrate,I] = sort(fieldmaxrate,'descend');
    fieldbounds_plot = fieldbounds_plot(I,:);
    fieldbounds_mat = fieldbounds_mat(I,:);
    linearobjbin = linearobjbin(I);
    
    for ii = 1:3 % Plot only for first 3 fields per map
        
        if ii > size(fieldbounds_plot,1)
            break;
        end
            
        % Plot base map and outline field
        fignum = str2double(['10' num2str(ii)]);
        h = figure(fignum);
        ax = gca;
        h.Name = ['BaseMap: Pixel ' num2str(ii)];
        switch objtype
            case 'place'
                plotplacemap(smoothplacemap);
                ax = gca;
                set(ax,'CLim',[0 max(smoothplacemap)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                    'XColor','none','YColor','none','ZColor','none',...
                    'FontSize',14,'GridLineStyle','none','Color','none');
                % Patch environment boundaries
                patchenvbounds('place');
                % Patch field of interest
                patch([fieldbounds_plot{ii,1}(1)-1 fieldbounds_plot{ii,1}(1)-1 fieldbounds_plot{ii,1}(end) fieldbounds_plot{ii,1}(end)],[fieldbounds_plot{ii,2}(1)-1 fieldbounds_plot{ii,2}(end) fieldbounds_plot{ii,2}(end) fieldbounds_plot{ii,2}(1)-1], [0 0 0 0] ,'EdgeColor','k','FaceColor','none');
            case 'spatialview'
                plotviewmap(smoothviewmap);
        end
        ax.Title.String = 'BaseMap';
        
        maxrate_px = 0;
        for yy = 1:size(fieldbounds_plot{ii,2},2) % Going left to right, bottom to top
            for xx = 1:size(fieldbounds_plot{ii,1},2)
                % Get base pixel
                if isnan(fieldbounds_plot{ii,1}(xx)) || isnan(fieldbounds_plot{ii,2}(yy)) || isnan(linearobjbin{ii}(size(linearobjbin{ii},1)-yy+1,xx))
                    continue;
                end
                
                % Get corresponding secondary pixels from pv object
                ind_pv = pv_base == linearobjbin{ii}(size(linearobjbin{ii},1)-yy+1,xx);
                sec = [];
                sec(:,1) = pv(ind_pv,3); % view px
                sec(:,2) = pv(ind_pv,4); % dur
                sec(:,3) = pv(ind_pv,5); % spikes
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
                if max(rate_components_px(:,4)) > maxrate_px 
                    maxrate_px = max(rate_components_px(:,4));
                end
                
            end
        end
        disp(ii);
        
        % Plot secondary maps by base pixel
        fignum = str2double(['10' num2str(ii) '1']);
        h = figure(fignum);
        h.Name = ['SecondaryMaps: Pixel ' num2str(ii)];
        for yy = 1:size(fieldbounds_plot{ii,2},2) % Going left to right, bottom to top
            for xx = 1:size(fieldbounds_plot{ii,1},2)
                if isempty(rate_components{3-yy+1,xx})
                    continue;
                end
                % Plot secondary map for each base pixel in field
                switch objtype
                    case 'place'
                        tempsecmap = nan(size(smoothviewmap));
                        tempsecmap(1,rate_components{3-yy+1,xx}(:,1)) = rate_components{3-yy+1,xx}(:,4);
                        ax = subplot('Position',[ 0.05+(xx-1)*(0.9/3) 0.05+(yy-1)*(0.9/3) 0.8/3 0.8/3 ]);
                        plotviewmap(tempsecmap);
                        ax.Title.String = horzcat('x',num2str(xx),'y',num2str(yy));
                        set(ax,'CLim',[0 max(maxrate_px)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
                            'XColor','none','YColor','none','ZColor','none',...
                            'FontSize',14,'GridLineStyle','none','Color','none');
                        % Patch environment boundaries
                        patchenvbounds('spatialview');
                        % Patch basemap pixel
                        patch([fieldbounds_plot{ii,1}(xx)-1 fieldbounds_plot{ii,1}(xx)-1 fieldbounds_plot{ii,1}(xx) fieldbounds_plot{ii,1}(xx)],[fieldbounds_plot{ii,2}(yy)-1 fieldbounds_plot{ii,2}(yy) fieldbounds_plot{ii,2}(yy) fieldbounds_plot{ii,2}(yy)-1], [0 0 0 0] ,'EdgeColor','k','FaceColor','r');
                        % Patch environment boundaries
                        patchenvbounds('spatialview');
                    case 'spatialview'
                        
                end
                disp(ii);
            end
        end
        
        
        
        
    end
    disp(ii);

end



% Plot place map
function [placemapG] = plotplacemap(placemapL)

% Set up surf frame for plotting
floor_x = repmat(0:40, 41, 1);
floor_y = flipud(repmat([0:40]', 1, 41));
floor_z = zeros(41,41);

placemapG = flipud(reshape(placemapL, 40, 40)');
surf(floor_x, floor_y, floor_z, placemapG);
alpha 1; shading flat;
colormap jet;
view(-35,20);

% Plot view map
function [viewmapG]= plotviewmap(viewmapL)

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

floor = flipud(reshape(viewmapL(3:3+1600-1), 40, 40)');
% ceiling follows floor mapping, top down view
ceiling = flipud(reshape(viewmapL(1603:1603+1600-1), 40, 40)');
% from top down, slit walls at bottom left corner, open outwards.
% start from row closest to ground, rightwards, then climb rows
walls = flipud(reshape(viewmapL(3203:3203+1280-1), 40*4, 8)');
% BL - bottom left, and so on, from top view, same slicing as walls
% pillar width 8, height 5
P1_BR = flipud(reshape(viewmapL(4483:4483+160-1), 8*4, 5)');
P1_BR = [P1_BR; nan(1,size(P1_BR,2))];
P1_BR = [P1_BR nan(size(P1_BR,1),1)];
P2_BL = flipud(reshape(viewmapL(4643:4643+160-1), 8*4, 5)');
P2_BL = [P2_BL; nan(1,size(P2_BL,2))];
P2_BL = [P2_BL nan(size(P2_BL,1),1)];        
P3_TR = flipud(reshape(viewmapL(4803:4803+160-1), 8*4, 5)');
P3_TR = [P3_TR; nan(1,size(P3_TR,2))];
P3_TR = [P3_TR nan(size(P3_TR,1),1)];                
P4_TL = flipud(reshape(viewmapL(4963:4963+160-1), 8*4, 5)');
P4_TL = [P4_TL; nan(1,size(P4_TL,2))];
P4_TL = [P4_TL nan(size(P4_TL,1),1)];

viewmapG = { NaN; NaN; floor; ceiling; walls; P1_BR; P2_BL; P3_TR; P4_TL };

% floor
surf(floor_x, floor_y, floor_z, floor);
alpha 1; shading flat;
colormap jet;
hold on;

% ceiling and walls
surf(ceiling_x, ceiling_y, ceiling_z, ceiling);
alpha 1; shading flat;
surf(walls_x, walls_y, walls_z, walls);      
alpha 1; shading flat;

disp(sum(sum(find(ceiling==0))) + sum(sum(find(floor==0))) + sum(sum(find(P4_TL==0))));

% pillars
surf(P1_x, P1_y, PX_z, P1_BR);
alpha 1; shading flat;
surf(P2_x, P2_y, PX_z, P2_BL);
alpha 1; shading flat;
surf(P3_x, P3_y, PX_z, P3_TR);
alpha 1; shading flat;
surf(P4_x, P4_y, PX_z, P4_TL);
alpha 1; shading flat; 
view(-35,20);

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
        % Pillar 1
        patch([24 24 32 32],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 24 24],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([26.88 26.88 29.12 29.12],[32 32 32 32],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Rabbit Poster on m_wall_25
        % Pillar 2
        patch([24 24 32 32],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 24 24],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([24 24 32 32],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([32 32 32 32],[10.88 10.88 13.12 13.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Cat poster on m_wall_10
        patch([24 24 24 24],[10.88 10.88 13.12 13.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Pig poster on m_wall_29
        % Pillar 3
        patch([8 8 16 16],[24 24 24 24],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 8 8],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[32 32 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[24 24 32 32],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[26.88 26.88 29.12 29.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Croc poster on m_wall_4
        patch([8 8 8 8],[26.88 26.88 29.12 29.12],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Donkey poster on m_wall_15
        % Pillar 4
        patch([8 8 16 16],[8 8 8 8],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 8 8],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([8 8 16 16],[16 16 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([16 16 16 16],[8 8 16 16],[16 21 21 16],[1 1 1 1],'FaceColor','none');
        patch([10.88 10.88 13.12 13.12],[8 8 8 8],[17.8 19.2 19.2 17.8],[1 1 1 1],'FaceColor','none'); % Camel poster on m_wall_20
        % Ceiling
        patch([0 0 40 40],[0 40 40 0],[40 40 40 40],[1 1 1 1],'FaceColor','none');
        
end
