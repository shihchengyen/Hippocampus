function glm_corrmap(filttype,pix,varargin)

%
%
%   Be in directory of cell, 
%   or input cell dir as varargin
%    
%

%%%% CHANGE THIS TO SOMETHING ELSE
% figdir = '/Volumes/User/huimin/Desktop/';
figdir = '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Figures/';

if nargin > 2
    cd(varargin{1});
end
cwd = pwd;
objdir = [cwd '/' filttype '/' num2str(pix) 'px/'];


%%% not part of placebyspatialview.m
cd ../../..
pv = load([num2str(pix) 'vmpv.mat']);
pv = pv.pv;
%%% over
cd(objdir);

% Load vmpc object
pc = load('vmpc.mat');
pc = pc.vmp.data;
% % Load vmsv object
% sv = load('vmsv.mat');
% sv = sv.vms.data;
cd(cwd);

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
% pvTfill = pvT;

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

    if ~isempty(subsample)
        
        % Get spikes and duration for place only
        place_durations(1,i) = sum(subsample(:,2));
        place_spikes(1,i) = sum(subsample(:,3));
        
        % back-filling spikes for view
        subsample(subsample(:,3)==0,3) = nan;
%         subsample(:,4) = circshift(subsample(:,2)~=0 ,-1);
        subsample(:,4) = circshift(subsample(:,2)~=0 ,0);
        subsample(isnan(subsample(:,3)) & subsample(:,4), 3) = 0;
        subsample(:,4) = [];
        subsample(:,3) = fillmissing(subsample(:,3), 'next');
        % back-filling time for view
        subsample(subsample(:,2)==0,2) = nan;
%         subsample(:,2) = fillmissing(subsample(:,2), 'previous');
        subsample(:,2) = fillmissing(subsample(:,2), 'next');
        % Put backfill into sessionTimeC array
%         pvTfill(inds,[3 4 5]) = subsample;

        % padding with 5122 bin
        if subsample(end,1) ~= 5122
            subsample = [subsample; [5122 NaN NaN]];
        end
%         subsample = [subsample; [5122 0 0]];
        % remove bad view spots
        subsample(isnan(subsample(:,1)),:) = [];
        % sum durations
        view_durations(:,i) = accumarray(subsample(:,1), subsample(:,2),[],[],NaN);
        % sum spikes
        view_spikes(:,i) = accumarray(subsample(:,1), subsample(:,3),[],[],NaN); 

    end
end

p_array_orig = place_spikes./place_durations;
sv_array_orig = nansum(view_spikes,2)./nansum(view_durations,2);
view_spikes_temp = view_spikes;
view_spikes_temp(isnan(view_spikes)) = 0; % to prevent error on factorial

% Maximum-likelihood maps
startorig = 1;
if startorig
    % Start with Original maps
    p_array = p_array_orig;
    sv_array = sv_array_orig;
else
    % Start with uniform array
    p_array = ones(1,size(place_durations,2));
    p_array(isnan(p_array_orig)) = nan;
    sv_array = ones(size(view_durations,1),1);
    sv_array(isnan(sv_array_orig)) = nan;
end
% Starting llh
llh = sum( nansum(view_spikes.*log(p_array.*view_durations.*sv_array)) - nansum(p_array.*view_durations.*sv_array) - nansum(log(factorial(view_spikes_temp))) );
llh_vec = llh;

p_array_set = [p_array];
sv_array_set = [sv_array];
p_spk_set = {view_spikes};
sv_spk_set = {view_spikes};
mlm = true;
count = 0;

while mlm == true

%     spk_adj = nan(size(view_spikes,1),size(view_spikes,2));
%     for ii = 1:size(view_spikes,2)
%         for jj = 1:size(view_spikes,1)
%             if ~isnan(view_spikes(jj,ii))
%                 if isnan(sv_array(jj))
%                     spk_adj(jj,ii) = nan;
%                 elseif sv_array(jj) == 0
%                     spk_adj(jj,ii) = 0; 
%                 else
%                     spk_adj(jj,ii) = view_spikes(jj,ii)/sv_array(jj);
%                 end
%             end
%         end
%     end
%     multfac = nansum(nansum(view_spikes)) / nansum(nansum(spk_adj));
%     spk_adj = multfac*spk_adj;
%     p_array = nansum(spk_adj,1)./nansum(view_durations,1);
%     p_array_set = [p_array_set; p_array];
%     p_spk_set{end+1} = spk_adj;
    
    dur_adj = nan(size(view_durations,1),size(view_durations,2));
    for ii = 1:size(view_durations,2)
        for jj = 1:size(view_durations,1)
            if ~isnan(view_durations(jj,ii))
                if isnan(sv_array(jj))
                    dur_adj(jj,ii) = 0;
                elseif sv_array(jj) == 0
                    dur_adj(jj,ii) = view_durations(jj,ii); 
                else
                    dur_adj(jj,ii) = view_durations(jj,ii)*sv_array(jj);
                end
            end
        end
    end
    p_array = nansum(view_spikes,1)./nansum(dur_adj,1);
    spk = p_array.*nansum(view_durations,1);
    totspk = nansum(spk);
    multfac = nansum(nansum(view_spikes)) / totspk;
    p_array = multfac * p_array;
    p_array_set = [p_array_set; p_array];
    
%     spk_adj = nan(size(view_spikes,1),size(view_spikes,2));
%     for jj = 1:size(view_spikes,1) 
%         for ii = 1:size(view_spikes,2)
%             if ~isnan(view_spikes(jj,ii))
%                 if isnan(p_array(ii))
%                     spk_adj(jj,ii) = nan; 
%                 elseif p_array(ii) == 0
%                     spk_adj(jj,ii) = 0; 
%                 else
%                     spk_adj(jj,ii) = view_spikes(jj,ii)/p_array(ii);
%                 end
%             end
%         end
%     end
%     multfac = nansum(nansum(view_spikes)) / nansum(nansum(spk_adj));
%     spk_adj = multfac*spk_adj;
%     sv_array = nansum(spk_adj,2)./nansum(view_durations,2);
%     sv_array_set = [sv_array_set sv_array];
%     sv_spk_set{end+1} = spk_adj;
    
    dur_adj = nan(size(view_durations,1),size(view_durations,2));
    for jj = 1:size(view_durations,1)
        for ii = 1:size(view_durations,2)
            if ~isnan(view_durations(jj,ii))
                if isnan(p_array(ii))
                    dur_adj(jj,ii) = 0;
                elseif p_array(ii) == 0
                    dur_adj(jj,ii) = view_durations(jj,ii); 
                else
                    dur_adj(jj,ii) = view_durations(jj,ii)*p_array(ii);
                end
            end
        end
    end
    sv_array = nansum(view_spikes,2)./nansum(dur_adj,2);
    spk = sv_array.*nansum(view_durations,2);
    totspk = nansum(spk);
    multfac = nansum(nansum(view_spikes)) / totspk;
    sv_array = multfac * sv_array;
    sv_array_set = [sv_array_set sv_array];

%     prev_llh = llh;
    llh = sum( nansum(view_spikes.*log(p_array.*view_durations.*sv_array)) - nansum(p_array.*view_durations.*sv_array) - nansum(log(factorial(view_spikes_temp))) );
    llh_vec(end+1,1) = llh;
    disp([num2str(size(llh_vec,1)) ':' num2str(llh)]);
    if startorig
        pick = find(llh_vec == max(llh_vec));
    else 
        pick = find(llh_vec == max(llh_vec(2:end)));
    end
    if size(pick,1) > 1
        disp('undifferentiable llh');
        pick = pick(1);
        picklabel = ['undiff' num2str(pick)];
        break;
    else 
        picklabel = num2str(pick);
    end

    if pick < size(llh_vec,1) && size(llh_vec,1) - pick > 1 
        disp('stop');
        count = count + 1;
    end
    if count > 10
        break;
    end
    if isinf(llh)
        mlm = false;
        disp('stop');
    end
end

corr_p_array = p_array_set(pick,:);
corr_sv_array = sv_array_set(:,pick);
% corr_spk = sv_spk_set{pick};
% corr_sv_spk = sv_spk_set{pick};
disp(pick);

%%% temporary %%%

% Collect date, session, array, channel, cell
s = regexp(cwd,'session');
identifiers = [str2double(cwd(s-9:s-2)) str2double(cwd(s+7:s+8)) ...
    str2double(cwd(s+15:s+16)) str2double(cwd(s+25:s+27)) str2double(cwd(s+33:s+34))];
ID = [num2str(identifiers(1,1)) 'ch' num2str(identifiers(1,4)) 'c' num2str(identifiers(1,5))];
figdir = [figdir filttype '/' num2str(pix) 'px/Corrmaps/' ID 'pick' picklabel '/'];
mkdir(figdir);
cd(figdir)
% mkdir(ID)
% cd('corrmaps')
%%%%%%%%%%%%%%%

save(['llh_history' '.mat'], 'llh_vec')

% Original place map
h = figure(1);
set(h,'Units','normalized','Position',[0 0 1 1]);
ax = subplot(1,2,1);
% h.Name = [ID 'vmpvRawPlaceMap'];
map = p_array_orig;
plotmap(map,'place');
maxC = nanmax(map);
ax.Title.String = ['vmpvRawPlaceMap'];
set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
clear map;

% Corrected place map 
ax = subplot(1,2,2);
map = corr_p_array;
plotmap(map,'place');
ax.Title.String = ['CorrPlaceMap' ' Pick ' picklabel ', MaxRate ' num2str(nanmax(map)) 'Hz'];
set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
clear map;
saveas(h,['placemap.fig']);
saveas(h,['placemap.png']);
close(h);

% Original view map
h = figure(2);
ax = subplot(1,2,1);
map = sv_array_orig;
plotmap(map,'spatialview');
maxC = nanmax(map);
ax.Title.String = ['vmpvRawViewMap'];
set(ax,'CLim',[0 nanmax(map)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
clear map;

% Corrected view map 
ax = subplot(1,2,2);
set(h,'Units','normalized','Position',[0 0 1 1]);
map = corr_sv_array;
plotmap(map,'spatialview');
ax.Title.String = ['CorrViewMap' ' Pick ' picklabel ', MaxRate ' num2str(nanmax(map)) 'Hz'];
set(ax,'CLim',[0 maxC],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
    'XColor','none','YColor','none','ZColor','none',...
    'FontSize',14,'GridLineStyle','none','Color','none');
clear map;
saveas(h,['viewmap.fig']);
saveas(h,['viewmap.png']);
close(h);

% LLH curve - short
h = figure(3);
if pick > 5
    start = pick - 5;
else 
    start = 1;
end
line(start:size(llh_vec,1), llh_vec(start:end));
saveas(h,['llh_short.fig']);
saveas(h,['llh_short.png']);
close(h);

% LLH curve - full
h = figure(4);
line(1:size(llh_vec,1), llh_vec(1:end));
saveas(h,['llh_full.fig']);
saveas(h,['llh_full.png']);
close(h);



% p_array(place_durations == 0) = NaN;
% sv_array(nansum(view_durations,2) == 0) = NaN;
% total_place_spikes = nansum(place_spikes);
% total_view_spikes = nansum(nansum(view_spikes));
% p_array(isnan(place_durations)) = NaN;
% sv_array(nansum(view_durations,2) == 0) = NaN;
% total_place_spikes = nansum(place_spikes);
% total_view_spikes = nansum(nansum(view_spikes));

% h = figure(11);
% ax = gca;
% % h.Name = [ID 'vmpvRawPlaceMap'];
% set(h,'Units','normalized','Position',[0 0 1 1]);
% % emptyplacegrids = all(isnan(full_rate),1);
% % rawplacemap1 = nansum(full_rate,1);
% % rawplacemap1(1,emptyplacegrids) = NaN;
% rawplacemap1 = place_spikes./place_durations;
% % rawplacemap1(place_durations == 0) = NaN;
% plotmap(rawplacemap1,'place');
% ax.Title.String = ['vmpvRawPlaceMap'];
% set(ax,'CLim',[0 max(rawplacemap1)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%     'XColor','none','YColor','none','ZColor','none',...
%     'FontSize',14,'GridLineStyle','none','Color','none');
% % saveas(h,['vmpvrawplacemap.fig']);
% % saveas(h,['vmpvrawplacemap.png']);
% % close(h)
% 
% h = figure(14);
% mlmplacefitmap = nansum(p_array.*view_durations.*sv_array,1) ./ nansum(view_durations,1);
% mlmplace_ssres = nansum((rawplacemap1 - mlmplacefitmap).^2);
% mlmplacerawavg = nanmean(rawplacemap1);
% mlmplace_sstot = nansum((rawplacemap1 - mlmplacerawavg).^2);
% mlmplace_goodness = 1 - (mlmplace_ssres / mlmplace_sstot);
% ax = gca;
% %h.Name = [ID 'MLMPlaceMap'];
% set(h,'Units','normalized','Position',[0 0 1 1]);
% [mlmplacefitmapG,~] = plotmap(mlmplacefitmap,'place');
% % [smoothplacemapG,placemapGdummy] = plotplacemap(smoothplacemap); % Do not use grid output as input for another plot as it will turn out rotated 90deg CCW
% ax.Title.String = ['MLMPlaceFitMap, Goodness of Fit R-squared: ' num2str(mlmplace_goodness,10)];
% set(ax,'CLim',[0 max(rawplacemap1)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%     'XColor','none','YColor','none','ZColor','none',...
%     'FontSize',14,'GridLineStyle','none','Color','none');
% saveas(h,['mlmfitplacemap' num2str(title_suffix) '.fig']);
% saveas(h,['mlmfitplacemap' num2str(title_suffix) '.png']);
% close(h)
% 
% h = figure(15);
% ax = gca;
% %h.Name = [ID 'MLMPlaceMap'];
% set(h,'Units','normalized','Position',[0 0 1 1]);
% mlmplacemap = nansum(p_array.*place_durations); % total spikes before scaling
% mlmplacemap = p_array * (total_place_spikes/mlmplacemap);
% [mlmplacemapG,~] = plotmap(mlmplacemap,'place');
% mlmplacemap_exceed = mlmplacemap(mlmplacemap > max(rawplacemap1));
% % [smoothplacemapG,placemapGdummy] = plotplacemap(smoothplacemap); % Do not use grid output as input for another plot as it will turn out rotated 90deg CCW
% ax.Title.String = {'MLMPlaceMap', ['values exceeding CLim: ' mat2str(mlmplacemap_exceed)]};
% set(ax,'CLim',[0 max(rawplacemap1)],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%     'XColor','none','YColor','none','ZColor','none',...
%     'FontSize',14,'GridLineStyle','none','Color','none');
% % patchenvbounds('place');
% % if savefig
% %     savefigure(h,h.Name,figdir);
% % end
% saveas(h,['mlmplacemap' num2str(title_suffix) '.fig']);
% saveas(h,['mlmplacemap' num2str(title_suffix) '.png']);
% close(h)
% 
% h = figure(21); 
% ax = gca;
% % h.Name = [ID 'vmpvRawViewMap'];
% set(h,'Units','normalized','Position',[0 0 1 1]);
% spikes = nansum(view_spikes,2);
% durations = nansum(view_durations,2);
% rawviewmap1 = spikes./durations;
% rawviewmap1(durations==0) = nan;
% % rawviewmap1 = emptyinsidepillar(rawviewmap1); % Temporary measure only! Remove data from inside of pillar where it should be empty
% plotmap(rawviewmap1,'spatialview');
% ax.Title.String = ['vmpvRawViewMap'];
% set(ax,'CLim',[0 max(rawviewmap1(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%     'XColor','none','YColor','none','ZColor','none',...
%     'FontSize',14,'GridLineStyle','none','Color','none');
% saveas(h,['vmpvrawviewmap.fig']);
% saveas(h,['vmpvrawviewmap.png']);
% close(h)
% 
% h = figure(24); 
% mlmviewfitmap = nansum(p_array.*view_durations.*sv_array,2) ./ nansum(view_durations,2);
% mlmview_ssres = nansum((rawviewmap1 - mlmviewfitmap).^2);
% mlmviewrawavg = nanmean(rawviewmap1);
% mlmview_sstot = nansum((rawviewmap1 - mlmviewrawavg).^2);
% mlmview_goodness = 1 - (mlmview_ssres / mlmview_sstot);
% ax = gca;
% %h.Name = [ID 'MLMViewMap'];
% set(h,'Units','normalized','Position',[0 0 1 1]);
% plotmap(mlmviewfitmap,'spatialview');
% ax.Title.String = ['MLMViewFitMap, Goodness of Fit R-squared: ' num2str(mlmview_goodness,10)];
% set(ax,'CLim',[0 max(rawviewmap1(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%     'XColor','none','YColor','none','ZColor','none',...
%     'FontSize',14,'GridLineStyle','none','Color','none');
% % patchenvbounds('spatialview');
% % if savefig
% %     savefigure(h,h.Name,figdir);
% % end
% saveas(h,['mlmfitviewmap' num2str(title_suffix) '.fig']);
% saveas(h,['mlmfitviewmap' num2str(title_suffix) '.png']);
% close(h)
% 
% h = figure(25); 
% ax = gca;
% %h.Name = [ID 'MLMViewMap'];
% set(h,'Units','normalized','Position',[0 0 1 1]);
% %mlmviewmap = nansum(p_array.*view_durations.*sv_array,2) ./ nansum(view_durations,2);
% mlmviewmap = nansum(sv_array.*nansum(view_durations,2)); % total spikes before scaling
% mlmviewmap = sv_array * (total_view_spikes/mlmviewmap);
% plotmap(mlmviewmap,'spatialview');
% mlmviewmap_exceed = mlmviewmap(mlmviewmap > max(rawviewmap1));
% ax.Title.String = {'MLMViewMap', ['values exceeding CLim: ' mat2str(mlmviewmap_exceed)]};
% set(ax,'CLim',[0 max(rawviewmap1(3:end))],'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1],...
%     'XColor','none','YColor','none','ZColor','none',...
%     'FontSize',14,'GridLineStyle','none','Color','none');
% % patchenvbounds('spatialview');
% % if savefig
% %     savefigure(h,h.Name,figdir);
% % end
% saveas(h,['mlmviewmap' num2str(title_suffix) '.fig']);
% saveas(h,['mlmviewmap' num2str(title_suffix) '.png']);
% close(h)

end

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

%disp(sum(sum(find(ceiling==0))) + sum(sum(find(floor==0))) + sum(sum(find(P4_TL==0))));

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
colorbar('Location','EastOutside','FontSize',36);

end