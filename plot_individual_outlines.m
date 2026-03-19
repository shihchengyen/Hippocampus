function plot_individual_outlines()
% PLOT_INDIVIDUAL_OUTLINES
% Creates TWO figures per cell:
%   1) Filled red place-field outlines + peak stars (plot.m style)
%   2) MapOnly full-color rate map (no outlines)
%
% Output folders:
%   individual_cell_outlines/
%   individual_cell_maps/

%% ================== CHECK DATA ==================
if ~evalin('base', 'exist(''data'', ''var'')')
    error('No ''data'' variable found in workspace. Run ProcessDirs first.');
end
data = evalin('base', 'data');

if isempty(data) || size(data,2) < 3
    error('Data variable has incorrect structure.');
end

%% ================== LOAD CELL LIST ==================
if ~exist('place_selective_cells.txt','file')
    error('Cannot find place_selective_cells.txt');
end

fid = fopen('place_selective_cells.txt','r');
cell_dirs = {};
while ~feof(fid)
    t = fgetl(fid);
    if ischar(t) && ~isempty(t)
        cell_dirs{end+1} = t;
    end
end
fclose(fid);

fprintf('Found %d cell directories\n', numel(cell_dirs));

%% ================== SAMPLE MAP ==================
vmp_file = fullfile(cell_dirs{1},'vmpc.mat');
tmp = load(vmp_file);
if isfield(tmp,'vmp')
    vmp = tmp.vmp.data.maps_adsm(1,:);
else
    vmp = tmp.vpc.data.maps_adsm(1,:);
end

grid_size = 40;
map_nan = reshape(vmp, grid_size, grid_size);

%% ================== EXTRACT CELL DATA ==================
cell_outlines = {};
valid_cells = 0;

for i = 1:size(data,1)
    rc = data{i,1};    % region_cells
    pi = data{i,2};    % peak_indices
    pv = data{i,3};    % peak_values
    if ~isempty(rc) && iscell(rc)
        valid_cells = valid_cells + 1;
        sic = NaN;
        if size(data,2) >= 4
            sic = data{i,4};
        end
        cell_outlines{valid_cells} = {rc, pi, pv, sic};
    end
end

fprintf('Valid cells with regions: %d\n', valid_cells);

%% ================== OUTPUT DIRS ==================
outline_dir = 'individual_cell_outlines';
map_dir     = 'individual_cell_maps';
if ~exist(outline_dir,'dir'), mkdir(outline_dir); end
if ~exist(map_dir,'dir'), mkdir(map_dir); end

%% ================== LOOP CELLS ==================
for c = 1:valid_cells

    cell_path = cell_dirs{c};
    parts = strsplit(cell_path, filesep);

    date_str=''; session_str=''; channel_str=''; cell_str='';
    for i = 1:numel(parts)
        if ~isempty(regexp(parts{i},'^\d{8}$','once')), date_str=parts{i}; end
        if ~isempty(regexp(parts{i},'^session\d+$','once')), session_str=parts{i}; end
        if ~isempty(regexp(parts{i},'^channel\d+$','once')), channel_str=parts{i}; end
        if ~isempty(regexp(parts{i},'^cell\d+$','once')), cell_str=parts{i}; end
    end

    title_str = sprintf('%s | %s | %s | %s', ...
        date_str, session_str, channel_str, cell_str);

    %% ================== LOAD MAP ==================
    tmp = load(fullfile(cell_path,'vmpc.mat'));
    if isfield(tmp,'vmp')
        map = tmp.vmp.data.maps_adsm(1,:);
    else
        map = tmp.vpc.data.maps_adsm(1,:);
    end
    map2d = reshape(map, grid_size, grid_size);

    %% ================== MAP ONLY FIGURE ==================
    fig1 = figure('Visible','off','Position',[100 100 700 700]);
    imagesc(map2d);
    axis equal tight;
    set(gca,'YDir','reverse');
    colormap(jet);
    colorbar;
    title([title_str ' | MapOnly'],'FontWeight','bold');
    saveas(fig1, fullfile(map_dir, ...
        sprintf('cell_map_%s_%s_%s_%s.png', ...
        date_str, session_str, channel_str, cell_str)));
    close(fig1);

    %% ================== OUTLINE FIGURE (CORRECT) ==================
    fig2 = figure('Visible','off','Position',[100 100 700 700]);

    region_cells = cell_outlines{c}{1};
    peak_idx     = cell_outlines{c}{2};
    sic          = cell_outlines{c}{4};

    % ---- Build RGB image like plot(vpc,'outline',1) ----
    rgb = ones(grid_size, grid_size, 3);   % white background

    for r = 1:grid_size
        for cc = 1:grid_size
            if ~isnan(map2d(r,cc))
                rgb(r,cc,:) = 0.7;        % gray maze
            end
        end
    end

    % ---- Fill regions in RED ----
    for reg = 1:numel(region_cells)
        inds = region_cells{reg};
        [rr,cc] = ind2sub([grid_size grid_size], inds);
        for k = 1:numel(rr)
            rgb(rr(k),cc(k),:) = [1 0 0];
        end
    end

    image(rgb);
    axis equal tight;
    set(gca,'YDir','reverse');
    hold on;

    % ---- Peak stars (WHITE) ----
    for p = 1:numel(peak_idx)
        [pr,pc] = ind2sub([grid_size grid_size], peak_idx(p));
        plot(pc, pr, 'w*','MarkerSize',14,'LineWidth',2);
    end

    % ---- Title ----
    pix = sum(cellfun(@numel, region_cells));
    if ~isnan(sic)
        title(sprintf('%s\n%d regions | %d pixels | SIC %.3f', ...
            title_str, numel(region_cells), pix, sic), ...
            'FontWeight','bold');
    else
        title(sprintf('%s\n%d regions | %d pixels', ...
            title_str, numel(region_cells), pix), ...
            'FontWeight','bold');
    end

    hold off;
    saveas(fig2, fullfile(outline_dir, ...
        sprintf('cell_outline_%s_%s_%s_%s.png', ...
        date_str, session_str, channel_str, cell_str)));
    close(fig2);

    fprintf('Saved cell %d\n', c);
end

fprintf('\nDONE:\n- Filled outlines + stars (correct style)\n- MapOnly figures\n');
end
