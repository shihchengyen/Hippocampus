function [gridmap] = lineartogrid(linmap,var,gridSize)

% This function converts linear array to 2D grid format
% Grid bins are numbered left to right, bottom to top
% e.g. 1600 x 1 linear place map to a single 40 x 40 grid map
% e.g. 5122 x 1 linear spatial view map to 9 2D grid maps
%  - Grid 1: Cue
%  - Grid 2: Hint
%  - Grid 3: Floor 40 x 40
%  - Grid 4: Ceiling 40 x 40 (top down view) 
%  - Grid 5: Walls 8 x 160, starting from bottom left corner
%  - Grid 6: Pillar 1 bottom right 8 x 32, starting from bottom left corner
%  - Grid 7: Pillar 2 bottom left 8 x 32
%  - Grid 8: Pillar 3 top right 8 x 32
%  - Grid 9: Pillar 4 top left 8 x 32

% Inputs:
% linmap: 1D linear map, must be 1600 x nshuff or 5122 x nshuff
% var: 'place' or 'spatialview'
% gridSize: specifications for dimensions of 2D map

switch var
    case 'place'
        if size(gridSize,1) ~= 1 || size(linmap,1) ~= sum(gridSize(:,1).*gridSize(:,2))
            error('Incorrect dimensions specified for grid size');
        end
        gridmap = cell(1,1);
    case 'view'
        if size(gridSize,1) ~= 9 || size(linmap,1) ~= sum(gridSize(:,1).*gridSize(:,2))
            error('Incorrect dimensions specified for grid size');
        end
        gridmap = cell(size(gridSize,1),1);
    case 'headdirection'
        if size(gridSize,1) ~= 1 || size(linmap,1) ~= sum(gridSize(:,1).*gridSize(:,2))
            error('Incorrect dimensions specified for grid size');
        end
        gridmap = cell(1,1);
end

for ii = 1:size(gridmap,1)
    % Initialize empty grids
    temp = nan(gridSize(ii,1),gridSize(ii,2),size(linmap,2));
    % Assign linear bin to grid bin. Note that this goes left
    % to right, bottom to top
    for mm = 1:gridSize(ii,1)*gridSize(ii,2) % For every point in linear map
        if mod(mm,gridSize(ii,2)) == 0
            y = gridSize(ii,2); % y dim of matrix
        else
            y = mod(mm,gridSize(ii,2));
        end
        x = gridSize(ii,1)-ceil(mm/gridSize(ii,2))+1; % x dim of matrix
        indbins_lin = mm + sum(gridSize(1:ii-1,1).*gridSize(1:ii-1,2));
        % Assign
        temp(x,y,:) = linmap(indbins_lin,:);
    end
    % Collect output
    gridmap{ii} = temp;
end


