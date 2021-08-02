function [linmap] = gridtolinear(gridmap,var,gridSize)

% This function reshapes a 2D grid map to a 1D linear map
% Opposite function to lineartogrid.m

switch var
    case 'place'
        if size(gridmap,1) ~= 1 % || sum(gridSize(:,1).*sum(gridSize(:,2))) ~= sum(gridmap(:,1).*sum(gridmap(:,2)))
            error('Incorrect dimensions specified for grid size');
        end
    case 'spatialview'
        if size(gridmap,1) ~= 9 % || sum(gridSize(:,1).*sum(gridSize(:,2))) ~= sum(gridmap(:,1).*sum(gridmap(:,2)))
            error('Incorrect dimensions specified for grid size');
        end
end
linmap = nan(sum(gridSize(:,1).*gridSize(:,2)),1);
for ii = 1:size(gridmap,1)
    temp = reshape(rot90(gridmap{ii},-1),size(gridmap{ii},1)*size(gridmap{ii},2),1);
    lin_inds = sum(gridSize(1:ii-1,1).*gridSize(1:ii-1,2))+1:sum(gridSize(1:ii,1).*gridSize(1:ii,2));
    linmap(lin_inds,1) = reshape(temp,1,gridSize(ii,1)*gridSize(ii,2));
end