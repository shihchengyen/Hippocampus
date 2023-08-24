function [linmap] = gridtolinear(gridmap,var,gridSize)

% This function reshapes a 2D x nshuff grid map to a 1D x nshuff linear map
% Opposite function to lineartogrid.m

switch var
    case 'place'
        if size(gridmap,1) ~= 1 % || sum(gridSize(:,1).*sum(gridSize(:,2))) ~= sum(gridmap(:,1).*sum(gridmap(:,2)))
            error('Incorrect dimensions specified for grid size');
        end
    case 'view'
        if size(gridmap,1) ~= 9 % || sum(gridSize(:,1).*sum(gridSize(:,2))) ~= sum(gridmap(:,1).*sum(gridmap(:,2)))
            error('Incorrect dimensions specified for grid size');
        end
    case 'headdirection'
        if size(gridmap,1) ~= 1 % || sum(gridSize(:,1).*sum(gridSize(:,2))) ~= sum(gridmap(:,1).*sum(gridmap(:,2)))
            error('Incorrect dimensions specified for grid size');
        end
end
linmap = nan(sum(gridSize(:,1).*gridSize(:,2)),size(gridmap{1},3));
for ii = 1:size(gridmap,1)
    if strcmp(var,'place') | strcmp(var,'view') % matrix and plor coords are opposite for place and view, but not hd
        temp = reshape(rot90(gridmap{ii},-1),size(gridmap{ii},1)*size(gridmap{ii},2),size(gridmap{1},3));
    else
        temp = reshape(gridmap{ii},size(gridmap{ii},1)*size(gridmap{ii},2),size(gridmap{1},3));
    end
    lin_inds = sum(gridSize(1:ii-1,1).*gridSize(1:ii-1,2))+1:sum(gridSize(1:ii,1).*gridSize(1:ii,2));
%     linmap(lin_inds,:) = reshape(temp,1,gridSize(ii,1)*gridSize(ii,2),size(gridmap{1},3));
    linmap(lin_inds,:) = temp;
end