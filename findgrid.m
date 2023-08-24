% For plotting: Find grid in 3D frame for pixel, get plotting coords
% instead of binning coords
function [gridnum_all,plotx_all,ploty_all,matx_all,maty_all] = findgrid(px_all,objtype)
% Inputs
% 1. px: bin number in linearised map
% 2. object type of linearised map: 'place','view','headdirection'
%
% Outputs:
% 1. grid number of 3D map. e.g. floor of place/spatial view map is grid 3
% 2 and 3. plot coords in 3D map (x goes left to right, y goes bottom
% to top)
% 4 and 5. coords in matrix convention (x goes top to bottom, y goes left
% to right).

px_all = reshape(px_all,length(px_all),1);
nanpx = find(isnan(px_all));
validpx = find(~isnan(px_all));
px = px_all(validpx);
% Set up predefined lookup table
switch objtype
    case 'place'
        mapLdummy = 1:1600;
        gridnum = repmat(3,size(px,1),size(px,2))';
        temp = flipud(reshape(mapLdummy, 40, 40)');
        for ii = 1:mapLdummy(end)
            [plotydummy(ii) plotxdummy(ii)] = find(temp == ii);
            [matxdummy(ii) matydummy(ii)] = find(temp == ii);
        end
        plotydummy = size(temp,1)-plotydummy+1;
    case 'view'
        mapLdummy = 1:5122;
        griddummy = mapLdummy;
        for gg = 1:9
            switch gg 
                case 1
                    griddummy(1) = gg;
                    temp{gg} = 1;
                case 2
                    griddummy(2) = gg;
                    temp{gg} = 2;
                case 3
                    griddummy(3:1602) = gg;
                    temp{gg} = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');
                case 4
                    griddummy(1603:3202) = gg;
                    temp{gg} = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');
                case 5
                    griddummy(3203:4482) = gg;
                    temp{gg} = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');
                case 6
                    griddummy(4483:4642) = gg;
                    temp{gg} = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
                case 7
                    griddummy(4643:4802) = gg;
                    temp{gg} = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
                case 8
                    griddummy(4803:4962) = gg;
                    temp{gg} = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
                case 9
                    griddummy(4963:5122) = gg;
                    temp{gg} = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
            end
        end
        for ii = 1:mapLdummy(end)
            [plotydummy(ii) plotxdummy(ii)] = find(temp{griddummy(ii)} == ii);
            [matxdummy(ii) matydummy(ii)] = find(temp{griddummy(ii)} == ii);

            plotydummy(ii) = size(temp{griddummy(ii)},1)-plotydummy(ii)+1;
        end

        gridnum = griddummy(px)';

        % if px == 1 % Cue
        %     gridnum = 1;
        %     x = 1;
        %     y = 1;
        % elseif px == 2 % Hint
        %     gridnum = 2;
        %     x = 1; 
        %     y = 1;
        % elseif px >= 3 && px <= 1602 % Floor
        %     gridnum = 3;
        %     temp = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');
        % elseif px >= 1603 && px <= 3202 % Ceiling
        %     gridnum = 4;
        %     temp = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');
        % elseif px >= 3203 && px <= 4482 % Walls
        %     gridnum = 5;
        %     temp = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');
        % elseif px >= 4483 && px <= 4642 % Pillar 1
        %     gridnum = 6;
        %     temp = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
        % elseif px >= 4643 && px <= 4802 % Pillar 2
        %     gridnum = 7;
        %     temp = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
        % elseif px >= 4803 && px <= 4962 % Pillar 3
        %     gridnum = 8;
        %     temp = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
        % elseif px >= 4963 && px <= 5122 % Pillar 4
        %     gridnum = 9;
        %     temp = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
        % end
    case 'headdirection'
        mapLdummy = (1:60)';
        gridnum = ones(size(px,1),size(px,2))';
        % temp = flipud(reshape(mapLdummy, 60, 1)');
        temp = mapLdummy;
        for ii = 1:mapLdummy(end)
            [plotydummy(ii) plotxdummy(ii)] = find(temp == ii);
            [matxdummy(ii) matydummy(ii)] = find(temp == ii);
        end
end
plotx= plotxdummy(px);
ploty = plotydummy(px);
matx = matxdummy(px);
maty = matydummy(px);

plotx_all = nan(size(px_all));
ploty_all = nan(size(px_all));
matx_all = nan(size(px_all));
maty_all = nan(size(px_all));
gridnum_all = nan(size(px_all));

plotx_all(validpx) = plotx;
ploty_all(validpx) = ploty;
matx_all(validpx) = matx;
maty_all(validpx) = maty;
gridnum_all(validpx) = gridnum;



% [y,x] = find(temp == px);
% y = size(temp,1)-y+1;