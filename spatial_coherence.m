function [c] = spatial_coherence(spatialvar,gridsize,linmap,n)

% Function to calculate spatial coherence of rate map
% Reflects the similarity of firing rate in adjacent bins of 2D map
% Correlates rate of bin, and mean rate of the surrounding
% bins, over all visited bins. Maps should be non-smoothed.
%
%       [c] = spatial_coherence('place',[40 40],map,1);
%
% Inputs:
%
% - spatialvar: 'place', 'view', 'headdirection'
% - gridsize: [40 40], [1 1;1 1;40 40;40 40; ...; 8 32], [60 1]
% - linmap: linear rate map 1600x1 or 5122x1 or 60x1
% - n: number of bins around central bin to consider in coherence measure. default is 1, so that 8 surrounding bins are considered.
%
% e.g. Kubie and Muller (1989).

%%%%%%%MAKE IT TAKE SHUFFLES. Currently only 1 map


% Make sure map is 2D
if any(size(linmap) == 1)
    linmap = reshape(linmap,length(linmap),1);
    map_grid = lineartogrid(reshape(linmap,length(linmap),1),spatialvar,gridsize);
end

% Pad maps for view and headdirection only 
switch spatialvar
    case 'place'
        map_pad = map_grid;
    case 'view'
        gazeSections = {'Cue', 'Hint', 'Ground', 'Ceiling', 'Walls', 'Pillar1', 'Pillar2', 'Pillar3', 'Pillar4'};
        binDepths = [1 1;
            1 1;
            40 40;
            40 40;
            8 160;
            5 32;
            5 32;
            5 32;
            5 32];
        padpillar = false;
        [emptyfloorref_pad,~] = padsvmap(n,map_grid,gazeSections,padpillar);
        padpillar = true;
        [map_pad,retrievemap] = padsvmap(n,map_grid,gazeSections,padpillar);
    case 'headdirection'
        map_g = map_grid{1};
        map_pad = [map_g(end-n+1:end,1);map_g;map_g(1:1+n-1)];
        map_pad = {map_pad};
end

% Average rates across nearest neighbors n bins away
rates_neighbor_pad = cell(size(gridsize,1),1);
for ii = 1:size(gridsize,1)
    map = map_pad{ii};
    
    rate_neighbor_pad = nan(size(map));
    for aa = 1:size(map,1)
        for bb = 1:size(map,2)
            index1 = aa-n:aa+n;
            index2 = bb-n:bb+n;
            index1(index1<1 | index1>size(map,1)) = [];
            index2(index2<1 | index2>size(map,2)) = [];
            sham = map;
            sham(aa,bb) = nan;
            map_neighbor = sham(index1,index2);
            rate_neighbor_pad(aa,bb) = mean(map_neighbor(:),'omitnan');
        end
    end
    rates_neighbor_pad{ii} = rate_neighbor_pad;
end
% Unpad
switch spatialvar
    case 'place'
        rate_neighbor = rates_neighbor_pad;
    case 'view'
        rate_neighbor = unpadsvmap(rates_neighbor_pad,retrievemap);
    case 'headdirection'
        rates_neighbor_pad = rates_neighbor_pad{1};
        rate_neighbor = rates_neighbor_pad(1+n:end-n);
        rate_neighbor = {rate_neighbor};
end

rate_neighbor = gridtolinear(rate_neighbor,spatialvar,gridsize);
vis = ~isnan(linmap);
rmat = corrcoef(linmap(vis),rate_neighbor(vis),'rows','complete');
r = rmat(2,1);
c = atanh(r);
