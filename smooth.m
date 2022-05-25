function [map_sm] = smooth(maps_raw,k,unvis,smoothtype)

% Smooth a place map. map=map, k=boxcar kernel, unvis=index of unvisited bins.
% inputs: double or logical, no cell
% 
% maps_raw: grid, not linear. b x b x n dimensions. b = num bins, n = num shuffles
% unvis: b x b x n array of unvisited bins
% smoothtype: 'boxcar' or 'disk'
% default k = 5

if max(size(k)) ~= 1 || any(rem(k,1)~=0)
    error('Input smoothing kernel as single integer');
end

% Define smoothing kernel
switch smoothtype
    case 'boxcar'
        k = ones(k); % k by k square
    case 'disk'
        k = fspecial('disk',k);
end
% 
% % Extract nan values in original map
% unvis = isnan(map_raw);

% Set nan values to zero to work with imfilter
maps_raw(unvis)=0;
map_filt=imfilter(maps_raw,k);

% Filter logical occupancy map to account for nan-zero edges
occ=ones(size(maps_raw));
occ(unvis)=0;
occ_filt=imfilter(occ,k);

% Normalize filtered rate map by filtered occupancy map
% warning('off', 'MATLAB:divideByZero');
map_sm=map_filt./occ_filt;
% warning('on', 'MATLAB:divideByZero');

% Restore nan values
map_sm(unvis)=nan;