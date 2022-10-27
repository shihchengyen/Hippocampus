function [map_sm]=smoothdir(map_raw,n,dirSteps)
% Sliding window average of n bins
% raw map input needs to be in column form. 
% M dir bins by N shuffles
dim1 = size(map_raw,1);
dim2 = size(map_raw,2);
flip = false;
if dim1 ~= dirSteps
    flip = true;
    map_raw = map_raw';
end
% Smooth a dir map.
if n==1; map_sm=map_raw; return; end
p = (n-1)/2;                                               % Pad for circular smooth
pad_map = [map_raw(end-p+1:end,:); map_raw; map_raw(1:p,:)];                    %  ..
map_sm = mean( im2col(pad_map, [n 1], 'sliding') );
% Reshape
map_sm = reshape(map_sm,dirSteps,size(map_sm,2)/dirSteps);
if flip
    map_sm = map_sm';
end