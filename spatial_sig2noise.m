function [sig2noise] = spatial_sig2noise(rate)

% Function to calculate the signal to noise ratio of a rate map

% Check inputs for correct dimensions
if min(size(rate)) > 1 || size(rate,2) ~= 1
    rate = reshape(rate,size(rate,1)*size(rate,2),1);
end

% Remove unvisited bins
rate(isnan(rate)) = [];

sig2noise = max(rate)/mean(rate);