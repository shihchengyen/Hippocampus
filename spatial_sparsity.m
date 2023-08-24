function [sparsity3] = spatial_sparsity(occ,rate)

% Function to calculate sparsity of a spatial map over the occupied bins
% Currently there are 3 implementations
%
% (Markus et al 1994) s = proportion of env in which cell fires
% (Wirth et al 2017) 1 for firing concentrated to 1 bin, 0 when homogeneous
% (O'Neill Csicsvari 2017/Jung McNaughton 1994) s = proportion of env in which cell fires, corrected for dwell time

% Check inputs for correct dimensions
if min(size(occ)) > 1 || size(occ,2) ~= 1
    occ = reshape(occ,size(occ,1)*size(occ,2),1);
end
if min(size(rate)) > 1 || size(rate,2) ~= 1
    rate = reshape(rate,size(rate,1)*size(rate,2),1);
end

% Remove unvisited bins
rate(isnan(rate)) = [];
occ(occ==0) = [];
if length(rate) ~= length(occ)
    error('occupied bins in rate and occ are different');
end

% Calculate sparsity 

% (Markus et al 1994)
p_i = occ./sum(occ);
numer = p_i.*(rate.*rate);
denom = mean(rate)*mean(rate);
sparsity = sum(numer/denom);

% (Wirth et al 2017) 
numbin = length(occ);
numer = numbin - (sum(rate)*sum(rate))/sum(rate.*rate);
denom = numbin - 1;
sparsity2 = numer/denom;

% (O'Neill Csicsvari 2017/Jung McNaughton 1994) 
numer = sum(p_i.*rate)*sum(p_i.*rate);
denom = sum(p_i.*(rate.*rate));
sparsity3 = numer/denom;
