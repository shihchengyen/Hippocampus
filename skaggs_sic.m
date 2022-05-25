function [sic_out] = skaggs_sic(ratemaps,posmaps)

% Inputs:
% ratemap - b x n linearised firing rate map (b = num bins, n = num shuffles)
% posmap - b x n linearised occupancy dur map (b = num bins, n = num shuffles)

Pi1 = posmaps./sum(posmaps,1); % consider nansum to play safe
lambda_i = ratemaps;
lambda_i(isnan(lambda_i)) = 0;
lambda_bar = sum(Pi1 .* lambda_i,1);
% divide firing for each position by the overall mean
FRratio = lambda_i./repmat(lambda_bar,size(ratemaps,1),1);
% compute first term in SIC
SIC1 = Pi1 .* lambda_i; 
SIC2 = log2(FRratio);
zeros_placing = SIC1==0;  

bits_per_sec = SIC1 .* SIC2 ./ lambda_bar;
bits_per_sec(zeros_placing) = NaN;
lambda_bar_ok = lambda_bar>0;
lambda_bar_bad = ~lambda_bar_ok;
sic_out = nansum(bits_per_sec, 1);
sic_out(lambda_bar_bad) = NaN;