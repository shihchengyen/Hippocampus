function [unpadgrid_sm] = unpadsvmap(padgrid_sm,retrievecoords,unpadgrid_raw)

% This functions removes the padding from smoothed sv maps, a counterpart
% of padsvmap.m
% 
% Inputs:
%  - padgrid_sm: padded grid map that has been smoothed
%  - retrievecoords: the x y limits of the original map within the padded
%  frame
%  - unpadgrid_raw: the unpadded raw grid map for reference of where the
%  empty pillar pixels should be
% 
% Spatial view map output will be in 9 sections
%  - Grid 1: Cue 1 x 1 x nshuff
%  - Grid 2: Hint1 x 1 x nshuff
%  - Grid 3: Floor 40 x 40 x nshuff
%  - Grid 4: Ceiling 40 x 40 x nshuff (top down view) 
%  - Grid 5: Walls 8 x 160 x nshuff, starting from bottom left corner
%  - Grid 6: Pillar 1 bottom right 8 x 32 x nshuff, starting from bottom left corner
%  - Grid 7: Pillar 2 bottom left 8 x 32 x nshuff
%  - Grid 8: Pillar 3 top right 8 x 32 x nshuff
%  - Grid 9: Pillar 4 top left 8 x 32 x nshuff

% Remove padding 
unpadgrid_sm = cell(size(padgrid_sm));
for jj = 1:size(padgrid_sm,1)
    if jj == 1 || jj == 2
        temp = padgrid_sm{jj};
    else
        temp = padgrid_sm{jj}(retrievecoords{jj}(1,1):retrievecoords{jj}(1,2),retrievecoords{jj}(2,1):retrievecoords{jj}(2,2),:);
    end
    if jj == 3 % remove the pillar fills if unpadding floor
        temp(unpadgrid_raw{jj} == 0) = nan;
    end
    unpadgrid_sm{jj} = temp;
end