function [out,retrievemap] = padsvmap(n,gridmap,gazeSections)

% This function pads each section of a spatial view map with bins from
% adjoining sections
% Spatial view map should be in 9 sections
%  - Grid 1: Cue 1 x 1 x nshuff
%  - Grid 2: Hint 1 x 1 x nshuff
%  - Grid 3: Floor 40 x 40 x nshuff
%  - Grid 4: Ceiling 40 x 40 x nshuff (top down view) 
%  - Grid 5: Walls 8 x 160 x nshuff, starting from bottom left corner
%  - Grid 6: Pillar 1 bottom right 8 x 32 x nshuff, starting from bottom left corner
%  - Grid 7: Pillar 2 bottom left 8 x 32 x nshuff
%  - Grid 8: Pillar 3 top right 8 x 32 x nshuff
%  - Grid 9: Pillar 4 top left 8 x 32 x nshuff

if size(gridmap,1) ~= 9
%     gazeSections = {'Cue', 'Hint', 'Ground', 'Ceiling', 'Walls', 'Pillar1', 'Pillar2', 'Pillar3', 'Pillar4'};
    error('Not enough sections of grids for spatial view map');
end

% Initialize output
retrievemap = cell(size(gridmap,1),1);
out = cell(size(gridmap,1),1);
for jj = 1:size(gridmap,1)
    switch gazeSections{jj}
        case 'Cue' % No need to pad
            map = gridmap{jj};
            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [1 1; ...
                           1 1];
            % Send vars for adaptive smoothing
            out{jj} = map;
        case 'Hint' % No need to pad
            map = gridmap{jj};
            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [1 1; ...
                           1 1];
            % Send vars for adaptive smoothing
            out{jj} = map;
        case 'Ground'
            map = gridmap{jj};
            wallsection_ind = strcmp(gazeSections,'Walls');
            wallmap = gridmap{wallsection_ind};
            pillar1_ind = strcmp(gazeSections,'Pillar1'); % bottom right
            pillar2_ind = strcmp(gazeSections,'Pillar2'); % bottom left
            pillar3_ind = strcmp(gazeSections,'Pillar3'); % top right
            pillar4_ind = strcmp(gazeSections,'Pillar4'); % top left
            pillar1map = gridmap{pillar1_ind};
            pillar2map = gridmap{pillar2_ind};
            pillar3map = gridmap{pillar3_ind};
            pillar4map = gridmap{pillar4_ind};

            % Move original map to middle
            maptemp = zeros(size(map,1)+2*n,size(map,2)+2*n,size(map,3));
            maptemp(n+1:n+size(map,1), n+1:n+size(map,2),:) = map;

            % Pad with wall data
            maptemp(1:n,n+1:n+size(map,1),:) = wallmap(size(wallmap,1)-n+1:end,1*size(map,1)+1:2*size(map,1),:); % top
            maptemp(n+1:n+size(map,1),size(map,1)+n+1:end,:) = rot90(wallmap(size(wallmap,1)-n+1:end,2*size(map,1)+1:3*size(map,1),:),-1); % right
            maptemp(size(map,1)+n+1:end,n+1:size(map,1)+n,:) = rot90(wallmap(size(wallmap,1)-n+1:end,3*size(map,1)+1:4*size(map,1),:),-2); % bottom
            maptemp(n+1:size(map,1)+n,1:n,:) = rot90(wallmap(size(wallmap,1)-n+1:end,0*size(map,1)+1:1*size(map,1),:),1); % left

            % Pad with pillar1 data - additively, since I can't figure out how not to overlap the pixels
            plpx = 4; % pad with 3, not 5 pixels, to leave a bit of space in the middle
            temp1 = zeros(8,8,size(map,3)); temp2 = temp1; temp3 = temp1; temp4 = temp1;
            temp1(:,1:plpx,:) = flipud(rot90(pillar1map(5-plpx+1:end,1:8,:),-1));
            temp2(1:plpx,:,:) = flipud(pillar1map(5-plpx+1:end,9:16,:));
            temp3(:,8-plpx+1:8,:) = fliplr(rot90(pillar1map(5-plpx+1:end,17:24,:),-1));
            temp4(8-plpx+1:8,:,:) = fliplr(pillar1map(5-plpx+1:end,25:end,:));
            temp = (temp1+temp2+temp3+temp4)/4;
            maptemp(n+25:n+32,n+25:n+32,:) = temp; % bottom right
            % Pad with pillar2 data
            temp1 = zeros(8,8,size(map,3)); temp2 = temp1; temp3 = temp1; temp4 = temp1;
            temp1(:,1:plpx,:) = flipud(rot90(pillar2map(5-plpx+1:end,1:8,:),-1));
            temp2(1:plpx,:,:) = flipud(pillar2map(5-plpx+1:end,9:16,:));
            temp3(:,8-plpx+1:8,:) = fliplr(rot90(pillar2map(5-plpx+1:end,17:24,:),-1));
            temp4(8-plpx+1:8,:,:) = fliplr(pillar2map(5-plpx+1:end,25:end,:));
            temp = (temp1+temp2+temp3+temp4)/4;
            maptemp(n+25:n+32,n+9:n+16,:) = temp; % bottom left
            % Pad with pillar 3 data
            temp1 = zeros(8,8,size(map,3)); temp2 = temp1; temp3 = temp1; temp4 = temp1;
            temp1(:,1:plpx,:) = flipud(rot90(pillar3map(5-plpx+1:end,1:8,:),-1));
            temp2(1:plpx,:,:) = flipud(pillar3map(5-plpx+1:end,9:16,:));
            temp3(:,8-plpx+1:8,:) = fliplr(rot90(pillar3map(5-plpx+1:end,17:24,:),-1));
            temp4(8-plpx+1:8,:,:) = fliplr(pillar3map(5-plpx+1:end,25:end,:));
            temp = (temp1+temp2+temp3+temp4)/4;
            maptemp(n+9:n+16,n+25:n+32,:) = temp; % top right
            % Pad with pillar 4 data
            temp1 = zeros(8,8,size(map,3)); temp2 = temp1; temp3 = temp1; temp4 = temp1;
            temp1(:,1:plpx,:) = flipud(rot90(pillar4map(5-plpx+1:end,1:8,:),-1));
            temp2(1:plpx,:,:) = flipud(pillar4map(5-plpx+1:end,9:16,:));
            temp3(:,8-plpx+1:8,:) = fliplr(rot90(pillar4map(5-plpx+1:end,17:24,:),-1));
            temp4(8-plpx+1:8,:,:) = fliplr(pillar4map(5-plpx+1:end,25:end,:));
            temp = (temp1+temp2+temp3+temp4)/4;
            maptemp(n+9:n+16,n+9:n+16,:) = temp; % top left

            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [n+1 n+size(map,1); ...
                           n+1 n+size(map,2)];
            % Send vars for adaptive smoothing
            out{jj} = maptemp;

        case 'Ceiling'
            map = gridmap{jj};
            wallsection_ind = strcmp(gazeSections,'Walls');
            wallmap = gridmap{wallsection_ind};

            % Flip walldata upside down
            wallmap = flipud(wallmap);

            % Move original map to middle
            maptemp = zeros(size(map,1)+2*n,size(map,2)+2*n,size(map,3));
            maptemp(n+1:n+size(map,1), n+1:n+size(map,2),:) = map;

            % Pad with wall data
            maptemp(1:n,n+1:n+size(map,1),:) = wallmap(size(wallmap,1)-n+1:end,1*size(map,1)+1:2*size(map,1),:); % top
            maptemp(n+1:n+size(map,1),size(map,1)+n+1:end,:) = rot90(wallmap(size(wallmap,1)-n+1:end,2*size(map,1)+1:3*size(map,1),:),-1); % right
            maptemp(size(map,1)+n+1:end,n+1:size(map,1)+n,:) = rot90(wallmap(size(wallmap,1)-n+1:end,3*size(map,1)+1:4*size(map,1),:),-2); % bottom
            maptemp(n+1:size(map,1)+n,1:n,:) = rot90(wallmap(size(wallmap,1)-n+1:end,0*size(map,1)+1:1*size(map,1),:),1); % left

            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [n+1 n+size(map,1); ...
                           n+1 n+size(map,2)];
            % Send vars for adaptive smoothing
            out{jj} = maptemp;

        case 'Walls'
            map = gridmap{jj};
            groundsection_ind = strcmp(gazeSections,'Ground');
            ground_map = gridmap{groundsection_ind};

            ceilingsection_ind = strcmp(gazeSections,'Ceiling');
            ceiling_map = gridmap{ceilingsection_ind};

            % Move original map to middle
            maptemp = zeros(size(map,1)+2*n,size(map,2)+2*n,size(map,3));
            maptemp(n+1:n+size(map,1), n+1:n+size(map,2),:) = map;

            % Pad with ground data
            maptemp(n+size(map,1)+1:end,n+1:size(ground_map,2)+n,:) = rot90(ground_map(:,1:n,:),-1);
            maptemp(n+size(map,1)+1:end,n+size(ground_map,2)+1:n+2*size(ground_map,2),:) = ground_map(1:n,:,:);
            maptemp(n+size(map,1)+1:end,n+2*size(ground_map,2)+1:n+3*size(ground_map,2),:) = rot90(ground_map(:,size(ground_map,1)-n+1:end,:),1);
            maptemp(n+size(map,1)+1:end,n+3*size(ground_map,1)+1:n+4*size(ground_map,1),:) = rot90(ground_map(size(ground_map,1)-n+1:end,:,:),2);
            
            % Pad with ceiling data
            maptemp(1:n,n+1:size(ceiling_map,1)+n,:) = fliplr(rot90(ceiling_map(:,1:n,:),1));
            maptemp(1:n,n+size(ceiling_map,1)+1:n+2*size(ceiling_map,1),:) = flipud(ceiling_map(1:n,:,:));
            maptemp(1:n,n+2*size(ceiling_map,1)+1:n+3*size(ceiling_map,1),:) = fliplr(rot90(ceiling_map(:,size(ceiling_map,1)-n+1:end,:),-1));
            maptemp(1:n,n+3*size(ceiling_map,1)+1:n+4*size(ceiling_map,1),:) = fliplr(ceiling_map(size(ceiling_map,1)-n+1:end,:,:));

            % Pad with wall data on either end
            maptemp(n+1:n+size(map,1),1:n,:) = map(:,size(map,2)-n+1:end,:);
            maptemp(n+1:n+size(map,1),size(maptemp,2)-n+1:end,:) = map(:,1:n,:);

            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [n+1 n+size(map,1); ...
                           n+1 n+size(map,2)];
            % Send vars for adaptive smoothing
            out{jj} = maptemp;

        case 'Pillar1' % Bottom right
            map = gridmap{jj};
            groundsection_ind = strcmp(gazeSections,'Ground');
            ground_map = gridmap{groundsection_ind};

            % Move original map to middle
            maptemp = zeros(size(map,1)+n,size(map,2)+2*n,size(map,3));
            maptemp(1:size(map,1), n+1:n+size(map,2),:) = map;

            % Pad with ground data
            maptemp(size(map,1)+1:end,n+1:(size(map,2)/4)+n,:) = flipud(rot90(ground_map(25:32,25-n:24,:),-1));
            maptemp(size(map,1)+1:end,n+(size(map,2)/4)+1:n+2*(size(map,2)/4),:) = flipud(ground_map(25-n:24,25:32,:));
            maptemp(size(map,1)+1:end,n+2*(size(map,2)/4)+1:n+3*(size(map,2)/4),:) = flipud(rot90(ground_map(25:32,33:32+n,:),1));
            maptemp(size(map,1)+1:end,n+3*(size(map,2)/4)+1:n+4*(size(map,2)/4),:) = flipud(rot90(ground_map(33:32+n,25:32,:),2));

            % Pad with pillar data on either end
            maptemp(1:size(map,1),1:n,:) = map(:,size(map,2)-n+1:end,:);
            maptemp(1:size(map,1),size(maptemp,2)-n+1:end,:) = map(:,1:n,:);

            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [1 size(map,1); ...
                           n+1 n+size(map,2)];
            % Send vars for adaptive smoothing
            out{jj} = maptemp;

        case 'Pillar2' % bottom left
            map = gridmap{jj};
            groundsection_ind = strcmp(gazeSections,'Ground');
            ground_map = gridmap{groundsection_ind};

            % Move original map to middle
            maptemp = zeros(size(map,1)+n,size(map,2)+2*n,size(map,3));
            maptemp(1:size(map,1), n+1:n+size(map,2),:) = map;

            % Pad with ground data
            maptemp(size(map,1)+1:end,n+1:(size(map,2)/4)+n,:) = flipud(rot90(ground_map(25:32,9-n:8,:),-1));
            maptemp(size(map,1)+1:end,n+(size(map,2)/4)+1:n+2*(size(map,2)/4),:) = flipud(ground_map(25-n:24,9:16,:));
            maptemp(size(map,1)+1:end,n+2*(size(map,2)/4)+1:n+3*(size(map,2)/4),:) = flipud(rot90(ground_map(25:32,17:16+n,:),1));
            maptemp(size(map,1)+1:end,n+3*(size(map,2)/4)+1:n+4*(size(map,2)/4),:) = flipud(rot90(ground_map(33:32+n,9:16,:),2));

            % Pad with pillar data on either end
            maptemp(1:size(map,1),1:n,:) = map(:,size(map,2)-n+1:end,:);
            maptemp(1:size(map,1),size(maptemp,2)-n+1:end,:) = map(:,1:n,:);

            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [1 size(map,1); ...
                           n+1 n+size(map,2)];
            % Send vars for adaptive smoothing
            out{jj} = maptemp;

        case 'Pillar3' % top right
            map = gridmap{jj};
            groundsection_ind = strcmp(gazeSections,'Ground');
            ground_map = gridmap{groundsection_ind};

            % Move original map to middle
            maptemp = zeros(size(map,1)+n,size(map,2)+2*n,size(map,3));
            maptemp(1:size(map,1), n+1:n+size(map,2),:) = map;

            % Pad with ground data
            maptemp(size(map,1)+1:end,n+1:(size(map,2)/4)+n,:) = flipud(rot90(ground_map(9:16,25-n:24,:),-1));
            maptemp(size(map,1)+1:end,n+(size(map,2)/4)+1:n+2*(size(map,2)/4),:) = flipud(ground_map(9-n:8,25:32,:));
            maptemp(size(map,1)+1:end,n+2*(size(map,2)/4)+1:n+3*(size(map,2)/4),:) = flipud(rot90(ground_map(9:16,33:32+n,:),1));
            maptemp(size(map,1)+1:end,n+3*(size(map,2)/4)+1:n+4*(size(map,2)/4),:) = flipud(rot90(ground_map(17:16+n,25:32,:),2));

            % Pad with pillar data on either end
            maptemp(1:size(map,1),1:n,:) = map(:,size(map,2)-n+1:end,:);
            maptemp(1:size(map,1),size(maptemp,2)-n+1:end,:) = map(:,1:n,:);

            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [1 size(map,1); ...
                           n+1 n+size(map,2)];
            % Send vars for adaptive smoothing
            out{jj} = maptemp;

        case 'Pillar4' % top left
            map = gridmap{jj};
            groundsection_ind = strcmp(gazeSections,'Ground');
            ground_map = gridmap{groundsection_ind};

            % Move original map to middle
            maptemp = zeros(size(map,1)+n,size(map,2)+2*n,size(map,3));
            maptemp(1:size(map,1), n+1:n+size(map,2),:) = map;

            % Pad with ground data
            maptemp(size(map,1)+1:end,n+1:(size(map,2)/4)+n,:) = flipud(rot90(ground_map(9:16,9-n:8,:),-1));
            maptemp(size(map,1)+1:end,n+(size(map,2)/4)+1:n+2*(size(map,2)/4),:) = flipud(ground_map(9-n:8,9:16,:));
            maptemp(size(map,1)+1:end,n+2*(size(map,2)/4)+1:n+3*(size(map,2)/4),:) = flipud(rot90(ground_map(9:16,17:16+n,:),1));
            maptemp(size(map,1)+1:end,n+3*(size(map,2)/4)+1:n+4*(size(map,2)/4),:) = flipud(rot90(ground_map(17:16+n,9:16,:),2));

            % Pad with pillar data on either end
            maptemp(1:size(map,1),1:n,:) = map(:,size(map,2)-n+1:end,:);
            maptemp(1:size(map,1),size(maptemp,2)-n+1:end,:) = map(:,1:n,:);

            % Save indices of original grid [from_x to_x; from_y to_y]
            retrievemap{jj} = [1 size(map,1); ...
                           n+1 n+size(map,2)];
            % Send vars for adaptive smoothing
            out{jj} = maptemp;
    end
end