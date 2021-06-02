% function [ise_out] = ise(actual_image, shuffled_images, dim1, dim2)
% 
%     % function to calculate image spatial entropy, which builds on shannon's entropy but
%     % uses quadrilateral markov random field to take into account spatial
%     % structure. see 'Computation of Image Spatial Entropy Using
%     % Quadrilateral Markov Random Field (2009)' and 'Fast computation
%     % methods for estimation of image spatial entropy (2009)'.
%     
%     % first input is 1x1600, second input is num_shufflesx1600
%     % third input is first reshaped dimension (larger, horizontal length
%     % for pillars and walls)
%     % fourth input is second reshaped dimension (smaller, vertical length
%     % for pillars and walls)
%     
%     % parameters to discretize maps
%     bin_resolution = 0.05; % 0.005 
% 
%     % binning each datapoint
%     actual_disc = floor(actual_image/bin_resolution)+1;
%     actual_disc(isnan(actual_disc)) = 0;
%     shuffled_disc = floor(shuffled_images/bin_resolution)+1;
%     shuffled_disc(isnan(shuffled_disc)) = 0;
%     combined_disc = [actual_disc; shuffled_disc];
%     
%     % reshaping to 2d structure, stacked by shuffles for 3d result
%     combined_disc = reshape(combined_disc, size(combined_disc,1), dim1, dim2);
%     
%     % H(X,Xu) computations
%     temp = cat(4, combined_disc(:,1:end-1,:), combined_disc(:,2:end,:));
%     temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*40 x 2
%     [vert_entropy] = shuffled_joint_entropy(temp);
%     
%     % H(Xl,Xu) computations
%     upper = combined_disc(:,1:end-1,:);
%     centre = combined_disc(:,2:end,:);
%     left = circshift(centre,1,3);
%     temp = cat(4, upper(:,:,2:end), left(:,:,2:end));
%     temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*39 x 2
%     [pos_angled_entropy] = shuffled_joint_entropy(temp);
% 
%     % H(Xr,Xu) computations
%     upper = combined_disc(:,1:end-1,:);
%     centre = combined_disc(:,2:end,:);
%     right = circshift(centre,-1,3);
%     temp = cat(4, upper(:,:,2:end), right(:,:,2:end));
%     temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*39 x 2
%     [neg_angled_entropy] = shuffled_joint_entropy(temp);    
%     
%     % H(Xr/X) computations ~ H(Xr,X) - H(X)
%     temp = cat(4, combined_disc(:,:,1:end-1), combined_disc(:,:,2:end));
%     temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*40 x 2
%     [hor_entropy] = shuffled_joint_entropy(temp);  
%     temp = cat(4, combined_disc(:,:,1:end-1), combined_disc(:,:,1:end-1)); % for H(X)
%     temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*40 x 2
%     [self_entropy] = shuffled_joint_entropy(temp);
%     hor_cond_entropy = hor_entropy - self_entropy;
%     
%     % final calculation
% %     mn = size(combined_disc,1)*size(combined_disc,2);
%     mn = size(combined_disc,2)*size(combined_disc,3);
%     ise_out = mn.*(vert_entropy + hor_cond_entropy);
%     ise_out = ise_out - (mn/2).*(pos_angled_entropy + neg_angled_entropy);
% %     ise_out = ise_out ./ mn;
%      
% end
% 
% 
% function [entropy_summed] = shuffled_joint_entropy(stacked_input)
% 
%     % expects (num_shuffles+1 x num_pairs_x*num_pairs_y x 2) array 
%     % outputs column vector of entropies for all shuffles
% 
%     saved_size2 = size(stacked_input,2);
%     saved_size1 = size(stacked_input,1);
%     temp = permute(stacked_input, [2 1 3]);
%     temp = reshape(temp, size(temp,1)*size(temp,2), size(temp,3)); % now (num_shuffles+1)*39*40 x 2
%     temp = [temp repelem(1:saved_size1, 1, saved_size2)'];
%     temp = sortrows(temp,[3 1 2]);
%     temp(find(sum(temp==0,2)>0),:) = [];
%     temp(:,4) = [0; (diff(temp(:,1))~=0 | diff(temp(:,2))~=0 + diff(temp(:,3))~=0)];
%     occur_count = diff([1; find(temp(:,4)>0); length(temp(:,4))+1]); % occurrences of the corresponding sequences below
%     unique_combi = temp([1; find(temp(:,4)>0)],1:3); % n x 3 (first 2 columns store binned firing rate, last column stores shuffle number)
%     total_count = accumarray(unique_combi(:,3), occur_count);
%     
%     occur_count = occur_count./total_count(unique_combi(:,3)); % converted to probability
%         
%     sum_edges = [diff(unique_combi(:,3)); 1];
%     entropy_base = -occur_count .* log2(occur_count);
%     entropy_summed = cumsum(entropy_base);
%     entropy_summed = entropy_summed(find(sum_edges>0));
%     entropy_summed = diff([0; entropy_summed]);
% 
% end

function [ise_out] = ise(actual_image, shuffled_images, dim1, dim2)

    % function to calculate image spatial entropy, which builds on shannon's entropy but
    % uses quadrilateral markov random field to take into account spatial
    % structure. see 'Computation of Image Spatial Entropy Using
    % Quadrilateral Markov Random Field (2009)' and 'Fast computation
    % methods for estimation of image spatial entropy (2009)'.
    
    % first input is 1x1600, second input is num_shufflesx1600
    % third input is first reshaped dimension (larger, horizontal length
    % for pillars and walls)
    % fourth input is second reshaped dimension (smaller, vertical length
    % for pillars and walls)
    
    % parameters to discretize maps
    bin_resolution = 0.05; % 0.005 

    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    actual_disc(isnan(actual_disc)) = 0;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    shuffled_disc(isnan(shuffled_disc)) = 0;
    combined_disc = [actual_disc; shuffled_disc];
    
    % reshaping to 2d structure, stacked by shuffles for 3d result
    combined_disc = reshape(combined_disc, size(combined_disc,1), dim1, dim2);
    
    % H(X,Xu) computations
    temp = cat(4, combined_disc(:,1:end-1,:), combined_disc(:,2:end,:));
    temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*40 x 2
    [vert_entropy] = shuffled_joint_entropy(temp);
    
    % H(Xl,Xu) computations
    upper = combined_disc(:,1:end-1,:);
    centre = combined_disc(:,2:end,:);
    left = circshift(centre,1,3);
    temp = cat(4, upper(:,:,2:end), left(:,:,2:end));
    temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*39 x 2
    [pos_angled_entropy] = shuffled_joint_entropy(temp);

    % H(Xr,Xu) computations
    upper = combined_disc(:,1:end-1,:);
    centre = combined_disc(:,2:end,:);
    right = circshift(centre,-1,3);
    temp = cat(4, upper(:,:,2:end), right(:,:,2:end));
    temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*39 x 2
    [neg_angled_entropy] = shuffled_joint_entropy(temp);    
    
    % H(Xr/X) computations ~ H(Xr,X) - H(X)
    temp = cat(4, combined_disc(:,:,1:end-1), combined_disc(:,:,2:end));
    temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*40 x 2
    [hor_entropy] = shuffled_joint_entropy(temp);  
    temp = cat(4, combined_disc(:,:,1:end-1), combined_disc(:,:,1:end-1)); % for H(X)
    temp = reshape(temp, size(temp,1), size(temp,2)*size(temp,3), 2); % now num_shuffles+1 x 39*40 x 2
    [self_entropy] = shuffled_joint_entropy(temp);
    hor_cond_entropy = hor_entropy - self_entropy;
    
    % final calculation
%     mn = size(combined_disc,1)*size(combined_disc,2);
    mn = size(combined_disc,3)*size(combined_disc,2);
    ise_out = mn.*(vert_entropy + hor_cond_entropy);
    ise_out = ise_out - (mn/2).*(pos_angled_entropy + neg_angled_entropy);
%     ise_out = ise_out ./ mn;
    
end


function [entropy_summed] = shuffled_joint_entropy(stacked_input)

    % expects (num_shuffles+1 x num_pairs_x*num_pairs_y x 2) array 
    % outputs column vector of entropies for all shuffles

    saved_size2 = size(stacked_input,2);
    saved_size1 = size(stacked_input,1);
    temp = permute(stacked_input, [2 1 3]);
    temp = reshape(temp, size(temp,1)*size(temp,2), size(temp,3)); % now (num_shuffles+1)*39*40 x 2
    temp = [temp repelem(1:saved_size1, 1, saved_size2)'];
    temp = sortrows(temp,[3 1 2]);
    temp(find(sum(temp==0,2)>0),:) = [];
    temp(:,4) = [0; (diff(temp(:,1))~=0 | diff(temp(:,2))~=0 + diff(temp(:,3))~=0)];
    occur_count = diff([1; find(temp(:,4)>0); length(temp(:,4))+1]); % occurrences of the corresponding sequences below
    unique_combi = temp([1; find(temp(:,4)>0)],1:3); % n x 3 (first 2 columns store binned firing rate, last column stores shuffle number)
    total_count = accumarray(unique_combi(:,3), occur_count);
    
    occur_count = occur_count./total_count(unique_combi(:,3)); % converted to probability
        
    sum_edges = [diff(unique_combi(:,3)); 1];
    entropy_base = -occur_count .* log2(occur_count);
    entropy_summed = cumsum(entropy_base);
    entropy_summed = entropy_summed(find(sum_edges>0));
    entropy_summed = diff([0; entropy_summed]);

end




