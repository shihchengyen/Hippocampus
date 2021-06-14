% new version of ise based on my understanding
% https://nptel.ac.in/content/storage2/courses/117104069/chapter_1/1_10.html


% Steps:
%   histogram count for each "pixel" [count,edge]=histcounts(image,bin#=20);
%   propability of each "pixel"
%       individual probability
%       joint probability
%   calculate each entropy
%       joint entropy H(X,Xu) and entropy H(X)
%   outputs a (column) vector of entropies for all shuffles
% 
% actual_image= [rand(1,8211)];
% shuffled_images= [rand(10,8211)];
% % [counts,binLocations] = imhist(actual_image)
% dim1=51;
% dim2=161;
% % [ise_out] = ise(actual_image, shuffled_images, dim1, dim2);
function [ise_out] = ise_ted(actual_image, shuffled_images, dim1, dim2)
 
    % parameters to discretize maps
    bin_resolution = 0.05; % 0.005 
%   input image
    actual_disc = floor(actual_image/bin_resolution)+1;
    actual_disc(isnan(actual_disc)) = 0;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    shuffled_disc(isnan(shuffled_disc)) = 0;
    combined = [actual_disc; shuffled_disc]; %1+shuffles x 1600
    length=size(combined,2); %total number of "pixel"= 1600
    
%    for each (image) row of combined calculate the probablity
    prob=[];
    bin=20;
    p=[];
    for row = 1:size(combined,1)
       image=combined(row,:);
       [count,edge]=histcounts(image,bin);
       count_p=count/length;
%        assign propability
%        use: ind=edge(1)<=image(c) && image(c)<edge(2);
%       image(ind)=count_p;
        for c=1:size(combined,2) %make more efficient later
           if edge(1)<=image(c) && image(c)<edge(2)
               p=[p count_p(1)];
           elseif edge(2)<=image(c) && image(c)<edge(3)
               p=[p count_p(2)];
           elseif edge(3)<=image(c) && image(c)<edge(4)
               p=[p count_p(3)];
           elseif edge(4)<=image(c) && image(c)<edge(5)
               p=[p count_p(4)];
           elseif edge(5)<=image(c) && image(c)<edge(6)
               p=[p count_p(5)];
           elseif edge(6)<=image(c) && image(c)<edge(7)
               p=[p count_p(6)];
           elseif edge(7)<=image(c) && image(c)<edge(8)
               p=[p count_p(7)];
           elseif edge(8)<=image(c) && image(c)<edge(9)
               p=[p count_p(8)];
           elseif edge(9)<=image(c) && image(c)<edge(10)
               p=[p count_p(9)];
           elseif edge(10)<=image(c) && image(c)<edge(11)
               p=[p count_p(10)];
           elseif edge(11)<=image(c) && image(c)<edge(12)
               p=[p count_p(11)];
           elseif edge(12)<=image(c) && image(c)<edge(13)
               p=[p count_p(12)];
           elseif edge(13)<=image(c) && image(c)<edge(14)
               p=[p count_p(13)];
           elseif edge(14)<=image(c) && image(c)<edge(15)
               p=[p count_p(14)];
           elseif edge(15)<=image(c) && image(c)<edge(16)
               p=[p count_p(15)];
           elseif edge(16)<=image(c) && image(c)<edge(17)
               p=[p count_p(16)];
           elseif edge(17)<=image(c) && image(c)<edge(18)
               p=[p count_p(17)];
           elseif edge(18)<=image(c) && image(c)<edge(19)
               p=[p count_p(18)];
           elseif edge(19)<=image(c) && image(c)<edge(20)
               p=[p count_p(19)];
           else
               p=[p count_p(20)];
           end    
        end
       prob=[prob(:); p'];
       p=[];
    end
    prob=reshape(prob,size(combined,1),size(combined,2));
%     reshape the prob
    combined = reshape(prob, size(combined,1), dim1, dim2); %now num_shuffles+1 x 51 x 161
    
% find Xu,X,Xl,Xr
    upper = combined(:,1:end-1,:); %Xu
    centre = combined(:,2:end,:); %X
    left = circshift(centre,1,3); %Xl
    right = circshift(centre,-1,3); %Xr
    
%   1) joint probability 2)entropy
    %   for  H(X,Xu) computations
    p_X_Xu=centre .* upper;
    [vert_entropy] = entropy(p_X_Xu);
    
    %   H(Xl,Xu) computations
    p_X1_Xu= left(:,:,2:end) .* upper(:,:,2:end); %idk why yet
    [pos_angled_entropy] = entropy(p_X1_Xu);
    
    %   H(Xr,Xu) computations
    p_Xr_Xu=right(:,:,2:end) .* upper(:,:,2:end);%idk why yet
    [neg_angled_entropy] = entropy(p_Xr_Xu);  
    
    %   H(Xr/X) computations ~ H(Xr,X) - H(X)
    %   H(Xr,X)
    p_Xr_X=combined(:,:,1:end-1).* combined(:,:,2:end); %idk why yet
    [hor_entropy] = entropy(p_Xr_X);
    p_X= combined(:,:,1:end-1); %according to Kian Wei's ise.m
    %   H(X)
    [self_entropy] = entropy(p_X);
%     %original
%     hor_cond_entropy = hor_entropy - self_entropy;
    %modify to avoid negative entropy
    hor_cond_entropy = -hor_entropy + self_entropy;
    
%   final calculation:
      mn = size(combined,2)*size(combined,3);
%       dim1*dim2
      ise_out = mn.*(vert_entropy + hor_cond_entropy);
      ise_out = ise_out - (mn/2).*(pos_angled_entropy + neg_angled_entropy);

end
function [entropy] = entropy(input)

    % expects (num_shuffles+1 x num_pairs_x x num_pairs_y) array 
    % outputs column vector of entropies for all shuffles
    temp=input;
    temp = reshape(input, size(input,1),size(input,2)*size(input,3)); % now (num_shuffles+1) x 39*40 
    temp1 = reshape(input, 1,[]);
    l=log2(temp);
    l2=log2(temp1);
%     should be:
    entropy= -temp .* log2(temp);
% %     without (-) version
%     entropy= temp .* log2(temp);
    entropy= cumsum(entropy,2);
    entropy=entropy(:,end);%only the sum value of each image
    
end
