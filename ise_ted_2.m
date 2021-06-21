function [ise_out] = ise_ted_2(actual_image, shuffled_images, dim1, dim2)
    tic;
    %discretize the map
    bin_resolution = 0.05; % 0.005 

    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    actual_disc(isnan(actual_disc)) = 0;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    shuffled_disc(isnan(shuffled_disc)) = 0;
    
combined = [actual_disc; shuffled_disc];
[counts,binLocations] = imhist(actual_disc);
%assign propability
%convert count to probability
counts=counts/sum(counts);
% empty index array
n=logical(zeros(size(combined,1),size(combined,2)));
%index
s=size(counts,1);
image=zeros(size(combined,1),size(combined,2));
for c=1:size(combined,1)
    i=1;
    for i=1:size(counts)
        if i ~= s
            temp=(binLocations(i)<=combined(c,:)) + (combined(c,:)<binLocations(i+1));
            ind=temp==2;
            n(c,:)=ind;
        else
            temp=(binLocations(i-1)<=combined(c,:)) + (combined(c,:)<=binLocations(i));
            ind=temp==2;
            n(c,:)=ind;
        end
        image(n(c,:))=counts(i);
    end
end




    ISE=[];
    for  i=1:size(combined,1)
        temp = reshape(image(i,:), 1, dim1, dim2);
        % find Xu,X,Xl,Xr
            upper = image(:,1:end-1,:); %Xu
            centre = image(:,2:end,:); %X
            left = circshift(centre,1,3); %Xl
            right = circshift(centre,-1,3); %Xr
        % H(X,Xu) computations
        %1) joint probability 2)entropy    
        j_X_Xu=joint_p(reshape(centre,[],1),reshape(upper,[],1));
        [vert_entropy] = entropy(j_X_Xu);
        %   H(Xl,Xu) computations
        j_X1_Xu= joint_p(reshape(left,[],1),reshape(upper,[],1));
        [pos_angled_entropy] = entropy(j_X1_Xu);
        %   H(Xr,Xu) computations
        j_Xr_Xu=joint_p(reshape(right,[],1),reshape(upper,[],1))
        [neg_angled_entropy] = entropy(j_Xr_Xu);  
        %   H(Xr/X) computations ~ H(Xr,X) - H(X)
        %   H(Xr,X)
        j_Xr_X=joint_p(reshape(right,[],1),reshape(centre,[],1)); 
        [hor_entropy] = entropy(j_Xr_X);
        p_X= centre; %according to Kian Wei's ise.m
        %   H(X)
        [self_entropy] = entropy(p_X);
        %H(Xr,X) - H(X)
        hor_cond_entropy = hor_entropy - self_entropy;
        
        %   final calculation:
        mn = dim2*dim3;
        ise_out = mn.*(vert_entropy + hor_cond_entropy);
        ise_out = ise_out - (mn/2).*(pos_angled_entropy + neg_angled_entropy);
      
        ISE= [ISE ; ise_out];
    end

    ise_out=ISE;
    disp(['time taken to calculate for ise_matlab: ' num2str(toc)]);
end
function [joint]=joint_p(arr1,arr2)
data=[arr1, arr2];

% [bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY)
[bandwidth,density,X,Y]=kde2d(data,8211);

    %need to devide the density with sum of density to normalize the joint
    %probability. 
    s=sum(sum(density));
    joint=density/s;
end
function [entropy] = entropy(input)

    % expects (num_shuffles+1 x num_pairs_x x num_pairs_y) array 
    % outputs column vector of entropies for all shuffles
    temp=input;
    temp = reshape(input, size(input,1),size(input,2)*size(input,3)); % now (num_shuffles+1) x 39*40 
   
    entropy= -temp .* log2(temp);
    entropy= cumsum(entropy,2);
    entropy=entropy(:,end);%only the sum value of each image
    
end
