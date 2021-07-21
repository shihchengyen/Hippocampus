% new version of ise based on my understanding
% https://nptel.ac.in/content/storage2/courses/117104069/chapter_1/1_10.html

function [ise_out] = ise_ted(actual_image, shuffled_images, dim1, dim2)
    tic;
    %discretize the map
    bin_resolution = 0.05; % 0.005 
%     bin_resolution = 0.01; 
%     bin_resolution = 0.005;
    
    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    actual_disc(isnan(actual_disc)) = 0;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    shuffled_disc(isnan(shuffled_disc)) = 0;
    
    combined = [actual_disc; shuffled_disc];


    ISE=[];
    for  i=1:size(combined,1)
        % empty index array
        n=logical(zeros(1,size(combined(i,:),2)));
        %index
        image=zeros(1,size(combined(i,:),2));
        if max(combined(i,:))==0
            disp('a row of NaN');%a row of shuffle is empty
            ise_out = 0;
            ISE= [ISE ; ise_out];
            continue
        end
        if max(combined(i,:))==1
            ise_out = 0;
            ISE= [ISE ; ise_out];
            disp('floor(combined(i,:))=0');%a row of shuffle is empty
            continue
        end
        binLocations = [1:1:max(combined(i,:))];
        [counts] = histcounts(combined(i,:),binLocations);
        %assign propability
        %convert count to probability
        counts=counts/sum(counts);
        s=size(counts,2);
        for ii=1:size(counts,2)
            if ii ~= s
                temp=(binLocations(ii)<=combined(i,:)) + (combined(i,:)<binLocations(ii+1));
                ind=temp==2;
            else
                temp=(binLocations(ii)<=combined(i,:)) + (combined(i,:)<=binLocations(ii+1));
                ind=temp==2;
            end
            image(ind)=counts(ii);
        end
        image = image/sum(image);
        image = reshape(image, dim1, dim2);
%         temp = reshape(image(i,:), dim1, dim2);
        temp = reshape(combined(i,:), dim1, dim2);
%         % find Xu,X,Xl,Xr
%             upper = [zeros(1,dim2); temp(1:end,:)]; %Xu
%             centre = [temp(1:end,:);zeros(1,dim2)]; %X
%             left = [zeros(dim1,1) temp(2:end,1:end-1)]; %Xl
%             right = [temp(2:end,2:end) zeros(dim1,1)]; %Xr
            
        % H(X,Xu) computations
        %1) joint probability 2)entropy    
        upper =[zeros(1,dim2); temp(1:end,:)]; %Xu
        centre = [temp(1:end,:);zeros(1,dim2)]; %X
        h=hist3([reshape(centre,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))});
        total=sum(sum(h))-h(1,1);
        j_X_Xu=h/total;
        j_X_Xu=reshape(j_X_Xu,[],1);
        %remove zero and 0 next to 0 probability
        j_X_Xu(j_X_Xu==0)=[];
        j_X_Xu(1)=[];
        [vert_entropy] = entropy(j_X_Xu)*dim2; %I added
        
        %   H(Xl,Xu) computations
        left =[zeros(dim1+1,1) [temp(1:end,1:end); zeros(1,dim2)]]; %Xl
        upper = [zeros(1,dim2+1); [temp(1:end,:) zeros(dim1,1)]]; %Xu
%         j_X1_Xu= joint_p(reshape(left,[],1),reshape(upper,[],1));
        h=hist3([reshape(left,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))});
        total=sum(sum(h))-h(1,1);
        j_X1_Xu=h/total;
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        %remove zero and 0 next to 0 probability
        j_X1_Xu(j_X1_Xu==0)=[];
        j_X1_Xu(1)=[];
        [pos_angled_entropy] = entropy(j_X1_Xu);
        
        %   H(Xr,Xu) computations
        right =[[temp(1:end,1:end); zeros(1,dim2)] zeros(dim1+1,1)]; %Xr
        upper = [zeros(1,dim2+1); [zeros(dim1,1) temp(1:end,:)]]; %Xu
%         j_Xr_Xu=joint_p(reshape(right,[],1),reshape(upper,[],1));
        h=hist3([reshape(right,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))});
        total=sum(sum(h))-h(1,1);
        j_Xr_Xu=h/total;
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        %remove zero and 0 next to 0 probability
        j_Xr_Xu(j_Xr_Xu==0)=[];
        j_Xr_Xu(1)=[];        
        [neg_angled_entropy] = entropy(j_Xr_Xu);  
        
        %   H(Xr/X) computations ~ H(Xr,X) - H(X)
        %   H(Xr,X)
        right =[temp(1:end,1:end) zeros(dim1,1)]; %Xr
        centre = [zeros(dim1,1) temp(1:end,:)]; %X
        h=hist3([reshape(right,[],1),reshape(centre,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))});
        total=sum(sum(h))-h(1,1);
        j_Xr_X=h/total;
        j_Xr_X=reshape(j_Xr_X,[],1);
        %remove zero and 0 next to 0 probability
        j_Xr_X(j_Xr_X==0)=[];
        j_Xr_X(1)=[];
        [hor_entropy] = entropy(j_Xr_X)*dim1;
        centre = image(2:end,:); %X
        p_X= reshape(centre,[],1); 
        p_X(p_X==0) = [];
%         p_X=p_X/(sum(p_X));
        %   H(X)
        [self_entropy] = entropy(p_X);
        %H(Xr,X) - H(X)
%         hor_cond_entropy = hor_entropy - self_entropy;
        hor_cond_entropy = abs(hor_entropy - self_entropy);
%         hor_cond_entropy = -hor_entropy + self_entropy; %I flipped it
        
        %   final calculation:
        mn = dim1*dim2;
        ise_out = mn.*(vert_entropy + hor_cond_entropy);
        ise_out = ise_out - (mn/2).*(pos_angled_entropy + neg_angled_entropy);
      
        ISE= [ISE ; ise_out];
    end

    ise_out=ISE;
    disp(['time taken to calculate for ise_ted: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects (num_shuffles+1 x num_pairs_x x num_pairs_y) array 
    % outputs column vector of entropies for all shuffles
    temp=input;
%     temp = reshape(input, size(input,1),size(input,2)*size(input,3)); % now (num_shuffles+1) x 39*40 
   
    entropy= -temp .* log2(temp);
    entropy= sum(entropy);%only the sum value of each image
    
end
