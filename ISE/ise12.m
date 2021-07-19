function [ise_out] = ise12(actual_image, shuffled_images, dim1, dim2)
    % Use QMRF
    tic; %show how long ise caculation will take
    
    % parameters to discretize maps
    bin_resolution = 0.005;
    
    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    actual_disc(isnan(actual_disc)) = 0;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    shuffled_disc(isnan(shuffled_disc)) = 0;
    
    combined = [actual_disc; shuffled_disc]; %(1+shuffle) x dim1*dim2
    
    %Image spatial entropy calculation
    ISE=[];
    for  i=1:size(combined,1)
        if max(combined(i,:))==0 %an empty image don't have ISE
            disp('an image consist of NaN');
            ise_out = 0;
            ISE= [ISE ; ise_out];
            continue
        end
        %an image being discrtized by not small enough bin_reolution is
        %excluded from ISE calculation
        if max(combined(i,:))==1 
            ise_out = 0;
            ISE= [ISE ; ise_out];
            disp('floor(combined(i,:))=0');%every intensity value rounded to zero
            continue
        end
%         % empty index array
%         n=logical(zeros(1,size(combined(i,:),2)));
%         %index
%         image=zeros(1,size(combined(i,:),2));
%         binLocations = [1:1:max(combined(i,:))];
%         [counts] = histcounts(combined(i,:),binLocations);
%         %convert count to probability
%         counts=counts/sum(counts);
%         s=size(counts,2);
%         %assign propability of how often does a intensity occur in this image
%         for ii=1:size(counts,2)
%             if ii ~= s
%                 temp=(binLocations(ii)<=combined(i,:)) + (combined(i,:)<binLocations(ii+1));
%                 ind=temp==2;
%             else
%                 temp=(binLocations(ii)<=combined(i,:)) + (combined(i,:)<=binLocations(ii+1));
%                 ind=temp==2;
%             end
%             image(ind)=counts(ii);
%         end
%         image = reshape(image, dim1, dim2); %image in the form of its pixel intensity probability
        temp = reshape(combined(i,:), dim1, dim2); %image in the form of its pixel intensity
        %image with two layer of zero border
        temp = [zeros(2,dim2); temp ; zeros(2,dim2)]; %add zeros above and below
        temp = [zeros(dim1+4,2) temp zeros(dim1+4,2)]; %add zeros the side
        
        %two pixel
        % H(X,Xu) computations   
        upper =temp(1:end-2,:); %Xu
        centre = temp(3:end,:); %X
        h=hist3([reshape(centre,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total;%convert histgram count to probability
%         j_X_Xu=h^2/total; 
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(j_X_Xu==0)=[]; %remove probability of zero
        j_X_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [vert_entropy] = entropy(j_X_Xu);
        
        %   H(Xl,Xu) computations
        left =temp(3:end,1:end-2); %Xl
        upper = temp(1:end-2,3:end); %Xu
        h=hist3([reshape(left,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X1_Xu=h/total; %convert histgram count to probability
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        j_X1_Xu(j_X1_Xu==0)=[]; %remove probability of zero
        j_X1_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [pos_angled_entropy] = entropy(j_X1_Xu);
        
        %   H(Xr,Xu) computations
        right =temp(3:end,3:end); %Xr
        upper = temp(1:end-2,1:end-2); %Xu
        h=hist3([reshape(right,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_Xu=h/total; %convert histgram count to probability
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        j_Xr_Xu(j_Xr_Xu==0)=[]; %remove probability of zero
        j_Xr_Xu(1)=[];  %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3       
        [neg_angled_entropy] = entropy(j_Xr_Xu);  
        
        %   H(Xr/X) computations = H(Xr,X) - H(X)
        %   H(Xr,X)
        right =temp(:,3:end); %Xr
        centre = temp(:,1:end-2); %X
        h=hist3([reshape(right,[],1),reshape(centre,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total; %convert histgram count to probability
%         j_Xr_X=h^2/total;
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(j_Xr_X==0)=[]; %remove probability of zero
        j_Xr_X(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [hor_entropy] = entropy(j_Xr_X);
        
%         %   H(X)
%         centre = image; %X
%         p_X= reshape(centre,[],1); 
%         p_X(p_X==0) = [];
%         p_X=unique(p_X);
%         [self_entropy] = entropy(p_X);
%         %H(Xr,X) - H(X)
%         hor_cond_entropy = hor_entropy - self_entropy;
%         
        %   final calculation:
        mn = dim1*dim2;
        ise_out = (vert_entropy + hor_entropy +pos_angled_entropy +neg_angled_entropy);
        ise_out = ise_out./mn;
        
        ISE= [ISE ; ise_out];
    end

    [ise_out]=ISE;
    disp(['time taken to calculate for ise11: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
    entropy= sum(entropy);%output sum of all entropy
    
end
