function [ise_out] = ise5_1(actual_image, shuffled_images, dim1, dim2)
    % Use one pixel neigher, with diagnal relationship
    tic; %show how long ise caculation will take
    
    % parameters to discretize maps
    bin_resolution = 0.05;
    
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
        % empty index array
        n=logical(zeros(1,size(combined(i,:),2)));
        %index
        image=zeros(1,size(combined(i,:),2));
        binLocations = [1:1:max(combined(i,:))];
        [counts] = histcounts(combined(i,:),binLocations);
        %convert count to probability
        counts=counts/sum(counts);
        s=size(counts,2);
        %assign propability of how often does a intensity occur in this image
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
        image = reshape(image, dim1, dim2); %image in the form of its pixel intensity probability
        temp = reshape(combined(i,:), dim1, dim2); %image in the form of its pixel intensity
        %image with zero border
        temp = [zeros(1,dim2); temp ; zeros(1,dim2)]; %add zeros above and below
        temp = [zeros(dim1+2,1) temp zeros(dim1+2,1)]; %add zeros the side
        
        % H(X,Xu) computations   
        upper = temp(1:end-1,:); %Xu
        centre = temp(2:end,:); %X
        h=hist3([reshape(centre,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability 
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total; %convert histgram count to probability
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(j_X_Xu==0)=[]; %remove probability of zero
        j_X_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [vert_entropy] = entropy(j_X_Xu);
        [vert_entropy] = entropy(j_X_Xu)*2; %for both H(X,Xu) and H(X,Xd) 
        
        %   H(Xr,X)
        right = temp(:,2:end); %Xr
        centre = temp(:,1:end-1); %X
        h=hist3([reshape(right,[],1),reshape(centre,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total; %convert histgram count to probability
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(j_Xr_X==0)=[]; %remove probability of zero
        j_Xr_X(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [hor_entropy] = entropy(j_Xr_X);
        [hor_entropy] = entropy(j_Xr_X)*2; %for both H(Xr,X) and H(Xl,X) 
        
        % H(X,Xul) computations compare upper left to the center  
        upper_left = temp(1:end-1,1:end-1); %Xul
        centre = temp(2:end,2:end); %X
        h=hist3([reshape(centre,[],1),reshape(upper_left,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability 
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xul=h/total; %convert histgram count to probability
        j_X_Xul=reshape(j_X_Xul,[],1);
        j_X_Xul(j_X_Xul==0)=[]; %remove probability of zero
        j_X_Xul(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [dia1_entropy] = entropy(j_X_Xul);
        [dia1_entropy] = entropy(j_X_Xul)*2; %for both H(X,Xul) and H(X,Xdr) 
        %because H(X,Xul)= H(X,Xdr) which is compare lower left to the center
        
        % H(X,Xur) computations compare upper right to the center  
        upper_right = temp(1:end-1,2:end); %Xur
        centre = temp(2:end, 1:end-1); %X
        h=hist3([reshape(centre,[],1),reshape(upper_right,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability 
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xur=h/total; %convert histgram count to probability
        j_X_Xur=reshape(j_X_Xur,[],1);
        j_X_Xur(j_X_Xur==0)=[]; %remove probability of zero
        j_X_Xur(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [dia2_entropy] = entropy(j_X_Xur);
        [dia2_entropy] = entropy(j_X_Xur)*2; %for both H(X,Xur) and H(X,Xdl) 
        %because H(X,Xur)= H(X,Xdl) which is compare lower right to the center
        
        %   H(X)
        centre = image; %X
        p_X= reshape(centre,[],1); 
        p_X(p_X==0) = [];
        [self_entropy] = entropy(p_X);
            
        %   final calculation:
        ise_out = vert_entropy + hor_entropy + dia1_entropy + dia2_entropy + self_entropy;
        
        ISE= [ISE ; ise_out];
    end

    ise_out=ISE;
    disp(['time taken to calculate for ise: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
    entropy= sum(entropy);%output sum of all entropy
    
end
