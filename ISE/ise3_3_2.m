function [ise_out] = ise3_3_2(actual_image, shuffled_images, dim1, dim2)
    % Use QMRF
    tic; %show how long ise caculation will take
    
    % parameters to discretize maps
    bin_resolution = 0.01;
    
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
        
        % H(X,Xu) computations   
        upper =[zeros(1,dim2); temp(1:end,:)]; %Xu
        centre = [temp(1:end,:);zeros(1,dim2)]; %X
        h=hist3([reshape(centre,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability 
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total; %convert histgram count to probability
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(j_X_Xu==0)=[]; %remove probability of zero
        j_X_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [vert_entropy] = entropy(j_X_Xu);
        
        %   H(Xl,Xu) computations
        left =[zeros(dim1+1,1) [temp(1:end,1:end); zeros(1,dim2)]]; %Xl
        upper = [zeros(1,dim2+1); [temp(1:end,:) zeros(dim1,1)]]; %Xu
        h=hist3([reshape(left,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X1_Xu=h/total; %convert histgram count to probability
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        j_X1_Xu(j_X1_Xu==0)=[]; %remove probability of zero
        j_X1_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [pos_angled_entropy] = entropy(j_X1_Xu);
        
        %   H(Xr,Xu) computations
        right =[[temp(1:end,1:end); zeros(1,dim2)] zeros(dim1+1,1)]; %Xr
        upper = [zeros(1,dim2+1); [zeros(dim1,1) temp(1:end,:)]]; %Xu
        h=hist3([reshape(right,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_Xu=h/total; %convert histgram count to probability
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        j_Xr_Xu(j_Xr_Xu==0)=[]; %remove probability of zero
        j_Xr_Xu(1)=[];  %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3       
        [neg_angled_entropy] = entropy(j_Xr_Xu);  
        
        %   H(Xr/X) computations = H(Xr,X) - H(X)
        %   H(Xr,X)
        right =[temp(1:end,1:end) zeros(dim1,1)]; %Xr
        centre = [zeros(dim1,1) temp(1:end,:)]; %X
        h=hist3([reshape(right,[],1),reshape(centre,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total; %convert histgram count to probability
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(j_Xr_X==0)=[]; %remove probability of zero
        j_Xr_X(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [hor_entropy] = entropy(j_Xr_X);
        centre = image; %X
        p_X= reshape(centre,[],1); 
        p_X(p_X==0) = [];
        %   H(X)
        [self_entropy] = entropy(p_X);
        %H(Xr,X) - H(X)
%         hor_cond_entropy = hor_entropy - self_entropy;
        hor_cond_entropy = abs(hor_entropy - self_entropy);% so won't have negative entropy
        
        %   final calculation:
        mn = dim1*dim2;
        ise_out = mn.*(vert_entropy + hor_cond_entropy);
        ise_out = ise_out - (mn/2).*(pos_angled_entropy + neg_angled_entropy);
        
        ISE= [ISE ; ise_out];
    end

    ise_out=ISE;
    disp(['time taken to calculate for ise3_1: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
    entropy= sum(entropy);%output sum of all entropy
    
end
