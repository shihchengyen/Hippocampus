function [vert_entropy, pos_angled_entropy,neg_angled_entropy, hor_entropy, self_entropy] = ise9(actual_image, shuffled_images, dim1, dim2)
    % Use QMRF
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
        image = reshape(image, dim1, dim2); %image in the form of its pixel intensity probability
        temp = reshape(combined(i,:), dim1, dim2); %image in the form of its pixel intensity
        %image with two layer of zero border
        temp = [zeros(2,dim2); temp ; zeros(2,dim2)]; %add zeros above and below
        temp = [zeros(dim1+4,2) temp zeros(dim1+4,2)]; %add zeros the side
        
        % H(X,Xu) computations   
        upper =temp(1:end-1,:); %Xu
        centre = temp(2:end,:); %X
        h=hist3([reshape(centre,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total;%convert histgram count to probability
        
        d1 = size(j_X_Xu,1);
        d2 = size(j_X_Xu,2);
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(j_X_Xu==0)=NaN;
        [vert_entropy] = entropy(j_X_Xu);
        j_X_Xu(isnan(j_X_Xu))=[];
         vert_entropy = reshape(vert_entropy,d1,d2);
         vert_entropy = vert_entropy.*h;
         vert_entropy(isnan(vert_entropy)) = [];
         vert_entropy(1) = [];
         ver=vert_entropy;
         vert_entropy = sum(vert_entropy);
         
        
        %   H(Xl,Xu) computations
        left =temp(2:end,1:end-1); %Xl
        upper = temp(1:end-1,2:end); %Xu
        h=hist3([reshape(left,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X1_Xu=h/total; %convert histgram count to probability
        
        d1 = size(j_X1_Xu,1);
        d2 = size(j_X1_Xu,2);
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        j_X1_Xu(j_X1_Xu==0)=NaN;
        [pos_angled_entropy] = entropy(j_X1_Xu);
        j_X1_Xu(isnan(j_X1_Xu))=[];
         pos_angled_entropy = reshape(pos_angled_entropy,d1,d2);
         pos_angled_entropy = pos_angled_entropy.*h;
         pos_angled_entropy(isnan(pos_angled_entropy)) = [];
         pos_angled_entropy(1) = [];
         pos=pos_angled_entropy;
         pos_angled_entropy = sum(pos_angled_entropy);
         
        
        %   H(Xr,Xu) computations
        right =temp(2:end,2:end); %Xr
        upper = temp(1:end-1,1:end-1); %Xu
        h=hist3([reshape(right,[],1),reshape(upper,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_Xu=h/total; %convert histgram count to probability
        
        d1 = size(j_Xr_Xu,1);
        d2 = size(j_Xr_Xu,2);
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        j_Xr_Xu(j_Xr_Xu==0)=NaN;
        [neg_angled_entropy] = entropy(j_Xr_Xu);
        j_Xr_Xu(isnan(j_Xr_Xu))=[];
         neg_angled_entropy = reshape(neg_angled_entropy,d1,d2);
         neg_angled_entropy = neg_angled_entropy.*h;
         neg_angled_entropy(isnan(neg_angled_entropy)) = [];
         neg_angled_entropy(1) = [];
         neg=neg_angled_entropy;
         neg_angled_entropy = sum(neg_angled_entropy);
          
        
        %   H(Xr/X) computations = H(Xr,X) - H(X)
        %   H(Xr,X)
        right =temp(:,2:end); %Xr
        centre = temp(:,1:end-1); %X
        h=hist3([reshape(right,[],1),reshape(centre,[],1)],{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total; %convert histgram count to probability
        
        d1 = size(j_Xr_X,1);
        d2 = size(j_Xr_X,2);
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(j_Xr_X==0)=NaN;
        [hor_entropy] = entropy(j_Xr_X);
        j_Xr_X(isnan(j_Xr_X))=[];
         hor_entropy = reshape(hor_entropy,d1,d2);
         hor_entropy = hor_entropy.*h;
         hor_entropy(isnan(hor_entropy)) = [];
         hor_entropy(1) = [];
         hor=hor_entropy;
         hor_entropy = sum(hor_entropy);
         
        
        %   H(X)
        centre = image; %X
        p_X= reshape(centre,[],1); 
        p_X(p_X==0) = [];
        [self_entropy] = sum(entropy(p_X));
        %H(Xr,X) - H(X)
        hor_cond_entropy = hor_entropy - self_entropy;
        
        %   final calculation:
        mn = dim1*dim2;
        ise_out = mn.*(vert_entropy + hor_cond_entropy);
        ise_out = ise_out - (mn/2).*(pos_angled_entropy + neg_angled_entropy);
        ise_out = ise_out./mn;
        
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
%     entropy= sum(entropy);%output sum of all entropy
    
end
