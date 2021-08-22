function [ise_out] = ise3pixels(actual_image, shuffled_images, dim1, dim2, bin_resolution)
    % three pixel neighbor
    tic; %show how long ise caculation will take
    
    % parameters to discretize maps
%     bin_resolution = 0.01;
    
    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    
    combined = [actual_disc; shuffled_disc]; %(1+shuffle) x dim1*dim2
    
    %Image spatial entropy calculation
    ISE=zeros(size(combined,1),2);
    for  i=1:size(combined,1)
        if max(combined(i,:))==0 %an empty image don't have ISE
            disp('an image consist of NaN');
            ise_out1 = 0;
            ise_out2 = 0;
            ISE(i,:) = [ ise_out1, ise_out2];
            continue
        end
        %an image being discrtized by not small enough bin_reolution is
        %excluded from ISE calculation
        if max(combined(i,:))==1 
            ise_out1 = 0;
            ise_out2 = 0;
            ISE(i,:) = [ ise_out1, ise_out2];
            disp('floor(combined(i,:))=0');%every intensity value rounded to zero
            continue
        end
        temp = reshape(combined(i,:), dim1, dim2); %image in the form of its pixel intensity
        
        %three pixel
        % H(X,Xu) computations   
        upper =temp(1:end-3,:); %Xu
        centre = temp(4:end,:); %X
        d=[reshape(centre,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total;%convert histgram count to probability
        %The joint probability of every pixel
        d1 = size(j_X_Xu,1);
        d2 = size(j_X_Xu,2);
        j_X_Xu2 = j_X_Xu;
        j_X_Xu2=reshape(j_X_Xu2,[],1);
        j_X_Xu2(j_X_Xu2==0)=NaN;
        [vert_entropy2] = entropy(j_X_Xu2);
        vert_entropy2(isnan(vert_entropy2))=0;
         vert_entropy2 = reshape(vert_entropy2,d1,d2);
         vert_entropy2 = vert_entropy2.*h;
         vert_entropy2(vert_entropy2==0)=[];
         vert_entropy2 = sum(vert_entropy2);
        %The joint probability of every intensity
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(j_X_Xu==0)=[]; %remove probability of zero
        j_X_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [vert_entropy] = sum(entropy(j_X_Xu));
        
        %   H(Xl,Xu) computations
        left =temp(4:end,1:end-3); %Xl
        upper = temp(1:end-3,4:end); %Xu
        d=[reshape(left,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X1_Xu=h/total; %convert histgram count to probability
        %The joint probability of every pixel
        d1 = size(j_X1_Xu,1);
        d2 = size(j_X1_Xu,2);
        j_X1_Xu2 = j_X1_Xu;
        j_X1_Xu2=reshape(j_X1_Xu2,[],1);
        j_X1_Xu2(j_X1_Xu2==0)=NaN;
        [pos_angled_entropy2] = entropy(j_X1_Xu2);
         pos_angled_entropy2(isnan(pos_angled_entropy2))=0;
         pos_angled_entropy2 = reshape(pos_angled_entropy2,d1,d2);
         pos_angled_entropy2 = pos_angled_entropy2.*h;
         pos_angled_entropy2(pos_angled_entropy2==0)=[];
         pos_angled_entropy2 = sum(pos_angled_entropy2);
        %The joint probability of every intensity
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        j_X1_Xu(j_X1_Xu==0)=[]; %remove probability of zero
        j_X1_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [pos_angled_entropy] = sum(entropy(j_X1_Xu));
        
        %   H(Xr,Xu) computations
        right =temp(4:end,4:end); %Xr
        upper = temp(1:end-3,1:end-3); %Xu
        d=[reshape(right,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_Xu=h/total; %convert histgram count to probability
        %The joint probability of every pixel
        d1 = size(j_Xr_Xu,1);
        d2 = size(j_Xr_Xu,2);
        j_Xr_Xu2 = j_Xr_Xu;
        j_Xr_Xu2=reshape(j_Xr_Xu2,[],1);
        j_Xr_Xu2(j_Xr_Xu2==0)=NaN;
        [neg_angled_entropy2] = entropy(j_Xr_Xu2);
        neg_angled_entropy2(isnan(neg_angled_entropy2))=0;
         neg_angled_entropy2 = reshape(neg_angled_entropy2,d1,d2);
         neg_angled_entropy2 = neg_angled_entropy2.*h;
         neg_angled_entropy2(neg_angled_entropy2==0)=[];
         neg_angled_entropy2 = sum(neg_angled_entropy2);
        %The joint probability of every intensity
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        j_Xr_Xu(j_Xr_Xu==0)=[]; %remove probability of zero
        j_Xr_Xu(1)=[];  %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3       
        [neg_angled_entropy] = sum(entropy(j_Xr_Xu));  
        
        %   H(Xr,X) computations 
        right =temp(:,4:end); %Xr
        centre = temp(:,1:end-3); %X
        d=[reshape(right,[],1),reshape(centre,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total; %convert histgram count to probability
        %The joint probability of every pixel
        d1 = size(j_Xr_X,1);
        d2 = size(j_Xr_X,2);
        j_Xr_X2 = j_Xr_X;
        j_Xr_X2=reshape(j_Xr_X2,[],1);
        j_Xr_X2(j_Xr_X2==0)=NaN;
        [hor_entropy2] = entropy(j_Xr_X2);
        hor_entropy2(isnan(hor_entropy2))=0;
         hor_entropy2 = reshape(hor_entropy2,d1,d2);
         hor_entropy2 = hor_entropy2.*h;
         hor_entropy2(hor_entropy2==0)=[];
         hor_entropy2 = sum(hor_entropy2);
        %The joint probability of every intensity
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(j_Xr_X==0)=[]; %remove probability of zero
        j_Xr_X(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [hor_entropy] = sum(entropy(j_Xr_X));
        
        %   final calculation:
        mn = dim1*dim2;
        %The joint probability of every intensity
        ise_out1 = (vert_entropy + hor_entropy +pos_angled_entropy +neg_angled_entropy);
        ise_out1 = ise_out1./mn;
        %The joint probability of every pixel
        ise_out2 = (vert_entropy2 + hor_entropy2 +pos_angled_entropy2 +neg_angled_entropy2);
        ise_out2 = ise_out2./mn;
        
        %The first column use every intensity joint probability.
        %The second column use every pixel joint probability.
        ISE(i,:) = [ ise_out1, ise_out2];
    end

    [ise_out]=ISE;
    disp(['time taken to calculate for ise3pixels: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);

end
