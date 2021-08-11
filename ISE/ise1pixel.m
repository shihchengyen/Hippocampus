function [ise_out] = ise1pixel(actual_image, shuffled_images, dim1, dim2, bin_resolution)
    % one pixel neighbor
    tic; %show how long ise caculation will take
    
    % parameters to discretize maps
%     bin_resolution = 0.005;
    
    % binning each datapoint
    actual_disc = floor(actual_image/bin_resolution)+1;
    shuffled_disc = floor(shuffled_images/bin_resolution)+1;
    
    combined = [actual_disc; shuffled_disc]; %(1+shuffle) x dim1*dim2
    
    %Image spatial entropy calculation
    ISE=zeros(size(combined,1),1);
    for  i=1:size(combined,1)
        if max(combined(i,:))==0 %an empty image don't have ISE
            disp('an image consist of NaN');
            ISE(i) = 0;
            continue
        end
        %an image being discrtized by not small enough bin_reolution is
        %excluded from ISE calculation
        if max(combined(i,:))==1 
            ISE(i) = 0;
            disp('floor(combined(i,:))=0');%every intensity value rounded to zero
            continue
        end
        if sum(sum(~isnan(combined(i,:))))==0 
            ISE(i) = 0;
            disp('All NaN');%every intensity value rounded to zero
            continue
        end
        temp = reshape(combined(i,:), dim1, dim2); %image in the form of its pixel intensity
        
        % H(X,Xu) computations   
        upper =temp(1:end-1,:); %Xu
        centre = temp(2:end,:); %X
        d=[reshape(centre,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        %For examination
        figure('Name','X,Xu joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        hist3(d,'CdataMode', 'auto','FaceColor','interp')
        colorbar
        view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total;%convert histgram count to probability
        %The joint probability of every intensity
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(j_X_Xu==0)=[]; %remove probability of zero
        j_X_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [vert_entropy] = sum(entropy(j_X_Xu));
        
        %   H(Xl,Xu) computations
        left =temp(2:end,1:end-1); %Xl
        upper = temp(1:end-1,2:end); %Xu
        d=[reshape(left,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        %For examination
        figure('Name','X,Xu joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        hist3(d,'CdataMode', 'auto','FaceColor','interp')
        colorbar
        view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X1_Xu=h/total; %convert histgram count to probability
        %The joint probability of every intensity
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        j_X1_Xu(j_X1_Xu==0)=[]; %remove probability of zero
        j_X1_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        [pos_angled_entropy] = sum(entropy(j_X1_Xu));
        
        %   H(Xr,Xu) computations
        right =temp(2:end,2:end); %Xr
        upper = temp(1:end-1,1:end-1); %Xu
        d=[reshape(right,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        %For examination
        figure('Name','X,Xu joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        hist3(d,'CdataMode', 'auto','FaceColor','interp')
        colorbar
        view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_Xu=h/total; %convert histgram count to probability
        %The joint probability of every intensity
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        j_Xr_Xu(j_Xr_Xu==0)=[]; %remove probability of zero
        j_Xr_Xu(1)=[];  %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3       
        [neg_angled_entropy] = sum(entropy(j_Xr_Xu));  
        
        %   H(Xr,X) computations
        right =temp(:,2:end); %Xr
        centre = temp(:,1:end-1); %X
        d=[reshape(right,[],1),reshape(centre,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
        %For examination
        figure('Name','X,Xu joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        hist3(d,'CdataMode', 'auto','FaceColor','interp')
        colorbar
        view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(combined(i,:)) 0:1:max(combined(i,:))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total; %convert histgram count to probability
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
                
        %The first column use every intensity joint probability.
        ISE(i) = ise_out1;
    end

    [ise_out]=ISE;
    disp(['time taken to calculate for ise1pixel: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
    
end
