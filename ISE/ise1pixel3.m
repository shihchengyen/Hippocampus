function [ise_out] = ise1pixel3(actual_image, shuffled_images, dim1, dim2)
    % one pixel neighbor
    tic; %show how long ise caculation will take
      
    combined = [actual_image; shuffled_images]; %(1+shuffle) x dim1*dim2
    
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
            disp('floor(combined(i,:))=0');%every intensity value rounded to zero
            continue
        end
        image = reshape(combined(i,:), dim1, dim2); %image in the form of its pixel intensity
        
        %define threshold
        threshold = prctile(combined(i,:),50);
       %pad a layer of NaN
        temp = [NaN(dim1,1) image NaN(dim1,1)];
        padded = [NaN(1,dim2+2); temp ;NaN(1,dim2+2)];
        dimr = size(padded,1);
        dimc = size(padded,2);
        pass = padded>threshold;
        if sum(sum(pass))<2 
            ISE(i) = 0;
            disp('floor(combined(i,:))=0');%every intensity value rounded to zero
            continue
        end
        %find the edge
        hor = pass + [pass(2:end,:); zeros(1,dimc)] +[zeros(1,dimc); pass(1:end-1,:) ];
        ver = pass + [pass(:,2:end) zeros(dimr,1)] +[zeros(dimr,1) pass(:,1:end-1) ];
        pos_ang = pass + [[zeros(dimr-1,1) pass(2:end,1:end-1)];zeros(1,dimc)]+[zeros(1,dimc);[pass(1:end-1,2:end) zeros(dimr-1,1)]];
        neg_ang = pass + [[pass(2:end,2:end) zeros(dimr-1,1)];zeros(1,dimc)]+[zeros(1,dimc);[zeros(dimr-1,1) pass(1:end-1,1:end-1)]];
        all= hor + ver + pos_ang+ neg_ang;
        %bin sieved
        bin_sieved = all>0;
        %bin sieved image
        process = NaN(dimr, dimc);
        process(bin_sieved) = padded(bin_sieved);
        ma = max(max(process)); %maximum intensity
        mi = min(min(process)); %minimum intensity 

        % parameters to discretize maps
        bin_resolution = (ma-mi)/15;
    
    % binning each datapoint
        temp = floor(process/bin_resolution)+1;
    % H(X,Xu) computations   
        upper =temp(1:end-1,:); %Xu
        centre = temp(2:end,:); %X
        d=[reshape(centre,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[]; %get rid of any pair of NaN
%         figure('Name','X,Xu joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total;%convert histgram count to probability
        
        %The joint probability of every intensity
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_X_Xu(j_X_Xu==0)=[]; %remove probability of zero
        [vert_entropy] = sum(entropy(j_X_Xu));
     
% H(Xl,Xu) computations   
        left =temp(2:end,1:end-1); %Xl
        upper = temp(1:end-1,2:end); %Xu
        d=[reshape(left,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[]; %get rid of any pair of NaN
%         figure('Name','Xl,Xu joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X1_Xu=h/total;%convert histgram count to probability
        
        %The joint probability of every intensity
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        j_X1_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_X1_Xu(j_X1_Xu==0)=[]; %remove probability of zero
        [pos_angled_entropy] = sum(entropy(j_X1_Xu));
        
% H(Xr,Xu) computations   
        right =temp(2:end,2:end); %Xr
        upper = temp(1:end-1,1:end-1); %Xu
        d=[reshape(right,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[]; %get rid of any pair of NaN
%         figure('Name','Xr,Xu joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_Xu=h/total;%convert histgram count to probability
        
        %The joint probability of every intensity
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        j_Xr_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_Xr_Xu(j_Xr_Xu==0)=[]; %remove probability of zero
        [neg_angled_entropy] = sum(entropy(j_Xr_Xu));
        
% H(Xr,X) computations   
        right =temp(:,2:end); %Xr
        centre = temp(:,1:end-1); %X
        d=[reshape(right,[],1),reshape(centre,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[]; %get rid of any pair of NaN
%         figure('Name','Xr,X joint histogram','NumberTitle','off');
%         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total;%convert histgram count to probability
        
        %The joint probability of every intensity
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_Xr_X(j_Xr_X==0)=[]; %remove probability of zero
        [hor_entropy] = sum(entropy(j_Xr_X));
        
%         results = zeros(4,1);
%         results(1) = vert_entropy;
%         results(2) = pos_angled_entropy;
%         results(3) = neg_angled_entropy;
%         results(4) = hor_entropy;      

        %   final calculation:
        mn = dim1*dim2;
        %The joint probability of every intensity
        ise_out1 = (vert_entropy + hor_entropy +pos_angled_entropy +neg_angled_entropy);
        ise_out1 = ise_out1./mn;
        
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
