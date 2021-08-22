%Image spatial entropy without using QMRF.
%Actual_image =actual map from data collected.
%Shuffled_images= randomly shuffle map base on the actual map.
%dim1, dim2 are the dimension of the input map. 

%bin_resolution is used to discretize the map intensity to multiple bin.
%Some are too big or too small for the cell's spike aplitude.
%When bin resolution is not 0, divide parameter won't be consider.

%When bin resolution is 0,divide is used to set how many bin you want to
%seperate the map intensity into.
%Thresh is used to calculating only parts of the map that is above the
%threshold.

%If the thresh=0, than change 1 and 0 have any effect. 
%since change 1 replace anything above the thershold with one.
%but it should change to maximum value pluse one. So it won't corrupt the
%original data.(I didn't change it yet).
%since change 0 replace anything below the thershold with zero.
%but it should change to minimum value minus one. So it won't corrupt the
%original data.(I didn't change it yet).

%pixel 1 means comparing value directly next to each other.
%pixel 2 means comparing value one across each other.

%Acceptable input value range.
%bin_resolution: 0~1
%divide>0
%thresh:0~99
%change: 2,1,0
%pixel: 1,2
function [ise_out] = ise_threshold(actual_image, shuffled_images, dim1, dim2,bin_resolution, divide, thresh, change, pixel)
    tic; %show how long ise caculation will take
    
    combined = [actual_image; shuffled_images]; %(1+shuffle) x dim1*dim2
    
    %Image spatial entropy calculation
    ISE=zeros(size(combined,1),1);
    for  i=1:size(combined,1)
        if max(combined(i,:))==1 %an image being discrtized by not small enough bin_reolution is excluded from ISE calculation
            ISE(i) = 0;
            disp('floor(combined(i,:))=0');%every intensity value rounded to zero
            continue
        elseif sum(sum(~isnan(combined(i,:))))==0 %an empty image don't have ISE
            ISE(i) = 0;
            disp('an image consist of NaN');%every intensity value rounded to zero
            continue
        end
        image = reshape(combined(i,:), dim1, dim2); %image in the form of its pixel intensity
        
        if thresh==0 %no threshold
            % parameters to discretize maps
            ma = max(max(image)); %maximum intensity
            mi = min(min(image)); %minimum intensity 
            if bin_resolution==0 %set to changeable bin_resolution
                bin_resolution = (ma-mi)/divide;
            end
            %else bin_resolution is the same as the input
            
            % binning each datapoint
            temp = floor(image/bin_resolution)+1;
        else %define threshold
            
            threshold = prctile(combined(i,:),thresh);
            %pad a layer of NaN
            temp = [NaN(dim1,1) image NaN(dim1,1)];
            padded = [NaN(1,dim2+2); temp ;NaN(1,dim2+2)];
            dimr = size(padded,1);
            dimc = size(padded,2);
            pass = padded>threshold;

            if sum(sum(pass))<2 
                ISE(i) = 0;
                disp('Only one intensity pass the threshold,image does not have enough data.');%every intensity value rounded to zero
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

            if change==0 %doesn't pass -->0 (consider background
                process = zeros(dimr, dimc);
                process(pass) = padded(pass);
            elseif change==1 %pass-->1 (exclude background
                process = NaN(dimr, dimc);
                process(bin_sieved) = padded(bin_sieved);
                o = ones(dimr, dimc);
                process(pass) = o(pass);
            elseif change==2 %no change
                process = NaN(dimr, dimc);
                process(bin_sieved) = padded(bin_sieved);
            else
                disp('Missing whether or not the user want changes to the intensity or not')
            end
            % parameters to discretize maps
            ma = max(max(process)); %maximum intensity
            mi = min(min(process)); %minimum intensity 
            if bin_resolution==0 %changeable bin_resolution
                    bin_resolution = (ma-mi)/divide;
            end
                %else bin_resolution is the same as the input

            % binning each datapoint
            temp = floor(process/bin_resolution)+1;
        end
                    
%         %For examination
%         figure('Name','binned intensity histogram','NumberTitle','off');
%         b=histogram(temp)
% %         b.BinWidth= 1;
% %         b.BinEdges= min(min(temp)):1:max(max(temp));
        % H(X,Xu) computations   
        if pixel == 1
            upper =temp(1:end-1,:); %Xu
            centre = temp(2:end,:); %X
        elseif pixel == 2
            upper =temp(1:end-2,:); %Xu
            centre = temp(3:end,:); %X
        end
        d=[reshape(centre,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
%         %For examination
%         figure('Name','X,Xu joint histogram','NumberTitle','off');
% %         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
%         hist3(d,'CdataMode', 'auto','FaceColor','interp')
%         colorbar
%         view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X_Xu=h/total;%convert histgram count to probability
        %The joint probability of every intensity
        j_X_Xu=reshape(j_X_Xu,[],1);
        j_X_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_X_Xu(j_X_Xu==0)=[]; %remove probability of zero
        [vert_entropy] = sum(entropy(j_X_Xu));
        
        %   H(Xl,Xu) computations
        if pixel == 1
            left =temp(2:end,1:end-1); %Xl
            upper = temp(1:end-1,2:end); %Xu
        elseif pixel == 2
            left =temp(3:end,1:end-2); %Xl
            upper = temp(1:end-2,3:end); %Xu
        end
        
        d=[reshape(left,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
%         %For examination
%         figure('Name','Xl,Xu joint histogram','NumberTitle','off');
% %         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
%         hist3(d,'CdataMode', 'auto','FaceColor','interp')
%         colorbar
%         view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_X1_Xu=h/total; %convert histgram count to probability
        %The joint probability of every intensity
        j_X1_Xu=reshape(j_X1_Xu,[],1);
        j_X1_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_X1_Xu(j_X1_Xu==0)=[]; %remove probability of zero
        [pos_angled_entropy] = sum(entropy(j_X1_Xu));
        
        %   H(Xr,Xu) computations
        if pixel == 1
            right =temp(2:end,2:end); %Xr
            upper = temp(1:end-1,1:end-1); %Xu
        elseif pixel == 2
            right =temp(3:end,3:end); %Xr
            upper = temp(1:end-2,1:end-2); %Xu
        end
        
        d=[reshape(right,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
%         %For examination
%         figure('Name','Xr,Xu joint histogram','NumberTitle','off');
% %         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
%         hist3(d,'CdataMode', 'auto','FaceColor','interp')
%         colorbar
%         view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_Xu=h/total; %convert histgram count to probability
        %The joint probability of every intensity
        j_Xr_Xu=reshape(j_Xr_Xu,[],1);
        j_Xr_Xu(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_Xr_Xu(j_Xr_Xu==0)=[]; %remove probability of zero
        [neg_angled_entropy] = sum(entropy(j_Xr_Xu));  
        
        %   H(Xr,X) computations
        if pixel == 1
            right =temp(:,2:end); %Xr
            centre = temp(:,1:end-1); %X
        elseif pixel == 2
            right =temp(:,3:end); %Xr
            centre = temp(:,1:end-2); %X
        end
        
        d=[reshape(right,[],1),reshape(centre,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[];
%         %For examination
%         figure('Name','Xr,X joint histogram','NumberTitle','off');
% %         hist3(d,{min(min(temp)):1:max(max(temp)) min(min(temp)):1:max(max(temp))})
%         hist3(d,'CdataMode', 'auto','FaceColor','interp')
%         colorbar
%         view(2) %what's view 3 for 3D
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total; %convert histgram count to probability
        %The joint probability of every intensity
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_Xr_X(j_Xr_X==0)=[]; %remove probability of zero
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
    disp(['time taken to calculate for ise_threshold: ' num2str(toc)]);
end

function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
    
end
