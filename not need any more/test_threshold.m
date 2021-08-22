close all
clear 
clc
%generate random 
image = rand(40,40);
%spike
% % %ones(1,4)
% image(87)=1;
% image(127)=1;
% image(167)=1;
% image(207)=1;
% %ones(4,1)
% image(87:90)=ones(4,1);
%%negative angle
% image(87)=1;
% image(128)=1;
% image(169)=1;
% image(210)=1;
%%positive angle
% image(90)=1;
% image(129)=1;
% image(168)=1;
% image(207)=1;
%%ones(2)
% image(4:5,4:5) = ones(2);
% %%ones(3)
% image(4:6,4:6) = ones(3);
% %%ones(4)
% image(4:7,4:7) = ones(4);
% %%ones(5)
% image(4:8,4:8) = ones(5);
%ones(6)
image(4:9,4:9) = ones(6);
% %%ones(7)
% image(4:10,4:10) = ones(7);

%histogram of intensity distribution
figure('Name','original','NumberTitle','off');
h=histogram(image)
dim1 = size(image,1);
dim2 = size(image,2);
%pad a layer of NaN
temp = [NaN(dim1,1) image NaN(dim1,1)];
padded = [NaN(1,dim2+2); temp ;NaN(1,dim2+2)];
dim1 = size(padded,1);
dim2 = size(padded,2);
%define threshold
threshold = 1;
pass = padded==threshold;
%find the edge
hor = pass + [pass(2:end,:); zeros(1,dim2)] +[zeros(1,dim2); pass(1:end-1,:) ];
ver = pass + [pass(:,2:end) zeros(dim1,1)] +[zeros(dim1,1) pass(:,1:end-1) ];
pos_ang = pass + [[zeros(dim1-1,1) pass(2:end,1:end-1)];zeros(1,dim2)]+[zeros(1,dim2);[pass(1:end-1,2:end) zeros(dim1-1,1)]];
neg_ang = pass + [[pass(2:end,2:end) zeros(dim1-1,1)];zeros(1,dim2)]+[zeros(1,dim2);[zeros(dim1-1,1) pass(1:end-1,1:end-1)]];
all= hor + ver + pos_ang+ neg_ang;
%bin sieved
bin_sieved = all>0;
%bin sieved image
process = NaN(dim1, dim2);
process(bin_sieved) = padded(bin_sieved)
figure('Name','Bin Sieved','NumberTitle','off');
h=histogram(process)


%ISE section
% parameters to discretize maps
    bin_resolution = 0.1;
    
    %histogram BinWidth=bin_resolution
    h.BinWidth=bin_resolution;
    
    % binning each datapoint
    actual_disc = floor(process/bin_resolution)+1;
    temp = actual_disc;
    
% H(X,Xu) computations   
        upper =temp(1:end-1,:); %Xu
        centre = temp(2:end,:); %X
        d=[reshape(centre,[],1),reshape(upper,[],1)];
        d(isnan(d(:,1)),:)=[];
        d(isnan(d(:,2)),:)=[]; %get rid of any pair of NaN
        figure('Name','X,Xu joint histogram','NumberTitle','off');
        hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))})
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
        figure('Name','Xl,Xu joint histogram','NumberTitle','off');
        hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))})
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
        figure('Name','Xr,Xu joint histogram','NumberTitle','off');
        hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))})
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
        figure('Name','Xr,X joint histogram','NumberTitle','off');
        hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))})
        h=hist3(d,{0:1:max(max(temp)) 0:1:max(max(temp))}); %generate joint probability
        total=sum(sum(h))-h(1,1); %total count = count of all intesity - count of NaN
        j_Xr_X=h/total;%convert histgram count to probability
        
        %The joint probability of every intensity
        j_Xr_X=reshape(j_Xr_X,[],1);
        j_Xr_X(1)=[]; %remove probability of zero happen to next to zero which is the first bin of the joint probability output from hist3
        j_Xr_X(j_Xr_X==0)=[]; %remove probability of zero
        [hor_entropy] = sum(entropy(j_Xr_X));
        
        results = zeros(4,1);
        results(1) = vert_entropy;
        results(2) = pos_angled_entropy;
        results(3) = neg_angled_entropy;
        results(4) = hor_entropy;
        
function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
    
end
