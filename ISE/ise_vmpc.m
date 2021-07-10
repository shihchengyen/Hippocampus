clear all
load('vmpc.mat');

%for ise_1_1.mat
origin = vmp.data.origin;
ise_out = ise1_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40);
ise = ise_out(1);
ise_sh = ise_out(2:end);
ise_sh(ise_sh==0) = []; %exclude zero ise
ise_2_5 = prctile(ise_sh, 2.5);
z =(ise-mean(ise_sh))/std(ise_sh);
%save new .mat file
save('ise_1_1.mat','origin','ise','ise_2_5','z'); 

%append ise_1_1.mat
load('ise_1_1.mat')
Origin = vmp.data.origin;
origin=[origin;Origin]; %file directory
ise_out = ise1_1(vmp.data.maps_adsm, vmp.data.maps_adsmsh,40,40); %ISE calculation
ISE = ise_out(1); 
ise = [ise; ISE];
ise_sh = ise_out(2:end); %ISE shuffle
ise_sh(ise_sh==0) = []; %exclude zero ise
ISE_2_5 = prctile(ise_sh, 2.5); %ISE 2.5 percentile
ise_2_5 = [ise_2_5;ISE_2_5];
Z =(ise-mean(ise_sh))/std(ise_sh); %ISE z score from 
z = [z;Z];
save('ise_1_1.mat', 'origin', 'ise','ise_2_5','z','-append'); 