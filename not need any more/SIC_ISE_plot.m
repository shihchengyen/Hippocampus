% clear
close all
load('vmpc.mat')

%SIC
sic = vmp.data.SIC;
sic_97 = prctile(vmp.data.SICsh, 97.5);
sic_z = (sic-mean(vmp.data.SICsh))/std(vmp.data.SICsh);
%ISE
ise1 = vmp.data.ISE;
ise_2 = prctile(vmp.data.ISEsh, 2.5);
ise_z = (ise1-mean(vmp.data.ISEsh))/std(vmp.data.ISEsh);

results = [sic, sic_97,sic_z, ise1, ise_2, ise_z];  
%map
actual = vmp.data.maps_adsm;
figure; plotmap(actual,'place');view(2);
l=isnan(actual);
v=~l;
map=NaN(1,1600);
map(v) = vmp.data.maps_adsmsh(1,v);
figure; plotmap(map,'place');view(2);
map(v) = vmp.data.maps_adsmsh(2,v);
figure; plotmap(map,'place');view(2);
map(v)=vmp.data.maps_adsmsh(3,v);
figure; plotmap(map,'place');view(2);
map(v)=vmp.data.maps_adsmsh(4,v);
figure; plotmap(map,'place');view(2);
map(v)=vmp.data.maps_adsmsh(10,v);
figure; plotmap(map,'place');view(2);
% figure; plotmap(vmp.data.maps_adsmsh(1,:),'place');view(2);
% figure; plotmap(vmp.data.maps_adsmsh(2,:),'place');view(2);
% figure; plotmap(vmp.data.maps_adsmsh(3,:),'place');view(2);
% figure; plotmap(vmp.data.maps_adsmsh(4,:),'place');view(2);
% figure; plotmap(vmp.data.maps_adsmsh(10,:),'place');view(2);

% %vmsv
% load('vmsv.mat')
% %SIC
% sic = vms.data.SIC;
% sic_97 = prctile(vms.data.SICsh, 97.5);
% sic_z = (sic-mean(vms.data.SICsh))/std(vms.data.SICsh);
% %ISE
% ise = vmp.data.ISE;
% ise_2 = prctile(vms.data.ISEsh, 2.5);
% ise_z = (ise-mean(vms.data.ISEsh))/std(vms.data.ISEsh);
% 
% results2 = [sic, sic_97,sic_z, ise, ise_2, ise_z];  
% %map
% figure; plotmap(vms.data.maps_adsm,'spatialview');
% figure; plotmap(vms.data.maps_adsmsh(1,:),'spatialview');
% figure; plotmap(vms.data.maps_adsmsh(2,:),'spatialview');
% figure; plotmap(vms.data.maps_adsmsh(3,:),'spatialview');
% figure; plotmap(vms.data.maps_adsmsh(4,:),'spatialview');
% figure; plotmap(vms.data.maps_adsmsh(10,:),'spatialview');