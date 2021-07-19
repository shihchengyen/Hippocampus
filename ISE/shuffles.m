close all
clear
clc
load('vmsv.mat');
% load('vmpc.mat');

bin_resolution=3:15;
% isesh=vmp.data.ISEsh(bin_resolution);
% ise2_5=prctile(vmp.data.ISEsh,2.5);
isesh=vms.data.ISEsh(bin_resolution);
ise2_5=prctile(vms.data.ISEsh,97.5);
figure;h=histogram(vms.data.ISEsh)

for i=1:size(bin_resolution,2)
% %shuffle image vmp
% im=vmp.data.maps_adsmsh(i,:);
% im(im==0)=NaN;
% figure('Name', ['Shuffle Image #',num2str(bin_resolution(i))],'NumberTitle','off');
% plotmap(im,'place') 


%vms
%shuffle image
im=vms.data.maps_adsmsh(i,:);
im(im==0)=NaN;
figure('Name', 'Shuffle Image 1','NumberTitle','off');
plotmap(im,'spatialview') 


end

% %check spiketrain data
% 
% load('vmpc.mat')
% figure('Name', 'Spike Image','NumberTitle','off');
% plotmap(vmp.data.maps_adsm,'place')
% figure
% histogram(vmp.data.maps_adsm)
% h=histogram(vmp.data.maps_adsm)
% h.BinWidth=0.01
% 
% load('vmsv.mat')
% figure('Name', 'Spike Image','NumberTitle','off');
% plotmap(vms.data.maps_adsm,'spatialview')
% figure; h=histogram(vms.data.maps_adsm)
% h.BinWidth=0.01