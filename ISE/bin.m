close all
clear
clc
load('vmsv.mat');
load('vmpc.mat');

bin_resolution=[0.05;0.01;0.005;0.001;0.0001];
t=zeros(12,4);
n=0;
for i=1:size(bin_resolution,1)
    map=vmp.data.maps_adsm;
    mi = min(map);
    r=range(map);
    ma = max(map);
    s=std(map);
    n=n+1;
    t(n,1) = mi;
    t(n,2) = r;
    t(n,3) = ma;
    t(n,4) = s;
    
    actual_disc = floor(vmp.data.maps_adsm/bin_resolution(i))+1;
    figure('Name', 'Actual Image','NumberTitle','off');
    plotmap(actual_disc,'place')
    actual_disc(isnan(actual_disc)) = 0;
    actual_disc(actual_disc==0)=[];
    mi = min(actual_disc);
    m=mod(mi,bin_resolution(i));
    r=range(actual_disc);
    ma = max(actual_disc);
    s=std(actual_disc);
    n=n+1;
    t(n,1) = mi;
    t(n,2) = r;
    t(n,3) = ma;
    t(n,4) = s;
    
figure('Name', num2str(bin_resolution(i)),'NumberTitle','off');
h=histogram(actual_disc,mi:1:ma) 
title(num2str(bin_resolution(i)));


%shuffle image
im=vmp.data.maps_adsmsh(1,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution(i))+1;
figure('Name', ['Shuffle Image 1',num2str(bin_resolution(i))],'NumberTitle','off');
plotmap(actual_disc,'place') 
im=vmp.data.maps_adsmsh(2,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution(i))+1;
figure('Name', 'Shuffle Image 2','NumberTitle','off');
plotmap(actual_disc,'place') 
im=vmp.data.maps_adsmsh(10,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution(i))+1;
figure('Name', 'Shuffle Image 10','NumberTitle','off');
plotmap(actual_disc,'place') 

%vms
map=vms.data.maps_adsm;
    mi = min(map);
    r=range(map);
    ma = max(map);
    s=std(map);
    n=n+1;
    t(n,1) = mi;
    t(n,2) = r;
    t(n,3) = ma;
    t(n,4) = s;
    
    actual_disc = floor(vms.data.maps_adsm/bin_resolution(i))+1;
    figure('Name', 'Actual Image','NumberTitle','off');
    plotmap(actual_disc,'spatialview')
    actual_disc(isnan(actual_disc)) = 0;
    actual_disc(actual_disc==0)=[];
    mi = min(actual_disc);
    m=mod(mi,bin_resolution(i));
    r=range(actual_disc);
    ma = max(actual_disc);
    s=std(actual_disc);
    n=n+1;
    t(n,1) = mi;
    t(n,2) = r;
    t(n,3) = ma;
    t(n,4) = s;
    
figure('Name', num2str(bin_resolution(i)),'NumberTitle','off');
h=histogram(actual_disc,mi:1:ma) 
title(num2str(bin_resolution(i)));


%shuffle image
im=vms.data.maps_adsmsh(1,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution(i))+1;
figure('Name', 'Shuffle Image 1','NumberTitle','off');
plotmap(actual_disc,'spatialview') 
im=vms.data.maps_adsmsh(2,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution(i))+1;
figure('Name', 'Shuffle Image 2','NumberTitle','off');
plotmap(actual_disc,'spatialview') 
im=vms.data.maps_adsmsh(10,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution(i))+1;
figure('Name', 'Shuffle Image 10','NumberTitle','off');
plotmap(actual_disc,'spatialview') 

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