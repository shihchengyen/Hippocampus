close all
clear
clc
load('vmsv.mat');
load('vmpc.mat');

bin_resolution=[0.05;0.01;0.005];
t=zeros(2,4);
t2=zeros(2,4);
for i=1:size(bin_resolution,1)
    map=vmp.data.maps_adsm;
    mi = min(map);
    r=range(map);
    ma = max(map);
    s=std(map);
    t(1,1) = mi;
    t(1,2) = r;
    t(1,3) = ma;
    t(1,4) = s;
    
    actual_disc = floor(vmp.data.maps_adsm/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
    figure('Name', 'Actual Image','NumberTitle','off');
    plotmap(actual_disc,'place')
    actual_disc(actual_disc==0)=[];
    mi = min(actual_disc);
    m=mod(mi,bin_resolution(i));
    r=range(actual_disc);
    ma = max(actual_disc);
    s=std(actual_disc);
    t(2,1) = mi;
    t(2,2) = r;
    t(2,3) = ma;
    t(2,4) = s;
    
figure('Name', num2str(bin_resolution(i)),'NumberTitle','off');
h=histogram(actual_disc,mi:1:ma) 
title(num2str(bin_resolution(i)));


%shuffle image
actual_disc = floor(vmp.data.maps_adsmsh(1,:)/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
figure('Name', 'Shuffle Image 1','NumberTitle','off');
plotmap(actual_disc,'place') 
actual_disc = floor(vmp.data.maps_adsmsh(2,:)/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
    
figure('Name', 'Shuffle Image 2','NumberTitle','off');
plotmap(actual_disc,'place') 
actual_disc = floor(vmp.data.maps_adsmsh(10,:)/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
    
figure('Name', 'Shuffle Image 10','NumberTitle','off');
plotmap(actual_disc,'place') 

%vms
map=vms.data.maps_adsm;
    mi = min(map);
    r=range(map);
    ma = max(map);
    s=std(map);
    t(1,1) = mi;
    t(1,2) = r;
    t(1,3) = ma;
    t(1,4) = s;
    
    actual_disc = floor(vms.data.maps_adsm/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
    figure('Name', 'Actual Image','NumberTitle','off');
    plotmap(actual_disc,'spatialview')
    actual_disc(actual_disc==0)=[];
    mi = min(actual_disc);
    m=mod(mi,bin_resolution(i));
    r=range(actual_disc);
    ma = max(actual_disc);
    s=std(actual_disc);
    t(2,1) = mi;
    t(2,2) = r;
    t(2,3) = ma;
    t(2,4) = s;
    
figure('Name', num2str(bin_resolution(i)),'NumberTitle','off');
h=histogram(actual_disc,mi:1:ma) 
title(num2str(bin_resolution(i)));


%shuffle image
actual_disc = floor(vms.data.maps_adsmsh(1,:)/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
   
figure('Name', 'Shuffle Image 1','NumberTitle','off');
plotmap(actual_disc,'spatialview') 
actual_disc = floor(vms.data.maps_adsmsh(2,:)/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
   
figure('Name', 'Shuffle Image 2','NumberTitle','off');
plotmap(actual_disc,'spatialview') 
actual_disc = floor(vms.data.maps_adsmsh(10,:)/bin_resolution(i))+1;
    actual_disc(isnan(actual_disc)) = 0;
    
figure('Name', 'Shuffle Image 10','NumberTitle','off');
plotmap(actual_disc,'spatialview') 
end