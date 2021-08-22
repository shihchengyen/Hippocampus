% shuffle image
close all
clear
clc
load('vmsv.mat');
load('vmpc.mat');

bin_resolution=0.005;

for i=1:10
%shuffle image
im=vmp.data.maps_adsmsh(i,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution)+1;
figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
plotmap(actual_disc,'place') 
% 
% im=vms.data.maps_adsmsh(i,:);
% im(im==0)=NaN;
% actual_disc = floor(im/bin_resolution)+1;
% figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
% plotmap(actual_disc,'spatialview') 
end
v<p;
i=5873;im=vmp.data.maps_adsmsh(i,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution)+1;
figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
plotmap(actual_disc,'place'); i=5874;im=vmp.data.maps_adsmsh(i,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution)+1;
figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
plotmap(actual_disc,'place'); i=5875;im=vmp.data.maps_adsmsh(i,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution)+1;
figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
plotmap(actual_disc,'place'); i=5876;im=vmp.data.maps_adsmsh(i,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution)+1;
figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
plotmap(actual_disc,'place'); i=5877;im=vmp.data.maps_adsmsh(i,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution)+1;
figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
plotmap(actual_disc,'place'); i=5878;im=vmp.data.maps_adsmsh(i,:);
im(im==0)=NaN;
actual_disc = floor(im/bin_resolution)+1;
figure('Name', ['Shuffle Image ',num2str(i)],'NumberTitle','off');
plotmap(actual_disc,'place'); 
vms.data.SICsh(1:10)

ans =

    0.1295
    0.1306
    0.1942
    0.0956
    0.1344
    0.2036
    0.1056
    0.1655
    0.1481
    0.1180

figure;histogram(vms.data.SICsh);figure;histogram(vms.data.ISEsh);
vms.data.ISEsh(1:10)

ans =

    2.9697
    2.9621
    2.9726
    2.9654
    2.9562
    2.9651
    2.9448
    2.9555
    2.9253
    2.9470

prctile(vms.data.ISEsh, 2.5);
p=prctile(vms.data.ISEsh, 2.5);
vms.data.ISEsh;
v=vms.data.ISEsh;