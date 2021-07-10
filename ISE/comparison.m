% clear all
% close all
% clc
% load('vmsv.mat');
% load('vmpc.mat');
% 
% % bin_resolution=[1; 0.5;0.1;  0.05;  0.01;0.005;0.001; 0.0005; 0.0001  ];
% bin_resolution=[0.1;  0.05;  0.01;0.005;0.001; 0.0005];
% for i=1:size(bin_resolution,1)
% %     actual_disc = floor(vms.data.maps_adsm/bin_resolution(i))+1;
% %     actual_disc(isnan(actual_disc)) = 0;
% % %     actual_disc = (actual_disc-mean(actual_disc))/std(actual_disc);
% % %     figure('Name', num2str(bin_resolution(i)),'NumberTitle','off');
% %     subplot(2,3,i)
% % h=histogram(actual_disc) 
% % % h=histogram(actual_disc, 0:(1/bin_resolution(i)/5):max(actual_disc))
% % title(num2str(bin_resolution(i)));
% % %     plotmap(actual_disc,'spatialview')
% 
% actual_disc = floor(vmp.data.maps_adsm/bin_resolution(i))+1;
%     actual_disc(isnan(actual_disc)) = 0;
% % %     actual_disc = (actual_disc-mean(actual_disc))/std(actual_disc);
% %     figure('Name',num2str(bin_resolution(i)),'NumberTitle','off');
% %     plotmap(actual_disc,'place')
% %     lim = caxis;
% %     caxis([lim(2)*2/3 lim(2)])
%     subplot(2,3,i)
% % h=histogram(actual_disc) 
% h=histogram(actual_disc, 0:(1/bin_resolution(i)/20):max(actual_disc))
% % low=(h.BinLimits(2)*8/10)/ h.BinWidth;
% low=(h.BinLimits(2)*0.02)/ h.BinWidth;
% h.BinLimits=[floor(low)*h.BinWidth  h.BinLimits(2)]
% title(num2str(bin_resolution(i)));
% end          
           
%%%%vmp
% load('actual_image.mat');ise_out = ise_ted_2(actual_image, shuffled_images, 51, 161);save('ise_ted.mat','ise_out');
clear
clc
vmp=vmpc('auto');save('vmpc_ted','vmp');load('vmpc_ted.mat'); vmp.data.ISEsh;
% load('vmpc_ted.mat')


% sic_vmp = vmp.data.SIC;
% sic_97 = prctile(vmp.data.SICsh, 97.5);
% z =(vmp.data.SIC-mean(vmp.data.SICsh))/std(vmp.data.SICsh);
% sic1 = vmp.data.SIC1;
% sic2 = vmp.data.SIC2;
% T1 = table(sic_vmp, sic_97, z,sic1, sic2)
% 
ise_vmp = vmp.data.ISE;
ise_2_5 = prctile(vmp.data.ISEsh, 2.5);
z =(vmp.data.ISE-mean(vmp.data.ISEsh))/std(vmp.data.ISEsh);
ise1 = vmp.data.ISE1;
ise2 = vmp.data.ISE2;
T2 = table(ise_vmp, ise_2_5, z, ise1, ise2)
figure
histogram(vmp.data.ISEsh)
% %%%%%vms
% sic_vs = vs.data.SIC;
% sic_95 = prctile(vs.data.SICsh, 95);
% z =(vs.data.SIC-mean(vs.data.SICsh))/std(vs.data.SICsh);
% sic1 = vs.data.SIC1;
% sic2 = vs.data.SIC2;
% T3 = table(sic_vs, sic_95, z,sic1, sic2)
% 
% ise_vs = vs.data.ISE;
% ise_2_5 = prctile(vs.data.ISEsh, 2.5);
% z =(vs.data.ISE-mean(vs.data.ISEsh))/std(vs.data.ISEsh);
% ise1 = vs.data.ISE1;
% ise2 = vs.data.ISE2;
% T4 = table(ise_vs, ise_2_5, z, ise1, ise2)


% %%%%%vms
% sic_vms = vms.data.SIC;
% sic_95 = prctile(vms.data.SICsh, 97.5);
% z =(vms.data.SIC-mean(vms.data.SICsh))/std(vms.data.SICsh);
% sic1 = vms.data.SIC1;
% sic2 = vms.data.SIC2;
% T3 = table(sic_vms, sic_95, z,sic1, sic2)
% 
% ise_vms = vms.data.ISE;
% ise_2_5 = prctile(vms.data.ISEsh, 2.5);
% z =(vms.data.ISE-mean(vms.data.ISEsh))/std(vms.data.ISEsh);
% ise1 = vms.data.ISE1;
% ise2 = vms.data.ISE2;
% T4 = table(ise_vms, ise_2_5, z, ise1, ise2)
% % %map out vmp
% %            figure('Name','maps_raw','NumberTitle','off');
% %            plotmap(vmp.data.maps_raw,'place')
% %            figure('Name','maps_adsm','NumberTitle','off');
% %            plotmap(vmp.data.maps_adsm,'place')
% %            figure('Name','maps_adsm1','NumberTitle','off');
% %            plotmap(vmp.data.maps_adsm1,'place')
% %            figure('Name','maps_adsm2','NumberTitle','off');
% %            plotmap(vmp.data.maps_adsm2,'place')
% %map out vms
%            figure('Name','maps_raw','NumberTitle','off');
%            plotmap(vms.data.maps_raw,'spatialview')
%            figure('Name','maps_adsm','NumberTitle','off');
%            plotmap(vms.data.maps_adsm,'spatialview')
%            figure('Name','maps_adsm1','NumberTitle','off');
%            plotmap(vms.data.maps_adsm1,'spatialview')
%            figure('Name','maps_adsm2','NumberTitle','off');
%            plotmap(vms.data.maps_adsm2,'spatialview')




%            figure('Name','maps_raw','NumberTitle','off');
%            plotmap(vms.data.maps_raw,'spatialview')
%            figure('Name','radii','NumberTitle','off');
%            plotmap(vms.data.radii,'spatialview')
%            figure('Name','dur_adsm','NumberTitle','off');
%            plotmap(vms.data.dur_adsm,'spatialview')
%            figure('Name','maps_adsm','NumberTitle','off');
%            plotmap(vms.data.maps_adsm,'spatialview')
%            figure('Name','maps_raw1','NumberTitle','off');
%            plotmap(vms.data.maps_raw1,'spatialview')
%            figure('Name','maps_adsm1','NumberTitle','off');
%            plotmap(vms.data.maps_adsm1,'spatialview')
%            figure('Name','maps_raw2','NumberTitle','off');
%            plotmap(vms.data.maps_raw2,'spatialview')
%            figure('Name','maps_adsm2','NumberTitle','off');
%            plotmap(vms.data.maps_adsm2,'spatialview')
% % Method = {'ise0.01';'ise0.05';'ise0.005';'ise_change';'ise_mat';'ise_mat0.05';'ise_ted'};
% % sv0=load('vmsv_ise0.01.mat');
% % sv1=load('vmsv_ise0.05.mat');
% % sv2=load('vmsv_ise0.005.mat');
% % sv3=load('vmsv_isechange.mat');
% % sv4=load('vmsv_mat.mat');
% % sv6=load('vmsv_mat0.05.mat');
% % sv5=load('vmsv_ted.mat');
% % 
% % ISE =[sv0.vs.data.ISE; sv1.vs.data.ISE; sv2.vs.data.ISE;  sv3.sv.data.ISE; sv4.vs.data.ISE; sv6.vs.data.ISE; sv5.vms.data.ISE];
% % ISEsh =[mean(sv0.vs.data.ISEsh);mean(sv1.vs.data.ISEsh); mean(sv2.vs.data.ISEsh);  mean(sv3.sv.data.ISEsh); mean(sv4.vs.data.ISEsh); mean(sv6.vs.data.ISEsh);mean(sv5.vms.data.ISEsh)];
% % ISE_ISEsh =ISE./ISEsh;
% % SIC =[sv0.vs.data.SIC; sv1.vs.data.SIC; sv2.vs.data.SIC;  sv3.sv.data.SIC; sv4.vs.data.SIC; sv6.vs.data.SIC; sv5.vms.data.SIC];
% % SICsh =[mean(sv0.vs.data.SICsh); mean(sv1.vs.data.SICsh); mean(sv2.vs.data.SICsh);  mean(sv3.sv.data.SICsh); mean(sv4.vs.data.SICsh); mean(sv6.vs.data.SICsh); mean(sv5.vms.data.SICsh)];
% % 
% % % Create a table, T, as a container for the workspace variables. The table function uses the workspace variable names as the names of the table variables in T. A table variable can have multiple columns. For example, the BloodPressure variable in T is a 5-by-2 array.
% % 
% % T = table(Method,ISE,ISEsh,ISE_ISEsh,SIC,SICsh)