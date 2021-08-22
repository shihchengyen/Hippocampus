bin_resolution=[1; 0.5 ; 0.1; 0.05; 0.01; 0.005; 0.001; 0.0005; 0.0001];

Bin = bin_resolution;

Ise = [];
Ise_sh = [];
SIC = [];
SICsh = [];
z = [];
% Std = [];
% pass = [];
% confidence_interval = [];
% min_ISEsh = [];
max_ISEsh = [];
mean_ISEsh = [];
ISE_meannISEsh = [];
% per = [];
per2_5 = [];
% per95 = [];
% for i=1:size(bin_resolution,1)
% name=['ise' num2str(bin_resolution(i)) '.mat'];
% load(name);
% 
% Ise =[Ise; ISE];
% 
% subplot(3,3,i)
% histogram((ISEsh-mean(ISEsh))/std(ISEsh)) 
% title(num2str(bin_resolution(i)));
% 
% % t=std(ISEsh);
% % Std =[Std; t];
% 
% % [h,p,ci,zval] = ztest(ISEsh ,ISE, std(ISEsh));
% zval = (ISE-mean(ISEsh))/std(ISEsh);
% z = [z ; zval];
% % pass = [pass; h]; %1 means not typical, statistically signigicant.
% % confidence_interval = [confidence_interval; ci'];
% % min_ISEsh = [min_ISEsh; min(ISEsh)];
% max_ISEsh = [max_ISEsh; max(ISEsh)];
% mean_ISEsh = [mean_ISEsh; mean(ISEsh)];
% % per = [per; p];
% ISE_meannISEsh = [ISE_meannISEsh; (ISE/mean(ISEsh))];
% per2_5 = [per2_5; prctile(ISEsh, 2.5)];
% % per95 = [per95;prctile(ISEsh, 95)];
% 
% % figure('Name',num2str(bin_resolution(i)),'NumberTitle','off');
% % name_fig = ['ise' num2str(bin_resolution(i)) '.mat'];
% end
load('ise_ted.mat');
Bin = 0.05;
ISE = ise_out(1);
Ise =[Ise; ISE];
ISEsh = ise_out(2:end);
figure
histogram(ISEsh) 
zval = (ISE-mean(ISEsh))/std(ISEsh);
z = [z ; zval];
% pass = [pass; h]; %1 means not typical, statistically signigicant.
% confidence_interval = [confidence_interval; ci'];
% min_ISEsh = [min_ISEsh; min(ISEsh)];
max_ISEsh = [max_ISEsh; max(ISEsh)];
mean_ISEsh = [mean_ISEsh; mean(ISEsh)];
% per = [per; p];
ISE_meannISEsh = [ISE_meannISEsh; (ISE/mean(ISEsh))];
per2_5 = [per2_5; prctile(ISEsh, 2.5)];


% smaller_min = Ise < min_ISEsh;
smaller_per2_5 = Ise < per2_5;
%Create a table, T, as a container for the workspace variables. The table function uses the workspace variable names as the names of the table variables in T. A table variable can have multiple columns. For example, the BloodPressure variable in T is a 5-by-2 array.
T = table(Bin, Ise, smaller_per2_5, per2_5 ,mean_ISEsh, max_ISEsh, z, ISE_meannISEsh )
% T = table(Bin, Ise, Ise_sh, min_ISEsh , max_ISEsh, mean_ISEsh, ISE_meannISEsh )

% multi={};
% multi_r={};
% for i=1:size(bin_resolution,1)
% % binning each datapoint
%     % not rounded actual disc
%     actual_disc = actual_image/bin_resolution(i);
%     actual_disc(isnan(actual_disc)) = 0;
%     actual_disc1 =actual_disc/sum(actual_disc);
%     %rounded actual disc
%     actual_disc = floor(actual_disc)+1;
%     actual_disc2 =actual_disc/sum(actual_disc);
%    
% % %     figure('Name',num2str(bin_resolution(i)),'NumberTitle','off');
% subplot(3,3,i)
% histogram(actual_disc) 
% title(num2str(bin_resolution(i)));
% 
% % histogram(actual_disc2,25);
% % hold on;
% % figure
% % histogram(actual_disc) 
% 
% 
% % %  map it out
% % mapL=actual_disc1;
% % [mapG,mapGdummy]=plotmap(mapL,'spatialview');
% % name1=['ise' num2str(bin_resolution(i)) '.fig'];
% % name1_2=['ise' num2str(bin_resolution(i)) '.png'];
% % savefig(name1);
% % saveas(gcf,name1_2);
% % 
% % mapL=actual_disc2;
% % [mapG,mapGdummy]=plotmap(mapL,'spatialview');
% % name2=['ise' num2str(bin_resolution(i)) '_r.fig'];
% % name2_2=['ise' num2str(bin_resolution(i)) '_r.png'];
% % savefig(name2);
% % saveas(gcf,name2_2);
% % % montage
% % img1 = imread(name1_2);
% % img2 = imread(name2_2);
% % % multi = cat(3, multi, img1);
% % % multi_r = cat(3, multi_r, img2);
% % multi{i} =img1;
% % multi_r{i} =img2;
% 
% end
% 
% % % Display a montage of the images in the multiframe image.
% % figure
% % montage(multi);
% % figure
% % montage(multi_r);
% % montage({imRGB, imGray, 'cameraman.tif'})
% % savefig('SineWave.fig')
% % % Close the figure, then reopen the saved figure using the openfig function.
% % 
% % close(gcf)
% % openfig('SineWave.fig')
% 
% % z*sig+avg=x
% histogram(ISEsh,4.9:0.01:5)
% max(ISEsh)
%     min(ISEsh)
%     histogram(shuffled_images(1,:))
% max(shuffled_images(1,:))
% 
% histogram(shuffled_images(100,:))
% openfig('ise0.1.fig')
% [mapG,mapGdummy] = plotmap(vmp.data.maps_adsm,'place');
% [mapG,mapGdummy] = plotmap(vms.data.maps_adsm,'spatialview');
