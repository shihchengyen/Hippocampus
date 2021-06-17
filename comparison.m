Method = {'ise0.01';'ise0.05';'ise0.005';'ise_change';'ise_mat';'ise_mat0.05';'ise_ted'};
sv0=load('vmsv_ise0.01.mat');
sv1=load('vmsv_ise0.05.mat');
sv2=load('vmsv_ise0.005.mat');
sv3=load('vmsv_isechange.mat');
sv4=load('vmsv_mat.mat');
sv6=load('vmsv_mat0.05.mat');
sv5=load('vmsv_ted.mat');

ISE =[sv0.vs.data.ISE; sv1.vs.data.ISE; sv2.vs.data.ISE;  sv3.sv.data.ISE; sv4.vs.data.ISE; sv6.vs.data.ISE; sv5.vms.data.ISE];
ISEsh =[mean(sv0.vs.data.ISEsh);mean(sv1.vs.data.ISEsh); mean(sv2.vs.data.ISEsh);  mean(sv3.sv.data.ISEsh); mean(sv4.vs.data.ISEsh); mean(sv6.vs.data.ISEsh);mean(sv5.vms.data.ISEsh)];
ISE_ISEsh =ISE./ISEsh;
SIC =[sv0.vs.data.SIC; sv1.vs.data.SIC; sv2.vs.data.SIC;  sv3.sv.data.SIC; sv4.vs.data.SIC; sv6.vs.data.SIC; sv5.vms.data.SIC];
SICsh =[mean(sv0.vs.data.SICsh); mean(sv1.vs.data.SICsh); mean(sv2.vs.data.SICsh);  mean(sv3.sv.data.SICsh); mean(sv4.vs.data.SICsh); mean(sv6.vs.data.SICsh); mean(sv5.vms.data.SICsh)];

% Create a table, T, as a container for the workspace variables. The table function uses the workspace variable names as the names of the table variables in T. A table variable can have multiple columns. For example, the BloodPressure variable in T is a 5-by-2 array.

T = table(Method,ISE,ISEsh,ISE_ISEsh,SIC,SICsh)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%change bin resolution from 0.5; 0.05; 0.005; 0.0005; 0.1; 0.01; 0.001; 0.0001; 1; 2; 4; 5; 10; 15; 20
bin_resolution=[0.5; 0.05; 0.005; 0.0005; 0.1; 0.01; 0.001; 0.0001; 1; 2; 4; 5; 10; 15; 20];

% for i=1:size(bin_resolution,1)
% [ise_out] = ise(actual_image, shuffled_images, dim1, dim2,bin_resolution(i))
% vs=vmsv('auto');
% name=['ise' num2str(bin_resolution(i))'.mat'];
% save(name,'vs');
% end

% comparison
Method = bin_resolution;
ISE = [];
ISEsh = [];
SIC = [];
SICsh = [];

for i=1:size(bin_resolution,1)
name=['ise' num2str(bin_resolution(i)) '.mat'];
sv=load(name);

ISE =[ISE; sv.vs.data.ISE];
ISEsh =[ISEsh; mean(sv.vs.data.ISEsh)];
SIC =[SIC; sv.vs.data.SIC];
SICsh =[SICsh; mean(sv.vs.data.SICsh)];

end
ISE_ISEsh =ISE./ISEsh;

% Create a table, T, as a container for the workspace variables. The table function uses the workspace variable names as the names of the table variables in T. A table variable can have multiple columns. For example, the BloodPressure variable in T is a 5-by-2 array.

T = table(Method,ISE,ISEsh,ISE_ISEsh,SIC,SICsh)