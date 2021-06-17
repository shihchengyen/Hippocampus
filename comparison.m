Method = {'ise0.05';'ise0.005';'ise_change';'ise_mat';'ise_ted'};
sv1=load('vmsvz_ise0.05.mat');
sv2=load('vmsvz_ise0.005.mat');
sv3=load('vmsvz_isechange.mat');
sv4=load('vmsvz_mat.mat');
sv5=load('vmsvz_ted.mat');

ISE =[sv1.vms.data.ISE; sv2.vms.data.ISE;  sv3.vms.data.ISE; sv4.vms.data.ISE; sv5.vms.data.ISE];
ISEsh =[mean(sv1.vms.data.ISEsh); mean(sv2.vms.data.ISEsh);  mean(sv3.vms.data.ISEsh); mean(sv4.vms.data.ISEsh); mean(sv5.vms.data.ISEsh)];
ISE_ISEsh =ISE./ISEsh;
SIC =[sv1.vms.data.SIC; sv2.vms.data.SIC;  sv3.vms.data.SIC; sv4.vms.data.SIC; sv5.vms.data.SIC];
SICsh =[mean(sv1.vms.data.SICsh); mean(sv2.vms.data.SICsh);  mean(sv3.vms.data.SICsh); mean(sv4.vms.data.SICsh); mean(sv5.vms.data.SICsh)];

% Create a table, T, as a container for the workspace variables. The table function uses the workspace variable names as the names of the table variables in T. A table variable can have multiple columns. For example, the BloodPressure variable in T is a 5-by-2 array.

T = table(Method,ISE,ISEsh,ISE_ISEsh,SIC,SICsh)