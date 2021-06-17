Method = {'ise0.01';'ise0.05';'ise0.005';'ise_change';'ise_mat';'ise_mat0.05';'ise_ted'};
sv0=load('vmsv_ise0.01.mat');
sv1=load('vmsv_ise0.05.mat');
sv2=load('vmsv_ise0.005.mat');
sv3=load('vmsv_isechange.mat');
sv4=load('vmsv_mat.mat');
sv6=load('vmsv_mat0.05.mat');
sv5=load('vmsv_ted.mat');

ISE =[sv0.vs.data.ISE; sv1.vs.data.ISE; sv2.vs.data.ISE;  sv3.sv.data.ISE; sv4.vs.data.ISE; sv6.vs.data.ISE; sv5.vms.data.ISE];
% ISE =[sv1.data.ISE; sv2.data.ISE;  sv3.data.ISE; sv4.data.ISE; sv5.data.ISE];
ISEsh =[mean(sv0.vs.data.ISEsh);mean(sv1.vs.data.ISEsh); mean(sv2.vs.data.ISEsh);  mean(sv3.sv.data.ISEsh); mean(sv4.vs.data.ISEsh); mean(sv6.vs.data.ISEsh);mean(sv5.vms.data.ISEsh)];
ISE_ISEsh =ISE./ISEsh;
SIC =[sv0.vs.data.SIC; sv1.vs.data.SIC; sv2.vs.data.SIC;  sv3.sv.data.SIC; sv4.vs.data.SIC; sv6.vs.data.SIC; sv5.vms.data.SIC];
SICsh =[mean(sv0.vs.data.SICsh); mean(sv1.vs.data.SICsh); mean(sv2.vs.data.SICsh);  mean(sv3.sv.data.SICsh); mean(sv4.vs.data.SICsh); mean(sv6.vs.data.SICsh); mean(sv5.vms.data.SICsh)];

% Create a table, T, as a container for the workspace variables. The table function uses the workspace variable names as the names of the table variables in T. A table variable can have multiple columns. For example, the BloodPressure variable in T is a 5-by-2 array.

T = table(Method,ISE,ISEsh,ISE_ISEsh,SIC,SICsh)