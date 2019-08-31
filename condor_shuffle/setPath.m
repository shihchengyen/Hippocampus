function setPath

cd ./DPV; 
nptAddPath; 
cd ..; 
cd ./newNpt; 
nptAddPath; 
cd ..; 
cd ./Hippocampus; 
addpath(pwd);
cd Compiler; 
addpath(genpath(pwd)); 
cd ../..;  
