function [] = plotStackedFFT()
%PLOTSTACKEDFFT Gather all the FFT png in the session01 and sessioneye in
%the current directory, plot the FFT png togther in dirrent figures showing
%each array, referred to certain geometry specified in function
%'assignChannelGeo'
% 
%   Detailed explanation goes here

numRow = 5;
numArray = 4;

channelGeo = assignChannelGeo(numRow,numArray); % a function to get the channel geometry

hpInfo = dir('**/*hp.png');

lfpInfo = dir('**/*lfp.png');
numLfp = length(lfpInfo);

%% plot stacked hp.png
close all

numCol = size(channelGeo,2);

offsetX = 1/numCol;
offsetY = 1/numRow;
setX = 0 : offsetX : 1;
setY = fliplr(0 : offsetY : 1-offsetY);

[pHp,pEHp] = insertFFT(hpInfo,numArray,channelGeo,setX,setY,offsetX,offsetY);
for i = 1:numArray
    saveas(pHp(i,1),['session01',filesep,'array0',num2str(i),filesep,'stackedhp.png']);
    delete(pHp(i,1));
    
    saveas(pEHp(i,1),['sessioneye',filesep,'array0',num2str(i),filesep,'stackedhp.png']);
    delete(pEHp(i,1));
end    


[pLfp,pELfp] = insertFFT(lfpInfo,numArray,channelGeo,setX,setY,offsetX,offsetY);
for i = 1:numArray
    saveas(pLfp(i,1),['session01',filesep,'array0',num2str(i),filesep,'stackedhp.png']);
    delete(pLfp(i,1));
    
    saveas(pELfp(i,1),['sessioneye',filesep,'array0',num2str(i),filesep,'stackedhp.png']);
    delete(pLfp(i,1));
end    


%% plot stacked lfp.png
for i = 1:numLfp
end

end

