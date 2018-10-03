function [p,pE] = insertFFT(info,numArray,channelGeo,setX,setY,offsetX,offsetY)
%INSERTFFT Insert the png into its corresponding figure
%   Detailed explanation goes here

for i = 1:numArray
    p(i,1) = figure; % for session01
    set(gcf, 'Position', get(0,'Screensize')-[0 0 0 80],'PaperPositionMode', 'auto');
    pE(i,1) = figure; % for sessioneye
    set(gcf, 'Position', get(0,'Screensize')-[0 0 0 80],'PaperPositionMode', 'auto');
end

numHp = length(info);

for i = 1:numHp
    indexArray = strfind(info(i).folder,'array');
    arrayNum = str2num(info(i).folder(indexArray+5:indexArray+6)); % array number
    
    indexChannel = strfind(info(i).folder,'channel');
    channelNum = str2num(info(i).folder(indexChannel+7:indexChannel+9)); % channel number
    
    [y,x] = find(channelGeo(:,:,arrayNum) == channelNum); % location of the channel in channel geometry    
    
    pTemp = fullfile(info(i,1).folder,info(i,1).name);

    if isempty(strfind(info(i,1).folder,'eye'))
        set(0,'CurrentFigure',p(arrayNum,1)); % change to the target figure of session01
    else
        set(0,'CurrentFigure',pE(arrayNum,1)); % change to the target figure of sessiontest
    end
    
    axes('pos',[setX(x), setY(y), offsetX, offsetY]);
    imshow(pTemp)
    
end

end

