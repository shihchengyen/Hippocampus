function channelGeo = assignChannelGeo(numRow,numArray)
%ASSIGNCHANNELGEO To assign the geometry of the channel into each array
% function channelGeo = assignChannelGeo(numRow,numArray)

channelTbl = [7,38; 31,70; 63,102; 95,118];

channelGeo = horzcat(nan, fliplr(1:6), nan); % first array first row of channel is different from others
channelGeo = vertcat(channelGeo, ...
    fliplr(reshape(transpose(reshape(channelTbl(1,1):channelTbl(1,2),[],numRow-1)),numRow-1,[])));

for i = 2 : numArray-1 % second and third array
    channelGeo(:,:,i) = fliplr(reshape(transpose(reshape(channelTbl(i,1):channelTbl(i,2),[],numRow)),numRow,[]));
end

channelGeo(:,:,4) = vertcat( ... % last array
    fliplr(reshape(transpose(reshape(channelTbl(4,1):channelTbl(4,2),[],numRow-2)),numRow-2,[])),...
    horzcat(fliplr(119:125), nan),...
    nan(1,size(channelGeo,2)));


end