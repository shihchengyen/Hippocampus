% Crop frames
for i = 0:1074; % edit accordingly   
    num = int2str(i);
    filename = strcat(num,'.png');
    pict = imread(filename);
    pictC = imcrop(pict,[1915 0 1913 1002]);
    imwrite(pictC, filename);
end

% USE IMAGEJ MACROS TO SUPERIMPOSE ON 1980x1020 BGD

