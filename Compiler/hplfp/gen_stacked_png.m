% Programme Description
% Generate stacked_hp.png and stacked_lfp.png from png.tar which contains information about hp.png and lfp.png using HPC

function [] = gen_stacked_png_hpc()
	% Set directory numbers
	array_dir = pwd;
	array_dir_no = array_dir(end);
    session_dir_no = array_dir(end - 9 : end - 8);
	
	% Plot stacked images
    plot_stacked_png(session_dir_no, array_dir_no);
end

% Gather all the lfp.png nad hp.png in the channels
% Plot the images togther in different figures referred to certain geometry specified in function 'assignChannelGeo'
% cf. ${GITHUB_MATLAB}/Hippocampus/Compiler/hplfp/plotStackedFFT
function [] = plot_stacked_png(session_dir_no, array_dir_no)
% session_dir_no : Session directory number
% array_dir_no : Array directory number
	
	% Set image geometry
    numRow = 5;
    numArray = 4;
    channelGeo = assignChannelGeo(numRow, numArray); % Get the channel geometry
    numCol = size(channelGeo, 2);
    resolutionValue = '-r1200'; % -r???, ??? is the resolution value
    offsetX = 1 / numCol;
    offsetY = 1 / numRow;
    setX = 0 : offsetX : 1;
    setY = fliplr(0 : offsetY : 1 - offsetY);
    
	% Check if png.tar exists
	if isfile('png.tar') % E
		untar('png.tar'); % Untar png.tar
	end
	
	% Get the information of images
    hpInfo = dir('**/hp*.png');
    lfpInfo = dir('**/lfp*.png');

	% Close all figures
    close all  
    
	% Check if lfp.png exists
    if isempty(lfpInfo) == 0 % E
		% Plot stacked_lfp.png
        [pLfp, pELfp] = insertFFT(lfpInfo, numArray, channelGeo, setX, setY, offsetX, offsetY);
		
		% Set the name of the image
        for i = 1 : numArray
			% Check if the array directory exists
            if any(i == str2num(array_dir_no)) == 1 % E
                if contains(session_dir_no, 'e') == 0 % session??
                    set(pLfp(i, 1), 'Visible', 'on'); 
                    fileNameTemp = ['stackedlfp', session_dir_no, '-a', num2str(i), '.png'];
                    print(fileNameTemp, '-dpng', resolutionValue); 
				else % sessioneye
                    set(pELfp(i, 1), 'Visible', 'on');
                    fileNameTemp = ['stackedlfpeye-a', num2str(i), '.png']; 
                    print(fileNameTemp, '-dpng', resolutionValue);
                    delete(pELfp(i, 1)); 
                end
            end
        end
	end
    
	% Check if hp.png exists
    if isempty(hpInfo) == 0 % E
		% Plot stacked_hp.png
        [pHp, pEHp] = insertFFT(hpInfo, numArray, channelGeo, setX, setY, offsetX, offsetY);
		
		% Set the name of the image
        for i = 1 : numArray
			% Check if the array directory exists
            if any(i == str2num(array_dir_no)) == 1 % E
                if contains(session_dir_no, 'e') == 0 % session??
                    set(pHp(i,1), 'Visible', 'on');
                    fileNameTemp = ['stackedhp' session_dir_no, '-a', num2str(i), '.png']; 
                    print(fileNameTemp, '-dpng', resolutionValue); 
                    delete(pHp(i, 1));
                else % sessioneye
                    set(pEHp(i, 1), 'Visible', 'on'); 
                    fileNameTemp = ['stackedhpeye-a',num2str(i),'.png'];
                    print(fileNameTemp, '-dpng', resolutionValue); 
                    delete(pEHp(i, 1));
                end
            end
        end    
    end
end

% Assign the geometry of the channel into each array
function channelGeo = assignChannelGeo(numRow, numArray)
% numRow : Row number
% numArray: Array number
    
    channelTbl = [7,38; 31,70; 63,102; 95,118];
    channelGeo = horzcat(nan, fliplr(1:6), nan); % First array first row of channel is different from others
    channelGeo = vertcat(channelGeo, ...
				 fliplr(reshape(transpose(reshape(channelTbl(1,1):channelTbl(1,2),[],numRow-1)),numRow-1,[])));

	% Second and third array
    for i = 2 : numArray - 1 
        channelGeo(:,:,i) = fliplr(reshape(transpose(reshape(channelTbl(i,1):channelTbl(i,2),[],numRow)),numRow,[]));
	end

	% Last array
    channelGeo(:,:,4) = vertcat( ... 
						fliplr(reshape(transpose(reshape(channelTbl(4,1):channelTbl(4,2),[],numRow-2)),numRow-2,[])),...
						horzcat(nan, fliplr(119:124), nan),...
						nan(1,size(channelGeo,2)));

end

function [p,pE] = insertFFT(info,numArray,channelGeo,setX,setY,offsetX,offsetY)

    p = gobjects(numArray,1);
    pE = gobjects(numArray,1);
    
    set(0, 'DefaultFigureVisible', 'off');

    for i = 1:numArray
        p(i,1) = figure; % session??
        pE(i,1) = figure; % sessioneye
    end

    numP = length(info);

    for i = 1:numP
        indexArray = strfind(info(i).folder,'array');
        arrayNum = str2double(info(i).folder(indexArray+5:indexArray+6)); 

        indexChannel = strfind(info(i).folder,'channel');
        channelNum = str2double(info(i).folder(indexChannel+7:indexChannel+9)); 

        if channelNum < 125
            [y,x] = find(channelGeo(:,:,arrayNum) == channelNum);
            pTemp = fullfile(info(i,1).folder,info(i,1).name);

            if ~isempty(strfind(info(i,1).folder,'eye'))
                set(0,'CurrentFigure',pE(arrayNum,1));
            else
                set(0,'CurrentFigure',p(arrayNum,1)); 
            end

            axes('pos',[setX(x), setY(y), offsetX, offsetY]);
            imshow(pTemp)
			
			% Set the title of each image (Unnecessary if lfp.png or hp.png already has its own title)
			path = strsplit(info(i).folder, '/');
			session_dir_no = regexprep(path(end - 2), 'session', '');
			array_dir_no = regexprep(path(end - 1), 'array', '');
			channel_dir_no = regexprep(path(end), 'channel', '');
			name = strcat({'p'}, path(end - 3), {'s'}, session_dir_no, {'a'}, array_dir_no, {'c'}, channel_dir_no);
			title(name, 'fontsize', 4);			
        end
	end
end
