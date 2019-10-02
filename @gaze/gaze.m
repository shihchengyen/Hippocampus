function [] = gaze(varargin)

% Generate spatial view object for each session. 

% @spatialview Constructor function for spatialview class
%   OBJ = spatialview(varargin)
%
%   OBJ = spatialview('auto') attempts to create a raycast object by 
%   extracting data from a csv file.
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on spatialview %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% example [as, Args] = spatialview('save','redo')
%
% dependencies: 


% Use same grid size as in place analysis
Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'RequiredFile','raycast.mat', 'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'FRSIC',0, 'UseAllTrials',0, 'UseMedian',0, ...
                'NumFRBins',4,'AdaptiveSmooth',1, 'FiltLowOcc',0);
Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'gaze';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'gz';

% To decide the method to create or load the object
[command,robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
    % l = load(Args.matname);
    % obj = eval(['l.' Args.matvarname]);
	obj = robj;
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

if(~isempty(dir(Args.RequiredFile)))
	% Load raycast data
    rcdata = load(Args.RequiredFile);
    rcdata = rcdata.el.data;

    gridSteps = Args.GridSteps;
    gridSize = 25 / gridSteps;

    %  Group gaze data into object types and get raw gaze position data
    %  relative to top left corner of fixated object type
    [fixatedObjNum,RelativeToFixdObjGazeAdj,gazeSections] = groupgazesections(rcdata.fixatedObj,rcdata.RelativeToFixdObjGaze);

    % Bin relative gaze position data into linear array for session
    [binnedRelGaze,binDepths] = bingazedata(fixatedObjNum, RelativeToFixdObjGazeAdj, gridSize);
    
    % Split binnedRelGaze into trials
    [binnedRelGazeTrial] = splittrial(binnedRelGaze,rcdata.trialTimestamps,rcdata.timestamps,binDepths);

    % Output data
    % New variables
    data.binnedRelGaze = binnedRelGaze;
    data.binnedRelGazeTrial = binnedRelGazeTrial;
    data.fixatedObjNum = fixatedObjNum;
    data.RelativeToFixdObjGazeAdj = RelativeToFixdObjGazeAdj;
    data.gridSize = gridSize;
    data.gridSteps = gridSteps;
    data.binDepths = binDepths;
    data.gazeSections = gazeSections;
    % Original raycast-object variables
    data.numTrials = rcdata.numSets;
    data.timestamps = rcdata.timestamps;
    data.trialTimestamps = rcdata.trialTimestamps;
%     data.fixatedObjLoc = rcdata.fixatedObjLoc;
%     data.RelativeToFixdObjGaze = rcdata.RelativeToFixdObjGaze;
%     data.fixatedObj = rcdata.fixatedObj;
%     data.RelativeToFixdObjGaze = rcdata.RelativeToFixdObjGaze;
%     data.index = rcdata.index;
%     data.fixIndex = rcdata.fixIndex;
%     data.rawGazeData = rcdata.rawGazeData;
%     data.playerLocation = rcdata.playerLocation;
%     data.playerGazeLocation = rcdata.playerGazeLocation;
    
	% create nptdata so we can inherit from it    
	data.numSets = 1;
	data.Args = Args;
	n = nptdata(1,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);

else
	obj = createEmptyObject(Args);
end

% save spatialview.mat data



function [fixatedObjNum,RelativeToFixdObjGazeNew,gazeSections] = groupgazesections(fixatedObj,RelativeToFixdObjGaze)
% Groups object of gaze into sections 1-9 below, and adjust gaze
% coordinates relative to top left corner instead of center of object
% 1. Cue
% 2. Hint
% 3. Maze floor
% 4. Maze ceiling
% 5. Maze walls
% 6-9. Maze pillars 1-4
fixatedObjNum = nan(size(fixatedObj,1),1);
RelativeToFixdObjGazeNew = nan(size(RelativeToFixdObjGaze,1),2);
gazeSections = {'Cue','Hint','Ground','Ceiling','Walls','Pillar1','Pillar2','Pillar3','Pillar4'};
for ii = 1:size(gazeSections,2)
    
    switch gazeSections{ii}
        case 'Cue'
            match = ~cellfun(@isempty,regexpi(fixatedObj,'CueImage'));
            
        case 'Hint'
            match = ~cellfun(@isempty,regexpi(fixatedObj,'HintImage'));
            
        case 'Ground'
            match = ~cellfun(@isempty,regexpi(fixatedObj,'Ground'));
            RelativeToFixdObjGazeNew(match,:) = RelativeToFixdObjGaze(match,:) + 12.5;
            
        case 'Ceiling'
            match = ~cellfun(@isempty,regexpi(fixatedObj,'Ceiling'));
            RelativeToFixdObjGazeNew(match,:) = RelativeToFixdObjGaze(match,:) + 12.5;
            %%%%%%% SORT CEILING LAMPS
        case 'Walls'
            match = ~cellfun(@isempty,regexpi(fixatedObj,'^wall'));
            
            adjustments = { '01'   0  % Adjust center of each m_wall to top left corner of extended wall
                            '02'   5
                            '03'   10
                            '04'   15
                            '05'   20
                            '06'   25
                            '07'   25
                            '08'   30
                            '09'   35
                            '10'   40
                            '11'   45
                            '12'   50
                            '13'   50
                            '14'   55
                            '15'   60
                            '16'   65
                            '17'   70
                            '18'   75
                            '19'   75
                            '20'   80
                            '21'   85
                            '22'   90
                            '23'   95
                            '24'   100
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('^wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + adjustments{jj,2};
                RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 2.5;
            end
            
        case 'Pillar1'
            match1 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_1'));
            match2 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_25'));
            match3 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_26'));
            match4 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_5'));
            match5 = ~cellfun(@isempty,regexpi(fixatedObj,'RabitPoster')); % m_wall_25
            match = match1 | match2 | match3 | match4 | match5;
            
            adjustments = { '5'     2.5
                            '26'    7.5
                            '25'    12.5
                            '1'     17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + adjustments{jj,2};
                RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 25 % RabitPoster
                    str = 'RabitPoster';
                    inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                    RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + 2.5 - 0.0387;
                    RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                end
            end 
            
        case 'Pillar2'
            match1 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_10'));
            match2 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_6'));
            match3 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_29'));
            match4 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_21'));
            match5 = ~cellfun(@isempty,regexpi(fixatedObj,'CatPoster')); % m_wall_10
            match6 = ~cellfun(@isempty,regexpi(fixatedObj,'PigPoster')); % m_wall_29
            match = match1 | match2 | match3 | match4 | match5 | match6;
            
            adjustments = { '21'    2.5
                            '29'    7.5
                            '6'     12.5
                            '10'    17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGazeNew(inds,1) + adjustments{jj,2};
                RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 10 % CatPoster
                    str = 'CatPoster';
                    inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                    RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + 2.5;
                    RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                elseif str2double(adjustments{jj,1}) == 29 % PigPoster
                    str = 'PigPoster';
                    inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                    RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + 2.5 - 0.01;
                    RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                end
            end
            
        case 'Pillar3'
            match1 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_4'));
            match2 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_24'));
            match3 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_15'));
            match4 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_3'));
            match5 = ~cellfun(@isempty,regexpi(fixatedObj,'DonkeyPoster')); % m_wall_15
            match6 = ~cellfun(@isempty,regexpi(fixatedObj,'CrocodilePoster')); % m_wall_4
            match = match1 | match2 | match3 | match4 | match5 | match6;
            
            adjustments = { '3'     2.5
                            '15'    7.5
                            '24'    12.5
                            '4'     17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGazeNew(inds,1) + adjustments{jj,2};
                RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 15 % DonkeyPoster
                    str = 'DonkeyPoster';
                    inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                    RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + 2.5 - 0.168;
                    RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                elseif str2double(adjustments{jj,1}) == 4 % PigPoster
                    str = 'CrocodilePoster';
                    inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                    RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + 2.5 + 0.021;
                    RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                end
            end
            
        case 'Pillar4'
            match1 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_7'));
            match2 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_8'));
            match3 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_12'));
            match4 = ~cellfun(@isempty,regexpi(fixatedObj,'m_wall_20'));
            match5 = ~cellfun(@isempty,regexpi(fixatedObj,'CamelPoster')); % m_wall_20
            match = match1 | match2 | match3 | match4 | match5;
            
            adjustments = { '20'   2.5
                            '12'   7.5
                            '8'    12.5
                            '7'    17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGazeNew(inds,1) + adjustments{jj,2};
                RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 20 % DonkeyPoster
                    str = 'CamelPoster';
                    inds = ~cellfun(@isempty,regexpi(fixatedObj,str));
                    RelativeToFixdObjGazeNew(inds,1) = RelativeToFixdObjGaze(inds,1) + 2.5 + 0.025;
                    RelativeToFixdObjGazeNew(inds,2) = RelativeToFixdObjGaze(inds,2) + 1.555;
                end
            end
            
    end
    % Check that there are no overlapping matches
    matched = unique(~isnan(fixatedObjNum(match)));
    if matched ~= 0 
        return;
    end
    fixatedObjNum(match) = ii;
end



function [binnedGaze,binDepths] = bingazedata(fixatedObjNum, RelativeToFixdObjGazeNew, gridSize)
% Linear bin format, if 40x40 grid on floor
% CueImage: bin 1
% HintImage: bin 2
% Ground: bin 3 - 1602 (25Wx25L in Unity units) left to right, top to bottom
% Ceiling: bin 1603 - 3202 bins (25Wx25L)
% Walls: bin 3203 - 5250 (5Hx100W) 
% Pillar1: bin 5251 - 5313 (3.11Hx20W) 
% Pillar2: bin 5314 - 5376 (3.11Hx20W) 
% Pillar3: bin 5377 - 5439 (3.11Hx20W) 
% Pillar4: bin 5440 - 5502 (3.11Hx20W) 

binDepths = [   1;  % CueImage
                1;  % HintImage
            	(25/gridSize)*(25/gridSize);  % Ground
            	(25/gridSize)*(25/gridSize);  % Ceiling
            	(100/gridSize)*(8/gridSize); % Walls
            	ceil((20/gridSize)*(3.11*gridSize)); % Pillar1
            	ceil((20/gridSize)*(3.11*gridSize)); % Pillar2
            	ceil((20/gridSize)*(3.11*gridSize)); % Pillar3
            	ceil((20/gridSize)*(3.11*gridSize)); % Pillar4
            ];


binnedGaze = nan(size(fixatedObjNum,1),1);
for ii = 1:size(fixatedObjNum,1) 
    
    if ~isnan(fixatedObjNum(ii))
    
        if fixatedObjNum(ii) == 1 % CueImage
            binnedGaze(ii) = 1;
        elseif fixatedObjNum(ii) == 2 % HintImage
            binnedGaze(ii) = sum(binDepths(1:fixatedObjNum(ii)-1)) + 1;
        else
            binnedGaze(ii) = sum(binDepths(1:fixatedObjNum(ii)-1)) + ceil(RelativeToFixdObjGazeNew(ii,1)/gridSize) + floor(RelativeToFixdObjGazeNew(ii,2)/gridSize)*40;
        end
        
    end
    
end


function [binnedRelGazeTrial] = splittrial(binnedRelGaze,trialTimestamps,gazetimestamps,binDepths)

binnedRelGazeTrial = nan(sum(binDepths,1),1);
for ii = 1:size(trialTimestamps,1)
    
    indTrial = ismember(gazetimestamps,trialTimestamps(ii,2):trialTimestamps(ii,3));
    trialbins = binnedRelGaze(indTrial);
    uniquebins = unique(trialbins(~isnan(trialbins)));
    for jj = 1:size(uniquebins,1)
        binnedRelGazeTrial(uniquebins(jj),ii) = sum(binnedRelGaze == uniquebins(jj));
    end
    
end



function obj = createEmptyObject(Args)

% create nptdata so we can inherit from it
% useful fields for most objects
data.numSets = 0;
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);


