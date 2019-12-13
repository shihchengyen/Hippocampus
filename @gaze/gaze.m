function [obj,varargout] = gaze(varargin)

% Generate gaze object for each session. 
% 1) Groups raycast data into sections of the environment for ease of
% smoothing in subsequent @spatialview object
% 2) Bins data using same grid size as @vmpc object
% 3) Moves (0,0) of each section to bottom left corner instead of center of
% object (Unity default)

% @gaze Constructor function for gaze class
%   OBJ = gaze(varargin)
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on spatialview %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% example [as, Args] = gaze('save','redo')
%
% dependencies: 


% Use same grid size as in place analysis
Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Session', 'RequiredFile','raycast.mat', 'GridSteps',40);
Args.flags = {'Auto','ArgsOnly'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'GridSteps'};                            

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
    disp('Loading raycast data...');
    rcdata = load(Args.RequiredFile);
    rcdata = rcdata.el.data;

    gridSteps = Args.GridSteps;
    gridSize = 25 / gridSteps;

    %  Group gaze data into object types and get raw gaze position data
    %  relative to bottom left corner of fixated object type
    disp('Re-grouping gaze data...');
    [fixObjNum,RelGazeRawAdj,gazeSections,playerLocAdj] = groupgazesections(rcdata.fixatedObj,rcdata.RelativeToFixdObjGaze,rcdata.playerLocation);

    % Bin relative gaze position data into linear array for session
    disp('Binning gaze data...');
    [binGazeLin,binLocLin,binDepths,binGridRef] = bingazedata(fixObjNum,RelGazeRawAdj,gridSize,playerLocAdj,gridSteps);
    
    % Split binnedRelGaze into trials
    disp('Splitting gaze data into trials...');
    trialTimestamps = rcdata.trialTimestamps;
    timestamps = rcdata.timestamps;
    [binGazeTrial,gpDurGaze,binLocTrial,gpDurLoc,binLocLin,trialInds,timestampsTrial,tTrial] = splittrial(binGazeLin,trialTimestamps,index,timestamps,binDepths,binLocLin,gridSteps);

    %%% Output data
    
    % New variables
    data.binGazeLin = binGazeLin;
    data.binGazeTrial = binGazeTrial;
    data.trialInds = trialInds;
    data.fixObjNum = fixObjNum;
    data.RelGazeRawAdj = RelGazeRawAdj;
    data.gridSize = gridSize;
    data.gridSteps = gridSteps;
    data.binDepths = binDepths;
    data.binGridRef = binGridRef;
    data.gazeSections = gazeSections;
    data.gpDurGaze = gpDurGaze;
    data.binLocLin = binLocLin;
    data.binLocTrial = binLocTrial;
    data.gpDurLoc = gpDurLoc;
    
    % Original raycast-object variables
    data.numTrials = rcdata.numSets;
    data.timestamps = timestamps;
    data.timestampsTrial = timestampsTrial;
    data.tTrial = tTrial;
%     data.RelativeToFixdObjGaze = rcdata.RelativeToFixdObjGaze;
    data.fixatedObj = rcdata.fixatedObj;
%     data.index = rcdata.index;
%     data.fixIndex = rcdata.fixIndex;
%     data.rawGazeData = rcdata.rawGazeData;
    data.playerLocation = rcdata.playerLocation;
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



function [fixObjNum,RelGazeRawAdj,gazeSections,locAdj] = groupgazesections(fixObj,RelGazeRaw,locRaw)
% Groups object of gaze into sections 1-9 below, and adjust gaze
% coordinates relative to bottom left corner instead of center of object
% 1. Cue
% 2. Hint
% 3. Maze floor
% 4. Maze ceiling
% 5. Maze walls
% 6-9. Maze pillars 1-4
fixObjNum = nan(size(fixObj,1),1);
RelGazeRawAdj = nan(size(RelGazeRaw,1),2);
gazeSections = {'Cue','Hint','Ground','Ceiling','Walls','Pillar1','Pillar2','Pillar3','Pillar4'};
for ii = 1:size(gazeSections,2)
    
    switch gazeSections{ii}
        case 'Cue'
            match = ~cellfun(@isempty,regexpi(fixObj,'CueImage'));
            
        case 'Hint'
            match = ~cellfun(@isempty,regexpi(fixObj,'HintImage'));
            
        case 'Ground'
            match = ~cellfun(@isempty,regexpi(fixObj,'Ground'));
            RelGazeRawAdj(match,1) = RelGazeRaw(match,1) + 12.5;
            RelGazeRawAdj(match,2) = RelGazeRaw(match,2) + 12.5;
            % Safeguard for any points exceeding the bounds of object
            RelGazeRawAdj((RelGazeRawAdj(:,1)>25 & match),1) = 25;
            RelGazeRawAdj((RelGazeRawAdj(:,1)<0 & match),1) = 0;
            RelGazeRawAdj((RelGazeRawAdj(:,2)>25 & match),2) = 25;
            RelGazeRawAdj((RelGazeRawAdj(:,2)<0 & match),2) = 0;
            
        case 'Ceiling'
            match1 = ~cellfun(@isempty,regexpi(fixObj,'Ceiling'));
            match2 = ~cellfun(@isempty,regexpi(fixObj,'ceil_'));
            match = match1 | match2;
            % Fix Ceiling
            RelGazeRawAdj(match1,1) = RelGazeRaw(match1,1) + 12.5;
            RelGazeRawAdj(match1,2) = RelGazeRaw(match1,2) + 12.5;
            % Fix ceiling lights
            adjustments = { '01'    9.28    9.47
                            '02'    9.28    5.34
                            '03'    9.28    0.34
                            '04'    9.28    -5.01
                            '05'    9.28    -8.94
                            '06'    0.06    -8.94
                            '07'    0.06    -5.01
                            '08'    0.06    0.34
                            '09'    0.06    5.34
                            '10'    0.06    9.47
                            '11'    -10.93  9.47
                            '12'    -10.93  5.34
                            '13'    -10.93  0.34
                            '14'    -10.93  -5.01
                            '15'    -10.93  -8.94
                            };
            for ll = 1:size(adjustments,1)
                str = horzcat('ceil_lamp',adjustments{ll,1});
                inds = ~cellfun(@isempty,regexpi(fixObj,str));
                RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + adjustments{ll,2} + 12.5;
                RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + adjustments{ll,3} + 12.5;
            end
            % Safeguard for any points exceeding the bounds of object
            RelGazeRawAdj((RelGazeRawAdj(:,1)>25 & match),1) = 25;
            RelGazeRawAdj((RelGazeRawAdj(:,1)<0 & match),1) = 0;
            RelGazeRawAdj((RelGazeRawAdj(:,2)>25 & match),2) = 25;
            RelGazeRawAdj((RelGazeRawAdj(:,2)<0 & match),2) = 0;
            
        case 'Walls'
            match = ~cellfun(@isempty,regexpi(fixObj,'^wall'));
            
            adjustments = { '01'   0  % Adjust center of each m_wall horizontally to top left corner of extended wall
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
                inds = ~cellfun(@isempty,regexpi(fixObj,str));
                RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + adjustments{jj,2};
                RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 2.5;
            end
            % Safeguard for any points exceeding the bounds of object
            RelGazeRawAdj((RelGazeRawAdj(:,1)>100 & match),1) = 100;
            RelGazeRawAdj((RelGazeRawAdj(:,1)<0 & match),1) = 0;
            RelGazeRawAdj((RelGazeRawAdj(:,2)>5 & match),2) = 5;
            RelGazeRawAdj((RelGazeRawAdj(:,2)<0 & match),2) = 0;
            
        case 'Pillar1'
            match1 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_1'));
            match2 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_25'));
            match3 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_26'));
            match4 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_5'));
            match5 = ~cellfun(@isempty,regexpi(fixObj,'RabitPoster')); % m_wall_25
            match = match1 | match2 | match3 | match4 | match5;
            
            adjustments = { '5'     2.5
                            '26'    7.5
                            '25'    12.5
                            '1'     17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixObj,str));
                RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + adjustments{jj,2};
                RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 25 % RabitPoster
                    str = 'RabitPoster';
                    inds = ~cellfun(@isempty,regexpi(fixObj,str));
                    RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + 2.5 - 0.0387;
                    RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                end
            end 
            % Safeguard for any points exceeding the bounds of object
            RelGazeRawAdj((RelGazeRawAdj(:,1)>20 & match),1) = 20;
            RelGazeRawAdj((RelGazeRawAdj(:,1)<0 & match),1) = 0;
            RelGazeRawAdj((RelGazeRawAdj(:,2)>3.11 & match),2) = 3.11;
            RelGazeRawAdj((RelGazeRawAdj(:,2)<0 & match),2) = 0;
            
        case 'Pillar2'
            match1 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_10'));
            match2 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_6'));
            match3 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_29'));
            match4 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_21'));
            match5 = ~cellfun(@isempty,regexpi(fixObj,'CatPoster')); % m_wall_10
            match6 = ~cellfun(@isempty,regexpi(fixObj,'PigPoster')); % m_wall_29
            match = match1 | match2 | match3 | match4 | match5 | match6;
            
            adjustments = { '21'    2.5
                            '29'    7.5
                            '6'     12.5
                            '10'    17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixObj,str));
                RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + adjustments{jj,2};
                RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 10 % CatPoster
                    str = 'CatPoster';
                    inds = ~cellfun(@isempty,regexpi(fixObj,str));
                    RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + 2.5;
                    RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                elseif str2double(adjustments{jj,1}) == 29 % PigPoster
                    str = 'PigPoster';
                    inds = ~cellfun(@isempty,regexpi(fixObj,str));
                    RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + 2.5 - 0.01;
                    RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                end
            end
            % Safeguard for any points exceeding the bounds of object
            RelGazeRawAdj((RelGazeRawAdj(:,1)>20 & match),1) = 20;
            RelGazeRawAdj((RelGazeRawAdj(:,1)<0 & match),1) = 0;
            RelGazeRawAdj((RelGazeRawAdj(:,2)>3.11 & match),2) = 3.11;
            RelGazeRawAdj((RelGazeRawAdj(:,2)<0 & match),2) = 0;
            
        case 'Pillar3'
            match1 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_4'));
            match2 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_24'));
            match3 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_15'));
            match4 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_3'));
            match5 = ~cellfun(@isempty,regexpi(fixObj,'DonkeyPoster')); % m_wall_15
            match6 = ~cellfun(@isempty,regexpi(fixObj,'CrocodilePoster')); % m_wall_4
            match = match1 | match2 | match3 | match4 | match5 | match6;
            
            adjustments = { '3'     2.5
                            '15'    7.5
                            '24'    12.5
                            '4'     17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixObj,str));
                RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + adjustments{jj,2};
                RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 15 % DonkeyPoster
                    str = 'DonkeyPoster';
                    inds = ~cellfun(@isempty,regexpi(fixObj,str));
                    RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + 2.5 - 0.168;
                    RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                elseif str2double(adjustments{jj,1}) == 4 % PigPoster
                    str = 'CrocodilePoster';
                    inds = ~cellfun(@isempty,regexpi(fixObj,str));
                    RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + 2.5 + 0.021;
                    RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                end
            end
            % Safeguard for any points exceeding the bounds of object
            RelGazeRawAdj((RelGazeRawAdj(:,1)>20 & match),1) = 20;
            RelGazeRawAdj((RelGazeRawAdj(:,1)<0 & match),1) = 0;
            RelGazeRawAdj((RelGazeRawAdj(:,2)>3.11 & match),2) = 3.11;
            RelGazeRawAdj((RelGazeRawAdj(:,2)<0 & match),2) = 0;
            
        case 'Pillar4'
            match1 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_7'));
            match2 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_8'));
            match3 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_12'));
            match4 = ~cellfun(@isempty,regexpi(fixObj,'m_wall_20'));
            match5 = ~cellfun(@isempty,regexpi(fixObj,'CamelPoster')); % m_wall_20
            match = match1 | match2 | match3 | match4 | match5;
            
            adjustments = { '20'   2.5
                            '12'   7.5
                            '8'    12.5
                            '7'    17.5
                            };
            for jj = 1:size(adjustments,1)
                str = horzcat('m_wall_',adjustments{jj,1});
                inds = ~cellfun(@isempty,regexpi(fixObj,str));
                RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + adjustments{jj,2};
                RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                if str2double(adjustments{jj,1}) == 20 % DonkeyPoster
                    str = 'CamelPoster';
                    inds = ~cellfun(@isempty,regexpi(fixObj,str));
                    RelGazeRawAdj(inds,1) = RelGazeRaw(inds,1) + 2.5 + 0.025;
                    RelGazeRawAdj(inds,2) = RelGazeRaw(inds,2) + 1.555;
                end
            end
            % Safeguard for any points exceeding the bounds of object
            RelGazeRawAdj((RelGazeRawAdj(:,1)>20 & match),1) = 20;
            RelGazeRawAdj((RelGazeRawAdj(:,1)<0 & match),1) = 0;
            RelGazeRawAdj((RelGazeRawAdj(:,2)>3.11 & match),2) = 3.11;
            RelGazeRawAdj((RelGazeRawAdj(:,2)<0 & match),2) = 0;
            
    end
    % Check that there are no overlapping matches
    matched = unique(~isnan(fixObjNum(match)));
    if matched ~= 0 
        return;
    end
    fixObjNum(match) = ii;
end
% Adjust player location to start from top left corner of floor
locAdj(:,1) = locRaw(:,1) + 12.5;
locAdj(:,3) = locRaw(:,3) + 12.5;
% Safeguard for any points exceeding the bounds of object
locAdj((locAdj(:,1)>25),1) = 25;
locAdj((locAdj(:,1)<0),1) = 0;
locAdj((locAdj(:,3)>25),2) = 25;
locAdj((locAdj(:,3)<0),2) = 0;


function [binGazeLin,binLocLin,binDepths,binGridRef] = bingazedata(fixObjNum,RelGazeAdj,gridSize,locAdj,gridSteps)
% GAZE BIN PARTITIONS
% Linear bin format, if 40x40 grid on floor
% CueImage: bin 1
% HintImage: bin 2
% Ground: bin 3 - 1602 (25Wx25L in Unity units) left to right, bottom to top
% Ceiling: bin 1603 - 3202 bins (25Wx25L)
% Walls: bin 3203 - 4482 (5Hx100W) 
% Pillar1: bin 4483 - 4642 (3.11Hx20W) 
% Pillar2: bin 4643 - 4802 (3.11Hx20W) 
% Pillar3: bin 4803 - 4962 (3.11Hx20W) 
% Pillar4: bin 4963 - 5122 (3.11Hx20W) 

% binDepths in matrix x-y form, but Unity axes are reversed
binDepths = [   1                       1       ;  % CueImage
                1                       1       ;  % HintImage
            	ceil(25/gridSize)       ceil(25/gridSize);  % Ground
            	ceil(25/gridSize)       ceil(25/gridSize);  % Ceiling
            	ceil(5/gridSize)        ceil(100/gridSize); % Walls
            	ceil(3.11/gridSize)     ceil(20/gridSize); % Pillar1
            	ceil(3.11/gridSize)     ceil(20/gridSize); % Pillar2
            	ceil(3.11/gridSize)     ceil(20/gridSize); % Pillar3
            	ceil(3.11/gridSize)     ceil(20/gridSize); % Pillar4
            ];
binGridRef = cell(size(binDepths,1),1);
for ii = 1:size(binGridRef,1)
    temp = nan(binDepths(ii,1),binDepths(ii,2));
    base = sum(binDepths(1:(ii-1),1).*binDepths(1:(ii-1),2));
    for jj = 1:binDepths(ii,1) % For each row
        temp(jj,:) = (jj-1)*binDepths(ii,2)+base+(1:binDepths(ii,2));
    end
    temp = flipud(temp);
    binGridRef{ii} = temp;
end


binGazeLin = nan(size(fixObjNum,1),1);
binLocLin = nan(size(fixObjNum,1),1);
% binH is the col number, binV is the row
% binGazeGrid = nan(size(fixObjNum,1),2);
% binLocGrid = nan(size(fixObjNum,1),2);
for ii = 1:size(fixObjNum,1) 
    
    if ~isnan(fixObjNum(ii))
        
        % GAZE
        if fixObjNum(ii) == 1 % CueImage
            binGazeLin(ii) = 1;
%             binGazeGrid(ii,:) = [1 1];
        elseif fixObjNum(ii) == 2 % HintImage
            binGazeLin(ii) = sum(binDepths(1:fixObjNum(ii)-1,1).*binDepths(1:fixObjNum(ii)-1,2)) + 1;
%             binGazeGrid(ii,:) = [1 1];
        else
            % Make sure values at 0 are correctly captured in grid bin 1
            for kk = 1:2
                if RelGazeAdj(ii,kk) == 0
                    RelGazeAdj(ii,kk) = RelGazeAdj(ii,kk) + gridSize/2;
                end
            end
            binGazeLin(ii) = sum(binDepths(1:fixObjNum(ii)-1,1).*binDepths(1:fixObjNum(ii)-1,2)) + ...
                floor(abs(RelGazeAdj(ii,2)/gridSize))*binDepths(fixObjNum(ii),2) + ...
                ceil(abs(RelGazeAdj(ii,1))/gridSize);
%             binGazeGrid(ii,:) = [ceil(abs(RelGazeAdj(ii,2)/gridSize)) ceil(abs(RelGazeAdj(ii,1))/gridSize)];
%             binGazeGrid(ii,:) = [binDepths(fixObjNum(ii),1)-ceil(abs(RelGazeAdj(ii,2)/gridSize))+1 ceil(abs(RelGazeAdj(ii,1))/gridSize)];
            % Sanity check to make sure binning is correct
            if floor(abs(RelGazeAdj(ii,2))/gridSize)*binDepths(fixObjNum(ii),2) < 0
                disp(ii);
            end
        end
        
        % POSITION
        if floor(abs(locAdj(ii,3)/gridSize)) == abs(locAdj(ii,3)/gridSize)
            binLocLin(ii) = (floor(abs(locAdj(ii,3)/gridSize))-1)*gridSteps + ...
                ceil(abs(locAdj(ii,1))/gridSize);
        else
            binLocLin(ii) = floor(abs(locAdj(ii,3)/gridSize))*gridSteps + ...
                ceil(abs(locAdj(ii,1))/gridSize);
        end
%         binLocGrid(ii,:) = [ceil(abs(locAdj(ii,3)/gridSize)) ceil(abs(locAdj(ii,1))/gridSize)];
    end
    
end


function [binGazeTrial,gpDurGaze,binLocTrial,gpDurLoc,binLocLin,trialInds,timestampsTrial,tTrial] = splittrial(binGazeLin,timestampsOrigEvent,index,timestampsOrigAll,binDepths,binLocLin,gridSteps)

ntrial = size(timestampsOrigEvent,1);
longestDur = max(timestampsOrigEvent(:,3) - timestampsOrigEvent(:,2))+2;
binGazeTrial = nan(longestDur+1,ntrial);
binLocTrial = nan(longestDur+1,ntrial);
timestampsTrial = nan(longestDur+1,ntrial);
tTrial = nan(longestDur+1,ntrial);
trialInds = nan(ntrial,2);
gpDurGaze = nan(sum(binDepths(:,1).*(binDepths(:,2))),ntrial);
gpDurLoc = nan(gridSteps*gridSteps,ntrial);

for ii = 1:ntrial
    
    % Find trial indices from cue offset
%     inds = find(ismember(timestampsOrigAll,timestampsOrigEvent(ii,2):0.001:timestampsOrigEvent(ii,3)));
    inds = index(ii,2):index(ii,3);
    % get indices for this trial 
    indsx = inds(2:end-1);
    % Adjust trial timestamps to seconds from ms
    timestampsWithinTrial = (timestampsOrigAll(inds))/1000;
    timestampsRelStartNav = timestampsWithinTrial(2:end-1)-timestampsWithinTrial(1);
    
    % get gaze grid positions for this trial - Remove 1st and last points which
    % correspond to event markers for cue offset and end of trial
    gaze = binGazeLin(indsx);
    binGazeTrial(1:size(indsx,1),ii) = gaze;
    trialInds(ii,:) = [indsx(1) indsx(end)];
%     timestampsTrial(1:size(indsx,1),ii) = timestampsWithinTrial(2:end-1);
    timestampsTrial(1:size(indsx,1),ii) = timestampsOrigAll(indsx);
    tTrial(1:size(indsx,1),ii) = timestampsRelStartNav;
    
    % get unique gaze positions
    ugaze = unique(gaze(~isnan(gaze)));
    if ugaze(end) > sum(binDepths(:,1).*(binDepths(:,2)))
        disp(ii)
    end
    % Get duration spend in each grid position for this trial
    for gg = 1:size(ugaze,1)
        tempgp = ugaze(gg);
        % find indices that have this grid position
        utgpidx = find(gaze==tempgp);
        gpDurGaze(tempgp,ii) = (size(utgpidx,1))/1000;
    end
    
    % get location grid positions for this trial - Remove 1st and last points which
    % correspond to event markers for cue offset and end of trial
    loc = binLocLin(indsx);
    binLocTrial(1:size(indsx,1),ii) = loc;
    
    % get unique location positions
    uloc = unique(loc(~isnan(loc)));
    % Get duration spend in each grid position for this trial
    for gg = 1:size(uloc,1)
        tempgp = uloc(gg);
        % find indices that have this grid position
        utgpidx = find(loc==tempgp);
        gpDurLoc(tempgp,ii) = (size(utgpidx,1))/1000;
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


