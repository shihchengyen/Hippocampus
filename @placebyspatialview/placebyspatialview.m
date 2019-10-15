function [obj, varargout] = placebyspatialview(varargin)
%@vmplacecell Constructor function for placebyspatialview class
%   OBJ = spatialview(varargin)
%
%   OBJ = spatialview('auto') attempts to create a spatialview object
%   using an unitymaze object, a rplparallel object, and the spiketrain.mat
%   in a Cell directory.
%   
%Example, to create a spatialview object in a cell directory:
%   sv = spatialview('auto');
%
%To look at the fields in the object:
%   sv.data
%
%To create a spatialview object, and save it in the current directory:
%   sv = spatialview('auto','save');
%
%To create a spatialview object with different GridSteps, which will be
%passed to the unitymaze object, and to save both the new unitymaze object
%and the spatialview objects:
%   sv = spatialview('auto','GridSteps',10,'SaveLevels',2);
%
%To create spatialview from a channel directory:
%   sv = ProcessLevel(vmplacecell,'Levels','Channel','save');
%
%To create spatialview from an array directory:
%   sv = ProcessLevel(spatialview,'Levels','Array','save');
%
%dependencies: unitymaze, rplparallel

Args = struct('RedoLevels',0, 'SaveLevels',1, 'Auto',0, 'ArgsOnly',0, ...
				'ObjectLevel','Cell', 'RequiredFile','spatialview.mat', ...
				'MaxTimeDiff',0.002,'MinTrials',5, 'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'UseAllTrials',1, 'AdaptiveSmooth',1, 'FiltLowOcc',0);
Args.flags = {'Auto','ArgsOnly','HPC','FRSIC','UseAllTrials','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'MinTrials','GridSteps','ShuffleLimits', ...
	'NumShuffles','UseAllTrials'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'placebyspatialview';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'psv';

% To decide the method to create or load the object
[command, robj] = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
	obj = robj;
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
    obj = createObject(Args,modvarargin{:});
end

function obj = createObject(Args,varargin)

if(~isempty(dir(Args.RequiredFile)))
	% load gaze object
	sv = spatialview('auto','NumShuffles',Args.NumShuffles,'GridSteps',40,'UseAllTrials',1,'FiltLowOcc',0);
    vp = vmplacecell('auto','NumShuffles',Args.NumShuffles,'GridSteps',40,'UseAllTrials',1,'FiltLowOcc',0);
    
    % Find place SI
    vpSI = vp.data.SIC;
    % Find place SI threshold
    vpSIthresh = prctile(vp.data.SICsh(2:end,1),95);
    % Find spatialview SI
    svSI = sv.data.SIC;
    % Find spatialview SI threshold
    siSIthresh = prctile(sv.data.SICsh(2:end,1),95);
    binDepths = sv.data.binDepths;
    % Find 
    % if place SI cross the SI threshold, analyse for spatial
    count = 1;
    count_i = [];
    if vpSI > vpSIthresh
        
        % Find fields with at least 1 pixel of > half peak rate
        peakrate = max(vp.data.maps_adsmooth);
        threshrate = peakrate/2;
        % Restructure linear map into grid map 
        mapadsmLin = vp.data.maps_adsmooth;
        maprawLin = vp.data.maps_raw;
        binNumLin = 1:Args.GridSteps*Args.GridSteps;
        mapadsmGrid = nan(Args.GridSteps);
        maprawGrid = nan(Args.GridSteps);
        binNumGrid = nan(Args.GridSteps);
        for ii = 1:Args.GridSteps
            mapadsmGrid(ii,:) = mapadsmLin( ((ii-1)*Args.GridSteps)+1:ii*Args.GridSteps );
            maprawGrid(ii,:) = maprawLin( ((ii-1)*Args.GridSteps)+1:ii*Args.GridSteps );
            binNumGrid(ii,:) = binNumLin( ((ii-1)*Args.GridSteps)+1:ii*Args.GridSteps );
        end
        
        
        ind_fields = mapadsmGrid > threshrate;
        % Find separate fields
        fieldlabel = bwlabel(ind_fields,8);
        fieldcount = max(max(fieldlabel));
        
        % For each field, look at spatial view for the 8 px around max px
        for ii = 1:fieldcount
            % Find local peak rate
            localpeakrate = max(max(mapadsmGrid(fieldlabel == ii)));
            ind_localpeak = mapadsmGrid == localpeakrate;
            xlocalpeak = find(sum(ind_localpeak,2) > 0);
            ylocalpeak = find(sum(ind_localpeak,1) > 0);
            % Find field boundary 
            xfieldpx = xlocalpeak-1:xlocalpeak+1;
            yfieldpx = ylocalpeak-1:ylocalpeak+1;
            % Initialize storage for view rate maps for each of 9 place px
            ratemaps = cell(3,3);
            % For each of 9 bins within boundary
            for xx = 1:size(xfieldpx,2)
                for yy = 1:size(yfieldpx,2)
                    placeBin = binNumGrid(xfieldpx(xx),yfieldpx(yy));
                    % Find which samples match this place bin
                    ts = sv.data.binLocLin == placeBin;
                    % Find viewed locations from this place bin
                    if sum(ts) > 0
                        viewedLin = sv.data.binGazeLin(ts);
                        viewedGrid = sv.data.binGazeGrid(ts,:);
                        [uviewedLin,i] = unique(viewedLin);
                        spikes = sv.data.spikepersample(ts);
                        fixObjNum = sv.data.fixObjNum(ts);
                        % Compute spike rates per view
                        rateview = nan(size(uviewedLin,1),1);
                        fixObj = nan(size(uviewedLin,1),1);
                        for kk = 1:size(uviewedLin,1)
                            ind = viewedLin == uviewedLin(kk); % s
                            rateview(kk) = sum(spikes(ind))/sum(ind);
                            fixObj(kk) = fixObjNum(find(ind,1));
                        end
                        % Convert lin rate map to grid
                        mapGrid = cell(size(binDepths,1),1);
                        for bb = 1:size(binDepths,1)
                            map = nan(binDepths(bb,1),binDepths(bb,2));
                            ind_obj = fixObj == bb;
                            if sum(ind_obj) > 0
                                binLin = uviewedLin(ind_obj);
                                binGrid = viewedGrid(i(ind_obj),:);
                                rates = rateview(ind_obj);
                                for cc = 1:size(rates,1)
                                   map(binGrid(cc,1),binGrid(cc,2)) = rates(cc);
                                end
                            end
                            mapGrid{bb} = map;
                        end
                    end
                    % Store for plotting later
                    ratemaps{xx,yy} = mapGrid;
                end
            end
            
            cwd = pwd;
            plotratemaps('place',0,5,0,'adaptive',cwd,ratemaps,binDepths,ii,[xlocalpeak ylocalpeak]);

            
        end
        
        
        
        
        
        
        
        
        
    end
    
    % Find max 
    

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

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.gridSteps = [];
data.meanFRs = [];
data.semFRs = [];
data.SIC = [];
data.SICsh = [];

% create nptdata so we can inherit from it
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);