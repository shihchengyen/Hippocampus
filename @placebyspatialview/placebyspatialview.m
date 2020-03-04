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
				'ObjectLevel','Cell', 'RequiredFile','vmsv.mat', ...
				'MaxTimeDiff',0.002,'MinTrials',5, 'GridSteps',40, ...
                'ShuffleLimits',[0.1 0.9], 'NumShuffles',10000, ...
                'AdaptiveSmooth',1, 'FiltLowOcc',1);
Args.flags = {'Auto','ArgsOnly','HPC','FRSIC','UseMedian'};
% Specify which arguments should be checked when comparing saved objects
% to objects that are being asked for. Only arguments that affect the data
% saved in objects should be listed here.
Args.DataCheckArgs = {'MinTrials','GridSteps','ShuffleLimits', ...
	'NumShuffles'};

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
    cwd = pwd;
	sv = vmsv('auto');
    vp = vmpc('auto');
    cd ..; cd ..; cd ..;
    rc = load('raycast.mat');
    cd(cwd);
%     st = load('spiketrain.mat');
%     ufile = unityfile('auto',varargin{:});
    um = umaze('auto','GridSteps',40,varargin{:});
    cd(cwd);
    
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
    if vpSI > 0 % vpSIthresh
        
%         sessionTime = um.data.sessionTime(:,1);
%         sessionTime(um.data.zero_indices,1) = [];
        [~,~,binrc] = histcounts((rc.el.data.timestamps/1000),um.data.sessionTime(:,1));
        uniquepos = unique(binrc(binrc~=0));
        pospersample = zeros(size(rc.el.data.timestamps,1),1);
        for ii = 1:size(uniquepos,1)
           ind = binrc == uniquepos(ii);
           pospersample(ind) = um.data.sessionTime(uniquepos(ii),2);
        end
        
%         % Bin spike trains by location
%         data.horGridBound = um.data.horGridBound;
%         data.vertGridBound = um.data.vertGridBound;
%         sTime = um.data.sessionTime;
%         sortedGPindinfo = um.data.sortedGPindinfo;
%         unityTime = ufile.data.unityTime;
% 
%         % load spike train
%         ts1 = (st.timestamps/1000)';
%         % remove any spike times before Unity started and after Unity ended
%         ts = ts1(ts1<sTime(end,1) & ts1>unityTime(1));
% 
%         % find spikes that occurred during navigation
%         % e.g.
%     %	row		sessionTime			ts		bins	zero_indices
%     %	1	0		0		7.8217	6.7649	1	1
%     %	2	7.8217	3.0000	3.6656	6.7917	1	6
%     %	3	11.4873	4.0000	1.5644	8.1015	2	14
%     %	4	13.0517	5.0000	0.3790	9.3976	2	20
%     %	5	13.4307	10.0000	1.1994	16.4411	6	26
%     %	6	14.6301	0		3.1062	16.4855	6	34
%     %	7	17.7363	10.0000	2.4343	16.5127	6	40
%     %	8	20.1706	15.0000	0.4862	16.5400	6	46
%     %	9	20.6568	14.0000	1.7571	16.8922	6	54
%     %	10	22.4139	13.0000	0.8260	16.9383	6	61
%         % value in bins match up with the values in zero_indices, so we can just find the values
%         % in bins that are not members of zero_indices
%         [~,~,bins] = histcounts(ts,sTime(:,1));
%         % compute the repetition number within each grid position each spike occurs in
%         % first find the row in sortedGPindices that corresponds to the occurence of the grid
%         % position that the spike occurred in
%         [~, spike_index] = ismember(bins,um.data.sortedGPindices);
%         % get indices that correspond to non-zero grid positions
%         n0bins = spike_index>sortedGPindinfo(1,2);
%         % get the grid position for each spike in non-zero grid positions
%         spike_gridposition = sTime(bins(n0bins),2);
%         % compute the repetition number for each spike by finding the row that corresponds to
%         % the start of the appropriate grid position
%         rep_num = mod(spike_index(n0bins),sortedGPindinfo(um.data.gp2ind(spike_gridposition),2)) + 1;
% 
%         % perform histogram on spike times that don't correspond to zero bins using unity time 
%         % bins in order to get the closest xy position
%         [~,~,ubin] = histcounts(ts(n0bins),ufile.data.unityTime);
%         spike_xy = ufile.data.unityData(ubin,3:4);
        
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
            binNumGrid(Args.GridSteps-(ii-1),:) = binNumLin( ((ii-1)*Args.GridSteps)+1:ii*Args.GridSteps );
        end
        
        ind_fields = mapadsmGrid > threshrate;
        % Find separate fields
        fieldlabel = bwlabel(ind_fields,8);
        fieldcount = max(max(fieldlabel));
        
%         % Sanity check for map orientation
%         binLocGrid = sv.data.binLocGrid;
%         locmapGrid = nan(vp.data.Args.GridSteps);
%         ratemapGrid = nan(vp.data.Args.GridSteps);
%         for pp = 1:vp.data.Args.GridSteps
%             for ppp = 1:vp.data.Args.GridSteps
%                 inds = binLocGrid(:,1)==pp & binLocGrid(:,2)==ppp;
%                 spikesGrid = sum(sv.data.spikepersample(inds));
%                 locmapGrid(pp,ppp) = sum(inds);
%                 ratemapGrid(pp,ppp) = spikesGrid/locmapGrid(pp,ppp);
%             end
%         end
%         locmapGrid = locmapGrid/(max(max(locmapGrid)));
%         locmapGrid(locmapGrid == 0) = nan;
%         figure(888);
%         ax = gca;
%         colormap(jet);
%         surfx = repmat((0:40)',1,41);
%         surfy = repmat(0:40,41,1);
%         surfz = zeros(41);
%         surf(surfx,surfy,surfz,ratemapGrid);
%         shading flat;
%         set(ax,'CLim',[0 max(max(ratemapGrid))]);
%         view(-60,90);
        
        % For each field, look at spatial view for the 8 px around max px
        
        for ii = 1:fieldcount
            % Find local peak rate
            localpeakrate = max(max(mapadsmGrid(fieldlabel == ii)));
            ind_localpeak = mapadsmGrid == localpeakrate;
            xlocalpeak = find(sum(ind_localpeak,2) > 0,1);
            ylocalpeak = find(sum(ind_localpeak,1) > 0,1);
            % Find field boundary 
            xfieldpx = xlocalpeak-1:xlocalpeak+1;
            yfieldpx = ylocalpeak-1:ylocalpeak+1;
            % Initialize storage for view rate maps for each of 9 place px
            ratemapsG = cell(3,3);
            ratemapsL = nan(sum(binDepths(:,1).*binDepths(:,2)),9);
            % For each of 9 bins within boundary
            for xx = 1:size(xfieldpx,2)
                for yy = 1:size(yfieldpx,2)
                    placeBin = binNumGrid(xfieldpx(xx),yfieldpx(yy));
                    % Find which samples match this place bin
%                     tt = sv.data.binLocLin == placeBin;
                    tt = pospersample == placeBin;
                    % Find viewed locations from this place bin
                    if sum(tt) > 0
%                         tts;
                        viewedLin = sv.data.binGazeLin(tt);
                        viewedGrid = sv.data.binGazeGrid(tt,:);
                        [uviewedLin,i] = unique(viewedLin);
                        uviewedLin = uviewedLin(~isnan(uviewedLin));
                        i = i(~isnan(uviewedLin));
                        spikes = sv.data.spikepersample(tt);
                        fixObjNum = sv.data.fixObjNum(tt);
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
                        mapL = nan(sum(binDepths(:,1).*binDepths(:,2)),1);
                        for bb = 1:size(binDepths,1)
                            mapG = nan(binDepths(bb,1),binDepths(bb,2));
                            ind_obj = fixObj == bb;
                            if sum(ind_obj) > 0
                                binLin = uviewedLin(ind_obj);
                                binGrid = viewedGrid(i(ind_obj),:);
                                rates = rateview(ind_obj);
                                for cc = 1:size(rates,1)
                                   mapG(binGrid(cc,1),binGrid(cc,2)) = rates(cc);
                                end
                                mapL(binLin) = rates; %%%??????
                            end
                            mapGrid{bb} = mapG;
                        end
                        
                    % Store for plotting later
                    ratemapsG{xx,yy} = mapGrid;
                    ratemapsL(:, yy+(xx-1)*yy) = mapL;
                    
                    h1 = figure(ii);
                    ax = gca;
                    % Plot main map
                    plotratemaps('Place',0,0,'adaptive',cwd,mapadsmGrid,ax,mapadsmLin,binDepths,ii,[40-xlocalpeak ylocalpeak]);
                    h2 = figure(ii+100);
                    % Plot spatialview maps per pixel
                    ax = subplot(3,3,(xx-1)*3+yy);
                    plotratemaps('Spatialview',0,0,'adaptive',cwd,mapGrid,ax,mapL,binDepths,ii,[40-xlocalpeak ylocalpeak]);
                    end
                end
            end
            % 
            
%             cwd = pwd;
%             h = figure(1);
%             for xx = 1:3
%                 for yy = 1:3
%                     ax = subplot(3,3,(xx-1)*3+yy);
%                     plotratemaps('Place',0,0,'adaptive',cwd,ratemapsG,ax,ratemapsL,binDepths,ii,[40-xlocalpeak ylocalpeak]);
%                     plotratemaps('Spatialview',0,0,'adaptive',cwd,ratemapsG{xx,yy},ax,ratemapsL((xx-1)*3+yy),binDepths,ii,[40-xlocalpeak ylocalpeak]);
%                 end
%             end

            
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