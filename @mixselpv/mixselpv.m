function [obj, varargout] = mixselpv(varargin)
%@dirfiles Constructor function for DIRFILES class
%   OBJ = mixselpv(varargin)
%
%   OBJ = mixselpv('auto') attempts to create a DIRFILES object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on mixselpv %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Example [as, Args] = mixselpv('save','redo')
%
%   Dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, 'ObjectLevel', 'Cell',...
                'RequiredFile','spiketrain.mat', 'GridSteps',40,...
                'AdaptiveSmooth',1, 'UseCorr',1);
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
Args.classname = 'mixselpv';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'mspv';

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

% example object
dlist = nptDir;
% get entries in directory
dnum = size(dlist,1);
% Variables in conjunction
objtype = {'place','spatialview'};

% check if the right conditions were met to create object
if(dnum>0)
	% these are object specific fields
	data.dlist = dlist;
	% set index to keep track of which data goes with which directory
	data.setIndex = [0; dnum];
    cwd = pwd;
    
    % Determine if cell is selective for both variables
    cd('FiltVel/1px');
    pc = load('vmpc.mat');
    pc = pc.vmp.data;
    sv = load('vmsv.mat');
    sv = sv.vms.data;
    corr = load('vmcorr.mat');
    corr = corr.vmcorr.data;
    if ~Args.UseCorr % If using original pc/sv objects to define fields
        pSI = pc.SIC;
        vSI = sv.SIC;
    else
        pSI = corr.SIC_corrp;
        vSI = corr.SIC_corrsv;
    end
    c_pc = load('/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Combined Objects/FiltVel/1px/c_vmpc.mat');
    c_pc = c_pc.vmp.data;
    c_sv = load('/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Combined Objects/FiltVel/1px/c_vmsv.mat');
    c_sv = c_sv.vms.data;
    pSIthr = max([prctile(c_pc.SICsh,95) prctile(pc.SICsh,95)]); % Both population and cell threshold
    vSIthr = max([prctile(c_sv.SICsh,95) prctile(sv.SICsh,95)]);
    cd(cwd);
    
    if pSI>pSIthr && vSI>vSIthr 
        data.mixsel = true;
        % Filter pv with same criteria used in vmpc/sv
        cd ..; cd ..; cd ..;
        pv = load('1vmpv.mat');
        pv = pv.pv;
        pvT = pv.data.sessionTimeC;
        pvT(:,4) = [diff(pvT(:,1)); 0]; % dwell time
        % Bin spikes
        cd(cwd);
        spiketrain = load('spiketrain.mat');
        spiketimes = spiketrain.timestamps ./ 1000; % in seconds
        binned = histcounts(spiketimes, pvT(:,1))';
        pvT(:,5) = [binned; 0];
        if pc.Args.ThresVel
            pvT(~get(pv,'SpeedLimit',pc.Args.ThresVel),:) = []; % pv here needs to be in object form, not structure
        end
        if pc.Args.UseMinObs
            pvT(~ismember(pvT(:,2),pv.data.place_good_bins),:) = [];
            pvT(~ismember(pvT(:,3),pv.data.view_good_bins),:) = [];
        end
        pvT(pvT(:,2)==0,:) = []; % Remove remaining rows without place bins
        
        % Backfill duration and spikes for view bins 
        pvTfill = pvT;
        view_durations = NaN(5122,1600);
        view_spikes = NaN(5122,1600);
        place_durations = NaN(1,1600);
        place_spikes = NaN(1,1600);
        for i = 1:1600
            inds = pvT(:,2)==i;
            subsample = [pvT(inds, [3 4 5])];
            if ~isempty(subsample)
                disp(i);
            end

            % Get spikes and duration for place only
            place_durations(1,i) = sum(subsample(:,2));
            place_spikes(1,i) = sum(subsample(:,3));

            % back-filling spikes for view
            subsample(subsample(:,3)==0,3) = nan;
            subsample(:,4) = circshift(subsample(:,2)~=0 ,-1);
            subsample(isnan(subsample(:,3)) & subsample(:,4), 3) = 0;
            subsample(:,4) = [];
            subsample(:,3) = fillmissing(subsample(:,3), 'next');
            % back-filling time for view
            subsample(subsample(:,2)==0,2) = nan;
            subsample(:,2) = fillmissing(subsample(:,2), 'previous');
            % Put backfill into sessionTimeC array
            pvTfill(inds,[3 4 5]) = subsample;

            % padding with 5122 bin
            subsample = [subsample; [5122 0 0]];
            % remove bad view spots
            subsample(isnan(subsample(:,1)),:) = [];
            % sum durations
            view_durations(:,i) = accumarray(subsample(:,1), subsample(:,2),[],[],NaN);
            % sum spikes
            view_spikes(:,i) = accumarray(subsample(:,1), subsample(:,3),[],[],NaN); 
        end
        % This is only for verifying similarity with pc/sv - Make maps from filtered pv object
        p_map = place_spikes./place_durations;
        v_map = nansum(view_spikes,2)./nansum(view_durations,2);
        v_map(nansum(view_durations,2)==0) = nan; % restore nans to unvisited bins
        data.place.p_map = p_map;
        data.spatialview.v_map = v_map;
        
        % Find fields in both place and view maps
        for oo = 1:size(objtype,2)
            switch objtype{oo}
                case 'place'
                    if Args.UseCorr
                       basemapLsm = corr.maps_adsm_corrp;
                       secmapLsm = corr.maps_adsm_corrsv;
                       basemapLrw = corr.maps_raw_corrp;
                    else
                        basemapLsm = pc.maps_adsm;
                        secmapLsm = sv.maps_adsm;
                        basemapLrw = pc.maps_raw;
                    end
                    SI = pSI;
                    SIthr = pSIthr;
                    [basemapGsm,dummygrid] = lineartogrid(basemapLsm,'place');
                    [basemapGrw,~] = lineartogrid(basemapLrw,'place');
                    peakrate_full = nanmax(basemapLsm);
                    peakrate_set = peakrate_full;
                    prI = 1;
                case 'spatialview'
                    if Args.UseCorr
                       basemapLsm = corr.maps_adsm_corrsv;
                       secmapLsm = corr.maps_adsm_corrp;
                       basemapLrw = corr.maps_raw_corrsv;
                    else
                        basemapLsm = sv.maps_adsm;
                        secmapLsm = pc.maps_adsm;
                        basemapLrw = sv.maps_raw;
                    end
                    SI = vSI;
                    SIthr = vSIthr;
                    [basemapGsm,dummygrid] = lineartogrid(basemapLsm,'spatialview');
                    [basemapGrw,~] = lineartogrid(basemapLrw,'spatialview');
                    for ii = 1:size(basemapGsm,1)
                        maxset(ii) = nanmax(reshape(basemapGsm{ii},size(basemapGsm{ii},1)*size(basemapGsm{ii},2),1));
                    end
                    [peakrate_full,prI] = max(maxset);
                    peakrate_set = max(maxset(3:end));
            end
            
            % Find 3 maxima of rate maps
            count = 0;
            for gg = 1:size(basemapGsm,1)
                % Skip cue and hint view maps 
                if size(basemapGsm{gg},1) == 1
                    continue;
                end
                % Find fields with at least 1 pixel of > 70% peak rate
                ind_fields = basemapGsm{gg} > 0.7*peakrate_full;
                % Find separate fields
                [fieldlabel,fieldcount] = bwlabel(ind_fields,4);
                % For walls and pillars, if there are fields that wrap around
                % split, merge them.
                if gg > 4
                   % Find possible split fields
                   if any(fieldlabel(:,1)) && any(fieldlabel(:,end))
                       boundary = [fieldlabel(:,1) fieldlabel(:,end)];
                       [blabel,bcount] = bwlabel(boundary,8);
                       for bb = 1:bcount
                           label = unique(fieldlabel(blabel(:,1)==bcount(bb),1));
                           label = label(1);
                           labelset = unique(boundary);
                           labelset = labelset(labelset>0);
                           labelset = setdiff(labelset,label);
                           for ll = 1:length(labelset)
                              fieldlabel(fieldlabel == labelset(ll)) = label; 
                              fieldcount = fieldcount - 1;
                           end
                       end
                   end
                end
                % Split fields that are too large i.e. > 1/4 of area
                for ii = 1:fieldcount
                    inds = fieldlabel == ii;
                    % Split fields that are too large
                    if sum(inds(:)) >= (size(basemapLsm,2)/4) 
                        subinds = inds & basemapGsm{gg} > 0.8*peakrate_full;
                        [sublabel,subcount] = bwlabel(subinds,4);
                        sublabel(subinds) = sublabel(subinds) + fieldcount;
                        fieldcount = fieldcount + subcount;
                        fieldlabel(inds) = 0;
                        fieldlabel(inds) = sublabel(inds);
                    end
                end
                % Find fields that are big enough (i.e. > 15 bins around peak, 5x5bins area)
                for ii = 1:fieldcount
                    inds = fieldlabel == ii;
                    % Make sure field size is > 15 bins 
                    if sum(inds(:)) >= 15
                        count = count + 1;
                        fieldmaxrate_rw(count,1) = max(basemapGrw{gg}(inds)); 
                        fieldmaxrate_sm(count,1) = max(basemapGsm{gg}(inds));
                        switch objtype{oo}
                            case 'place'
                                gridnum(count,1) = 3;
                            case 'spatialview'
                                gridnum(count,1) = gg;
                        end
                        [fieldcoordx_mat, fieldcoordy_mat] = find(inds);
                        fieldcoord{count,1} = [fieldcoordy_mat size(basemapGrw{gg},1)-fieldcoordx_mat+1]; % In plot coords. x left to right, y bottom to top
                        seclinbin{count,1} = dummygrid{gg}(inds);
                    end
                end
            end

            % If no significant fields, skip to next cell
            if count == 0
                % create nptdata so we can inherit from it
                data.(objtype{oo}).sigfields = 0;
                data.numSets = 1;    
                data.Args = Args;
                n = nptdata(data.numSets,0,pwd);
                d.data = data;
                obj = class(d,Args.classname,n);
                saveObject(obj,'ArgsC',Args);
                return;
            else 
                data.(objtype{oo}).sigfields = count;
            end

            % Sort fields
            [fieldmaxrate_sm,I] = sort(fieldmaxrate_sm,'descend');
            fieldmaxrate_rw = fieldmaxrate_rw(I);
            gridnum = gridnum(I);
            fieldcoord = fieldcoord(I);
            seclinbin = seclinbin(I);

            % Get secondary pixels for each base field 
            for ii = 1:3 % Limit to first 3 fields per base map

                if ii > size(fieldcoord,1)
                    break;
                end
                % Initialise variables
                rate_components_field = cell(size(fieldcoord{ii},1),1);
                maxrate_px = nan(size(fieldcoord{ii},1),1);
                secs = cell(size(fieldcoord{ii},1),1);
                usecpxs = [];
                % Find pixels of secondary map
                for pp = 1:size(fieldcoord{ii},1)
                    %Get corresponding secondary pixels from pv object
                    sec = [];
                    switch objtype{oo}
                        case 'place'
                            ind_pv = pvT(:,2) == seclinbin{ii}(pp);
                            sec(:,1) = pvT(ind_pv,3); % view px
                            sec(:,2) = pvT(ind_pv,4); % dur
                            sec(:,3) = pvT(ind_pv,5); % spikes
                        case 'spatialview'
                            ind_pv = pvTfill(:,3) == seclinbin{ii}(pp);
                            sec(:,1) = pvTfill(ind_pv,2); % place px
                            sec(:,2) = pvTfill(ind_pv,4); % dur
                            sec(:,3) = pvTfill(ind_pv,5); % spikes
                    end
                    % Get firing rates 
                    usecpx = unique(sec(:,1));
                    usecpx(isnan(usecpx)) = [];
                    rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
                    rate_components_px(:,1) = usecpx;
                    ind_timechange = find(sec(:,2) > 0);
                    for cc = 1:size(ind_timechange,1) % For each instance of being in this base pixel
                        linbintemp = nan(size(usecpx,1),2); % Temp dur and spikes for this secondary pixel(s)
                        if cc == size(ind_timechange,1)
                            newind = ind_timechange(cc):size(sec,1);
                        else
                            newind = ind_timechange(cc):ind_timechange(cc+1)-1; % index into secondary pixel(s) for this instance
                        end
                        newview = sec(newind,1); 
                        newset = ismember(usecpx,newview);
                        linbintemp(newset,1) = sec(newind(1),2); % duration for this instance listed with first secondary pixel
                        linbintemp(newset,2) = sec(newind(end),3); % spikes for this instance listed with last secondary pixel
                        rate_components_px(:,2) = nansum( [rate_components_px(:,2) linbintemp(:,1)] ,2); % Sum duration for this sec pixel across instances
                        rate_components_px(:,3) = nansum( [rate_components_px(:,3) linbintemp(:,2)] ,2); % Sum spikes for this sec pixel across instances
                    end
                    rightfulnans = rate_components_px(:,2) == 0;
                    rate_components_px(rightfulnans,:) = nan;
                    rate_components_px(:,4) = rate_components_px(:,3)./rate_components_px(:,2); % Compute firing rates

                    % Collect all sec rate information for the base pixels
                    rate_components_field{pp} = ( rate_components_px );
                    if ~isempty(rate_components_px)
                        maxrate_px(pp) = max(rate_components_px(:,4));
                    else 
                        maxrate_px(pp) = NaN;
                    end
                    secs{pp} = sec;
                    usecpxs = union(usecpxs,usecpx); % set of secondary pixels covered in whole session

                end
                disp(ii);

                session_seclinbin = usecpxs; 
                session_seclinbin(:,2:3) = NaN;
                % Consolidate occupancy and spikes across base field
                for pp = 1:size(fieldcoord{ii},1)
                    tempbins = nan(size(session_seclinbin,1),2);
                    tempbins(ismember(session_seclinbin(:,1),rate_components_field{pp}(:,1)),1) = rate_components_field{pp}(:,2); % Duration 
                    tempbins(ismember(session_seclinbin(:,1),rate_components_field{pp}(:,1)),2) = rate_components_field{pp}(:,3); % Spikes 
                    session_seclinbin(:,2) = nansum( [session_seclinbin(:,2) tempbins(:,1)] ,2); % Sum duration across session
                    session_seclinbin(:,3) = nansum( [session_seclinbin(:,3) tempbins(:,2)] ,2); % Sum spikes across session
                end
                session_seclinbin(:,4) = session_seclinbin(:,3)./session_seclinbin(:,2);

                rate_components_full{ii,1} = rate_components_field;
                seclinbin_full{ii,1} = session_seclinbin;
            end

            % Store data
            data.(objtype{oo}).SI = SI;
            data.(objtype{oo}).SIthr = SIthr;
            data.(objtype{oo}).basemapLsm = basemapLsm;
            data.(objtype{oo}).secmapLsm = secmapLsm;
            data.(objtype{oo}).basemapGsm = basemapGsm;
            data.(objtype{oo}).basemapGrw = basemapGrw;
            data.(objtype{oo}).dummygrid = dummygrid;
            data.(objtype{oo}).peakrate_full = peakrate_full;
            data.(objtype{oo}).peakrate_set = peakrate_set;
            data.(objtype{oo}).fieldmaxrate_sm = fieldmaxrate_sm;
            data.(objtype{oo}).fieldmaxrate_rw = fieldmaxrate_rw;
            data.(objtype{oo}).gridnum = gridnum;
            data.(objtype{oo}).fieldcoord = fieldcoord;
            data.(objtype{oo}).linbin = seclinbin;
            data.(objtype{oo}).set_sec_linbin = seclinbin_full;
            data.(objtype{oo}).rate_components = rate_components_full;

            clear fieldmaxrate_sm; clear fieldmaxrate_rw; clear gridnum; clear fieldcoord; clear linbin; clear session_seclinbin;
            clear seclinbin_full; % clear secgridmap_full; 
            clear rate_components_full; clear dummygrid; clear SI; clear SIthr;
        
        end
        
        % Stats
        for oo = 1:size(objtype,2)
            % Output variables
            secfieldrates = cell(size(data.(objtype{oo}).rate_components,1),1);
            secfieldrates_sh = cell(size(data.(objtype{oo}).rate_components,1),1);
            for ii = 1:size(data.(objtype{oo}).rate_components,1) % For each base field
                % Output variables
                secfieldrate = nan(size(data.(objtype{2-oo+1}).rate_components,1),1);
                secfieldrate_sh = cell(size(data.(objtype{2-oo+1}).rate_components,1),1);
                % Secondary pixels/map sampled for the entire session
                session_seclinbin = data.(objtype{oo}).set_sec_linbin{ii};
                dummygridsec = data.(objtype{2-oo+1}).dummygrid;
                session_seclinmap = nan(size(data.(objtype{2-oo+1}).basemapLsm));
                session_seclinmap(1,session_seclinbin(:,1)) = session_seclinbin(:,4);
                [session_secgridmap,~] = lineartogrid(session_seclinmap,objtype{2-oo+1});
                
                % Test if secondary pixels sampled from this base field are more likely to fall within any of the secondary fields than outside
                for jj = 1:size(data.(objtype{2-oo+1}).rate_components,1) % For each secondary field

                    % Get mean firing rate within secondary field
                    seclinbin = data.(objtype{2-oo+1}).linbin{jj}; % pixels that make up the sec field. Not all of these will be sampled from this base field
                    seclinbin_sampled = seclinbin(ismember(seclinbin,session_seclinbin(:,1)));
                    inds_infield = ismember(session_seclinbin(:,1),seclinbin); % find pixels of sec field that are sampled from this base field
                    if sum(inds_infield) == 0 
                        continue;
                    end
                    meanrate_infield = sum(session_seclinbin(inds_infield,3)) / sum(session_seclinbin(inds_infield,2));
                    
                    % Get mean firing rates for 10000 pseudorandom same-size fields outside of secondary field
                    meanrate_outfield = nan(10000,1);
                    session_seclinbin_out = session_seclinbin(~inds_infield,:);
                    if ~isempty(session_seclinbin_out) % If all pixels are within sec field ###### FIXXXXXXXXXX ????????
                        % Generate a psuedopopulation of outfields same size as sec field 
                        for ff = 1:10000
                            % Start from a random pixel that is outside of sec field
                            startpx = randsample(1:size(session_seclinbin_out,1),1);
                            startpx = session_seclinbin_out(startpx,1);
                            if ismember(startpx,seclinbin) % If random px overlaps with sec field, repeat
                                ff = ff - 1;
                                continue;
                            end
                            % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
                            switch objtype{2-oo+1}
                                case 'place'
                                    gnum = 1;
                                case 'spatialview'
                                    if startpx == 1 || startpx == 2 % Make sure not cue or hint
                                       ff = ff - 1;
                                       continue;
                                    end
                                    [gnum,~,~] = findgrid(startpx,objtype{2-oo+1});
                            end
                            % Get the grid coords of starting px
                            [startindx,startindy] = find(dummygridsec{gnum} == startpx);
                            tempmap = session_secgridmap{gnum}; % the actual sampled sec grid map
                            % Expand radius around starting px until hit the requisite number of px 
                            % while sum(sum(~isnan(tempmap(startindx,startindy)))) < sum(inds_infield)
                            while length(startindx)*length(startindy) < size(dummygridsec{gnum},1)*size(dummygridsec{gnum},2)
                                if startindx(1) > 1 && startindx(end) < size(dummygridsec{gnum},1)
                                    startindx = [startindx(1)-1 startindx startindx(end)+1];
                                elseif startindx(1) == 1 && startindx(end) < size(dummygridsec{gnum},1)
                                    startindx = [startindx startindx(end)+1];
                                elseif startindx(1) > 1 && startindx(end) == size(dummygridsec{gnum},1)
                                    startindx = [startindx(1)-1 startindx];
                                end
                                % If reach required num of px
                                if sum(sum(~isnan(tempmap(startindx,startindy)))) > sum(inds_infield)% length(startindx)*length(startindy) > size(seclinbin,1)
                                    break;
                                end
                                if startindy(1) > 1 && startindy(end) < size(dummygridsec{gnum},2)
                                    startindy = [startindy(1)-1 startindy startindy(end)+1];
                                elseif startindy(1) == 1 && startindy(end) < size(dummygridsec{gnum},2)
                                    startindy = [startindy startindy(end)+1];
                                elseif startindy(1) > 1 && startindy(end) == size(dummygridsec{gnum},2)
                                    startindy = [startindy(1)-1 startindy];
                                end
                                % If exceed map bounds
                                if startindx(1) == 1 && startindx(end) == size(dummygridsec{gnum},1) && startindy(1) == 1 && startindy(end) == size(dummygridsec{gnum},2)
                                    break;
                                end
                            end
                            % get another starting pixel if cannot sample enough px 
                            if length(startindx)*length(startindy) < sum(inds_infield) || sum(sum(~isnan(tempmap(startindx,startindy)))) < sum(inds_infield)
                                ff = ff - 1;
                                continue;
                            end
                            % Remove empty pixels from sampled field
                            pxsub = dummygridsec{gnum}(startindx,startindy);
                            inds_sampled = ~isnan(tempmap(startindx,startindy));
                            pxsub = pxsub(inds_sampled);
                            if sum(inds_infield)>length(pxsub)
                                error('Not enough sample pixels to draw pseudopopulation from');;
                            end
                            inds_keep = sort(randsample(1:length(pxsub),sum(inds_infield)))';
                            pxsub = pxsub(inds_keep);
                            % Start over if exactly the same px as field of interest
                            if isempty(setdiff(pxsub,seclinbin))
                                ff = ff - 1;
                                continue;
                            end
                            % Get sec bins actually sampled 
                            outfieldlinbin = pxsub;
                            inds_outfield = ismember(session_seclinbin(:,1),outfieldlinbin);
                            meanrate_outfield(ff,1) = sum(session_seclinbin(inds_outfield,3)) / sum(session_seclinbin(inds_outfield,2));

                        end
%                         meanrate_thr = prctile(meanrate_outfield(:),95);
                    end
                    % Store data
                    secfieldrate(jj,1) = meanrate_infield;
                    secfieldrate_sh{jj,1} = meanrate_outfield;
                end
                % Store data
                secfieldrates{ii,1} = secfieldrate;
                secfieldrates_sh{ii,1} = secfieldrate_sh;
            end
            % Store data
            data.(objtype{oo}).secfieldrates = secfieldrates;
            data.(objtype{oo}).secfieldrates_sh = secfieldrates_sh;
        end
    else % If cell is not selective for both place and view
        data.mixsel = false;
        data.place.sigfields = NaN;
        data.place.SI = pSI;
        data.place.SIthr = pSIthr;
        data.spatialview.sigfields = NaN;
        data.spatialview.SI = vSI;
        data.spatialview.SIthr = vSIthr;
        
    end
    
	% create nptdata so we can inherit from it
	data.numSets = 1;    
    data.Args = Args;
	n = nptdata(data.numSets,0,pwd);
	d.data = data;
	obj = class(d,Args.classname,n);
	saveObject(obj,'ArgsC',Args);
else
	% create empty object
	obj = createEmptyObject(Args);
end

function obj = createEmptyObject(Args)

% these are object specific fields
data.dlist = [];
data.setIndex = [];

% create nptdata so we can inherit from it
% useful fields for most objects
data.numSets = 0;
data.Args = Args;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);

function [mapG,mapGdummy]= lineartogrid(mapL,objtype)

% Insert floor place map into larger 3D view setting
if strcmp(objtype,'place')
    mapLtemp = mapL;
    mapL = nan(1,5122);
    mapL(3:3+1600-1) = mapLtemp;
    mapG = {flipud(reshape(mapLtemp, 40, 40)')};
    mapGdummy = {flipud(reshape(1:1600, 40, 40)')};
end

mapLdummy = 1:length(mapL);
% Set up surf frame for plotting
floor_x = repmat(0:40, 41, 1);
floor_y = flipud(repmat([0:40]', 1, 41));
floor_z = zeros(41,41);

ceiling_x = floor_x;
ceiling_y = floor_y;
ceiling_z = 40.*ones(41,41);

walls_x = repmat([0.*ones(1,40) 0:39 40.*ones(1,40) 40:-1:0], 9, 1);
walls_y = repmat([0:39 40.*ones(1,40) 40:-1:1 0.*ones(1,41)], 9, 1);
walls_z = repmat([24:-1:16]', 1, 40*4 + 1);

P1_x = repmat([24.*ones(1,8) 24:31 32.*ones(1,8) 32:-1:24], 6, 1);
P1_y = repmat([8:15 16.*ones(1,8) 16:-1:9 8.*ones(1,9)], 6, 1);
PX_z = repmat([21:-1:16]', 1, 8*4 + 1);

P2_x = repmat([8.*ones(1,8) 8:15 16.*ones(1,8) 16:-1:8], 6, 1);
P2_y = P1_y;

P3_x = P1_x;
P3_y = repmat([24:31 32.*ones(1,8) 32:-1:25 24.*ones(1,9)], 6, 1);

P4_x = P2_x;
P4_y = P3_y;

floor = flipud(reshape(mapL(3:3+1600-1), 40, 40)');
floordum = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');

% ceiling follows floor mapping, top down view
ceiling = flipud(reshape(mapL(1603:1603+1600-1), 40, 40)');
ceilingdum = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');

% from top down, slit walls at bottom left corner, open outwards.
% start from row closest to ground, rightwards, then climb rows
walls = flipud(reshape(mapL(3203:3203+1280-1), 40*4, 8)');
wallsdum = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');

% BL - bottom left, and so on, from top view, same slicing as walls
% pillar width 8, height 5
P1_BR = flipud(reshape(mapL(4483:4483+160-1), 8*4, 5)');
P1_BRdum = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
P2_BL = flipud(reshape(mapL(4643:4643+160-1), 8*4, 5)');
P2_BLdum = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
P3_TR = flipud(reshape(mapL(4803:4803+160-1), 8*4, 5)');
P3_TRdum = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
P4_TL = flipud(reshape(mapL(4963:4963+160-1), 8*4, 5)');
P4_TLdum = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
if strcmp(objtype,'spatialview')
    mapG = { NaN; NaN; floor; ceiling; walls; P1_BR; P2_BL; P3_TR; P4_TL };
    mapGdummy = { NaN; NaN; floordum; ceilingdum; wallsdum; P1_BRdum; P2_BLdum; P3_TRdum; P4_TLdum };
end

% Find grid in spatial view frame for 1 pixel
function [gridnum,x,y] = findgrid(px,objtype)
% returns grid number and plot coords (x goes left to right, y goes bottom
% to top)

switch objtype
    case 'place'
        mapLdummy = 1:1600;
        gridnum = 3;
        temp = flipud(reshape(mapLdummy, 40, 40)');
    case 'spatialview'
        mapLdummy = 1:5122;
        if px == 1 % Cue
            gridnum = 1;
            x = 1;
            y = 1;
        elseif px == 2 % Hint
            gridnum = 2;
            x = 1; 
            y = 1;
        elseif px >= 3 && px <= 1602 % Floor
            gridnum = 3;
            temp = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');
        elseif px >= 1603 && px <= 3202 % Ceiling
            gridnum = 4;
            temp = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');
        elseif px >= 3203 && px <= 4482 % Walls
            gridnum = 5;
            temp = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');
        elseif px >= 4483 && px <= 4642 % Pillar 1
            gridnum = 6;
            temp = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
        elseif px >= 4643 && px <= 4802 % Pillar 2
            gridnum = 7;
            temp = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
        elseif px >= 4803 && px <= 4962 % Pillar 3
            gridnum = 8;
            temp = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
        elseif px >= 4963 && px <= 5122 % Pillar 4
            gridnum = 9;
            temp = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
        end
end
[y,x] = find(temp == px);
y = size(temp,1)-y+1;
