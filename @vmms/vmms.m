function [obj, varargout] = vmms(varargin)

% Mixed selectivity
%   OBJ = vmms(varargin)
%
%   OBJ = vmms('auto') attempts to create a DIRFILES object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on vmms %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Example [as, Args] = vmms('save','redo')
%
%   Dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, 'ObjectLevel', 'Cell',...
                'RequiredFile','spiketrain.mat', 'GridSteps',40,...
                'AdaptiveSmooth',1, 'UseCorr',1,'FieldThr',0.7,'FieldSplitThr',0.8);
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
Args.classname = 'vmms';
Args.matname = [Args.classname num2str(Args.UseCorr) '.mat'];
Args.matvarname = 'vmms';

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
msobj = {'place','spatialview'};
Args.msobj = msobj;

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
    pSIthr = prctile(c_pc.SICsh,95); % Population threshold only
    vSIthr = prctile(c_sv.SICsh,95); % Population threshold only
%     pSIthr = max([prctile(c_pc.SICsh,95) prctile(pc.SICsh,95)]); % Both population and cell threshold
%     vSIthr = max([prctile(c_sv.SICsh,95) prctile(sv.SICsh,95)]); % Both population and cell threshold
    cd(cwd);
    
    data.mixsel = false;
    data.placesel = false;
    data.spatialviewsel = false;
    if pSI>pSIthr || vSI>vSIthr 
        % Record selectivity
        if pSI>pSIthr && vSI>vSIthr 
            data.mixsel = true; % Mixed selective
            data.placesel = true;
            data.spatialviewsel = true;
        elseif pSI>pSIthr && vSI<vSIthr
            data.placesel = true; % Place selective only
        elseif pSI<pSIthr && vSI>vSIthr
            data.spatialviewsel = true; % View selective only
        end
        
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
        
        % Filter
        conditions = ones(size(pvT,1),1);
        if pc.Args.UseAllTrials == 0
            conditions = conditions & pv.data.good_trial_markers;
        end
        if pc.Args.ThresVel > 0
            conditions = conditions & get(pv,'SpeedLimit',pc.Args.ThresVel); % pv here needs to be in object form, not structure
%             pvT(~get(pv,'SpeedLimit',pc.Args.ThresVel),:) = []; % pv here needs to be in object form, not structure
        end
        if pc.Args.UseMinObs
            bins_sieved_p = pv.data.place_good_bins;
            bins_removed_p = setdiff(1:size(pv.data.place_intervals_count,1),bins_sieved_p);
            bins_sieved_sv = pv.data.view_good_bins;
            bins_removed_sv = setdiff(1:size(pv.data.view_intervals_count,1),bins_sieved_sv);
            conditions = conditions & (pv.data.pv_good_rows); % Make sure maps take into account both place and view filters
%             pvT(~ismember(pvT(:,2),pv.data.place_good_bins),:) = [];
%             pvT(~ismember(pvT(:,3),pv.data.view_good_bins),:) = [];
        else
            bins_sieved_p = 1:(pc.Args.GridSteps * pc.Args.GridSteps);
            bins_removed_p = [];
            bins_sieved_sv = 1:size(pv.data.view_intervals_count,1);
            bins_removed_sv = [];
        end
        pvT(conditions ~= 1,:) = []; % Remove filtered rows
        pvT(pvT(:,2)==0,:) = []; % Remove remaining rows without place bins
        
        % Backfill duration and spikes for view bins - following code from
        % vmcorr
        pvTfill = nan(size(pvT));
        view_durations = zeros(5122,1600);
        view_spikes = zeros(5122,1600);
        place_durations = zeros(1,1600);
        place_spikes = zeros(1,1600);
        for i = 1:1600
            inds = pvT(:,2)==i;
            indnums = find(inds);
            subsample = [pvT(inds,:)]; % [time place view dur spk]
            
            % Consider only samples where both place and view are sampled
%             if any(isnan(subsample(:,3)))
%                 disp(find(isnan(subsample(:,3))));
%             end
            indnums(isnan(subsample(:,3)),:) = [];
            subsample(isnan(subsample(:,3)),:) = [];
            if isempty(subsample)
                continue;
            end

            % Get spikes and duration for place only
            place_durations(1,i) = sum(subsample(:,4));
            place_spikes(1,i) = sum(subsample(:,5));

            % back-filling spikes for view
            subsample(subsample(:,5)==0,5) = nan;
            subsample(:,6) = subsample(:,4)~=0;
            subsample(isnan(subsample(:,5)) & subsample(:,6), 5) = 0;
            subsample(:,6) = [];
            subsample(:,5) = fillmissing(subsample(:,5), 'next');
            % back-filling time for view
            subsample(subsample(:,4)==0,4) = nan;
            subsample(:,4) = fillmissing(subsample(:,4), 'next');
            
            % remove bad view spots
            indnums(isnan(subsample(:,3)),:) = [];
            subsample(isnan(subsample(:,3)),:) = [];
            
            % Put backfill into sessionTimeC array
            pvTfill(indnums,:) = subsample;

            % padding with 5122 bin
            subsample = [subsample; [0 0 5122 0 0]];
            
            % sum durations
            view_durations(:,i) = accumarray(subsample(:,3), subsample(:,4),[],[],NaN);
            % sum spikes
            view_spikes(:,i) = accumarray(subsample(:,3), subsample(:,5),[],[],NaN); 
        end
        pvTfill(isnan(pvTfill(:,1)),:) = []; % Remove all rows without both place and view data
        % Remove low obs bins
        pvTfill(~ismember(pvTfill(:,2),bins_sieved_p) | ~ismember(pvTfill(:,3),bins_sieved_sv),:) = [];
        place_durations(isnan(place_durations)) = 0;  
        view_durations(isnan(view_durations)) = 0; 
        place_durations(bins_removed_p) = 0;
        view_durations(bins_removed_sv,:) = 0;
        view_durations(:,bins_removed_p) = 0;
        place_spikes(bins_removed_p) = 0;
        view_spikes(bins_removed_sv,:) = 0;
        view_spikes(:,bins_removed_p) = 0;
        
        % This is only for verifying similarity with pc/sv - Make maps from filtered pv object
        p_map = place_spikes./place_durations;
        v_map = nansum(view_spikes,2)./nansum(view_durations,2);
        v_map(nansum(view_durations,2)==0) = nan; % restore nans to unvisited bins
        data.place.pvmap = p_map;
        data.spatialview.pvmap = v_map';
        
        % Find fields in both place and view maps
        for oo = 1:size(msobj,2)
            switch msobj{oo}
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
                    gridSize{oo} = [pc.Args.GridSteps pc.Args.GridSteps];
                    basemapGsm = lineartogrid(basemapLsm','place',gridSize{oo});
                    dummygrid = lineartogrid((1:size(basemapLsm,2))','place',gridSize{oo});
                    basemapGrw = lineartogrid(basemapLrw','place',gridSize{oo});
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
                    gridSize{oo} = sv.binDepths;
                    basemapGsm = lineartogrid(basemapLsm','spatialview',gridSize{oo});
                    dummygrid = lineartogrid((1:size(basemapLsm,2))','spatialview',gridSize{oo});
                    basemapGrw = lineartogrid(basemapLrw','spatialview',gridSize{oo});
                    for ii = 1:size(basemapGsm,1)
                        maxset(ii) = nanmax(reshape(basemapGsm{ii},size(basemapGsm{ii},1)*size(basemapGsm{ii},2),1));
                    end
                    [peakrate_full,prI] = max(maxset); % Max including cue/hint
                    peakrate_set = max(maxset(3:end)); % Max excluding cue/hint
            end
            
            % Find 3 maxima of rate maps
            count = 0;
            maxfieldcount = 3;
            discardfieldnum = 0;
            discardfieldreason = [];
            if data.([msobj{oo} 'sel'])
                for gg = 1:size(basemapGsm,1)
                    % Skip cue and hint view maps 
                    if size(basemapGsm{gg},1) == 1
                        continue;
                    end
                    % Find fields with at least 1 pixel of > 70% peak rate
                    ind_fields = basemapGsm{gg} > Args.FieldThr*peakrate_set;
                    % Find separate fields
                    [fieldlabel,fieldcount] = bwlabel(ind_fields,4); % Only count adjacent pixels if they share edge, not if they share corners
                    % For walls and pillars, if there are fields that wrap around split, merge them.
                    if gg > 4
                       % Find possible split fields
                       if any(fieldlabel(:,1)) && any(fieldlabel(:,end))
                           boundary = [fieldlabel(:,1) fieldlabel(:,end)];
                           [blabel,bcount] = bwlabel(boundary,4);
                           for bb = 1:bcount
                               label = blabel == bb;
                               label = unique(boundary(label));
                               labeltochange = label(2:end);
                               if ~isempty(labeltochange)
                                   for ll = 1:size(labeltochange,1)
                                       % Merge fields around corners
                                       fieldlabel(fieldlabel == labeltochange(ll)) = label(1); 
%                                        fieldcount = fieldcount - 1;
                                   end
                               end
%                                label = blabel(:,1) == bb;
%                                label = unique(boundary(label,1));
%                                labeltochange = unique(fieldlabel(blabel(:,end)==bb,end));
%                                for ll = 1:size(labeltochange,1)
%                                    fieldlabel(fieldlabel == labeltochange(ll)) = label; 
%                                    fieldcount = fieldcount - 1;
%                                end
                           end
                       end
                    end
                    % Split fields that are too large i.e. > 1/4 of usable area
                    for ii = 1:fieldcount
                        inds = fieldlabel == ii;
                        % Split fields that are too large
                        if sum(inds(:)) >= ((size(pc.maps_raw,2)-4*8*8)/4) 
                            subinds = inds & basemapGsm{gg} > Args.FieldSplitThr*peakrate_set;
                            [sublabel,subcount] = bwlabel(subinds,4);
                            sublabel(subinds) = sublabel(subinds) + fieldcount;
                            fieldcount = fieldcount + subcount;
                            fieldlabel(inds) = 0; % There will be some px that are not part of new fields cos of different rate thresholds in subfields
                            fieldlabel(inds) = sublabel(inds);
                        end
                    end
                    % Save fields that are big enough (i.e. > 15 bins around peak, 5x5bins area)
                    for ii = 1:fieldcount
                        inds = fieldlabel == ii;
                        % Make sure field size is > 15 bins 
                        if sum(inds(:)) >= 15
                            % Discard field if there are not at least ONE active pixel within this field in the raw map
                            if sum(basemapGrw{gg}(inds)>0) >= 1
                                count = count + 1;
                                fieldmaxrate_rw(count,1) = max(basemapGrw{gg}(inds)); 
                                fieldmaxrate_sm(count,1) = max(basemapGsm{gg}(inds));
                                switch msobj{oo}
                                    case 'place'
                                        gridnum(count,1) = 3;
                                    case 'spatialview'
                                        gridnum(count,1) = gg;
                                end
                                [fieldcoordx_mat, fieldcoordy_mat] = find(inds);
                                fieldcoord{count,1} = [fieldcoordy_mat size(basemapGrw{gg},1)-fieldcoordx_mat+1]; % In plot coords. x left to right, y bottom to top
                                linbin{count,1} = dummygrid{gg}(inds);
                            else 
                                disp(['discarding field: no spikes, ' num2str(num2str(sum(inds(:)))) 'px']);
                                discardfieldnum = discardfieldnum + 1;
                                discardfieldreason(1,end+1) = sum(inds(:));
                            end
                        elseif sum(inds(:)) == 0 % Field either merged or split before
                            continue;
                        else
                            disp(['discarding field: too small, ' num2str(sum(inds(:))) 'px']);
                            discardfieldnum = discardfieldnum + 1;
                            discardfieldreason(1,end+1) = sum(inds(:));
                        end
                    end
                end
            end

            % If no significant fields, skip to next var
            if count == 0
                % create nptdata so we can inherit from it
                data.(msobj{oo}).sigfields = 0; % Selective but no sig fields (if cell non-selective, sigfields = nan, see below)
                data.(msobj{oo}).discardfieldnum = discardfieldnum;
                data.(msobj{oo}).discardfieldreason = discardfieldreason;
                data.(msobj{oo}).SI = SI;
                data.(msobj{oo}).SIthr = SIthr;
                data.(msobj{oo}).basemapLsm = basemapLsm;
                data.(msobj{oo}).secmapLsm = secmapLsm;
                data.(msobj{oo}).basemapGsm = basemapGsm;
                data.(msobj{oo}).basemapGrw = basemapGrw;
                data.(msobj{oo}).dummygrid = dummygrid;
                data.(msobj{oo}).peakrate_full = [];
                data.(msobj{oo}).peakrate_set = [];
                data.(msobj{oo}).fieldmaxrate_sm = [];
                data.(msobj{oo}).fieldmaxrate_rw = [];
                data.(msobj{oo}).gridnum = [];
                data.(msobj{oo}).fieldcoord = {};
                data.(msobj{oo}).linbin = {};
                data.(msobj{oo}).set_sec_linbin = {};
                data.(msobj{oo}).rate_components = {};
                data.(msobj{oo}).secfieldrates = {};
                data.(msobj{oo}).secfieldrates_sh = {};
                continue;
            else 
                data.(msobj{oo}).sigfields = min([count maxfieldcount]);
                data.(msobj{oo}).discardfieldnum = discardfieldnum;
                data.(msobj{oo}).discardfieldreason = discardfieldreason;
            end

            % Sort fields
            [fieldmaxrate_sm,I] = sort(fieldmaxrate_sm,'descend');
            fieldmaxrate_rw = fieldmaxrate_rw(I);
            gridnum = gridnum(I);
            fieldcoord = fieldcoord(I);
            linbin = linbin(I);

            % Get secondary pixels for each base field 
            for ii = 1:maxfieldcount % Limit to first 3 fields per base map

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
                    switch msobj{oo}
                        case 'place'
                            ind_pv = pvTfill(:,2) == linbin{ii}(pp);
                            sec(:,1) = pvTfill(ind_pv,3); % view px
                            sec(:,2) = pvTfill(ind_pv,4); % dur
                            sec(:,3) = pvTfill(ind_pv,5); % spikes
                        case 'spatialview'
                            ind_pv = pvTfill(:,3) == linbin{ii}(pp);
                            sec(:,1) = pvTfill(ind_pv,2); % place px
                            sec(:,2) = pvTfill(ind_pv,4); % dur
                            sec(:,3) = pvTfill(ind_pv,5); % spikes
                    end
                    % Get firing rates 
                    usecpx = unique(sec(:,1));
                    usecpx(isnan(usecpx)) = [];
                    rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
                    rate_components_px(:,1) = usecpx;
                    if any(isnan(usecpx))
                        disp(nan);
                    end
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
                    rate_components_px(rightfulnans,2) = nan;
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
%                 disp(ii);

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
            data.(msobj{oo}).SI = SI;
            data.(msobj{oo}).SIthr = SIthr;
            data.(msobj{oo}).basemapLsm = basemapLsm;
            data.(msobj{oo}).secmapLsm = secmapLsm;
            data.(msobj{oo}).basemapGsm = basemapGsm;
            data.(msobj{oo}).basemapGrw = basemapGrw;
            data.(msobj{oo}).dummygrid = dummygrid;
            data.(msobj{oo}).peakrate_full = peakrate_full;
            data.(msobj{oo}).peakrate_set = peakrate_set;
            data.(msobj{oo}).fieldmaxrate_sm = fieldmaxrate_sm;
            data.(msobj{oo}).fieldmaxrate_rw = fieldmaxrate_rw;
            data.(msobj{oo}).gridnum = gridnum;
            data.(msobj{oo}).fieldcoord = fieldcoord;
            data.(msobj{oo}).linbin = linbin;
            data.(msobj{oo}).set_sec_linbin = seclinbin_full;
            data.(msobj{oo}).rate_components = rate_components_full;

            clear fieldmaxrate_sm; clear fieldmaxrate_rw; clear gridnum; clear fieldcoord; clear linbin; clear session_seclinbin;
            clear seclinbin_full; % clear secgridmap_full; 
            clear rate_components_full; clear dummygrid; clear SI; clear SIthr;
        
        end
        
        % Stats
        for oo = 1:size(msobj,2)
            % Output variables
            secfieldrates = cell(size(data.(msobj{oo}).rate_components,1),1);
            secfieldrates_sh = cell(size(data.(msobj{oo}).rate_components,1),1);
            for ii = 1:size(data.(msobj{oo}).rate_components,1) % For each base field
                % If there are secondary fields
                if ~isempty(data.(msobj{2-oo+1}).rate_components)
                    % Output variables
                    secfieldrate = nan(size(data.(msobj{2-oo+1}).rate_components,1),1);
                    secfieldrate_sh = cell(size(data.(msobj{2-oo+1}).rate_components,1),1);
                    % Secondary pixels/map sampled for the entire session
                    session_seclinbin = data.(msobj{oo}).set_sec_linbin{ii};
                    dummygridsec = data.(msobj{2-oo+1}).dummygrid;
                    session_seclinmap = nan(size(data.(msobj{2-oo+1}).basemapLsm));
                    session_seclinmap(1,session_seclinbin(:,1)) = session_seclinbin(:,4);
                    session_secgridmap = lineartogrid(session_seclinmap',msobj{2-oo+1},gridSize{2-oo+1});

                    % Test if secondary pixels sampled from this base field are more likely to fall within any of the secondary fields than outside
                    for jj = 1:size(data.(msobj{2-oo+1}).rate_components,1) % For each secondary field

                        % Get mean firing rate within secondary field
                        linbin = data.(msobj{2-oo+1}).linbin{jj}; % pixels that make up the sec field. Not all of these will be sampled from this base field
                        seclinbin_sampled = linbin(ismember(linbin,session_seclinbin(:,1)));
                        inds_infield = ismember(session_seclinbin(:,1),seclinbin_sampled); % find pixels of sec field that are sampled from this base field
                        if sum(inds_infield) == 0 
                            continue;
                        end
                        meanrate_infield = sum(session_seclinbin(inds_infield,3)) / sum(session_seclinbin(inds_infield,2));

                        % Get mean firing rates for 10000 pseudorandom same-size fields outside of secondary field
                        meanrate_outfield = nan(10000,1);
                        session_seclinbin_out = session_seclinbin(~inds_infield,:);
                        if size(session_seclinbin_out,1) > sum(inds_infield) % Filter 1 of 2: Make sure there are enough outfield px to generate stats (This catches situations where there are just too few outfield px to begin with) 
                            % Generate a psuedopopulation of outfields same size as sec field 
                            ff = 1;
                            attempt = 0; % Filter 2 of 2: Make sure there are enough outfield px to generate stats (This catches situations where outfield px are enough in number but scattered so that can't form coherent shuffled field without overlapping with original sec field)
                            abandon = false;
                            while ff <= 10000 && ~abandon
                                attempt = attempt + 1;
                                % Start from a random pixel that is outside of sec field
                                startpx = randsample(1:size(session_seclinbin_out,1),1);
                                startpx = session_seclinbin_out(startpx,1);
                                if ismember(startpx,linbin) % If random px overlaps with sec field, repeat
    %                                 reset = true;
                                    continue;
                                end
                                % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
                                switch msobj{2-oo+1}
                                    case 'place'
                                        gnum = 1;
                                    case 'spatialview'
                                        if startpx == 1 || startpx == 2 % Make sure not cue or hint
    %                                        reset = true;
                                           continue;
                                        end
                                        [gnum,~,~] = findgrid(startpx,msobj{2-oo+1});
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
                                    if sum(sum(~isnan(tempmap(startindx,startindy)))) > sum(inds_infield)% length(startindx)*length(startindy) > size(linbin,1)
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
                                % sampled pixels
                                sampledpx = dummygridsec{gnum}(startindx,startindy);
                                sampledpx = sampledpx(~isnan(tempmap(startindx,startindy)));
                                % get another starting pixel if cannot sample enough px / if sampled field overlaps too much with sec field 
                                if length(startindx)*length(startindy) < sum(inds_infield) || sum(sum(~isnan(tempmap(startindx,startindy)))) < sum(inds_infield)
                                    continue;
                                elseif size(intersect(sampledpx,linbin),1) > 0.25*size(linbin,1) % If sampled field overlaps with more than half of sec field
                                    if attempt == 10000 && ff < 10
                                        abandon = true;
                                    end
                                    continue;
                                end
                                % Remove empty pixels from sampled field
                                pxsub = dummygridsec{gnum}(startindx,startindy);
                                inds_sampled = ~isnan(tempmap(startindx,startindy));
                                pxsub = pxsub(inds_sampled);
                                if sum(inds_infield)>length(pxsub)
                                    error('Not enough sample pixels to draw pseudopopulation from');
                                end
                                inds_keep = sort(randsample(1:length(pxsub),sum(inds_infield)))';
                                pxsub = pxsub(inds_keep);
                                % Start over if exactly the same px as field of interest
                                if isempty(setdiff(pxsub,linbin))
    %                                 reset = true;
                                    continue;
                                end
                                % Get sec bins actually sampled 
                                outfieldlinbin = pxsub;
                                inds_outfield = ismember(session_seclinbin(:,1),outfieldlinbin);
                                meanrate_outfield(ff,1) = sum(session_seclinbin(inds_outfield,3)) / sum(session_seclinbin(inds_outfield,2));
%                                 if mod(ff,1000) == 0
%                                     disp(ff);
%                                 end
                                ff = ff+1;
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
%                 else
%                     % Store data
%                     secfieldrates{ii,1} = NaN;
%                     secfieldrates_sh{ii,1} = nan(10000,1);
                else 
                    if data.mixsel
                        disp('no sec fields');
                    end
                end
            end
            % Store data
            data.(msobj{oo}).secfieldrates = secfieldrates;
            data.(msobj{oo}).secfieldrates_sh = secfieldrates_sh;
        end
    else % If cell is not selective for both place or view
        for oo = 1:size(msobj,2)
            data.(msobj{oo}).pvmap = [];
            data.(msobj{oo}).sigfields = NaN; % Not selective
            data.(msobj{oo}).discardfieldnum = 0;
            data.(msobj{oo}).discardfieldreason = {};
            data.(msobj{oo}).SI = [];
            data.(msobj{oo}).SIthr = [];
            data.(msobj{oo}).basemapLsm = [];
            data.(msobj{oo}).secmapLsm = [];
            data.(msobj{oo}).basemapGsm = {};
            data.(msobj{oo}).basemapGrw = {};
            data.(msobj{oo}).dummygrid = {};
            data.(msobj{oo}).peakrate_full = [];
            data.(msobj{oo}).peakrate_set = [];
            data.(msobj{oo}).fieldmaxrate_sm = [];
            data.(msobj{oo}).fieldmaxrate_rw = [];
            data.(msobj{oo}).gridnum = [];
            data.(msobj{oo}).fieldcoord = {};
            data.(msobj{oo}).linbin = {};
            data.(msobj{oo}).set_sec_linbin = {};
            data.(msobj{oo}).rate_components = {};
            data.(msobj{oo}).secfieldrates = {};
            data.(msobj{oo}).secfieldrates_sh = {};
        end
    end
    
	% create nptdata so we can inherit from it
	data.numSets = 1;    
    data.Args = Args;
    data.origin = {cwd};
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
