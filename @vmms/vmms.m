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
                'RequiredFile','spiketrain.mat', 'GridSteps',40, 'pix', 1, ...
                'UseCorr',1,'FieldThr',0.6,'FieldSplitThr',0.7,'NumShuffles',1000, ...
                'FieldThrPseudo',0.6,'FieldSplitThrPseudo',0.7,'UseAllTrials',1,'ThresVel',1,'UseMinObs',0);
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

% check if the right conditions were met to create object
if(dnum>0)
    
    data.dlist = dlist;
    % set index to keep track of which data goes with which directory
    data.setIndex = [0; dnum];
    cwd = pwd;

    % Variables in conjunction
    spatialvarpairs = {{'place','view'},{'place','headdirection'},{'headdirection','view'}};
    msobjpairs = {{'pc','sv'},{'pc','hd'},{'hd','sv'}};
    Args.spatialvarpairs = spatialvarpairs;
    
    % Load corrected object
    cd(['FiltVel/' num2str(Args.pix) 'px']);
    cr = load('vmcorr.mat');
    cr = cr.vmcorr.data;
    pc = load('vmpc.mat');
    pc = pc.vmp.data;
    sv = load('vmsv.mat');
    sv = sv.vms.data;
    hd = load('vmhd.mat');
    hd = hd.vmd.data;
    cd ..; cd ..;

    if ~pc.discard && ~sv.discard && ~hd.discard
        data.discard = false;
    else 
        data.discard = true;
    end
    
    % Set up output structure
    for pair = 1:size(spatialvarpairs,2)

        msvar = spatialvarpairs{pair}; % For now, do in pairs only
        msobj = msobjpairs{pair};
        msvar_short = {msvar{1}(1),msvar{2}(1)};
        pairname_short = [msvar{1}(1) msvar{2}(1)];
    
        %% Set up output variables
        for oo = 1:size(msvar,2)

            %% Set up all output variables

            % Base var 
            data.(pairname_short).(msvar{oo}).varname = msvar{oo};
            data.(pairname_short).(msvar{oo}).varname_short = msvar_short{oo};
            data.(pairname_short).(msvar{oo}).pvmap = [];
            data.(pairname_short).(msvar{oo}).sigfields = NaN; % Not selective
            data.(pairname_short).(msvar{oo}).discardfieldnum = 0;
            data.(pairname_short).(msvar{oo}).discardfieldreason = {};
            data.(pairname_short).(msvar{oo}).SI = [];
            data.(pairname_short).(msvar{oo}).SIthr = [];
            data.(pairname_short).(msvar{oo}).basemapLsm = [];
            data.(pairname_short).(msvar{oo}).basemapLrw = [];
            data.(pairname_short).(msvar{oo}).basemapGsm = {};
            data.(pairname_short).(msvar{oo}).basemapGrw = {};
            data.(pairname_short).(msvar{oo}).dummygrid = {};
            data.(pairname_short).(msvar{oo}).peakrate_full = [];
            data.(pairname_short).(msvar{oo}).peakrate_set = [];
            data.(pairname_short).(msvar{oo}).fieldmaxrate_sm = [];
            data.(pairname_short).(msvar{oo}).fieldmaxrate_rw = [];
            data.(pairname_short).(msvar{oo}).gridnum = [];
            data.(pairname_short).(msvar{oo}).fieldcoord = {};
            data.(pairname_short).(msvar{oo}).fieldlinbin = {};
            data.(pairname_short).(msvar{oo}).basemaps_dist = {};
            data.(pairname_short).(msvar{oo}).base_distratio = [];
            data.(pairname_short).(msvar{oo}).base_distcorr = [];
            % Sec vars
            data.(pairname_short).(msvar{oo}).secmapLsm = [];
            data.(pairname_short).(msvar{oo}).condbase_map_rw = {}; 
            data.(pairname_short).(msvar{oo}).condbase_componentsperpx = {};
            data.(pairname_short).(msvar{oo}).secmaps_dist = {};
            data.(pairname_short).(msvar{oo}).sec_distratio = [];
            data.(pairname_short).(msvar{oo}).sec_distcorr = [];
            data.(pairname_short).(msvar{oo}).condbase_insecfieldrates = {};
            data.(pairname_short).(msvar{oo}).condbase_outsecfieldrates = {};
            data.(pairname_short).(msvar{oo}).condbase_linkedfield = {};
            data.(pairname_short).(msvar{oo}).condbase_map_sm = {};
        %         data.(pairname_short).(msvar{oo}).condbase_map_adsm = {};
        %         data.(pairname_short).(msvar{oo}).condbase_map_bcsm = {};
        %         data.(pairname_short).(msvar{oo}).condbase_map_dksm = {};
            data.(pairname_short).(msvar{oo}).condbase_SIC_sm = [];
        %         data.(pairname_short).(msvar{oo}).condbase_SIC_adsm = [];
        %         data.(pairname_short).(msvar{oo}).condbase_SIC_bcsm = [];
        %         data.(pairname_short).(msvar{oo}).condbase_SIC_dksm = [];
            data.(pairname_short).(msvar{oo}).condbase_sigfields = []; 
            data.(pairname_short).(msvar{oo}).condbase_gridnum = [];
            data.(pairname_short).(msvar{oo}).condbase_fieldcoord = {};
            data.(pairname_short).(msvar{oo}).condbase_fieldlinbin = {};
            data.(pairname_short).(msvar{oo}).condbase_fieldsizepercent = {};
            data.(pairname_short).(msvar{oo}).condbase_fieldoverlapind = {};
            data.(pairname_short).(msvar{oo}).condbase_fieldoverlap = {};
            data.(pairname_short).(msvar{oo}).pseudosecmaps_sm = {};
            data.(pairname_short).(msvar{oo}).pseudosecSIC_adsm = {};
            data.(pairname_short).(msvar{oo}).condpseudo_secdataperfield = {};
            data.(pairname_short).(msvar{oo}).secfieldnumbinset = {};
            data.(pairname_short).(msvar{oo}).condfield_inbinsset = {};
            data.(pairname_short).(msvar{oo}).condfield_nonoverlapbinset = {};
            % data.(pairname_short).(msvar{oo}).tertrate_orig = {};
            % data.(pairname_short).(msvar{oo}).tertrate_pseudo = {};
        end
    end
        
    % For all data, regardless of selectivity
    if ~data.discard
        
        %% Work out mixed selective properties wuth each spatialvar as base
        for pair = 1:size(spatialvarpairs,2)
            
            msvar = spatialvarpairs{pair}; % For now, do in pairs only
            msobj = msobjpairs{pair};
            msvar_short = {msvar{1}(1),msvar{2}(1)};
            pairname_short = [msvar{1}(1) msvar{2}(1)];
            for oo = 1:size(msvar,2)

                    %% Find base fields 
                    stcfilt = eval(msobj{oo}).stcfilt;
                    stcvars = {'timestamp','place','headdirection','view','duration','spikes'};
                    xbase = strcmp(stcvars,msvar{oo});
                    xsec = strcmp(stcvars,msvar{2-oo+1});
                    disp(['Delineating ' msvar{oo} ' fields ...']);
                    
                    baseobj = eval(msobj{oo});
                    secobj = eval(msobj{2-oo+1});
                    if Args.UseCorr
                        basemapLsm = cr.([msvar_short{1} msvar_short{2}]).(['maps_sm_corr' msvar_short{oo}]);
                        secmapLsm = cr.([msvar_short{1} msvar_short{2}]).(['maps_sm_corr' msvar_short{2-oo+1}]);
                        basemapLrw = cr.([msvar_short{1} msvar_short{2}]).(['maps_raw_corr' msvar_short{oo}]);
                    else
                        basemapLsm = baseobj.(['maps_sm']);
                        secmapLsm = secobj.(['maps_sm']);
                        basemapLrw = baseobj.maps_raw;
                    end
                    crit = eval(msobj{oo}).crit_sm;    
                    gridSize{oo} = cr.([msvar_short{1} msvar_short{2}]).([msvar{oo} 'binDepths']);
                    basemapGsm = lineartogrid(basemapLsm',msvar{oo},gridSize{oo});
                    dummygrid = lineartogrid((1:size(basemapLsm,2))',msvar{oo},gridSize{oo});
                    basemapGrw = lineartogrid(basemapLrw',msvar{oo},gridSize{oo});

                    switch msvar{oo}
                        case 'place'
                            peakrate_full = nanmax(basemapLsm);
                            peakrate_subset = peakrate_full;
                            prI = 1;
                            fieldsizethreshold = 9;
                            % fieldsizethreshold = 15;
                        case 'view'
                            for ii = 1:size(basemapGsm,1)
                                maxset(ii) = nanmax(reshape(basemapGsm{ii},size(basemapGsm{ii},1)*size(basemapGsm{ii},2),1));
                            end
                            [peakrate_full,prI] = max(maxset); % Max including cue/hint
                            peakrate_subset = max(maxset(3:end)); % Max excluding cue/hint
                            fieldsizethreshold = 9;
                            % fieldsizethreshold = 15;
                        case 'headdirection'
                            peakrate_full = nanmax(basemapLsm);
                            peakrate_subset = peakrate_full;
                            prI = 1;
                            fieldsizethreshold = 5;
                    end

                    % Find 3 maxima of rate maps
                    count = 0;
                    maxfieldcount = 3;
                    discardfieldnum = 0;
                    discardfieldreason = [];
                    fieldlinbin = {};
                    fieldcoord = {};
                    fieldmaxrate_rw = [];
                    fieldmaxrate_sm = [];
                    gridnum = [];

                    for gg = 1:size(basemapGsm,1)
                        % Skip cue and hint view maps 
                        if size(basemapGsm{gg},1) == 1 && size(basemapGsm{gg},2) == 1
                            continue;
                        end
                        % if there is too little variance in firing rates, reject
                        if Args.FieldThr*peakrate_subset < nanmean(basemapLsm) % + nanstd(basemapLsm) % nanmax(basemapLsm) < nanmean(basemapLsm) + 2*nanstd(basemapLsm)
                            continue;
                        end

                        % Find fields with at least 1 pixel of > 70% peak
                        % rate and > 1 std dev from mean
                        ind_fields = basemapGsm{gg} > Args.FieldThr*peakrate_subset & ...
                            basemapGsm{gg} > (nanmean(basemapLsm)+nanstd(basemapLsm));
                        % Find separate fields
%                             [fieldlabel,fieldcount] = bwlabel(ind_fields,4); % Only count adjacent pixels if they share edge, not if they share corners
                        [fieldlabel,fieldcount] = bwlabel(ind_fields,8); % Count adjacent pixels even if just they share corners
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
                                       end
                                   end
                               end
                           end
                        end
                        % Split fields that are too large i.e. > 1/4 of usable area
                        for ii = 1:fieldcount
                            inds = fieldlabel == ii;
                            % Split fields that are too large
                            switch msvar{oo}
                                case 'place'
                                    if sum(inds(:)) >= ((size(pc.maps_raw,2)-4*8*8)/4) % 1/4 of occupiable space
                                        subinds = inds & basemapGsm{gg} > Args.FieldSplitThr*peakrate_subset;
                                        [sublabel,subcount] = bwlabel(subinds,4);
                                        sublabel(subinds) = sublabel(subinds) + fieldcount;
                                        fieldcount = fieldcount + subcount;
                                        fieldlabel(inds) = 0; % There will be some px that are not part of new fields cos of different rate thresholds in subfields
                                        fieldlabel(inds) = sublabel(inds);
                                    end
                                case 'view' 
                                     if sum(inds(:)) >= ((size(sv.maps_raw,2)-4*8*8)/4) % 1/4 of occupiable space
                                        subinds = inds & basemapGsm{gg} > Args.FieldSplitThr*peakrate_subset;
                                        [sublabel,subcount] = bwlabel(subinds,4);
                                        sublabel(subinds) = sublabel(subinds) + fieldcount;
                                        fieldcount = fieldcount + subcount;
                                        fieldlabel(inds) = 0; % There will be some px that are not part of new fields cos of different rate thresholds in subfields
                                        fieldlabel(inds) = sublabel(inds);
                                     end
                                case 'headdirection'
                                     if sum(inds(:)) >= ((size(hd.maps_raw,2))/4) % 1/4 of occupiable space
                                        subinds = inds & basemapGsm{gg} > Args.FieldSplitThr*peakrate_subset;
                                        [sublabel,subcount] = bwlabel(subinds,4);
                                        sublabel(subinds) = sublabel(subinds) + fieldcount;
                                        fieldcount = fieldcount + subcount;
                                        fieldlabel(inds) = 0; % There will be some px that are not part of new fields cos of different rate thresholds in subfields
                                        fieldlabel(inds) = sublabel(inds);
                                     end
                            end

                        end
                        % Save fields that are big enough (i.e. > 15 bins around peak, 5x5bins area)
                        for ii = 1:fieldcount
                            inds = fieldlabel == ii;
                            % Make sure field size is big enough 
                            if sum(inds(:)) >= fieldsizethreshold
                                % Discard field if there are not at least ONE active pixel within this field in the raw map
                                if sum(basemapGrw{gg}(inds)>0) >= 1
                                    count = count + 1;
                                    fieldmaxrate_rw(count,1) = max(basemapGrw{gg}(inds)); 
                                    fieldmaxrate_sm(count,1) = max(basemapGsm{gg}(inds));
                                    switch msvar{oo}
                                        case 'place'
                                            gridnum(count,1) = 1;
                                        case 'view'
                                            gridnum(count,1) = gg;
                                        case 'headdirection'
                                            gridnum(count,1) = 1;
                                    end
                                    [fieldcoordx_mat, fieldcoordy_mat] = find(inds);
                                    fieldcoord{count,1} = [fieldcoordy_mat size(basemapGrw{gg},1)-fieldcoordx_mat+1]; % In plot coords. x left to right, y bottom to top
                                    fieldlinbin{count,1} = dummygrid{gg}(inds);
                                else 
%                                     disp(['discarding field: no spikes, ' num2str(num2str(sum(inds(:)))) 'px']);
                                    discardfieldnum = discardfieldnum + 1;
                                    discardfieldreason(1,end+1) = sum(inds(:));
                                end
                            elseif sum(inds(:)) == 0 % Field either merged or split before
                                continue;
                            else
%                                 disp(['discarding field: too small, ' num2str(sum(inds(:))) 'px']);
                                discardfieldnum = discardfieldnum + 1;
                                discardfieldreason(1,end+1) = sum(inds(:));
                            end
                        end
                    end

                    % If no significant fields, skip to next var
                    if count == 0
                        % create nptdata so we can inherit from it
                        data.(pairname_short).(msvar{oo}).sigfields = 0; % Selective but no sig fields (if cell non-selective, sigfields = nan, see below)
                        data.(pairname_short).(msvar{oo}).discardfieldnum = discardfieldnum;
                        data.(pairname_short).(msvar{oo}).discardfieldreason = discardfieldreason;
                        data.(pairname_short).(msvar{oo}).SI = crit;
                        % data.(pairname_short).(msvar{oo}).SIthr = SIthr;
                        data.(pairname_short).(msvar{oo}).basemapLsm = basemapLsm;
                        data.(pairname_short).(msvar{oo}).basemapLrw = basemapLrw;
                        data.(pairname_short).(msvar{oo}).secmapLsm = secmapLsm;
                        data.(pairname_short).(msvar{oo}).basemapGsm = basemapGsm;
                        data.(pairname_short).(msvar{oo}).basemapGrw = basemapGrw;
                        data.(pairname_short).(msvar{oo}).dummygrid = dummygrid;
                        continue;
                    else 
                        data.(pairname_short).(msvar{oo}).sigfields = min([count maxfieldcount]);
                        data.(pairname_short).(msvar{oo}).discardfieldnum = discardfieldnum;
                        data.(pairname_short).(msvar{oo}).discardfieldreason = discardfieldreason;
                    end

                    % Sort fields by firing rate
                    [fieldmaxrate_sm,I] = sort(fieldmaxrate_sm,'descend');
                    fieldmaxrate_rw = fieldmaxrate_rw(I);
                    gridnum = gridnum(I);
                    fieldcoord = fieldcoord(I);
                    fieldlinbin = fieldlinbin(I);

                    %% Get conditioned maps for each base field 
                    disp(['Parsing conditioned pixels for each base ' msvar{oo} ' field ...']);

                    condbase_componentsperpx = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                    condbase_map_rw = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                    basemaps_dist = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                    secmaps_dist = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                    base_distratio = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
                    sec_distratio = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
                    base_distcorr = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
                    sec_distcorr = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
                    for ii = 1:maxfieldcount % Limit to first 3 fields per base map

                        if ii > size(fieldcoord,1)
                            break;
                        end
                        % Initialise variables
                        rate_components_perbasepx = cell(size(fieldcoord{ii},1),1);
                        maxrate = nan(size(fieldcoord{ii},1),1);
                        secpxs = cell(size(fieldcoord{ii},1),1);
                        usecpxs = [];

                        % Find pixels of secondary map
                        ind_pvh = ismember(stcfilt(:,xbase),fieldlinbin{ii});
                        secpx = stcfilt(ind_pvh,[find(strcmp(stcvars,msvar{2-oo+1})) 5 6]); % [bin, dur, spk]
                        usecpx = unique(secpx(:,1));
                        secpx(end+1,:) = [length(data.(pairname_short).(msvar{2-oo+1}).basemapLrw) 0 0]; % without max bin, accumarray doesn't work right
                        secdata = nan(length(data.(pairname_short).(msvar{2-oo+1}).basemapLrw),3);
                        secdata(usecpx,1) = usecpx;
                        secdata(:,2) = accumarray(secpx(:,1),secpx(:,2));
                        secdata(:,3) = accumarray(secpx(:,1),secpx(:,3));
                        secdata(:,4) = secdata(:,3)./secdata(:,2);


                        % Find pixels of secondary map
                        for pp = 1:size(fieldcoord{ii},1)
                            %Get corresponding secondary pixels from pv object
                            secpx = [];
                            ind_pvh = stcfilt(:,xbase) == fieldlinbin{ii}(pp); % base var px
                            secpx(:,1) = stcfilt(ind_pvh,xsec); % sec var px
                            secpx(:,2) = stcfilt(ind_pvh,5); % dur
                            secpx(:,3) = stcfilt(ind_pvh,6); % spikes
                            % Get firing rates 
                            usecpx = unique(secpx(:,1));
                            usecpx(isnan(usecpx)) = [];
                            rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
                            rate_components_px(:,1) = usecpx;
                            if any(isnan(usecpx))
                                disp(nan);
                            end
                            ind_timechange = find(secpx(:,2) > 0);
                            for cc = 1:size(ind_timechange,1) % For each sample of base pixel
                                linbintemp = nan(size(usecpx,1),2); % Temp dur and spikes for this secondary pixel(s)
                                if cc == size(ind_timechange,1)
                                    sampind = ind_timechange(cc):size(secpx,1);
                                else
                                    sampind = ind_timechange(cc):ind_timechange(cc+1)-1; % index into secondary pixel(s) for this instance
                                end
                                secsamp = secpx(sampind,1); % sec px from this base px instance
                                setind = ismember(usecpx,secsamp);
                                linbintemp(setind,1) = secpx(sampind(1),2); % duration for this instance listed with first secondary pixel
                                linbintemp(setind,2) = secpx(sampind(end),3); % spikes for this instance listed with last secondary pixel
                                rate_components_px(:,2) = nansum( [rate_components_px(:,2) linbintemp(:,1)] ,2); % Sum duration for this sec pixel across instances
                                rate_components_px(:,3) = nansum( [rate_components_px(:,3) linbintemp(:,2)] ,2); % Sum spikes for this sec pixel across instances
                            end
                            rightfulnans = rate_components_px(:,2) == 0;
                            rate_components_px(rightfulnans,2) = nan;
                            rate_components_px(:,4) = rate_components_px(:,3)./rate_components_px(:,2); % Compute firing rates

                            % Collect all sec rate information for the base pixels
                            rate_components_perbasepx{pp} = ( rate_components_px );
                            if ~isempty(rate_components_px)
                                maxrate(pp) = max(rate_components_px(:,4));
                            else 
                                maxrate(pp) = NaN;
                            end
                            secpxs{pp} = secpx;
                            usecpxs = union(usecpxs,usecpx); % set of secondary pixels covered in whole session

                        end

                        usecpxs = reshape(usecpxs,length(usecpxs),1);
                        condbase_rawdata = usecpxs; 
                        condbase_rawdata(:,2:3) = NaN;
                        full_dur = zeros(size(basemapLrw,2),size(secmapLsm,2));
                        full_spk = full_dur;
                        % Consolidate occupancy and spikes across base field
                        for pp = 1:size(fieldcoord{ii},1)
                            tempbins = nan(size(condbase_rawdata,1),2);
                            % debug
                            if fieldlinbin{ii}(pp) == 1300
                                disp(fieldlinbin{ii}(pp));
                            end
                            tempbins(ismember(condbase_rawdata(:,1),rate_components_perbasepx{pp}(:,1)),1) = rate_components_perbasepx{pp}(:,2); % Duration 
                            tempbins(ismember(condbase_rawdata(:,1),rate_components_perbasepx{pp}(:,1)),2) = rate_components_perbasepx{pp}(:,3); % Spikes 
                            full_dur(fieldlinbin{ii}(pp),rate_components_perbasepx{pp}(:,1)) = rate_components_perbasepx{pp}(:,2);
                            full_spk(fieldlinbin{ii}(pp),rate_components_perbasepx{pp}(:,1)) = rate_components_perbasepx{pp}(:,3);
                            condbase_rawdata(:,2) = nansum( [condbase_rawdata(:,2) tempbins(:,1)] ,2); % Sum duration across session
                            condbase_rawdata(:,3) = nansum( [condbase_rawdata(:,3) tempbins(:,2)] ,2); % Sum spikes across session
                        end
                        condbase_rawdata(:,4) = condbase_rawdata(:,3)./condbase_rawdata(:,2);
                        if sum(condbase_rawdata(:,4)>0) == 0
                            disp('empty sec map');
                        end

                        condbase_componentsperpx{ii,1} = rate_components_perbasepx;
                        condbase_map_rw{ii,1} = condbase_rawdata;

                        %% Distributive hypothesis testing

                        % Null hypothesis: no influence of view other than that caused by a
                        % place effect. Vice versa, no influence of place other than
                        % that caused by a view effect. 
                        % e.g. For a place cell, assuming that firing is ideally location-specific, 
                        % firing rate as a function of view can be calculated
                        % knowing only rate as a function of position and dwell
                        % time as a function of position and view. 
                        % The adequacy of the distributive hypothesis can then be
                        % shown by comparing the expected and observed view firing
                        % distributions. 
                        % Outcome is exactly the same whether or not full array
                        % is used. 

                        % Set up
                        disp('Distributive hypothesis testing');
                        base_array_orig = basemapLrw(1,sort(fieldlinbin{ii}));
                        sec_array_orig = condbase_rawdata(:,4);
                        full_dur1 = full_dur;
                        full_spk1 = full_spk;
                        full_dur1(~ismember(1:size(basemapLrw,2),fieldlinbin{ii}),:) = [];
                        full_dur1(:,~ismember(1:size(secmapLsm,2),usecpxs)) = [];
                        full_spk1(~ismember(1:size(basemapLrw,2),fieldlinbin{ii}),:) = [];
                        full_spk1(:,~ismember(1:size(secmapLsm,2),usecpxs)) = [];
                        if strcmp(msvar{oo},'place')
                            sec_array_orig(usecpxs<3,:) = []; % removing cue and hint for distributive hypothesis testing
                            full_dur1(:,usecpxs<3) = []; % removing cue and hint for distributive hypothesis testing
                            full_spk1(:,usecpxs<3) = []; % removing cue and hint for distributive hypothesis testing
                        end

                        % Predicted rate as a function of sec variable
                        topterm = nansum(base_array_orig'.*full_dur1,1);
                        bottomterm = nansum(full_dur1,1);
                        sec_array_pred = topterm./bottomterm; % raw map
                        % If cell firing is only mod by base var, and sec var
                        % influence is attributable only to inhomogenous sampling,
                        % then dr = 0;
                        % If dr is significant, cell is modulated by sec var. 
                        ratio = log((1+sec_array_orig')./(1+sec_array_pred)); % Muller 1994
            %                 ratio = log(1+sec_array_orig')./(1+sec_array_pred); % Cacucci 2004
                        dr_sec = nansum(abs(ratio))/sum(bottomterm>0); 
                        % Correlation of orig with predicted sec maps
                        vis = ~isnan(sec_array_orig) & ~isnan(sec_array_pred');
                        dcorr_sec = corr2(sec_array_orig(vis),sec_array_pred(vis)');

                        % Predicted rate as a function of base variable
                        topterm = nansum(sec_array_orig'.*full_dur1,2);
                        bottomterm = nansum(full_dur1,2);
                        base_array_pred = topterm./bottomterm; % raw map
                        ratio = log((1+base_array_orig')./(1+base_array_pred)); % Muller 1994
                        dr_base = nansum(abs(ratio))/sum(bottomterm>0); 
                        % Correlation of orig with predicted base maps
                        vis = ~isnan(base_array_orig') & ~isnan(base_array_pred);
                        dcorr_base = corr2(base_array_orig(vis)',base_array_pred(vis));

                        % Store data
                        temp1 = nan(size(basemapLsm));
                        temp2 = nan(size(secmapLsm));
                        temp1(sort(fieldlinbin{ii})) = base_array_pred;
                        if strcmp(msvar{oo},'place')
                            temp2(usecpxs(usecpxs>=3)) = sec_array_pred;
                        else
                            temp2(usecpxs) = sec_array_pred;
                        end
                        basemaps_dist{ii} = temp1;
                        secmaps_dist{ii} = temp2;
                        base_distratio(ii) = dr_base;
                        sec_distratio(ii) = dr_sec;
                        base_distcorr(ii) = dcorr_base;
                        sec_distcorr(ii) = dcorr_sec;

                    end

                    % Store data
                    data.(pairname_short).(msvar{oo}).SI = crit;
                    % data.(pairname_short).(msvar{oo}).SIthr = SIthr;
                    data.(pairname_short).(msvar{oo}).basemapLsm = basemapLsm;
                    data.(pairname_short).(msvar{oo}).basemapLrw = basemapLrw;
                    data.(pairname_short).(msvar{oo}).secmapLsm = secmapLsm;
                    data.(pairname_short).(msvar{oo}).basemapGsm = basemapGsm;
                    data.(pairname_short).(msvar{oo}).basemapGrw = basemapGrw;
                    data.(pairname_short).(msvar{oo}).dummygrid = dummygrid;
                    data.(pairname_short).(msvar{oo}).peakrate_full = peakrate_full;
                    data.(pairname_short).(msvar{oo}).peakrate_set = peakrate_subset;
                    data.(pairname_short).(msvar{oo}).fieldmaxrate_sm = fieldmaxrate_sm;
                    data.(pairname_short).(msvar{oo}).fieldmaxrate_rw = fieldmaxrate_rw;
                    data.(pairname_short).(msvar{oo}).gridnum = gridnum;
                    data.(pairname_short).(msvar{oo}).fieldcoord = fieldcoord;
                    data.(pairname_short).(msvar{oo}).fieldlinbin = fieldlinbin;
                    data.(pairname_short).(msvar{oo}).condbase_map_rw = condbase_map_rw;
                    data.(pairname_short).(msvar{oo}).condbase_componentsperpx = condbase_componentsperpx;
                    data.(pairname_short).(msvar{oo}).basemaps_dist = basemaps_dist;
                    data.(pairname_short).(msvar{oo}).secmaps_dist = secmaps_dist;
                    data.(pairname_short).(msvar{oo}).base_distratio = base_distratio;
                    data.(pairname_short).(msvar{oo}).sec_distratio = sec_distratio;
                    data.(pairname_short).(msvar{oo}).base_distcorr = base_distcorr;
                    data.(pairname_short).(msvar{oo}).sec_distcorr = sec_distcorr;

                    clear fieldmaxrate_sm; clear fieldmaxrate_rw; clear gridnum; clear fieldcoord; clear fieldlinbin; clear condbase_map_rw;
                    clear condbasemapLrw; clear rate_components_perbasepx; clear condbase_rawdata; clear secpx;
                    clear condbase_componentsperpx; clear dummygrid; clear crit; clear SIthr; clear sec_distratio; clear secmaps_dist;
                    clear basemaps_dist; clear base_distratio; clear base_distcorr; clear sec_distcorr;

            end



            %% Smooth conditioned maps
            disp('Smoothing conditioned maps and identifying conditioned fields ...');

            for oo = 1:size(msvar,2)
                condbase_map_sm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            %         condbase_map_adsm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            %         condbase_map_bcsm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            %         condbase_map_dksm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbase_SIC_sm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
            %         condbase_SIC_adsm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
            %         condbase_SIC_bcsm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
            %         condbase_SIC_dksm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);

                condbase_sigfields = zeros(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbase_gridnum = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbase_fieldcoord = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbasefieldlinbin = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbase_fieldsizepercent = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbase_fieldoverlapind = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbase_fieldoverlap = cell(data.(pairname_short).(msvar{oo}).sigfields,1);

                for ii = 1:data.(pairname_short).(msvar{oo}).sigfields
                    % Get conditioned pixel map
                    condbase_rawdata = data.(pairname_short).(msvar{oo}).condbase_map_rw{ii};
                    condbasemapLrw = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                    condbasemapLrw(1,condbase_rawdata(:,1)) = condbase_rawdata(:,4);
            %                 secgridmap = lineartogrid(seclinmap',msvar{2-oo+1},gridSize{2-oo+1});
            %                 dummygridsec = data.(pairname_short).(msvar{2-oo+1}).dummygrid;
                    seclindur = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                    seclinspk = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                    seclindur(1,condbase_rawdata(:,1)) = condbase_rawdata(:,2); 
                    seclinspk(1,condbase_rawdata(:,1)) = condbase_rawdata(:,3); 
                    %% Smooth
                    switch msvar{2-oo+1}
                        case 'place'
                            gridSteps = gridSize{2-oo+1};
                            durG = cell2mat(lineartogrid(seclindur','place',gridSteps));
                            spkG = cell2mat(lineartogrid(seclinspk','place',gridSteps));
                            rateG = cell2mat(lineartogrid(condbasemapLrw','place',gridSteps));

                            % Adaptive smoothing
                            [maps_smG,~,dur_smG] = adsmooth(durG,spkG,1e2);
                            secmap_sm = gridtolinear({maps_smG},'place',gridSteps);
                            secdur_sm = gridtolinear({dur_smG},'place',gridSteps);
                            secdur_sm(isnan(secdur_sm)) = 0;

            %                     % Boxcar smoothing
            %                     unvis = ~(durG>0);
            %                     maps_bcsmG=smooth(rateG,5,unvis,'boxcar');
            %                     dur_bcsmG = smooth(durG,5,unvis,'boxcar');
            %                     secmapbcsm = gridtolinear({maps_bcsmG},'place',gridSteps);
            %                     secdurbcsm = gridtolinear({dur_bcsmG},'place',gridSteps);
            %                     secdurbcsm(isnan(secdurbcsm)) = 0;
            %                     % Disk smoothing
            %                     maps_dksmG=smooth(rateG,5,unvis,'disk');
            %                     dur_dksmG=smooth(durG,5,unvis,'disk');
            %                     secmapdksm = gridtolinear({maps_dksmG},'place',gridSteps);
            %                     secdurdksm = gridtolinear({dur_dksmG},'place',gridSteps);
            %                     secdurdksm(isnan(secdurdksm)) = 0;

                        case 'view'

                            gazeSections = sv.gazeSections;
                            binDepths = gridSize{2-oo+1};

                            % Assign linear bin to grid bin - left to right, bottom to top
                            durG = lineartogrid(seclindur','view',binDepths);
                            spkG = lineartogrid(seclinspk','view',binDepths);
                            rateG = lineartogrid(condbasemapLrw','view',binDepths);

                            % Pad sv map with 5 extra rows
                            n = 5;
                            padpillar = false;
                            [emptyfloorref_pad,~] = padsvmap(n,durG,gazeSections,padpillar);
                            padpillar = true;
                            [durGpad,retrievemap] = padsvmap(n,durG,gazeSections,padpillar);
                            [spkGpad,~] = padsvmap(n,spkG,gazeSections,padpillar);
                            [rateGpad,~] = padsvmap(n,rateG,gazeSections,padpillar);

                            % Adaptive smooth
                            maps_adsmGpad = cell(size(durGpad));
                            dur_adsmGpad = cell(size(durGpad));
                            for jj = 1:size(binDepths,1)
                                if jj == 1 || jj == 2
                                    maps_adsmGpad{jj} = spkGpad{jj}/durGpad{jj};
                                    dur_adsmGpad{jj} = durGpad{jj};
                                else
                                    [maps_adsmGpad{jj},spk_adsmG,dur_adsmGpad{jj}] = adsmooth(durGpad{jj},spkGpad{jj},1e2);
                                end
                            end

            %                     % Boxcar/Disk smooth
            %                     maps_bcsmGpad = cell(size(binDepths,1),1);
            %                     dur_bcsmGpad = cell(size(binDepths,1),1);
            %                     maps_dksmGpad = cell(size(binDepths,1),1);
            %                     dur_dksmGpad = cell(size(binDepths,1),1);
            %                     for jj = 1:size(binDepths,1)
            %                         if jj == 1 || jj == 2
            %                             maps_bcsmGpad{jj} = rateGpad{jj};
            %                             maps_dksmGpad{jj} = rateGpad{jj};
            %                             dur_bcsmGpad{jj} = durGpad{jj};
            %                             dur_dksmGpad{jj} = durGpad{jj};
            %                         else
            %                             unvis = ~(emptyfloorref_pad{jj}>0) | isnan(rateGpad{jj});
            %                             % Boxcar smoothing
            %                             maps_bcsmGpad{jj}=smooth(rateGpad{jj},5,unvis,'boxcar');
            %                             dur_bcsmGpad{jj}=smooth(durGpad{jj},5,unvis,'boxcar');
            %                             % Disk smoothing
            %                             maps_dksmGpad{jj}=smooth(rateGpad{jj},5,unvis,'disk');
            %                             dur_dksmGpad{jj}=smooth(durGpad{jj},5,unvis,'disk');
            %                         end
            %                     end

                            % Unpad smoothed map
                            maps_smG = unpadsvmap(maps_adsmGpad,retrievemap,rateG);
                            dur_smG = unpadsvmap(dur_adsmGpad,retrievemap,rateG);
            %                     maps_bcsmG = unpadsvmap(maps_bcsmGpad,retrievemap,rateG);
            %                     dur_bcsmG = unpadsvmap(dur_bcsmGpad,retrievemap,rateG);
            %                     maps_dksmG = unpadsvmap(maps_dksmGpad,retrievemap,rateG);
            %                     dur_dksmG = unpadsvmap(dur_dksmGpad,retrievemap,rateG);
                            % Convert grid map back to linear sv map
                            secmap_sm = gridtolinear(maps_smG,'view',binDepths);
                            secdur_sm = gridtolinear(dur_smG,'view',binDepths);
                            secdur_sm(isnan(secdur_sm)) = 0;
            %                     secmapbcsm = gridtolinear(maps_bcsmG,'view',binDepths);
            %                     secdurbcsm = gridtolinear(dur_bcsmG,'view',binDepths);
            %                     secdurbcsm(isnan(secdurbcsm)) = 0;
            %                     secmapdksm = gridtolinear(maps_dksmG,'view',binDepths);
            %                     secdurdksm = gridtolinear(dur_dksmG,'view',binDepths);
            %                     secdurdksm(isnan(secdurdksm)) = 0;
                        case 'headdirection'
                            n = 5;
            %                     headdirectionsmooth = 'boxcar';
                            % Smooth
                            rateG = condbasemapLrw'; % raw
                            secmap_sm = smoothdir(condbasemapLrw',n,cr.headdirectionbins);
                            maps_smG = secmap_sm; % no difference between linear and grid map for HD
                            secdur_sm = smoothdir(seclindur',n,cr.headdirectionbins);
                    end
                    %% Calculate SIC/RV
                    if ~strcmp(msvar{2-oo+1},'headdirection')
                        % Calculate SIC from adaptively smoothed map
                        sic_sm = skaggs_sic(secmap_sm,secdur_sm);
            %                 sic_bcsm = skaggs_sic(secmapLsm,secdurLsm);
            %                 sic_dksm = skaggs_sic(secmapdksm,secdurdksm);
                    else
                        % Rayleigh vector
                        map = secmap_sm';
                        binSize=(pi*2)/length(map(1,:));
                        binAngles=(0:binSize:( (359.5/360)*2*pi )) + binSize/2;
                        binWeights=map./(max(map,[],2));
                        S=nansum( sin(binAngles).*binWeights , 2);
                        C=nansum( cos(binAngles).*binWeights , 2);
                        R=sqrt(S.^2+C.^2);
                        meanR=R./nansum(binWeights,2);
                        sic_sm = meanR';
                    end
                    % Output vars
                    condbase_map_sm{ii,1} = secmap_sm';
            %             condbase_map_bcsm{ii,1} = secmapLsm';
            %             condbase_map_dksm{ii,1} = secmapdksm';
                    condbase_SIC_sm(ii,1) = sic_sm;
            %             condbase_SIC_bcsm(ii,1) = sic_bcsm;
            %             condbase_SIC_dksm(ii,1) = sic_dksm;

                    switch msvar{2-oo+1}
                        case 'place'
                            peakrate_sec = nanmax(secmap_sm);
                            meanrate_sec = nanmean(secmap_sm);
                            stdrate_sec = nanstd(secmap_sm);
                            dummygrid = lineartogrid((1:size(secmap_sm,1))',msvar{2-oo+1},gridSize{2-oo+1});
                            rateG = {rateG};
                            maps_smG = {maps_smG};
                        case 'view'
                            peakrate_sec = nanmax(secmap_sm(3:end));
                            meanrate_sec = nanmean(secmap_sm(3:end));
                            stdrate_sec = nanstd(secmap_sm(3:end));
                            dummygrid = lineartogrid((1:size(secmap_sm,1))',msvar{2-oo+1},gridSize{2-oo+1});
                        case 'headdirection'
                            peakrate_sec = nanmax(secmap_sm);
                            meanrate_sec = nanmean(secmap_sm);
                            stdrate_sec = nanstd(secmap_sm);
                            dummygrid = lineartogrid((1:size(secmap_sm,1))',msvar{2-oo+1},gridSize{2-oo+1});
                            rateG = {rateG};
                            maps_smG = {maps_smG};
                    end

                    %% Find max 3 conditioned fields
                    count = 0;
                    maxfieldcount = 3;
                    gridnum = [];
                    fieldcoord = {};
                    fieldlinbin = {};
                    fieldsizepercent = [];
                    fieldoverlapind = {};
                    fieldoverlap = [];
                    fieldmaxrate_rw = [];
                    fieldmaxrate_sm = [];
                    if peakrate_sec >= 0.7 % && sic_adsm>data.(pairname_short).(msvar{2-oo+1}).SIthr
                        for gg = 1:size(maps_smG,1)
                            % Skip cue and hint view maps 
                            if size(maps_smG{gg},1) == 1 && size(maps_smG{gg},2) == 1
                                continue;
                            end
                            % Find fields with at least 1 pixel of > 70% peak rate
            %                         ind_fields = maps_adsmG{gg} > Args.FieldThrPseudo*peakrate_sec;
    %                         ind_fields = maps_smG{gg} > (meanrate_sec+(2*stdrate_sec));
                            ind_fields = maps_smG{gg} > Args.FieldThr*peakrate_subset & ...
                                maps_smG{gg} > (nanmean(secmap_sm)+nanstd(secmap_sm));
                            % Find separate fields
            %                         [fieldlabel,fieldcount] = bwlabel(ind_fields,4); % Only count adjacent pixels if they share edge, not if they share corners
                            [fieldlabel,fieldcount] = bwlabel(ind_fields,8);
                            % For walls and pillars, if there are fields that wrap around split, merge them.
                            if gg > 4
                               % Find possible split fields
                               if any(fieldlabel(:,1)) && any(fieldlabel(:,end))
                                   boundary = [fieldlabel(:,1) fieldlabel(:,end)];
            %                                [blabel,bcount] = bwlabel(boundary,4);
                                   [blabel,bcount] = bwlabel(boundary,8);
                                   for bb = 1:bcount
                                       label = blabel == bb;
                                       label = unique(boundary(label));
                                       labeltochange = label(2:end);
                                       if ~isempty(labeltochange)
                                           for ll = 1:size(labeltochange,1)
                                               % Merge fields around corners
                                               fieldlabel(fieldlabel == labeltochange(ll)) = label(1); 
                                           end
                                       end
                                   end
                               end
                            end
                            % Split fields that are too large i.e. > 1/4 of usable area
                            for ff = 1:fieldcount
                                inds = fieldlabel == ff;
                                % Split fields that are too large
                                if sum(inds(:)) >= ((size(pc.maps_raw,2)-4*8*8)/4) 
            %                                 subinds = inds & maps_adsmG{gg} > Args.FieldSplitThrPseudo*peakrate_sec;
                                    subinds = inds & maps_smG{gg} > meanrate_sec+(3*stdrate_sec);
            %                                 [sublabel,subcount] = bwlabel(subinds,4);
                                    [sublabel,subcount] = bwlabel(subinds,8); % Include bins in field if just corners are touching
                                    sublabel(subinds) = sublabel(subinds) + fieldcount;
                                    fieldcount = fieldcount + subcount;
                                    fieldlabel(inds) = 0; % There will be some px that are not part of new fields cos of different rate thresholds in subfields
                                    fieldlabel(inds) = sublabel(inds);
                                end
                            end
                            % Save fields that are big enough (i.e. >=4 bins around peak)
                            for ff = 1:fieldcount
                                inds = fieldlabel == ff;
                                % Make sure field size is > 15 bins 
            %                             if sum(inds(:)) >= 15
                                if sum(inds(:)) >= 4
            %                                 % Discard field if there are not at least ONE active pixel within this field in the raw map
                                    if sum(rateG{gg}(inds)>0) >= 1
                                        count = count + 1;
                                        fieldmaxrate_rw(count,1) = max(rateG{gg}(inds)); 
                                        fieldmaxrate_sm(count,1) = max(maps_smG{gg}(inds));
                                        switch msvar{2-oo+1}
                                            case 'place'
                                                gridnum(count,1) = 1;
                                            case 'view'
                                                gridnum(count,1) = gg;
                                            case 'headdirection'
                                                gridnum(count,1) = 1;
                                        end
                                        [fieldcoordx_mat, fieldcoordy_mat] = find(inds);
                                        fieldcoord{count,1} = [fieldcoordy_mat size(rateG{gg},1)-fieldcoordx_mat+1]; % In plot coords. x left to right, y bottom to top
                                        fieldlinbin{count,1} = dummygrid{gg}(inds);
                                        fieldsizepercent(count,1) = (sum(sum(inds))/sum(~isnan(secmap_sm))) * 100;

                                        % Check overlap of this conditioned field with original sec base field 
                                        overlap = false(size(data.(pairname_short).(msvar{2-oo+1}).fieldlinbin,1),1);
                                        for mm = 1:size(data.(pairname_short).(msvar{2-oo+1}).fieldlinbin,1)
                                            if intersect(fieldlinbin{count,1},data.(pairname_short).(msvar{2-oo+1}).fieldlinbin{mm})...
                                                    >= 0.5*size(data.(pairname_short).(msvar{2-oo+1}).fieldlinbin{mm},1)
                                                overlap(mm,1) = true;
                                            end
                                        end
                                        fieldoverlapind{count,1} = overlap;
                                        fieldoverlap(count,1) = any(overlap);

            %                                 else 
            %                                     disp(['discarding field: no spikes, ' num2str(num2str(sum(inds(:)))) 'px']);
            %                                     discardfieldnum = discardfieldnum + 1;
            %                                     discardfieldreason(1,end+1) = sum(inds(:));
                                    end
                                elseif sum(inds(:)) == 0 % Field either merged or split before
                                    continue;
            %                             else
            %                                 disp(['discarding field: too small, ' num2str(sum(inds(:))) 'px']);
            %                                 discardfieldnum = discardfieldnum + 1;
            %                                 discardfieldreason(1,end+1) = sum(inds(:));
                                end
                            end
                        end
                    end

                    % If no significant fields, skip to next var
                    if count == 0
                        condbase_sigfields(ii) = 0;
                        continue;
                    else 
                        condbase_sigfields(ii) = min([count maxfieldcount]);
                    end

                    % Sort fields
                    [fieldmaxrate_sm,I] = sort(fieldmaxrate_sm,'descend');
                    fieldmaxrate_rw = fieldmaxrate_rw(I);
                    gridnum = gridnum(I);
                    fieldcoord = fieldcoord(I);
                    fieldlinbin = fieldlinbin(I);
                    fieldsizepercent = fieldsizepercent(I);
                    fieldoverlapind = fieldoverlapind(I);
                    fieldoverlap = fieldoverlap(I);

                    condbase_gridnum{ii,1} = gridnum;
                    condbase_fieldcoord{ii,1} = fieldcoord;
                    condbasefieldlinbin{ii,1} = fieldlinbin;
                    condbase_fieldsizepercent{ii,1} = fieldsizepercent;
                    condbase_fieldoverlapind{ii,1} = fieldoverlapind;
                    condbase_fieldoverlap{ii,1} = fieldoverlap;
                end

                data.(pairname_short).(msvar{oo}).condbase_sigfields  = condbase_sigfields;
                data.(pairname_short).(msvar{oo}).condbase_gridnum = condbase_gridnum;
                data.(pairname_short).(msvar{oo}).condbase_fieldcoord = condbase_fieldcoord;
                data.(pairname_short).(msvar{oo}).condbase_fieldlinbin = condbasefieldlinbin;
                data.(pairname_short).(msvar{oo}).condbase_fieldsizepercent = condbase_fieldsizepercent;
                data.(pairname_short).(msvar{oo}).condbase_fieldoverlapind = condbase_fieldoverlapind;
                data.(pairname_short).(msvar{oo}).condbase_fieldoverlap = condbase_fieldoverlap;
                data.(pairname_short).(msvar{oo}).condbase_map_sm = condbase_map_sm;
            %         data.(pairname_short).(msvar{oo}).condbase_map_bcsm = condbase_map_bcsm;
            %         data.(pairname_short).(msvar{oo}).condbase_map_dksm = condbase_map_dksm;
                data.(pairname_short).(msvar{oo}).condbase_SIC_sm = condbase_SIC_sm;
            %         data.(pairname_short).(msvar{oo}).condbase_SIC_bcsm = condbase_SIC_bcsm;
            %         data.(pairname_short).(msvar{oo}).condbase_SIC_dksm = condbase_SIC_dksm;

                clear count; clear gridnum; clear fieldcoord; clear fieldlinbin; clear fieldsizepercent;
                clear fieldoverlapind; clear fieldoverlap; clear fieldmaxrate_rw; clear fieldmaxrate_sm;

            end

            %% Is activity in conditioned map concentrated in same spot as original sec field? (are fields linked?)
            %   Compare
            %   1. Mean rate of map conditioned on base field, limited to bins within sec fields (meanrate_infield)
            %   2. nshuff samples of Mean rate of pseudopopulation of fields outside of sec fields (meanrate_outfield)
            %      pseudofield needs to be similar size to original sec field, but no spike requirement
            for oo = 1:size(msvar,2)
                disp(['Comparing secondary infield and outfield firing for base ' msvar{oo} ' field ...'])
                % Output variables

                condbase_insecfieldrates = cell(size(data.(pairname_short).(msvar{oo}).condbase_componentsperpx,1),1);
                condbase_outsecfieldrates = cell(size(data.(pairname_short).(msvar{oo}).condbase_componentsperpx,1),1);
                condbase_linkedfield = cell(size(data.(pairname_short).(msvar{oo}).condbase_componentsperpx,1),1);
                for ii = 1:size(data.(pairname_short).(msvar{oo}).condbase_componentsperpx,1) % For each base field

                    userawmap = 0;
                    if userawmap
                        condbase_rawdata = data.(pairname_short).(msvar{oo}).condbase_map_rw{ii};
                        condbase_secbins = condbase_rawdata(:,1);
                        condbase_mapL = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                        condbase_mapL(1,condbase_rawdata(:,1)) = condbase_rawdata(:,4);
                        condbase_mapG = lineartogrid(condbase_mapL',msvar{2-oo+1},gridSize{2-oo+1});
                    else
                        condbase_mapL = data.(pairname_short).(msvar{oo}).condbase_map_sm{ii};
                        condbase_secbins = find(~isnan(condbase_mapL));
                        condbase_mapG = lineartogrid(condbase_mapL',msvar{2-oo+1},gridSize{2-oo+1});
                    end
                    dummygridsec = data.(pairname_short).(msvar{2-oo+1}).dummygrid;

                    % If there are secondary fields
                    if ~isempty(data.(pairname_short).(msvar{2-oo+1}).fieldlinbin)
                        % Output variables
                        insecfieldrates = nan(size(data.(pairname_short).(msvar{2-oo+1}).fieldlinbin,1),1);
                        outsecfieldrates = cell(size(data.(pairname_short).(msvar{2-oo+1}).fieldlinbin,1),1);
                        linkedfield = false(size(data.(pairname_short).(msvar{2-oo+1}).fieldlinbin,1),1);

                        % Test if secondary pixels sampled from this base field are more likely to fall within any of the secondary fields than outside
                        for jj = 1:size(data.(pairname_short).(msvar{2-oo+1}).condbase_fieldlinbin,1) % For each secondary field

                            % Get mean firing rate within secondary field 
                            secfieldlinbin = data.(pairname_short).(msvar{2-oo+1}).fieldlinbin{jj}; % pixels that make up the sec field. Not all of these will be sampled from this base field
                            bins_infield = condbase_secbins(ismember(condbase_secbins,secfieldlinbin)); % find pixels of sec field that are sampled from this base field
                            if isempty(bins_infield) 
                                disp(['Cond map of base field ' num2str(ii) ' does not overlap with sec field ' num2str(jj)]);
                                continue;
                            end
                            meanrate_infield = mean(condbase_mapL(bins_infield),'omitnan');
                            % debug
                            seclinbin_sampled = secfieldlinbin(ismember(secfieldlinbin,condbase_secbins));
                            if ~isempty(setdiff(seclinbin_sampled,bins_infield))
                                error('error with getting overlap of sec bins');
                            end
                            
                            tic;
                            % Get mean firing rates (raw) for nshuff pseudorandom same-size fields outside of secondary field
                            meanrate_outfield = nan(Args.NumShuffles,1);
                            bins_outfield = condbase_secbins(~ismember(condbase_secbins,secfieldlinbin));
                            % Generate a psuedopopulation of outfields same size as sec field 
                            ff = 1;
                            attempt = 0; % Filter 2 of 2: Make sure there are enough outfield px to generate stats (This catches situations where outfield px are enough in number but scattered so that can't form coherent shuffled field without overlapping with original sec field)
                            abandon = false;
                            if length(bins_outfield) > length(bins_infield) % Filter 1 of 2: Make sure there are enough outfield px to generate stats (This catches situations where there are just too few outfield px to begin with) 
                                
                                while ff <= Args.NumShuffles && ~abandon
                                    attempt = attempt + 1;
                                    % Start from a random pixel that is outside of base field and has a spike
                                    startpx = randsample(1:length(bins_outfield),1);
                                    startpx = bins_outfield(1,startpx);
                                    if ismember(startpx,secfieldlinbin) % If random px overlaps with base field, repeat
                    %                                 reset = true;
                                        continue;
                                    end
                                    % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
                                    [gnum,~,~,startindx,startindy] = findgrid(startpx,msvar{2-oo+1});
                                    if ~strcmp(msvar{2-oo+1},'view')
                                        gnum = 1;
                                    end
                                    
                                    % Get the grid coords of starting px
                                    tempmap = condbase_mapG{gnum}; % the actual sampled sec grid map
                                    % Expand radius around starting px until hit the requisite number of px 
                                    while length(startindx)*length(startindy) < 0.5 * size(dummygridsec{gnum},1)*size(dummygridsec{gnum},2)
    
                                        startindx = [startindx(1)-1 startindx startindx(end)+1];
                                        startindy = [startindy(1)-1 startindy startindy(end)+1];
                                        % Keep within env bounds
                                        startindx(startindx < 1 | startindx >size(dummygridsec{gnum},1)) = [];
                                        startindy(startindy < 1 | startindy >size(dummygridsec{gnum},2)) = [];
            
                                        % If reach num of sec bins in original sec map
                                        lin_inds = dummygridsec{gnum}(startindx,startindy);
                                        lin_inds = lin_inds(~isnan(tempmap(startindx,startindy))); % base
                                        ind_pvh = ismember(stcfilt(:,strcmp(stcvars,msvar{oo})),lin_inds);
                                        growingpx = unique(stcfilt(ind_pvh,strcmp(stcvars,msvar{2-oo+1}))); % sec
                                        growingspikes = sum(stcfilt(ind_pvh,6),[],'omitnan');
            
                                        % If fulfill criteria
                                        if sum(sum(~isnan(tempmap(startindx,startindy)))) > length(bins_infield)
                                            break;
                                        end
                                    end

            %                         attempt = attempt + 1;
            %                         % Start from a random pixel that is outside of sec field
            %                         startpx = randsample(1:length(bins_outfield),1);
            %                         startpx = bins_outfield(startpx);
            %                         if ismember(startpx,secfieldlinbin) % If random px overlaps with sec field, repeat
            % %                                 reset = true;
            %                             continue;
            %                         end
            %                         % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
            %                         switch msvar{2-oo+1}
            %                             case 'place'
            %                                 gnum = 1;
            %                             case 'view'
            %                                 if startpx == 1 || startpx == 2 % Make sure not cue or hint
            % %                                        reset = true;
            %                                    continue;
            %                                 end
            %                                 [gnum,~,~] = findgrid(startpx,msvar{2-oo+1});
            %                             case 'headdirection'
            %                                 gnum = 1;
            %                         end
            %                         % Get the grid coords of starting px
            %                         [startindx,startindy] = find(dummygridsec{gnum} == startpx);
            %                         tempmap = condbase_mapG{gnum}; % the actual sampled sec grid map
            %                         % Expand radius around starting px until hit the requisite number of px 
            %                         % while sum(sum(~isnan(tempmap(startindx,startindy)))) < sum(inds_infield)
            %                         while length(startindx)*length(startindy) < size(dummygridsec{gnum},1)*size(dummygridsec{gnum},2)
            %                             if startindx(1) > 1 && startindx(end) < size(dummygridsec{gnum},1)
            %                                 startindx = [startindx(1)-1 startindx startindx(end)+1];
            %                             elseif startindx(1) == 1 && startindx(end) < size(dummygridsec{gnum},1)
            %                                 startindx = [startindx startindx(end)+1];
            %                             elseif startindx(1) > 1 && startindx(end) == size(dummygridsec{gnum},1)
            %                                 startindx = [startindx(1)-1 startindx];
            %                             end
            %                             % If reach required num of px
            %                             if sum(sum(~isnan(tempmap(startindx,startindy)))) > length(bins_infield)% length(startindx)*length(startindy) > size(linbin,1)
            %                                 break;
            %                             end
            %                             if startindy(1) > 1 && startindy(end) < size(dummygridsec{gnum},2)
            %                                 startindy = [startindy(1)-1 startindy startindy(end)+1];
            %                             elseif startindy(1) == 1 && startindy(end) < size(dummygridsec{gnum},2)
            %                                 startindy = [startindy startindy(end)+1];
            %                             elseif startindy(1) > 1 && startindy(end) == size(dummygridsec{gnum},2)
            %                                 startindy = [startindy(1)-1 startindy];
            %                             end
            %                             % If exceed map bounds
            %                             if startindx(1) == 1 && startindx(end) == size(dummygridsec{gnum},1) && startindy(1) == 1 && startindy(end) == size(dummygridsec{gnum},2)
            %                                 break;
            %                             end
            %                         end

                                    % sampled pixels
                                    sampledpx = dummygridsec{gnum}(startindx,startindy);
                                    sampledpx = sampledpx(~isnan(tempmap(startindx,startindy)));
                                    % get another starting pixel if sampled field overlaps too much with sec field 
                                    if size(intersect(sampledpx,secfieldlinbin),1) > 0.25*length(secfieldlinbin) % If sampled field overlaps with more than half of sec field
                                        if attempt == Args.NumShuffles && ff < 10
                                            abandon = true;
                                            disp(['Abandoning finding pseudo sec fields for ' msvar{2-oo+1} ' field ' num2str(jj)]);
                                        end
                                        continue;
                                    end
                                    % Remove empty pixels from sampled field
                                    pxsub = dummygridsec{gnum}(startindx,startindy);
                                    inds_sampled = ~isnan(tempmap(startindx,startindy));
                                    pxsub = pxsub(inds_sampled);
                                    if length(bins_infield)>length(pxsub)
                                        error('Not enough sample pixels to draw pseudopopulation from');
                                    end
                                    inds_keep = sort(randsample(1:length(pxsub),length(bins_infield)))';
                                    pxsub = pxsub(inds_keep);
                                    % Start over if exactly the same px as field of interest
                                    if isempty(setdiff(pxsub,secfieldlinbin))
            %                                 reset = true;
                                        continue;
                                    end
                                    % Get sec bins actually sampled 
                                    outfieldlinbin = pxsub;
                                    bins_outfield = condbase_secbins(ismember(condbase_secbins,outfieldlinbin));
                                    meanrate_outfield(ff,1) = mean(condbase_mapL(bins_outfield),'omitnan');
                                    % meanrate_outfield(ff,1) = sum(condbase_rawdata(bins_outfield,3),[],'omitnan') / sum(condbase_rawdata(bins_outfield,2),[],'omitnan');

                                    ff = ff+1;
                                end
                            end
                            disp([num2str(Args.NumShuffles) ' outfield rates for base ' msvar{oo} ' field ' num2str(ii) ' sec ' msvar{2-oo+1} ' field ' num2str(jj) ' took ' num2str(toc) 's']);
                            % Store data
                            insecfieldrates(jj,1) = meanrate_infield;
                            outsecfieldrates{jj,1} = meanrate_outfield;
                            if meanrate_infield > prctile(meanrate_outfield,95)
                                linkedfield(jj,1) = true;
                            end
                        end
                        % Store data
                        condbase_insecfieldrates{ii,1} = insecfieldrates;
                        condbase_outsecfieldrates{ii,1} = outsecfieldrates;
                        condbase_linkedfield{ii,1} = linkedfield;
                    else 
                        % if data.(pairname_short).mixsel
                            disp(['no sec ' msvar{2-oo+1} ' fields']);
                        % end
                    end
                end
                % Store data
                data.(pairname_short).(msvar{oo}).condbase_insecfieldrates = condbase_insecfieldrates;
                data.(pairname_short).(msvar{oo}).condbase_outsecfieldrates = condbase_outsecfieldrates;
                data.(pairname_short).(msvar{oo}).condbase_linkedfield = condbase_linkedfield;

            end

            %% Is the SIC of the conditioned map spuriously high (because of sparse map) or valid?
            %   Compare:
            %   1. SIC of original conditioned map
            %   2. SIC of conditioned map from pseudopopulation of base maps with 
            %       - comparable numbers of spikes
            %       - at least min num of base pixels
            %       - within +- 25% size of conditioned map
            %       - no need to have overlap with conditioned/sec fields


            %%  Method 1. Use cellfun to find nshuff pseudo base fields instead of looping
            % for oo = 1:size(msvar,2) 
            %     disp(['Creating pseudopopulation of base maps for ' msvar{oo} ' field ...']);
            % 
            %     for ii = 1:size(data.(pairname_short).(msvar{oo}).condbase_componentsperpx,1) % For each base field
            %         % Load base data
            %         userawmap = 1;
            %         if ~userawmap
            %             baselinmap = data.(pairname_short).(msvar{oo}).basemapLrw; % whole base map - linear
            %             basegridmap = data.(pairname_short).(msvar{oo}).basemapGrw; % whole base map - grid
            %             condbasemap = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
            %             condbasemap(data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,1)) = data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,4);
            %             condbaseallbins = data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,1);
            %         else
            %             baselinmap = data.(pairname_short).(msvar{oo}).basemapLsm; % whole base map - linear
            %             basegridmap = data.(pairname_short).(msvar{oo}).basemapGsm; % whole base map - grid
            %             condbasemap = data.(pairname_short).(msvar{oo}).condbase_map_sm{ii};
            %             condbaseallbins = find(~isnan(data.(pairname_short).(msvar{oo}).condbase_map_sm{ii}));
            %         end
            %         dummygridbase = data.(pairname_short).(msvar{oo}).dummygrid;
            %         dummygridsec = data.(pairname_short).(msvar{2-oo+1}).dummygrid;
            %         % Get constrain parameters
            %         basefieldnumspk = sum(data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,3)); % number of spikes within base field
            %         basefieldbins = data.(pairname_short).(msvar{oo}).fieldlinbin{ii}; % bins of this base field
            %         basedrawpx = find(~isnan(baselinmap)); % base bins sampled, full map, spike or no spike (restricting to starting on px with spikes will force resampling of same fields for low firing count cells
            %         basedrawpx = setdiff(basedrawpx,basefieldbins); % base bins sampled outside of this base field
            %         if strcmp(msvar{oo},'view')
            %             basedrawpx(basedrawpx == 1 | basedrawpx == 2) = [];
            %         end
            % 
            %         tic;
            % 
            %         stc_set = {};
            %         lin_inds_set = {};
            %         while size(stc_set,2) < Args.NumShuffles
            %             % Start from a random pixel that is outside of base field and has a spike
            %             startpx = randsample(1:length(basedrawpx),2*Args.NumShuffles,true);
            %             startpx = basedrawpx(1,startpx);
            %             startpx = sort(startpx);
            %             % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
            %             if strcmp(msvar{oo},'view')
            %                 [gnum,~,~,startindx,startindy] = findgrid(startpx,msvar{oo});
            %             else
            %                 [gnum,~,~,startindx,startindy] = findgrid(startpx,msvar{oo});
            %                 gnum(:) = 1;
            %             end
            %             tempmap = basegridmap(gnum); % the actual sampled sec grid map
            % 
            %             % Initialize variables
            %             indx = mat2cell(startindx,1,[repmat(1,1,length(startindx))]);
            %             indy = mat2cell(startindy,1,[repmat(1,1,length(startindy))]);
            %             lin_inds = cell(size(startindx));
            % 
            %             for gg = min(gnum):max(gnum)
            %                 disp(['gg = ' num2str(gg)]);
            %                 countx = 1;
            %                 county = 1;
            %                 indg = gnum == gg;
            %                 indtodo = indg;
            %                 grid = dummygridbase{gg};
            %                 while (countx < size(grid,1)/2 || county < size(grid,2)/2) && sum(indtodo) ~= 0
            %                     % if gg == 5
            %                     %     disp('test');
            %                     % end
            %                     % Expand radius around starting px until hit the requisite number of px 
            %                     indx(indtodo) = cellfun(@(x) [x(1)-1 x x(end)+1],indx(indtodo),'UniformOutput',false);
            %                     indy(indtodo) = cellfun(@(x) [x(1)-1 x x(end)+1],indy(indtodo),'UniformOutput',false);
            %                     [limx limy] = size(tempmap{find(indtodo,1)});
            %                     indx(indtodo) = cellfun(@(x) x(x>0 & x<=limx),indx(indtodo),'UniformOutput',false);
            %                     indy(indtodo) = cellfun(@(x) x(x>0 & x<=limy),indy(indtodo),'UniformOutput',false);
            % 
            %                     % Find sec bins conditioned on these pseudo base bins and corresponding num of spikes
            %                     lin_inds(indtodo) = cellfun(@(x,y) grid(x,y),indx(indtodo),indy(indtodo),'UniformOutput',false);
            %                     lin_inds = cellfun(@(x) x(:),lin_inds,'UniformOutput',false);
            %                     lin_inds = cellfun(@(x) x(~isnan(baselinmap(x))),lin_inds,'UniformOutput',false);
            %                     lin_inds_num = cellfun(@length,lin_inds,'UniformOutput',true);
            %                     ind_pvh = cellfun(@(x) ismember(stcfilt(:,strcmp(stcvars,msvar{oo})),x),lin_inds,'UniformOutput',false);
            %                     growingpx = cellfun(@(x) unique(stcfilt(x,strcmp(stcvars,msvar{2-oo+1}))),ind_pvh,'UniformOutput',false); % sec
            %                     growingpxnum = cellfun(@length,growingpx,'UniformOutput',true);
            %                     growingspikes = cellfun(@(x) sum(stcfilt(x,6),[],'omitnan'),ind_pvh,'UniformOutput',true);
            %                     % If reach num of bins in orig conditioned map and comparable num of spikes, stop expanding
            %                     sufficient = indg' & lin_inds_num > 0.8*length(basefieldbins) & growingpxnum > 0.8*length(condbaseallbins) & ...
            %                         growingspikes > 0.7*basefieldnumspk;
            %                     indtodo = indtodo & ~sufficient';
            % 
            %                     if countx == 1 || mod(countx,10) == 5 || sum(indtodo) == 0
            %                         disp(['Sum indg = ' num2str(sum(indtodo)) ', round = ' num2str(countx)]);
            %                     end
            %                     countx = countx + 1;
            %                     county = county + 1;
            %                 end
            %                 % done = ;
            %                 disp(['Sum indg = ' num2str(sum(indtodo)) ', round = ' num2str(countx)]);
            %                 stc_set(end+1:end+sum(sufficient)) = ind_pvh(sufficient);
            %                 lin_inds_set(end+1:end+sum(sufficient)) = lin_inds(sufficient);
            %             end
            %         end
            %         toc;
            %         % Reduce to nshuff
            %         indkeep = randsample(1:size(stc_set,2),Args.NumShuffles);
            %         stc_setkeep = stc_set(indkeep);
            %         lin_indskeep = lin_inds(indkeep);
            % 
            % 
            %     end
            % end


            %% Method 2: Find nshuff pseudo base maps by looping
            secfieldnumbinset = cell(size(msvar,2),1);
            condfield_inbinsset = cell(size(msvar,2),1);
            condfield_nonoverlapbinset = cell(size(msvar,2),1);
            condfield_tertlinbin = cell(size(msvar,2),1);
            for oo = 1:size(msvar,2) 
                disp(['Creating pseudopopulation of base maps for ' msvar{oo} ' field ...'])
                
                % Output variables
                condpseudo_secdataperfield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                secfieldnumbinperfield = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
                condbase_fieldlinbin = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                condfield_nonoverlapbinsperfield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                cond_tertlinbinsperfield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                for ii = 1:size(data.(pairname_short).(msvar{oo}).condbase_componentsperpx,1) % For each base field
                    
                    userawmap = 0;
                    % Load up maps
                    if ~userawmap
                        baselinmap = data.(pairname_short).(msvar{oo}).basemapLrw; % whole base map - linear
                        basegridmap = data.(pairname_short).(msvar{oo}).basemapGrw; % whole base map - grid
                        condbasemap = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
                        condbasemap(data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,1)) = data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,4);
                        condbaseallbins = data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,1);
                    else
                        baselinmap = data.(pairname_short).(msvar{oo}).basemapLsm; % whole base map - linear
                        basegridmap = data.(pairname_short).(msvar{oo}).basemapGsm; % whole base map - grid
                        condbasemap = data.(pairname_short).(msvar{oo}).condbase_map_sm{ii};
                        condbaseallbins = find(~isnan(data.(pairname_short).(msvar{oo}).condbase_map_sm{ii}));
                    end
                    dummygridbase = data.(pairname_short).(msvar{oo}).dummygrid;
                    dummygridsec = data.(pairname_short).(msvar{2-oo+1}).dummygrid;

                    % Get constrain parameters
                    basefieldnumspk = sum(data.(pairname_short).(msvar{oo}).condbase_map_rw{ii}(:,3)); % number of spikes within base field
                    basefieldbins = data.(pairname_short).(msvar{oo}).fieldlinbin{ii}; % bins of this base field
                    basedrawpx = find(~isnan(baselinmap)); % base bins sampled, full map, spike or no spike
                    basedrawpx = setdiff(basedrawpx,basefieldbins); % base bins sampled outside of this base field
                    if strcmp(msvar{oo},'view')
                        basedrawpx(basedrawpx == 1 | basedrawpx == 2) = [];
                    end
                    
                    condbasenumbin = length(condbaseallbins); % number of bins occupied in conditioned map
                    condbasefieldlinbin = find(condbasemap > mean(condbasemap,'omitnan'))'; % simplistically, fields in sec/conditioned maps are defined as pixels exceeding mean. no contiguity requirement
                    condbasenonfieldlinbin = setdiff(condbaseallbins,condbasefieldlinbin);

            %             % Get 1000 samples of field same size as conditioned fields
            %             ff = 1;
            %             shuff = Args.NumShuffles;
            %             attempt = 0;
            %             pseudosecdur = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
            %             pseudosecspk = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
            %             pseudosecmap_raw = nan(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
            %             pseudobasebins = cell(shuff,1);
            %             pseudononoverlap = cell(shuff,1);
            %             pseudotertbins = cell(shuff,1);
            %             tic;
            %                 while ff <= shuff
            %                     if rem(ff,100)==0
            %                         disp(['ff=' num2str(ff)]);
            %                     end
            %                     % Start from a random sec pixel on base map
            %                     startpxind = randsample(1:length(basedrawpx),1);
            %                     startpx = basedrawpx(startpxind);
            %                     attempt = attempt + 1;
            %                     
            %                     % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
            %                     switch msvar{oo}
            %                         case 'place'
            %                             gnum = 1;
            %                         case 'view'
            %                             if startpx == 1 || startpx == 2 % Make sure not cue or hint
            % %                                        reset = true;
            %                                continue;
            %                             end
            %                             [gnum,~,~] = findgrid(startpx,msvar{oo});
            %                         case 'headdirection'
            %                             gnum = 1;
            %                     end
            %                     
            %                     % Get the grid coords of starting px
            %                     [startindx,startindy] = find(dummygridbase{gnum} == startpx);
            %                     tempmap = basegridmap{gnum}; % the actual sampled sec grid map
            %                     
            %                     % Expand radius around starting px collected sec map includes a sec field that is same size as conditioned field
            %                     while length(startindx)*length(startindy) < size(dummygridbase{gnum},1)*size(dummygridbase{gnum},2)
            %                         if startindx(1) > 1 && startindx(end) < size(dummygridbase{gnum},1)
            %                             startindx = [startindx(1)-1 startindx startindx(end)+1];
            %                         elseif startindx(1) == 1 && startindx(end) < size(dummygridbase{gnum},1)
            %                             startindx = [startindx startindx(end)+1];
            %                         elseif startindx(1) > 1 && startindx(end) == size(dummygridbase{gnum},1)
            %                             startindx = [startindx(1)-1 startindx];
            %                         end
            %                         
            %                         % Get sec map pixels
            %                         lin_inds = dummygridbase{gnum}(startindx,startindy);
            %                         lin_inds = lin_inds(~isnan(tempmap(startindx,startindy))); % base
            %                         ind_pv = ismember(stcfill(:,strcmp(stcvars,msvar{oo})),lin_inds);
            %                         growingpx = unique(stcfill(ind_pv,strcmp(stcvars,msvar{2-oo+1}))); % sec
            % %                         switch msvar{oo}
            % %                             case 'place'
            % %                                 ind_pv = ismember(stcfill(:,2),lin_inds);
            % %                                 growingpx = unique(stcfill(ind_pv,4));
            % %                             case 'view'
            % %                                 ind_pv = ismember(stcfill(:,4),lin_inds);
            % %                                 growingpx = unique(stcfill(ind_pv,2));
            % %                         end
            %                         growingspikes = nansum(stcfill(ind_pv,6));
            %                         
            %                         % If captured sec field same size as conditioned fields but no overlap
            %                         nonoverlap = setdiff(growingpx,condfield_inbins);
            %                         if size(nonoverlap,1) > size(condfield_inbins,1) && size(lin_inds,1)>0.9*size(basefieldbins,1)
            %                             break;
            %                         end
            %                         if startindy(1) > 1 && startindy(end) < size(dummygridbase{gnum},2)
            %                             startindy = [startindy(1)-1 startindy startindy(end)+1];
            %                         elseif startindy(1) == 1 && startindy(end) < size(dummygridbase{gnum},2)
            %                             startindy = [startindy startindy(end)+1];
            %                         elseif startindy(1) > 1 && startindy(end) == size(dummygridbase{gnum},2)
            %                             startindy = [startindy(1)-1 startindy];
            %                         end
            %                         
            %                         % If exceed map bounds
            %                         if startindx(1) == 1 && startindx(end) == size(dummygridbase{gnum},1) && startindy(1) == 1 && startindy(end) == size(dummygridbase{gnum},2)
            %                             break;
            %                         end
            %                     end
            %                     lin_inds = dummygridbase{gnum}(startindx,startindy);
            %                     lin_inds = lin_inds(~isnan(tempmap(startindx,startindy)));
            %                     ind_pv = ismember(stcfill(:,strcmp(stcvars,msvar{oo})),lin_inds);
            %                     growingpx = unique(stcfill(ind_pv,strcmp(stcvars,msvar{2-oo+1})));
            % %                     switch msvar{oo}
            % %                         case 'place'
            % %                             ind_pv = ismember(stcfill(:,2),lin_inds);
            % %                             growingpx = unique(stcfill(ind_pv,4));
            % %                         case 'view'
            % %                             ind_pv = ismember(stcfill(:,4),lin_inds);
            % %                             growingpx = unique(stcfill(ind_pv,2));
            % %                     end
            %                     growingspikes = nansum(stcfill(ind_pv,6));
            % 
            %                     % If captured sec field same size as conditioned fields but no overlap
            %                     nonoverlap = setdiff(growingpx,condfield_inbins);
            %                     if size(nonoverlap,1) < size(condfield_inbins,1) || size(lin_inds,1) < 0.9*size(basefieldbins,1) ...
            %                             || growingspikes < 0.3*basefieldnumspk
            %                         continue;
            %                     elseif size(intersect(lin_inds,basefieldbins),1) > 0 % 0.5*size(basefieldbins,1) % If sampled field overlaps with more than half of sec field
            % %                         if attempt == shuff && ff < 10
            %                         if strcmp(cwd,'/Volumes/Hippocampus/Data/picasso-misc/20180828/session01/array01/channel021/cell03') % attempt > 100000
            % %                             abandon = true;
            %                             if size(intersect(lin_inds,basefieldbins),1) > 0.3*size(basefieldbins,1)
            %                                 continue;
            %                             end
            %                         else
            %                             continue;
            %                         end
            %                     end
            %                     
            %                     % Get tertiary field
            % %                     for kk = 1:size(data.(pairname_short).(msvar{oo}).condbase_fieldlinbin{ii},1)
            %                     % Make sure size of the shared sec field between base and pseudo base field is same or larger than size of base field
            %                         tertlinbin = intersect(condfield_outbins,nonoverlap);
            %                         if size(tertlinbin,1) < 0.3*size(basefieldbins,1) % Tried 0.5 but couldn't get any ff for 1102ch19c2
            %                             if attempt > 100000 && strcmp(cwd,'/Volumes/Hippocampus/Data/picasso-misc/20180905/session01/array01/channel020/cell02')
            %                                 if size(tertlinbin,1) < 0.1*size(basefieldbins,1)
            %                                     continue;
            %                                 end
            %                             else
            %                                 continue;
            %                             end
            %                         end
            % %                     end
            %                     
            %                     % Find pseudo sec maps
            %                     usecpxs = [];
            %                     for pp = 1:size(lin_inds,1)
            %                         secpx = [];
            %                         ind_pv = stcfill(:,strcmp(stcvars,msvar{oo})) == lin_inds(pp);
            %                         secpx(:,1) = stcfill(ind_pv,strcmp(stcvars,msvar{2-oo+1})); % sec px
            % %                         switch msvar{oo}
            % %                             case 'place'
            % %                                 ind_pv = stcfill(:,2) == lin_inds(pp);
            % %                                 secpx(:,1) = stcfill(ind_pv,4); % sec px
            % %                             case 'view'
            % %                                 ind_pv = stcfill(:,4) == lin_inds(pp);
            % %                                 secpx(:,1) = stcfill(ind_pv,2); % sec px
            % %                                 
            % %                         end
            %                         secpx(:,2) = stcfill(ind_pv,5); % dur
            %                         secpx(:,3) = stcfill(ind_pv,6); % spk
            %                         
            %                         % Get firing rates 
            %                         usecpx = unique(secpx(:,1));
            %                         usecpx(isnan(usecpx)) = [];
            %                         rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
            %                         rate_components_px(:,1) = usecpx;
            %                         if any(isnan(usecpx))
            %                             disp(nan);
            %                         end
            %                         ind_timechange = find(secpx(:,2) > 0);
            %                         for cc = 1:size(ind_timechange,1) % For each instance of being in this base pixel
            %                             linbintemp = nan(size(usecpx,1),2); % Temp dur and spikes for this secondary pixel(s)
            %                             if cc == size(ind_timechange,1)
            %                                 sampind = ind_timechange(cc):size(secpx,1);
            %                             else
            %                                 sampind = ind_timechange(cc):ind_timechange(cc+1)-1; % index into secondary pixel(s) for this instance
            %                             end
            %                             secsamp = secpx(sampind,1); 
            %                             setind = ismember(usecpx,secsamp);
            %                             linbintemp(setind,1) = secpx(sampind(1),2); % duration for this instance listed with first secondary pixel
            %                             linbintemp(setind,2) = secpx(sampind(end),3); % spikes for this instance listed with last secondary pixel
            %                             rate_components_px(:,2) = nansum( [rate_components_px(:,2) linbintemp(:,1)] ,2); % Sum duration for this sec pixel across instances
            %                             rate_components_px(:,3) = nansum( [rate_components_px(:,3) linbintemp(:,2)] ,2); % Sum spikes for this sec pixel across instances
            %                         end
            %                         rightfulnans = rate_components_px(:,2) == 0;
            %                         rate_components_px(rightfulnans,2) = nan;
            %                         rate_components_px(:,4) = rate_components_px(:,3)./rate_components_px(:,2); % Compute firing rates
            % 
            %                         % Collect all sec rate information for the base pixels
            %                         rate_components_perbasepx{pp} = ( rate_components_px );
            %                         if ~isempty(rate_components_px)
            %                             maxrate(pp) = max(rate_components_px(:,4));
            %                         else 
            %                             maxrate(pp) = NaN;
            %                         end
            %                         secpxs{pp} = secpx;
            %                         usecpxs = union(usecpxs,usecpx); % set of secondary pixels covered in whole session
            %                     end
            %                     
            %                     secdata = usecpxs; 
            %                     secdata(:,2:3) = NaN;
            %                     % Consolidate occupancy and spikes across base field
            %                     for pp = 1:size(lin_inds,1)
            %                         tempbins = nan(size(secdata,1),2);
            %                         tempbins(ismember(secdata(:,1),rate_components_perbasepx{pp}(:,1)),1) = rate_components_perbasepx{pp}(:,2); % Duration 
            %                         tempbins(ismember(secdata(:,1),rate_components_perbasepx{pp}(:,1)),2) = rate_components_perbasepx{pp}(:,3); % Spikes 
            %                         secdata(:,2) = nansum( [secdata(:,2) tempbins(:,1)] ,2); % Sum duration across session
            %                         secdata(:,3) = nansum( [secdata(:,3) tempbins(:,2)] ,2); % Sum spikes across session
            %                     end
            %                     secdata(:,4) = secdata(:,3)./secdata(:,2);
            %                     
            %                     pseudodur = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
            %                     pseudodur(secdata(:,1)) = secdata(:,2);
            %                     pseudospk = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
            %                     pseudospk(secdata(:,1)) = secdata(:,3);
            %                     pseudomap = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
            %                     pseudomap(secdata(:,1)) = secdata(:,4);
            % 
            %                     pseudosecdur(ff,:) = pseudodur;
            %                     pseudosecspk(ff,:) = pseudospk;
            %                     pseudosecmap_raw(ff,:) = pseudomap;
            %                     pseudobasebins{ff,:} = lin_inds;
            %                     pseudononoverlap{ff,:} = nonoverlap;
            %                     pseudotertbins{ff,:} = tertlinbin;
            %                     
            %                     ff = ff + 1;
            %                     
            %                 end

                    % Generate a psuedopopulation of 1000 base fields with same num of spikes as original base field 
                    ff = 1;
                    shuff = Args.NumShuffles;
                    attempt = 0; % Filter 2 of 2: Make sure there are enough outfield px to generate stats (This catches situations where outfield px are enough in number but scattered so that can't form coherent shuffled field without overlapping with original sec field)
                    abandon = false;
                    pseudosecmap_sm = nan(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudosecSIC_sm = nan(shuff,1);
                    pseudosecdur = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudosecspk = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudosecmap_raw = nan(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudobasebins = cell(shuff,1);
                    pseudononoverlap = cell(shuff,1);
                    pseudotertbins = cell(shuff,1);
                    % tic;
                   
                    while ff <= shuff && ~abandon

                        attempt = attempt + 1;
                        % Start from a random pixel that is outside of base field and has a spike
                        startpx = randsample(1:length(basedrawpx),1);
                        startpx = basedrawpx(1,startpx);
                        if ismember(startpx,basefieldbins) % If random px overlaps with base field, repeat
        %                                 reset = true;
                            continue;
                        end
                        % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
                        [gnum,~,~,startindx,startindy] = findgrid(startpx,msvar{oo});
                        if ~strcmp(msvar{oo},'view')
                            gnum = 1;
                        end
                        
                        % Get the grid coords of starting px
                        % [startindx,startindy] = find(dummygridbase{gnum} == startpx);
                        tempmap = basegridmap{gnum}; % the actual sampled sec grid map
                        % Expand radius around starting px until hit the requisite number of px 
                        while length(startindx)*length(startindy) < 0.5 * size(dummygridbase{gnum},1)*size(dummygridbase{gnum},2)

                            startindx = [startindx(1)-1 startindx startindx(end)+1];
                            startindy = [startindy(1)-1 startindy startindy(end)+1];
                            % Keep within env bounds
                            startindx(startindx < 1 | startindx >size(dummygridbase{gnum},1)) = [];
                            startindy(startindy < 1 | startindy >size(dummygridbase{gnum},2)) = [];

                            % If reach num of sec bins in original sec map
                            lin_inds = dummygridbase{gnum}(startindx,startindy);
                            lin_inds = lin_inds(~isnan(tempmap(startindx,startindy))); % base
                            ind_pvh = ismember(stcfilt(:,strcmp(stcvars,msvar{oo})),lin_inds);
                            growingpx = unique(stcfilt(ind_pvh,strcmp(stcvars,msvar{2-oo+1}))); % sec
                            growingspikes = sum(stcfilt(ind_pvh,6),[],'omitnan');

                            % If fulfill criteria, include this pseudo base field
                            if length(growingpx) > 0.8*condbasenumbin && length(lin_inds) > 0.8*length(basefieldbins) && ...
                                    growingspikes > 0.7*basefieldnumspk
                                break;
                            end
                        end
                        % disp(['Finding pseudo base field: ' num2str(toc) 's']);

                        nonoverlap = setdiff(growingpx,condbasefieldlinbin); % sec pixels resultant from pseudo base field, excluding fields from orig conditioned maps
                        if size(intersect(lin_inds,basefieldbins),1) > 0.5*size(basefieldbins,1) % If pseudo base field overlaps with more than half of orig base field
                            if attempt == shuff && ff < 10
                                abandon = true; % run too many rounds and failing to find pseudomaps of appropriate criteria. give up
                            end
                            continue;
                        end
                        % tertlinbin = intersect(condbasenonfieldlinbin,nonoverlap); % pix overlap between orig conditioned map and sec pixels of pseudo base field , minus the orig conditioned field

                        % Finding sec maps conditioned on pseudo base maps
                        secpx = stcfilt(ind_pvh,[find(strcmp(stcvars,msvar{2-oo+1})) 5 6]); % [bin, dur, spk]
                        usecpx = unique(secpx(:,1));
                        secpx(end+1,:) = [length(data.(pairname_short).(msvar{2-oo+1}).basemapLrw) 0 0]; % without max bin, accumarray doesn't work right
                        secdata = nan(length(data.(pairname_short).(msvar{2-oo+1}).basemapLrw),3);
                        secdata(usecpx,1) = usecpx;
                        secdata(:,2) = accumarray(secpx(:,1),secpx(:,2));
                        secdata(:,3) = accumarray(secpx(:,1),secpx(:,3));
                        secdata(:,4) = secdata(:,3)./secdata(:,2);
                        
                        pseudosecdur(ff,:) = secdata(:,2);
                        pseudosecspk(ff,:) = secdata(:,3);
                        pseudosecmap_raw(ff,:) = secdata(:,4);
                        pseudobasebins{ff} = lin_inds;
                        pseudononoverlap{ff,:} = nonoverlap;

                        % tic;
                        % % Find pseudo sec maps
                        % usecpxs = [];
                        % for pp = 1:size(lin_inds,1)
                        %     secpx = [];
                        %     ind_pvh = stcfilt(:,strcmp(stcvars,msvar{oo})) == lin_inds(pp);
                        %     secpx(:,1) = stcfilt(ind_pvh,strcmp(stcvars,msvar{2-oo+1})); % sec px
                        %     secpx(:,2) = stcfilt(ind_pvh,5); % dur
                        %     secpx(:,3) = stcfilt(ind_pvh,6); % spk
                        % 
                        %     % Get firing rates 
                        %     usecpx = unique(secpx(:,1));
                        %     usecpx(isnan(usecpx)) = [];
                        %     rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
                        %     rate_components_px(:,1) = usecpx;
                        %     if any(isnan(usecpx))
                        %         disp(nan);
                        %     end
                        %     ind_timechange = find(secpx(:,2) > 0);
                        %     for cc = 1:size(ind_timechange,1) % For each instance of being in this base pixel
                        %         linbintemp = nan(size(usecpx,1),2); % Temp dur and spikes for this secondary pixel(s)
                        %         if cc == size(ind_timechange,1)
                        %             newind = ind_timechange(cc):size(secpx,1);
                        %         else
                        %             newind = ind_timechange(cc):ind_timechange(cc+1)-1; % index into secondary pixel(s) for this instance
                        %         end
                        %         newview = secpx(newind,1); 
                        %         newset = ismember(usecpx,newview);
                        %         linbintemp(newset,1) = secpx(newind(1),2); % duration for this instance listed with first secondary pixel
                        %         linbintemp(newset,2) = secpx(newind(end),3); % spikes for this instance listed with last secondary pixel
                        %         rate_components_px(:,2) = nansum( [rate_components_px(:,2) linbintemp(:,1)] ,2); % Sum duration for this sec pixel across instances
                        %         rate_components_px(:,3) = nansum( [rate_components_px(:,3) linbintemp(:,2)] ,2); % Sum spikes for this sec pixel across instances
                        %     end
                        %     rightfulnans = rate_components_px(:,2) == 0;
                        %     rate_components_px(rightfulnans,2) = nan;
                        %     rate_components_px(:,4) = rate_components_px(:,3)./rate_components_px(:,2); % Compute firing rates
                        % 
                        %     % Collect all sec rate information for the base pixels
                        %     rate_components_perbasepx{pp} = ( rate_components_px );
                        %     if ~isempty(rate_components_px)
                        %         maxrate(pp) = max(rate_components_px(:,4));
                        %     else 
                        %         maxrate(pp) = NaN;
                        %     end
                        %     secpxs{pp} = secpx;
                        %     usecpxs = union(usecpxs,usecpx); % set of secondary pixels covered in whole session
                        % end
                        % disp(['Finding pseudo sec rate components: ' num2str(toc) 's']);
                        % 
                        % tic;
                        % condpseudo_rawdata = usecpxs; 
                        % condpseudo_rawdata(:,2:3) = NaN;
                        % % Consolidate occupancy and spikes across base field
                        % for pp = 1:size(lin_inds,1)
                        %     tempbins = nan(size(condpseudo_rawdata,1),2);
                        %     tempbins(ismember(condpseudo_rawdata(:,1),rate_components_perbasepx{pp}(:,1)),1) = rate_components_perbasepx{pp}(:,2); % Duration 
                        %     tempbins(ismember(condpseudo_rawdata(:,1),rate_components_perbasepx{pp}(:,1)),2) = rate_components_perbasepx{pp}(:,3); % Spikes 
                        %     condpseudo_rawdata(:,2) = nansum( [condpseudo_rawdata(:,2) tempbins(:,1)] ,2); % Sum duration across session
                        %     condpseudo_rawdata(:,3) = nansum( [condpseudo_rawdata(:,3) tempbins(:,2)] ,2); % Sum spikes across session
                        % end
                        % condpseudo_rawdata(:,4) = condpseudo_rawdata(:,3)./condpseudo_rawdata(:,2);
                        % disp(['Finding pseudo sec map: ' num2str(toc) 's']);
                        % 
                        % pseudodur = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
                        % pseudodur(condpseudo_rawdata(:,1)) = condpseudo_rawdata(:,2);
                        % pseudospk = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
                        % pseudospk(condpseudo_rawdata(:,1)) = condpseudo_rawdata(:,3);
                        % pseudomap = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
                        % pseudomap(condpseudo_rawdata(:,1)) = condpseudo_rawdata(:,4);
                        % 
                        % pseudosecdur(ff,:) = pseudodur;
                        % pseudosecspk(ff,:) = pseudospk;
                        % pseudosecmap_raw(ff,:) = pseudomap;
                        % pseudobasebins{ff} = lin_inds;
                        % pseudononoverlap{ff,:} = nonoverlap;
                        % % pseudotertbins{ff,:} = tertlinbin;

                        ff = ff + 1;
                    end
                    % disp([num2str(Args.NumShuffles) ' pseudo base fields for ' msvar{oo} ' field ' num2str(ii) ' took ' num2str(toc) 's']);
                    condpseudo_secdataperfield{ii} = {pseudobasebins pseudosecdur pseudosecspk pseudosecmap_raw};
                    secfieldnumbinperfield(ii) = condbasenumbin;
                    condbase_fieldlinbin{ii} = condbasefieldlinbin;
                    condfield_nonoverlapbinsperfield{ii} = pseudononoverlap;
                    % cond_tertlinbinsperfield{ii} = pseudotertbins;
                end
                data.(pairname_short).(msvar{oo}).condpseudo_secdataperfield = condpseudo_secdataperfield;
                data.(pairname_short).(msvar{oo}).secfieldnumbinset = secfieldnumbinperfield;
                data.(pairname_short).(msvar{oo}).condfield_inbinsset = condbase_fieldlinbin;
                data.(pairname_short).(msvar{oo}).condfield_nonoverlapbinset = condfield_nonoverlapbinsperfield;
                % condfield_tertlinbin{oo} = cond_tertlinbinsperfield;
            end

            % Smooth
            disp('Smoothing pseudopopulation of secondary maps ...');
            shuff = Args.NumShuffles;
            for oo = 1:size(msvar,2)
                pseudosecmaps_sm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                pseudosecSIC_sm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                % tertrate_orig = nan(shuff,data.(pairname_short).(msvar{oo}).sigfields);
                % tertrate_pseudo = nan(shuff,data.(pairname_short).(msvar{oo}).sigfields);
                for ii = 1:size(data.(pairname_short).(msvar{oo}).condbase_componentsperpx,1) % For each base field

                    % tic;
                    temp = data.(pairname_short).(msvar{oo}).condpseudo_secdataperfield{ii};
                    seclindur = temp{2};
                    seclinspk = temp{3};
                    seclinmap = temp{4};

                    switch msvar{2-oo+1}
                        case 'place'

                            gridSteps = [pc.Args.GridSteps pc.Args.GridSteps];
                            durG = cell2mat(lineartogrid(seclindur','place',gridSteps));
                            spkG = cell2mat(lineartogrid(seclinspk','place',gridSteps));

                            % Adaptive smoothing
                            [maps_smG,~,dur_smG] = adsmooth(durG,spkG,1e2);
                            secmap_sm = gridtolinear({maps_smG},'place',gridSteps);
                            secdur_sm = gridtolinear({dur_smG},'place',gridSteps);
                            secdur_sm(isnan(secdur_sm)) = 0;

                        case 'view'

                            gazeSections = sv.gazeSections;
                            binDepths = sv.binDepths;

                            % Assign linear bin to grid
                            [durG] = lineartogrid(seclindur','view',binDepths);
                            [spkG] = lineartogrid(seclinspk','view',binDepths);

                            % Pad sv map with 5 extra rows
                            n = 5;
                            padpillar = false;
                            [emptyfloorref_pad,~] = padsvmap(n,durG,gazeSections,padpillar);
                            padpillar = true;
                            [durGpad,retrievemap] = padsvmap(n,durG,gazeSections,padpillar);
                            [spkGpad,~] = padsvmap(n,spkG,gazeSections,padpillar);

                            % Adaptive smoothing of padded grids
                            alpha = 1e2; % Args.Alpha;
                            grid_smoothed_Gaze = cell(size(binDepths,1),1);
                            grid_smoothed_Gaze{1} = spkGpad{1}./durGpad{1}; % No need to smooth cue
                            grid_smoothed_Gaze{2} = spkGpad{2}./durGpad{2}; % No need to smooth hint
                            grid_smoothed_dur = cell(size(durGpad,1),1);
                            grid_smoothed_dur{1} = durGpad{1};
                            grid_smoothed_dur{2} = durGpad{2};
                            grid_ad_size = cell(size(durGpad,1),1);
                            grid_ad_size{1} = nan(size(durGpad{1}));
                            grid_ad_size{2} = nan(size(durGpad{1}));
                            for jj = 3:size(durGpad,1) % for each grid, floor onwards
                                
                %                 disp(['    ...Smoothing grid ' num2str(jj)]);
                                wip = ones(Args.NumShuffles,1);
                                gpdur1 = durGpad{jj};
                                preset_to_zeros = gpdur1(:,:,1);
                                preset_to_zeros(find(preset_to_zeros>0)) = 1;
                                preset_to_zeros(find(preset_to_zeros~=1)) = 0;
                                preset_to_zeros = ~preset_to_zeros;
                                preset_to_zeros = repmat(preset_to_zeros, [1,1,size(gpdur1,3)]);
                                
                                firing_counts_full1 = spkGpad{jj};
                                gpdur1(isnan(gpdur1)) = 0;
                                firing_counts_full1(isnan(firing_counts_full1)) = 0;
                                
                %                 to_compute = 1:0.5:Args.GridSteps/2; % unit bin is actually fspecial(...0.5)
                                to_compute = 1:0.5:(max(size(durGpad{jj}(:,:,1))))/2;
                                
                                possible = NaN(2,size(firing_counts_full1,1),size(firing_counts_full1,2),Args.NumShuffles + 1);
                                to_fill = NaN(size(possible,2), size(possible,3), size(possible,4));
                                to_fill(preset_to_zeros) = 0;
                                to_fill_smoothed_duration = NaN(size(possible,2), size(possible,3), size(possible,4));
                                to_fill_smoothed_duration(preset_to_zeros) = 0;
                                to_fill_size = NaN(size(possible,2), size(possible,3), size(possible,4));
                                to_fill_size(preset_to_zeros) = 0;
                                
                                for idx = 1:length(to_compute)
                                    
                                    f=fspecial('disk',to_compute(idx));
                                    f(f>=(max(max(f))/3))=1;
                                    f(f~=1)=0;
                                    
                                    possible(1,:,:,:) = repmat(imfilter(gpdur1(:,:,1), f, 'conv'), 1,1,Args.NumShuffles+1);
                                    possible(2,:,:,find(wip)) = imfilter(firing_counts_full1(:,:,find(wip)), f, 'conv');
                                    
                                    logic1 = squeeze(alpha./(possible(1,:,:,:).*sqrt(possible(2,:,:,:))) <= to_compute(idx));
                                    
                                    %debug
                                    %                         logic1(~logic1) = 1;
                                    
                                    slice1 = squeeze(possible(1,:,:,:));
                                    slice2 = squeeze(possible(2,:,:,:));
                                    
                                    to_fill(logic1 & isnan(to_fill)) = slice2(logic1 & isnan(to_fill))./slice1(logic1 & isnan(to_fill));
                                    to_fill_smoothed_duration(logic1 & isnan(to_fill_smoothed_duration)) = slice1(logic1 & isnan(to_fill_smoothed_duration));
                                    to_fill_size(logic1 & isnan(to_fill_size)) = to_compute(idx);
                                    
                                    
                                    remaining = sum(sum(sum(isnan(to_fill(:,:,:)))));
                %                     disp(['smoothed grid ' num2str(jj) ' with kernel size ' num2str(to_compute(idx)) ', leaving ' num2str(remaining) ' grids undone']);
                                    
                                    check = squeeze(sum(sum(isnan(to_fill),2),1));
                                    wip(check==0) = 0;
                                    
                                    if remaining == 0
                %                         disp('done');
                                        break;
                                    end
                                end
                                
                                to_fill(preset_to_zeros) = nan;
                                to_fill_size(preset_to_zeros) = nan;
                                grid_smoothed_Gaze{jj} = to_fill;
                                grid_smoothed_dur{jj} = to_fill_smoothed_duration;
                                grid_ad_size{jj} = to_fill_size;
                                
                            end
                            % disp(['Quick adaptive smoothing for ' msvar{2-oo+1} ' field ' num2str(ii) ' took ' num2str(toc) 's']);

                            % %%% Smoothing loop
                            % tic;
                            % secmap_sm = nan(size(seclindur'));
                            % secdur_sm = nan(size(seclindur'));
                            % for kk = 1:size(seclindur,1)
                            %     seclindursingle = seclindur(kk,:);
                            %     seclinspksingle = seclinspk(kk,:);
                            %     % Assign linear bin to grid bin - left to right, bottom to top
                            %     durG = lineartogrid(seclindursingle','view',binDepths);
                            %     spkG = lineartogrid(seclinspksingle','view',binDepths);
                            % 
                            %     % Pad sv map with 5 extra rows
                            %     n = 5;
                            %     padpillar = false;
                            %     [emptyfloorref_pad,~] = padsvmap(n,durG,gazeSections,padpillar);
                            %     padpillar = true;
                            %     [durGpad,retrievemap] = padsvmap(n,durG,gazeSections,padpillar);
                            %     [spkGpad,~] = padsvmap(n,spkG,gazeSections,padpillar);
                            % 
                            %     % Adaptive smooth
                            %     maps_adsmGpad = cell(size(durGpad));
                            %     dur_adsmGpad = cell(size(durGpad));
                            %     for jj = 1:size(binDepths,1)
                            %         if jj == 1 || jj == 2
                            %             maps_adsmGpad{jj} = spkGpad{jj}./durGpad{jj};
                            %             dur_adsmGpad{jj} = durGpad{jj};
                            %         else
                            %             [maps_adsmGpad{jj},spk_adsmG,dur_adsmGpad{jj}] = adsmooth(durGpad{jj},spkGpad{jj},1e2);
                            %         end
                            %     end
                            % 
                            %     % Unpad smoothed map
                            %     maps_smG = unpadsvmap(maps_adsmGpad,retrievemap,durG);
                            %     dur_smG = unpadsvmap(dur_adsmGpad,retrievemap,durG);
                            %     % Convert grid map back to linear sv map
                            %     secmap_sm(:,kk) = gridtolinear(maps_smG,'view',binDepths);
                            %     secdur_sm(:,kk) = gridtolinear(dur_smG,'view',binDepths);
                            %     secdur_sm(isnan(secdur_sm(:,kk)),kk) = 0;
                            % 
                            % end
                        case 'headdirection'
                            n = 5;
                            % Smooth
            %                     rateG = seclinmap'; % raw
                            secmap_sm = smoothdir(seclinmap',n,cr.headdirectionbins);
            %                     maps_smG = secmap_sm; % no difference between linear and grid map for HD
                            secdur_sm = smoothdir(seclindur',n,cr.headdirectionbins);
                    end
                    % disp(['Smoothing for ' msvar{2-oo+1} ' field ' num2str(ii) ' took ' num2str(toc) 's']);

                    if ~strcmp(msvar{2-oo+1},'headdirection')
                        % Calculate SIC from adaptively smoothed map
                        sic_sm = skaggs_sic(secmap_sm,secdur_sm);
                    else
                        % Rayleigh vector
                        binSize=(pi*2)/length(secmap_sm(:,1));
                        binAngles=(0:binSize:( (359.5/360)*2*pi )) + binSize/2;
                        binWeights=secmap_sm'./(max(secmap_sm',[],2));
                        S=nansum( sin(binAngles).*binWeights , 2);
                        C=nansum( cos(binAngles).*binWeights , 2);
                        R=sqrt(S.^2+C.^2);
                        sic_sm=R./nansum(binWeights,2);
                        sic_sm = sic_sm';
                    end

                    pseudosecmaps_sm{ii,1} = secmap_sm';
                    pseudosecSIC_sm{ii,1} = sic_sm';

                    % % Get mean firing rate of bins outside of conditioned
                    % % fields of original map
                    % rate_orig = nan(shuff,1);
                    % rate_pseudo = nan(shuff,1);
                    % for kk = 1:shuff
                    %     px = condfield_tertlinbin{oo}{ii}{kk};
                    %     rate_orig(kk) = nanmean(data.(pairname_short).(msvar{oo}).condbase_map_sm{ii}(px));
                    %     rate_pseudo(kk) = nanmean(secmap_sm(px,kk));
                    % end
                    % 
                    % tertrate_orig(:,ii) = rate_orig;
                    % tertrate_pseudo(:,ii) = rate_pseudo;

  
                end

                % Store data
                data.(pairname_short).(msvar{oo}).pseudosecmaps_sm = pseudosecmaps_sm;
                data.(pairname_short).(msvar{oo}).pseudosecSIC_adsm = pseudosecSIC_sm;
                % data.(pairname_short).(msvar{oo}).tertrate_orig = tertrate_orig;
                % data.(pairname_short).(msvar{oo}).tertrate_pseudo = tertrate_pseudo;
            end

       
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


% % Find grid in spatial view frame for 1 pixel
% function [gridnum,x,y] = findgrid(px,objtype)
% % returns grid number and plot coords (x goes left to right, y goes bottom
% % to top)
% 
% switch objtype
%     case 'place'
%         mapLdummy = 1:1600;
%         gridnum = 3;
%         temp = flipud(reshape(mapLdummy, 40, 40)');
%     case 'view'
%         mapLdummy = 1:5122;
%         if px == 1 % Cue
%             gridnum = 1;
%             x = 1;
%             y = 1;
%         elseif px == 2 % Hint
%             gridnum = 2;
%             x = 1; 
%             y = 1;
%         elseif px >= 3 && px <= 1602 % Floor
%             gridnum = 3;
%             temp = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');
%         elseif px >= 1603 && px <= 3202 % Ceiling
%             gridnum = 4;
%             temp = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');
%         elseif px >= 3203 && px <= 4482 % Walls
%             gridnum = 5;
%             temp = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');
%         elseif px >= 4483 && px <= 4642 % Pillar 1
%             gridnum = 6;
%             temp = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
%         elseif px >= 4643 && px <= 4802 % Pillar 2
%             gridnum = 7;
%             temp = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
%         elseif px >= 4803 && px <= 4962 % Pillar 3
%             gridnum = 8;
%             temp = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
%         elseif px >= 4963 && px <= 5122 % Pillar 4
%             gridnum = 9;
%             temp = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
%         end
%     case 'headdirection'
%         mapLdummy = 1:60;
%         gridnum = 1;
%         temp = flipud(reshape(mapLdummy, 60, 1)');
% end
% [y,x] = find(temp == px);
% y = size(temp,1)-y+1;
