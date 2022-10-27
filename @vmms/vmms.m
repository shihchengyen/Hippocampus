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
                'RequiredFile','spiketrain.mat', 'GridSteps',40, ...
                'UseCorr',1,'FieldThr',0.7,'FieldSplitThr',0.75,'NumShuffles',10000, ...
                'FieldThrPseudo',0.6,'FieldSplitThrPseudo',0.7);
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
    Args.spatialvarpairs = spatialvarpairs;
    
    %% Selectivity of whole cell
    
    data.placesel = false;
    data.viewsel = false;
    data.headdirectionsel = false;
%     data.mixsel = false;
    data.discard = false;
    
    cd('FiltVel/1px');
    pc = load('vmpc.mat');
    pc = pc.vmp.data;
    sv = load('vmsv.mat');
    sv = sv.vms.data;
    hd = load('vmhd.mat');
    hd = hd.vmd.data;
    cr = load('vmcorr.mat');
    cr = cr.vmcorr.data;
    % Get cell SI
    if ~Args.UseCorr % If using original pc/sv objects to define fields
        pSI = pc.SIC_adsm;
        vSI = sv.SIC_adsm;
        hSI = hd.crit_sm;
    else
        %%% PATCH (we have two instances each of pair corrections, which to pick?)
        pSI = cr.pv.SIC_sm_corrp;
        vSI = cr.pv.SIC_sm_corrv;
        hSI = cr.ph.SIC_sm_corrh;
    end
    % Get number of spikes retained after filtering
    numspk = sum(pc.spk_raw);
    data.filtspkcount = numspk; 
    if numspk < 100
        data.discard = true;
    end

    % Get population SI threshold without the cells that have < 100 spikes. Note: population is defined by size of combined object, not cell list
    c_pc = load('/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Combined Objects/FiltVel/1px/c_vmpc.mat');
    c_pc = c_pc.vmp.data;
    c_sv = load('/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Combined Objects/FiltVel/1px/c_vmsv.mat');
    c_sv = c_sv.vms.data;
    c_hd = load('/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/Combined Objects/FiltVel/1px/c_vmhd.mat');
    c_hd = c_hd.vmd.data;
    cell_indDiscardP = ismember(c_pc.origin,c_pc.origin(sum(c_pc.spk_raw,2)<100)); %%% Change in the future to a fixed var that is same in pc/sv
    cell_indDiscardV = ismember(c_sv.origin,c_sv.origin(sum(c_pc.spk_raw,2)<100));
    pSIset = c_pc.SICsh_adsm;
    vSIset = c_sv.SICsh_adsm;
    hSIset = c_hd.critsh_sm;
    cell_numDiscardP = find(cell_indDiscardP);
    cell_numDiscardV = find(cell_indDiscardV);
    for dd = 1:size(cell_numDiscardP,1)
        ind = cell_numDiscardP(dd);
        pSIset((ind-1)*pc.Args.NumShuffles+1:ind*pc.Args.NumShuffles) = nan;
    end
    for dd = 1:size(cell_numDiscardV,1)
        ind = cell_numDiscardV(dd);
        vSIset((ind-1)*sv.Args.NumShuffles+1:ind*sv.Args.NumShuffles) = nan;
    end
    pSIthr = prctile([pSI; pSIset],95); % Population threshold only
    vSIthr = prctile([vSI; vSIset],95); % Population threshold only
    hSIthr = prctile([hSI; hSIset],95);
%     pSIthr = max([prctile(c_pc.SICsh,95) prctile(pc.SICsh,95)]); % Both population and cell threshold
%     vSIthr = max([prctile(c_sv.SICsh,95) prctile(sv.SICsh,95)]); % Both population and cell threshold
    cd(cwd);
    
    % Get peak rate of adaptive-smooted rate map
    peakrate_pc = nanmax(pc.maps_adsm);
    peakrate_sv = nanmax(sv.maps_adsm); % exclude cue and hint???
    peakrate_hd = nanmax(hd.maps_sm);
    
    % Get selectivity of this cell
    if pSI>pSIthr && numspk>=100 && peakrate_pc>=0.7
        data.placesel = true;
    end
    if vSI>vSIthr && numspk>=100 && peakrate_sv>=0.7
        data.viewsel = true;
    end
    if hSI>hSIthr && numspk>=100 && peakrate_hd>=0.7
        data.headdirectionsel = true;
    end
    
    for pair = 1:size(spatialvarpairs,2)
        
        if pair == 2
            disp('stop');
        end
        msvar = spatialvarpairs{pair}; % For now, do in pairs only
        msvar_short = {msvar{1}(1),msvar{2}(1)};
        pairname_short = [msvar{1}(1) msvar{2}(1)];
        if data.([msvar{1} 'sel']) && data.([msvar{2} 'sel'])
            data.(pairname_short).mixsel = true;
        else 
            data.(pairname_short).mixsel = false;
        end
    
        % Set up associated standard variables
        msobj = cell(1,2);
        crit = cell(1,2);
        maptype = cell(1,2);
        gridSize = cell(1,2);
        for oo = 1:2 % for each var in pair
            switch msvar{oo}
                case 'place'
                    msobj{oo} = 'pc';
                    crit{oo} = 'crit';
                    maptype{oo} = '_sm';
                case 'view'
                    msobj{oo} = 'sv';
                    crit{oo} = 'crit';
                    maptype{oo} = '_sm';

                case 'headdirection'
                    msobj{oo} = 'hd';
                    crit{oo} = 'crit';
                    maptype{oo} = '_sm';
            end
        end


        %% Work out mixed selective properties wuth each spatialvar as base
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
            data.(pairname_short).(msvar{oo}).secmaps_raw = {}; % Fix name to align with base
            data.(pairname_short).(msvar{oo}).seccomponents_perbasepx = {};
            data.(pairname_short).(msvar{oo}).secmaps_dist = {};
            data.(pairname_short).(msvar{oo}).sec_distratio = [];
            data.(pairname_short).(msvar{oo}).sec_distcorr = [];
            data.(pairname_short).(msvar{oo}).sec_infieldrates = {};
            data.(pairname_short).(msvar{oo}).sec_outfieldrates = {};
            data.(pairname_short).(msvar{oo}).secmaps_sm = {};
        %         data.(pairname_short).(msvar{oo}).secmaps_adsm = {};
        %         data.(pairname_short).(msvar{oo}).secmaps_bcsm = {};
        %         data.(pairname_short).(msvar{oo}).secmaps_dksm = {};
            data.(pairname_short).(msvar{oo}).secSIC_sm = [];
        %         data.(pairname_short).(msvar{oo}).secSIC_adsm = [];
        %         data.(pairname_short).(msvar{oo}).secSIC_bcsm = [];
        %         data.(pairname_short).(msvar{oo}).secSIC_dksm = [];
            data.(pairname_short).(msvar{oo}).secsigcondfields = []; 
            data.(pairname_short).(msvar{oo}).seccond_gridnum = [];
            data.(pairname_short).(msvar{oo}).seccond_fieldcoord = {};
            data.(pairname_short).(msvar{oo}).seccond_fieldlinbin = {};
            data.(pairname_short).(msvar{oo}).seccond_fieldsizepercent = {};
            data.(pairname_short).(msvar{oo}).seccond_fieldoverlapind = {};
            data.(pairname_short).(msvar{oo}).seccond_fieldoverlap = {};
            data.(pairname_short).(msvar{oo}).pseudosecmaps_adsm = {};
            data.(pairname_short).(msvar{oo}).pseudosecSIC_adsm = {};
            data.(pairname_short).(msvar{oo}).pseudosecdataperfield = {};
            data.(pairname_short).(msvar{oo}).tertrate_orig = {};
            data.(pairname_short).(msvar{oo}).tertrate_pseudo = {};

        %         for vv = 1:size(msvar,2)-1 % This is for when there is more than 1 secondary variable
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'mapLsm']) = [];
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'maps_raw']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'components_perbasepx']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'maps_dist']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) '_distratio']) = [];
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) '_distcorr']) = [];
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) '_infieldrates']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) '_outfieldrates']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'maps_adsm']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'maps_bcsm']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'maps_dksm']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'SIC_adsm']) = [];
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'SIC_bcsm']) = [];
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'SIC_dksm']) = [];
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'sigcondfields']) = []; 
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'cond_gridnum']) = [];
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'cond_fieldcoord']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'cond_fieldlinbin']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'cond_fieldsizepercent']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'cond_fieldoverlapind']) = {};
        %             data.(pairname_short).(msvar{oo}).(['sec' num2str(vv) 'cond_fieldoverlap']) = {};
        %             data.(pairname_short).(msvar{oo}).(['pseudosec' num2str(vv) 'maps_adsm']) = {};
        %             data.(pairname_short).(msvar{oo}).(['pseudosec' num2str(vv) 'SIC_adsm']) = {};
        %             data.(pairname_short).(msvar{oo}).(['pseudosec' num2str(vv) 'dataperfield']) = {};
        %             data.(pairname_short).(msvar{oo}).(['tert' num2str(vv) 'rate_orig']) = {};
        %             data.(pairname_short).(msvar{oo}).(['tert' num2str(vv) 'rate_pseudo']) = {};
        %         end

        %         if data.([msvar{oo} 'sel']) 

                % Load pv object
                cd ..; cd ..; cd ..;
                pv = load('1vmpv.mat');
                pv = pv.pv;
                stc = pv.data.sessionTimeC;
                stc(:,5) = [diff(stc(:,1)); 0]; % dwell time

                % Bin spikes
                disp('Binning spikes ...');
                cd(cwd);
                spiketrain = load('spiketrain.mat');
                spiketimes = spiketrain.timestamps ./ 1000; % in seconds
                spiketimes(spiketimes(:,1) < stc(1,1),:) = []; 
                binned = histcounts(spiketimes, stc(:,1))';
                stc(:,6) = [binned; 0];

                %% Filter pv with same criteria used in vmpc/sv
                disp('Filtering ...');
                conditions = ones(size(stc,1),1);
                if pc.Args.UseAllTrials == 0
                    conditions = conditions & pv.data.good_trial_markers;
                end
                if pc.Args.ThresVel > 0
                    conditions = conditions & get(pv,'SpeedLimit',pc.Args.ThresVel); % pv here needs to be in object form, not structure
                end
                if pc.Args.UseMinObs
                    bins_sieved_p = pv.data.place_good_bins;
                    bins_removed_p = setdiff(1:size(pv.data.place_intervals_count,1),bins_sieved_p);
                    bins_sieved_v = pv.data.view_good_bins;
                    bins_removed_v = setdiff(1:size(pv.data.view_intervals_count,1),bins_sieved_v);
                    conditions = conditions & (pv.data.pv_good_rows); % Make sure maps take into account both place and view filters
                else
                    bins_sieved_p = 1:(pc.Args.GridSteps * pc.Args.GridSteps);
                    bins_removed_p = [];
                    bins_sieved_v = 1:size(pv.data.view_intervals_count,1);
                    bins_removed_v = [];
                end
                stc(conditions ~= 1,:) = []; % Filter stc. This will be kept for drawing pixels from in field analysis

                % Backfill duration and spikes for view bins - following code from
                % vmcorr
                stcvars = {'time','place','headdirection','view','dur','spike'};
        %             for ss = 2:size(msvar,2)
                stcfill = nan(size(stc)); stcfill(:,5:6) = 0; % To avoid issues where last dur and spk count are 0 and missing next data point to back fill. 
                sec_durations = zeros(cr.([msvar{2-oo+1} 'bins']),cr.([msvar{oo} 'bins']));
                sec_spikes = zeros(cr.([msvar{2-oo+1} 'bins']),cr.([msvar{oo} 'bins']));
                base_durations = zeros(1,cr.([msvar{oo} 'bins']));
                base_spikes = zeros(1,cr.([msvar{oo} 'bins']));
                xbase = strcmp(stcvars,msvar{oo});
                xsec = strcmp(stcvars,msvar{2-oo+1});
                for i = 1:cr.([msvar{oo} 'bins'])

                    inds = stc(:,xbase)==i;
                    indnums = find(inds);
                    subsample = [stc(inds,:)]; % [time place hd view dur spk]

                    % Consider only samples where both place and view are sampled
                    if strcmp(msvar{2-oo+1},'view')
                        indnums(isnan(subsample(:,xsec)),:) = [];
                        subsample(isnan(subsample(:,xsec)),:) = [];
                    else
                        indnums(subsample(:,xsec)<1,:) = [];
                        subsample(subsample(:,xsec)<1,:) = [];
                    end
                    if isempty(subsample)
                        continue;
                    end

                    % Get spikes and duration for place only
                    base_durations(1,i) = sum(subsample(:,5));
                    base_spikes(1,i) = sum(subsample(:,6));

                    % back-filling spikes for view
                    subsample(subsample(:,6)==0,6) = nan;
                    subsample(:,7) = subsample(:,5)~=0;
                    subsample(isnan(subsample(:,6)) & subsample(:,7), 6) = 0;
                    subsample(:,7) = [];
                    subsample(:,6) = fillmissing(subsample(:,6), 'next');
                    % back-filling time for view
                    subsample(subsample(:,5)==0,5) = nan;
                    subsample(:,5) = fillmissing(subsample(:,5), 'next');

                    % remove bad view spots
                    indnums(isnan(subsample(:,xsec)),:) = [];
                    subsample(isnan(subsample(:,xsec)),:) = [];

                    % Put backfill into sessionTimeC array
                    stcfill(indnums,:) = subsample;

                    % padding with 5122 bin
                    subsample = [subsample; [0 1600 60 5122 0 0]];

                    % sum durations
                    sec_durations(:,i) = accumarray(subsample(:,xsec), subsample(:,5),[],[],NaN);
                    % sum spikes
                    sec_spikes(:,i) = accumarray(subsample(:,xsec), subsample(:,6),[],[],NaN); 
                end
                stcfill(isnan(stcfill(:,1)),:) = []; % Remove all rows without both place and view data
                % Remove low obs bins
                stcfill(~ismember(stcfill(:,2),bins_sieved_p) | ~ismember(stcfill(:,3),bins_sieved_v),:) = [];

                %%%%%% FIXXXXXXX
                if isnan(stcfill(end,5) )
                    stcfill(end,5) = 0;
                end
                if isnan(stcfill(end,6))
                    stcfill(end,6) = 0;
                end

                % Replace nan duration values with zero
                base_durations(isnan(base_durations)) = 0; % Necessary because NaNs seem to mess up the smoothing. 
                sec_durations(isnan(sec_durations)) = 0; 
                % Remove low obs bins
                if ~strcmp(msvar{oo},'headdirection')
                    base_durations(eval(['bins_removed_' msvar_short{oo}])) = 0;
                    sec_durations(:,eval(['bins_removed_' msvar_short{oo}])) = 0;
                    base_spikes(eval(['bins_removed_' msvar_short{oo}])) = 0;
                    sec_spikes(:,eval(['bins_removed_' msvar_short{oo}])) = 0;
                end
                if ~strcmp(msvar{2-oo+1},'headdirection')
                    sec_durations(eval(['bins_removed_' msvar_short{2-oo+1}]),:) = 0;
                    sec_spikes(eval(['bins_removed_' msvar_short{2-oo+1}]),:) = 0;
                end

                % Make maps from filtered pv object
                base_map = base_spikes./base_durations;
                sec_map = nansum(sec_spikes,2)./nansum(sec_durations,2);
                sec_map(nansum(sec_durations,2)==0) = nan; % restore nans to unvisited bins
                %%% Store pv maps
                data.(pairname_short).(msvar{oo}).pvmap = base_map;
                data.(pairname_short).(msvar{2-oo+1}).pvmap = sec_map';
                if numspk ~= sum(base_spikes)
                    disp('filtered spikes different between this object and vmpc');
                end


                %% Find base fields 
        %             for oo = 1:size(msobj,2)
                disp(['Delineating ' msvar{oo} ' fields ...']);
                baseobj = eval(msobj{oo});
                secobj = eval(msobj{2-oo+1});
                if Args.UseCorr
                    basemapLsm = cr.([msvar_short{1} msvar_short{2}]).(['maps_sm_corr' msvar_short{oo}]);
                    secmap_sm = cr.([msvar_short{1} msvar_short{2}]).(['maps_sm_corr' msvar_short{2-oo+1}]);
                    basemapLrw = cr.([msvar_short{1} msvar_short{2}]).(['maps_raw_corr' msvar_short{oo}]);
                else
                    basemapLsm = baseobj.(['maps' maptype{oo}]);
                    secmap_sm = secobj.(['maps' maptype{2-oo+1}]);
                    basemapLrw = baseobj.maps_raw;
                end
                SI = eval([msvar_short{oo} 'SI']);
                SIthr = eval([msvar_short{oo} 'SIthr']);
                gridSize{oo} = cr.([msvar_short{1} msvar_short{2}]).([msvar{oo} 'binDepths']);
                basemapGsm = lineartogrid(basemapLsm',msvar{oo},gridSize{oo});
                dummygrid = lineartogrid((1:size(basemapLsm,2))',msvar{oo},gridSize{oo});
                basemapGrw = lineartogrid(basemapLrw',msvar{oo},gridSize{oo});

                switch msvar{oo}
                    case 'place'
                        peakrate_full = nanmax(basemapLsm);
                        peakrate_subset = peakrate_full;
                        prI = 1;
                        fieldsizethreshold = 15;
                    case 'view'
                        for ii = 1:size(basemapGsm,1)
                            maxset(ii) = nanmax(reshape(basemapGsm{ii},size(basemapGsm{ii},1)*size(basemapGsm{ii},2),1));
                        end
                        [peakrate_full,prI] = max(maxset); % Max including cue/hint
                        peakrate_subset = max(maxset(3:end)); % Max excluding cue/hint
                        fieldsizethreshold = 15;
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
                if strcmp(msvar{oo},'headdirection')
                    disp('stop')
                end
%                 if data.([msvar{oo} 'sel']) % For use when only considering selective cells. If commented out, will consider any and every field regardless of overall selectivity
                    for gg = 1:size(basemapGsm,1)
                        % Skip cue and hint view maps 
                        if size(basemapGsm{gg},1) == 1 && size(basemapGsm{gg},2) == 1
                            continue;
                        end
                        % if there is too little variance in firing rates, reject
                        if nanstd(basemapLsm)/nanmean(basemapLsm) < 0.25
                            continue;
                        end
                        % Find fields with at least 1 pixel of > 70% peak
                        % rate and > 1 std dev from mean
                        ind_fields = basemapGsm{gg} > Args.FieldThr*peakrate_subset & ...
                            basemapGsm{gg} > (nanmean(basemapLsm)+nanstd(basemapLsm));
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
                                    if sum(inds(:)) >= ((size(pc.maps_raw,2)-4*8*8)/4) 
                %                         if sum(inds(:)) >= sum(sum(~isnan(basemapGsm{gg})))/4 % 
                                        subinds = inds & basemapGsm{gg} > Args.FieldSplitThr*peakrate_subset;
                                        [sublabel,subcount] = bwlabel(subinds,4);
                                        sublabel(subinds) = sublabel(subinds) + fieldcount;
                                        fieldcount = fieldcount + subcount;
                                        fieldlabel(inds) = 0; % There will be some px that are not part of new fields cos of different rate thresholds in subfields
                                        fieldlabel(inds) = sublabel(inds);
                                    end
                                case 'view'
                                     if sum(inds(:)) >= ((size(sv.maps_raw,2)-4*8*8)/4) 
                %                         if sum(inds(:)) >= sum(sum(~isnan(basemapGsm{gg})))/4 % 
                                        subinds = inds & basemapGsm{gg} > Args.FieldSplitThr*peakrate_subset;
                                        [sublabel,subcount] = bwlabel(subinds,4);
                                        sublabel(subinds) = sublabel(subinds) + fieldcount;
                                        fieldcount = fieldcount + subcount;
                                        fieldlabel(inds) = 0; % There will be some px that are not part of new fields cos of different rate thresholds in subfields
                                        fieldlabel(inds) = sublabel(inds);
                                     end
                                case 'headdirection'
                                     if sum(inds(:)) >= ((size(hd.maps_raw,2))/2) 
                %                         if sum(inds(:)) >= sum(sum(~isnan(basemapGsm{gg})))/4 % 
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
                            % Make sure field size is > 15 bins 
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
%                 end

                % If no significant fields, skip to next var
                if count == 0
                    % create nptdata so we can inherit from it
                    data.(pairname_short).(msvar{oo}).sigfields = 0; % Selective but no sig fields (if cell non-selective, sigfields = nan, see below)
                    data.(pairname_short).(msvar{oo}).discardfieldnum = discardfieldnum;
                    data.(pairname_short).(msvar{oo}).discardfieldreason = discardfieldreason;
                    data.(pairname_short).(msvar{oo}).SI = SI;
                    data.(pairname_short).(msvar{oo}).SIthr = SIthr;
                    data.(pairname_short).(msvar{oo}).basemapLsm = basemapLsm;
                    data.(pairname_short).(msvar{oo}).basemapLrw = basemapLrw;
                    data.(pairname_short).(msvar{oo}).secmapLsm = secmap_sm;
                    data.(pairname_short).(msvar{oo}).basemapGsm = basemapGsm;
                    data.(pairname_short).(msvar{oo}).basemapGrw = basemapGrw;
                    data.(pairname_short).(msvar{oo}).dummygrid = dummygrid;
                    continue;
                else 
                    data.(pairname_short).(msvar{oo}).sigfields = min([count maxfieldcount]);
                    data.(pairname_short).(msvar{oo}).discardfieldnum = discardfieldnum;
                    data.(pairname_short).(msvar{oo}).discardfieldreason = discardfieldreason;
                end

                % Sort fields
                [fieldmaxrate_sm,I] = sort(fieldmaxrate_sm,'descend');
                fieldmaxrate_rw = fieldmaxrate_rw(I);
                gridnum = gridnum(I);
                fieldcoord = fieldcoord(I);
                fieldlinbin = fieldlinbin(I);

                %% Get secondary pixels for each base field 
                disp(['Parsing secondary pixels for each base ' msvar{oo} ' field ...']);

                seccomponents_perbasepx = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
                secmap_perbasefield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
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
                    for pp = 1:size(fieldcoord{ii},1)
                        %Get corresponding secondary pixels from pv object
                        secpx = [];
                        ind_pv = stcfill(:,xbase) == fieldlinbin{ii}(pp); % base var px
                        secpx(:,1) = stcfill(ind_pv,xsec); % sec var px
                        secpx(:,2) = stcfill(ind_pv,5); % dur
                        secpx(:,3) = stcfill(ind_pv,6); % spikes
        %                         switch msvar{oo}
        %                             case 'place'
        %                                 ind_pv = stcfill(:,2) == fieldlinbin{ii}(pp);
        %                                 secpx(:,1) = stcfill(ind_pv,4); % view px
        %                                 secpx(:,2) = stcfill(ind_pv,5); % dur
        %                                 secpx(:,3) = stcfill(ind_pv,6); % spikes
        %                             case 'view'
        %                                 ind_pv = stcfill(:,4) == fieldlinbin{ii}(pp);
        %                                 secpx(:,1) = stcfill(ind_pv,2); % place px
        %                                 secpx(:,2) = stcfill(ind_pv,5); % dur
        %                                 secpx(:,3) = stcfill(ind_pv,6); % spikes
        %                         end
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
                    secdata = usecpxs; 
                    secdata(:,2:3) = NaN;
                    full_dur = zeros(size(basemapLrw,2),size(secmap_sm,2));
                    full_spk = full_dur;
                    % Consolidate occupancy and spikes across base field
                    for pp = 1:size(fieldcoord{ii},1)
                        tempbins = nan(size(secdata,1),2);
                        % debug
                        if fieldlinbin{ii}(pp) == 1300
                            disp(fieldlinbin{ii}(pp));
                        end
                        tempbins(ismember(secdata(:,1),rate_components_perbasepx{pp}(:,1)),1) = rate_components_perbasepx{pp}(:,2); % Duration 
                        tempbins(ismember(secdata(:,1),rate_components_perbasepx{pp}(:,1)),2) = rate_components_perbasepx{pp}(:,3); % Spikes 
                        full_dur(fieldlinbin{ii}(pp),rate_components_perbasepx{pp}(:,1)) = rate_components_perbasepx{pp}(:,2);
                        full_spk(fieldlinbin{ii}(pp),rate_components_perbasepx{pp}(:,1)) = rate_components_perbasepx{pp}(:,3);
                        secdata(:,2) = nansum( [secdata(:,2) tempbins(:,1)] ,2); % Sum duration across session
                        secdata(:,3) = nansum( [secdata(:,3) tempbins(:,2)] ,2); % Sum spikes across session
                    end
                    secdata(:,4) = secdata(:,3)./secdata(:,2);


                    seccomponents_perbasepx{ii,1} = rate_components_perbasepx;
                    secmap_perbasefield{ii,1} = secdata;

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
                    base_array_orig = basemapLrw(1,sort(fieldlinbin{ii}));
                    sec_array_orig = secdata(:,4);
                    full_dur1 = full_dur;
                    full_spk1 = full_spk;
                    full_dur1(~ismember(1:size(basemapLrw,2),fieldlinbin{ii}),:) = [];
                    full_dur1(:,~ismember(1:size(secmap_sm,2),usecpxs)) = [];
                    full_spk1(~ismember(1:size(basemapLrw,2),fieldlinbin{ii}),:) = [];
                    full_spk1(:,~ismember(1:size(secmap_sm,2),usecpxs)) = [];
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
                    temp2 = nan(size(secmap_sm));
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
                data.(pairname_short).(msvar{oo}).SI = SI;
                data.(pairname_short).(msvar{oo}).SIthr = SIthr;
                data.(pairname_short).(msvar{oo}).basemapLsm = basemapLsm;
                data.(pairname_short).(msvar{oo}).basemapLrw = basemapLrw;
                data.(pairname_short).(msvar{oo}).secmapLsm = secmap_sm;
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
                data.(pairname_short).(msvar{oo}).secmaps_raw = secmap_perbasefield;
                data.(pairname_short).(msvar{oo}).seccomponents_perbasepx = seccomponents_perbasepx;
                data.(pairname_short).(msvar{oo}).basemaps_dist = basemaps_dist;
                data.(pairname_short).(msvar{oo}).secmaps_dist = secmaps_dist;
                data.(pairname_short).(msvar{oo}).base_distratio = base_distratio;
                data.(pairname_short).(msvar{oo}).sec_distratio = sec_distratio;
                data.(pairname_short).(msvar{oo}).base_distcorr = base_distcorr;
                data.(pairname_short).(msvar{oo}).sec_distcorr = sec_distcorr;

                clear fieldmaxrate_sm; clear fieldmaxrate_rw; clear gridnum; clear fieldcoord; clear fieldlinbin; clear secmap_perbasefield;
                clear seclinmap; clear rate_components_perbasepx; clear secdata; clear secpx;
                clear seccomponents_perbasepx; clear dummygrid; clear SI; clear SIthr; clear sec_distratio; clear secmaps_dist;
                clear basemaps_dist; clear base_distratio; clear base_distcorr; clear sec_distcorr;

        %             end
        %             end
        end



        %% Smooth secondary maps
        disp('Smoothing secondary maps and identifying secondary fields ...');

        for oo = 1:size(msvar,2)
            secmaps_sm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
        %         secmaps_adsm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
        %         secmaps_bcsm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
        %         secmaps_dksm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            secSIC_sm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
        %         secSIC_adsm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
        %         secSIC_bcsm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
        %         secSIC_dksm = nan(data.(pairname_short).(msvar{oo}).sigfields,1);

            secsigconfdields = zeros(data.(pairname_short).(msvar{oo}).sigfields,1);
            seccond_gridnum = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            seccond_fieldcoord = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            seccond_fieldlinbin = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            seccond_fieldsizepercent = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            seccond_fieldoverlapind = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            seccond_fieldoverlap = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            for ii = 1:data.(pairname_short).(msvar{oo}).sigfields
                % Get secondary pixel map
                secdata = data.(pairname_short).(msvar{oo}).secmaps_raw{ii};
                seclinmap = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                seclinmap(1,secdata(:,1)) = secdata(:,4);
        %                 secgridmap = lineartogrid(seclinmap',msvar{2-oo+1},gridSize{2-oo+1});
        %                 dummygridsec = data.(pairname_short).(msvar{2-oo+1}).dummygrid;
                seclindur = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                seclinspk = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                seclindur(1,secdata(:,1)) = secdata(:,2); 
                seclinspk(1,secdata(:,1)) = secdata(:,3); 
                %% Smooth
                switch msvar{2-oo+1}
                    case 'place'

        %                     gridSteps = [pc.Args.GridSteps pc.Args.GridSteps];
                        gridSteps = gridSize{2-oo+1};
                        durG = cell2mat(lineartogrid(seclindur','place',gridSteps));
                        spkG = cell2mat(lineartogrid(seclinspk','place',gridSteps));
                        rateG = cell2mat(lineartogrid(seclinmap','place',gridSteps));

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
        %                     binDepths = sv.binDepths;
                        binDepths = gridSize{2-oo+1};

                        % Assign linear bin to grid bin - left to right, bottom to top
                        durG = lineartogrid(seclindur','view',binDepths);
                        spkG = lineartogrid(seclinspk','view',binDepths);
                        rateG = lineartogrid(seclinmap','view',binDepths);

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
        %                             [maps_adsmGpad{jj},spk_adsmG,dur_adsmG] = adsmooth(durGpad{jj},spkGpad{jj},sv.Args.Alpha);
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
                        rateG = seclinmap'; % raw
                        secmap_sm = smoothdir(seclinmap',n,cr.headdirectionbins);
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
                secmaps_sm{ii,1} = secmap_sm';
        %             secmaps_bcsm{ii,1} = secmapLsm';
        %             secmaps_dksm{ii,1} = secmapdksm';
                secSIC_sm(ii,1) = sic_sm;
        %             secSIC_bcsm(ii,1) = sic_bcsm;
        %             secSIC_dksm(ii,1) = sic_dksm;

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

                %% Find max 3 secondary fields
                count = 0;
                maxfieldcount = 3;
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
                        if size(maps_smG{gg},1) == 1
                            continue;
                        end
                        % Find fields with at least 1 pixel of > 70% peak rate
        %                         ind_fields = maps_adsmG{gg} > Args.FieldThrPseudo*peakrate_sec;
                        ind_fields = maps_smG{gg} > (meanrate_sec+(2*stdrate_sec));
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
                    % create nptdata so we can inherit from it
        %                     data.(pairname_short).(msvar{oo}).secsigconfdields = 0; 
                    secsigconfdields(ii) = 0;
                    continue;
                else 
                    secsigconfdields(ii) = min([count maxfieldcount]);
        %                     data.(pairname_short).(msvar{oo}).secsigconfdields = min([count maxfieldcount]);
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

                seccond_gridnum{ii,1} = gridnum;
                seccond_fieldcoord{ii,1} = fieldcoord;
                seccond_fieldlinbin{ii,1} = fieldlinbin;
                seccond_fieldsizepercent{ii,1} = fieldsizepercent;
                seccond_fieldoverlapind{ii,1} = fieldoverlapind;
                seccond_fieldoverlap{ii,1} = fieldoverlap;
            end

            data.(pairname_short).(msvar{oo}).secsigconfdields  = secsigconfdields;
            data.(pairname_short).(msvar{oo}).seccond_gridnum = seccond_gridnum;
            data.(pairname_short).(msvar{oo}).seccond_fieldcoord = seccond_fieldcoord;
            data.(pairname_short).(msvar{oo}).seccond_fieldlinbin = seccond_fieldlinbin;
            data.(pairname_short).(msvar{oo}).seccond_fieldsizepercent = seccond_fieldsizepercent;
            data.(pairname_short).(msvar{oo}).seccond_fieldoverlapind = seccond_fieldoverlapind;
            data.(pairname_short).(msvar{oo}).seccond_fieldoverlap = seccond_fieldoverlap;
            data.(pairname_short).(msvar{oo}).secmaps_sm = secmaps_sm;
        %         data.(pairname_short).(msvar{oo}).secmaps_bcsm = secmaps_bcsm;
        %         data.(pairname_short).(msvar{oo}).secmaps_dksm = secmaps_dksm;
            data.(pairname_short).(msvar{oo}).secSIC_sm = secSIC_sm;
        %         data.(pairname_short).(msvar{oo}).secSIC_bcsm = secSIC_bcsm;
        %         data.(pairname_short).(msvar{oo}).secSIC_dksm = secSIC_dksm;

        end

                %% Stats for infield and outfield firing of sec field
        for oo = 1:size(msvar,2)
            disp(['Comparing secondary infield and outfield firing for base ' msvar{oo} ' field ...'])
            % Output variables

            sec_infieldrates = cell(size(data.(pairname_short).(msvar{oo}).seccomponents_perbasepx,1),1);
            sec_outfieldrates = cell(size(data.(pairname_short).(msvar{oo}).seccomponents_perbasepx,1),1);
            for ii = 1:size(data.(pairname_short).(msvar{oo}).seccomponents_perbasepx,1) % For each base field

                secdata = data.(pairname_short).(msvar{oo}).secmaps_raw{ii};
                seclinmap = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLsm));
                seclinmap(1,secdata(:,1)) = secdata(:,4);
                secgridmap = lineartogrid(seclinmap',msvar{2-oo+1},gridSize{2-oo+1});
                dummygridsec = data.(pairname_short).(msvar{2-oo+1}).dummygrid;

                % If there are secondary fields
                if ~isempty(data.(pairname_short).(msvar{2-oo+1}).seccomponents_perbasepx)
                    % Output variables
                    sec_infieldrate = nan(size(data.(pairname_short).(msvar{2-oo+1}).seccomponents_perbasepx,1),1);
                    sec_outfieldrate = cell(size(data.(pairname_short).(msvar{2-oo+1}).seccomponents_perbasepx,1),1);

                    % Test if secondary pixels sampled from this base field are more likely to fall within any of the secondary fields than outside
                    for jj = 1:size(data.(pairname_short).(msvar{2-oo+1}).seccomponents_perbasepx,1) % For each secondary field

                        % Get mean firing rate within secondary field
                        fieldlinbin = data.(pairname_short).(msvar{2-oo+1}).fieldlinbin{jj}; % pixels that make up the sec field. Not all of these will be sampled from this base field
                        seclinbin_sampled = fieldlinbin(ismember(fieldlinbin,secdata(:,1)));
                        inds_infield = ismember(secdata(:,1),seclinbin_sampled); % find pixels of sec field that are sampled from this base field
                        if sum(inds_infield) == 0 
                            continue;
                        end
                        meanrate_infield = sum(secdata(inds_infield,3)) / sum(secdata(inds_infield,2));

                        % Get mean firing rates for 10000 pseudorandom same-size fields outside of secondary field
                        meanrate_outfield = nan(Args.NumShuffles,1);
                        session_seclinbin_out = secdata(~inds_infield,:);
                        if size(session_seclinbin_out,1) > sum(inds_infield) % Filter 1 of 2: Make sure there are enough outfield px to generate stats (This catches situations where there are just too few outfield px to begin with) 
                            % Generate a psuedopopulation of outfields same size as sec field 
                            ff = 1;
                            attempt = 0; % Filter 2 of 2: Make sure there are enough outfield px to generate stats (This catches situations where outfield px are enough in number but scattered so that can't form coherent shuffled field without overlapping with original sec field)
                            abandon = false;
                            while ff <= Args.NumShuffles && ~abandon
                                attempt = attempt + 1;
                                % Start from a random pixel that is outside of sec field
                                startpx = randsample(1:size(session_seclinbin_out,1),1);
                                startpx = session_seclinbin_out(startpx,1);
                                if ismember(startpx,fieldlinbin) % If random px overlaps with sec field, repeat
        %                                 reset = true;
                                    continue;
                                end
                                % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
                                switch msvar{2-oo+1}
                                    case 'place'
                                        gnum = 1;
                                    case 'view'
                                        if startpx == 1 || startpx == 2 % Make sure not cue or hint
        %                                        reset = true;
                                           continue;
                                        end
                                        [gnum,~,~] = findgrid(startpx,msvar{2-oo+1});
                                    case 'headdirection'
                                        gnum = 1;
                                end
                                % Get the grid coords of starting px
                                [startindx,startindy] = find(dummygridsec{gnum} == startpx);
                                tempmap = secgridmap{gnum}; % the actual sampled sec grid map
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
                                elseif size(intersect(sampledpx,fieldlinbin),1) > 0.25*size(fieldlinbin,1) % If sampled field overlaps with more than half of sec field
                                    if attempt == Args.NumShuffles && ff < 10
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
                                if isempty(setdiff(pxsub,fieldlinbin))
        %                                 reset = true;
                                    continue;
                                end
                                % Get sec bins actually sampled 
                                outfieldlinbin = pxsub;
                                inds_outfield = ismember(secdata(:,1),outfieldlinbin);
                                meanrate_outfield(ff,1) = sum(secdata(inds_outfield,3)) / sum(secdata(inds_outfield,2));
        %                                 if mod(ff,1000) == 0
        %                                     disp(ff);
        %                                 end
                                ff = ff+1;
                            end
        %                         meanrate_thr = prctile(meanrate_outfield(:),95);
                        end
                        % Store data
                        sec_infieldrate(jj,1) = meanrate_infield;
                        sec_outfieldrate{jj,1} = meanrate_outfield;
                    end
                    % Store data
                    sec_infieldrates{ii,1} = sec_infieldrate;
                    sec_outfieldrates{ii,1} = sec_outfieldrate;
                else 
                    if data.(pairname_short).mixsel
                        disp(['no sec ' msvar{2-oo+1} ' fields']);
                    end
                end
            end
            % Store data
            data.(pairname_short).(msvar{oo}).sec_infieldrates = sec_infieldrates;
            data.(pairname_short).(msvar{oo}).sec_outfieldrates = sec_outfieldrates;

        end

        %% Get pseudopopulation of smoothed base maps and compare to secondary smoothed maps
        secfieldnumbinset = cell(size(msvar,2),1);
        condfield_inbinsset = cell(size(msvar,2),1);
        condfield_nonoverlapbinset = cell(size(msvar,2),1);
        condfield_tertlinbin = cell(size(msvar,2),1);
        for oo = 1:size(msvar,2)
            disp(['Creating pseudopopulation of base maps for ' msvar{oo} ' field ...'])
            test = [msvar_short{1} msvar_short{2}];
            if strcmp(test,'hv') 
                disp('stop')
            end
                
            % Output variables
            pseudosecdataperfield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            secfieldnumbinperfield = nan(data.(pairname_short).(msvar{oo}).sigfields,1);
            condfield_inbinsperfield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            condfield_nonoverlapbinsperfield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            cond_tertlinbinsperfield = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            for ii = 1:size(data.(pairname_short).(msvar{oo}).seccomponents_perbasepx,1) % For each base field

                basefieldnumspk = sum(data.(pairname_short).(msvar{oo}).secmaps_raw{ii}(:,3));
                basefieldbins = data.(pairname_short).(msvar{oo}).fieldlinbin{ii};
        %                 secfieldnumbin = size(data.(pairname_short).(msvar{oo}).secmaps_raw{ii},1);
                secfieldnumbin = sum(~isnan(data.(pairname_short).(msvar{oo}).secmaps_sm{ii}));
                baselinmap = data.(pairname_short).(msvar{oo}).basemapLrw;
                switch msvar{oo}
                    case 'place'
                        baselinspk = pc.spk_raw;
                    case 'view'
                        baselinspk = sv.spk_raw;
                    case 'headdirection'
                        baselinspk = hd.spk_raw;
                end
        %                 basedrawpx = find(~isnan(baselinmap) & baselinspk>0); % Sampled and has a spike
                basedrawpx = find(~isnan(baselinmap)); % Sampled, spike or no spike
                basedrawpx = setdiff(basedrawpx,basefieldbins);
                basegridmap = data.(pairname_short).(msvar{oo}).basemapGrw;
                dummygridbase = data.(pairname_short).(msvar{oo}).dummygrid;
                dummygridsec = data.(pairname_short).(msvar{2-oo+1}).dummygrid;

                % Generate a pseudopopulation of 500 base fields drawn from coverage of conditioned maps of each of the sec fields 
        %                 condfields = data.(pairname_short).(msvar{oo}).seccond_fieldlinbin{ii}; % bins of each of the conditioned fields from secmap
        %                 condfield_inbins = [];
        %                 for jj = 1:size(condfields,1)
        %                     condfield_inbins = [condfield_inbins; condfields{jj}];
        %                 end
        %                 condfield_inbins = unique(condfield_inbins);
                allsecbins = data.(pairname_short).(msvar{oo}).secmaps_raw{ii}(:,1);
                condfield_inbins = find(data.(pairname_short).(msvar{oo}).secmaps_sm{ii} > nanmean(data.(pairname_short).(msvar{oo}).secmaps_sm{ii}))'; % simplistically, fields in sec/conditioned maps are defined as pixels exceeding mean. no contiguity requirement
                condfield_outbins = setdiff(allsecbins,condfield_inbins);

        %             % Get 1000 samples of field same size as conditioned fields
        %             ff = 1;
        %             shuff = Args.NumShuffles;
        %             attempt = 0;
        %             pseudosecdur = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
        %             pseudosecspk = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
        %             pseudosecmap_raw = nan(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
        %             pseudobasebins = cell(shuff,1);
        %             pseudocondbins_out = cell(shuff,1);
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
        % %                     for kk = 1:size(data.(pairname_short).(msvar{oo}).seccond_fieldlinbin{ii},1)
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
        %                     pseudocondbins_out{ff,:} = nonoverlap;
        %                     pseudotertbins{ff,:} = tertlinbin;
        %                     
        %                     ff = ff + 1;
        %                     
        %                 end



                    % Generate a psuedopopulation of 1000 base fields with same num of spikes as original base field 
                    ff = 1;
                    shuff = 1000;
%                     if ~strcmp(test,'hv')
%                         shuff = 10;
%                     else
%                         shuff = 1000; % Args.NumShuffles;
%                     end
                    attempt = 0; % Filter 2 of 2: Make sure there are enough outfield px to generate stats (This catches situations where outfield px are enough in number but scattered so that can't form coherent shuffled field without overlapping with original sec field)
                    abandon = false;
                    pseudosecmap_adsm = nan(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudosecSIC_sm = nan(shuff,1);
                    pseudosecdur = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudosecspk = zeros(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudosecmap_raw = nan(shuff,size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw,2));
                    pseudobasebins = cell(shuff,1);
                    pseudocondbins_out = cell(shuff,1);
                    pseudotertbins = cell(shuff,1);
                    tic;
                    if strcmp(msvar{oo},'headdirection')
                        disp('stop');
                    end
                    while ff <= shuff && ~abandon

                        if mod(ff,100) == 0 % ff == shuff/2 % mod(ff,100) == 0
                            disp(['ff = ' num2str(ff)]);
                        end
                        attempt = attempt + 1;
                        % Start from a random pixel that is outside of base field and has a spike
                        startpx = randsample(1:length(basedrawpx),1);
                        startpx = basedrawpx(1,startpx);
                        if ismember(startpx,basefieldbins) % If random px overlaps with sec field, repeat
        %                                 reset = true;
                            continue;
                        end
                        % Constrain the pseudorandom population to same grid number (for spatial view) e.g. pillar only
                        switch msvar{oo}
                            case 'place'
                                gnum = 1;
                            case 'view'
                                if startpx == 1 || startpx == 2 % Make sure not cue or hint
        %                                        reset = true;
                                   continue;
                                end
                                [gnum,~,~] = findgrid(startpx,msvar{oo});
                            case 'headdirection'
                                gnum = 1;
                        end
                        % Get the grid coords of starting px
                        [startindx,startindy] = find(dummygridbase{gnum} == startpx);
                        tempmap = basegridmap{gnum}; % the actual sampled sec grid map
                        % Expand radius around starting px until hit the requisite number of px 
                        % while sum(sum(~isnan(tempmap(startindx,startindy)))) < sum(inds_infield)
                        while length(startindx)*length(startindy) < size(dummygridbase{gnum},1)*size(dummygridbase{gnum},2)
                            if startindx(1) > 1 && startindx(end) < size(dummygridbase{gnum},1)
                                startindx = [startindx(1)-1 startindx startindx(end)+1];
                            elseif startindx(1) == 1 && startindx(end) < size(dummygridbase{gnum},1)
                                startindx = [startindx startindx(end)+1];
                            elseif startindx(1) > 1 && startindx(end) == size(dummygridbase{gnum},1)
                                startindx = [startindx(1)-1 startindx];
                            end

                            % If reach num of sec bins in original sec map
                            lin_inds = dummygridbase{gnum}(startindx,startindy);
                            lin_inds = lin_inds(~isnan(tempmap(startindx,startindy))); % base
                            ind_pv = ismember(stcfill(:,strcmp(stcvars,msvar{oo})),lin_inds);
                            growingpx = unique(stcfill(ind_pv,strcmp(stcvars,msvar{2-oo+1}))); % sec
        %                         switch msvar{oo}
        %                             case 'place'
        %                                 ind_pv = ismember(stcfill(:,2),lin_inds);
        %                                 growingpx = unique(stcfill(ind_pv,3));
        %                             case 'view'
        %                                 ind_pv = ismember(stcfill(:,3),lin_inds);
        %                                 growingpx = unique(stcfill(ind_pv,2));
        %                         end
                            growingspikes = nansum(stcfill(ind_pv,6));

                            % If num of pseudo sec px matches num of original sec px
                            if size(growingpx,1) > 0.9*secfieldnumbin % 0.9*secfieldnumbin %|| growingspikes > 0.5*basefieldnumspk
                                break;
                            end
                            if startindy(1) > 1 && startindy(end) < size(dummygridbase{gnum},2)
                                startindy = [startindy(1)-1 startindy startindy(end)+1];
                            elseif startindy(1) == 1 && startindy(end) < size(dummygridbase{gnum},2)
                                startindy = [startindy startindy(end)+1];
                            elseif startindy(1) > 1 && startindy(end) == size(dummygridbase{gnum},2)
                                startindy = [startindy(1)-1 startindy];
                            end
                            % If exceed map bounds
                            if startindx(1) == 1 && startindx(end) == size(dummygridbase{gnum},1) && startindy(1) == 1 && startindy(end) == size(dummygridbase{gnum},2)
                                break;
                            end
                        end
                        % pixels of pseudo base field
                        lin_inds = dummygridbase{gnum}(startindx,startindy);
                        lin_inds = lin_inds(~isnan(tempmap(startindx,startindy))); % base
                        ind_pv = ismember(stcfill(:,strcmp(stcvars,msvar{oo})),lin_inds);
                        growingpx = unique(stcfill(ind_pv,strcmp(stcvars,msvar{2-oo+1}))); % sec
        %                     switch msvar{oo}
        %                         case 'place'
        %                             ind_pv = ismember(stcfill(:,2),lin_inds);
        %                             growingpx = unique(stcfill(ind_pv,3));
        %                         case 'view'
        %                             ind_pv = ismember(stcfill(:,3),lin_inds);
        %                             growingpx = unique(stcfill(ind_pv,2));
        %                     end
                        growingspikes = nansum(stcfill(ind_pv,6));

                        nonoverlap = setdiff(growingpx,condfield_inbins);
                        if size(growingpx,1) < 0.9*secfieldnumbin % || growingspikes < 0.5*basefieldnumspk % commented out because for very directional cells, this will never pass
        %                         disp(growingspikes);
        %                         if attempt == 100
        %                             disp('stop')
        %                         end
                            continue;
                        elseif size(intersect(lin_inds,basefieldbins),1) > 0.5*size(basefieldbins,1) % If sampled field overlaps with more than half of sec field
                            if attempt == shuff && ff < 10
                                abandon = true;
                            end
                            continue;
                        end
                        tertlinbin = intersect(condfield_outbins,nonoverlap);

                        % Find pseudo sec maps
                        usecpxs = [];
                        for pp = 1:size(lin_inds,1)
                            secpx = [];
                            ind_pv = stcfill(:,strcmp(stcvars,msvar{oo})) == lin_inds(pp);
                            secpx(:,1) = stcfill(ind_pv,strcmp(stcvars,msvar{2-oo+1})); % sec px
        %                         switch msvar{oo}
        %                             case 'place'
        %                                 ind_pv = stcfill(:,2) == lin_inds(pp);
        %                                 secpx(:,1) = stcfill(ind_pv,3); % sec px
        %                             case 'view'
        %                                 ind_pv = stcfill(:,3) == lin_inds(pp);
        %                                 secpx(:,1) = stcfill(ind_pv,2); % sec px
        %                         end
                            secpx(:,2) = stcfill(ind_pv,5); % dur
                            secpx(:,3) = stcfill(ind_pv,6); % spk

                            % Get firing rates 
                            usecpx = unique(secpx(:,1));
                            usecpx(isnan(usecpx)) = [];
                            rate_components_px = nan(length(usecpx),4); % Collect dur and spikes for secondary pixels for calculating firing rates
                            rate_components_px(:,1) = usecpx;
                            if any(isnan(usecpx))
                                disp(nan);
                            end
                            ind_timechange = find(secpx(:,2) > 0);
                            for cc = 1:size(ind_timechange,1) % For each instance of being in this base pixel
                                linbintemp = nan(size(usecpx,1),2); % Temp dur and spikes for this secondary pixel(s)
                                if cc == size(ind_timechange,1)
                                    newind = ind_timechange(cc):size(secpx,1);
                                else
                                    newind = ind_timechange(cc):ind_timechange(cc+1)-1; % index into secondary pixel(s) for this instance
                                end
                                newview = secpx(newind,1); 
                                newset = ismember(usecpx,newview);
                                linbintemp(newset,1) = secpx(newind(1),2); % duration for this instance listed with first secondary pixel
                                linbintemp(newset,2) = secpx(newind(end),3); % spikes for this instance listed with last secondary pixel
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

                        secdata = usecpxs; 
                        secdata(:,2:3) = NaN;
                        % Consolidate occupancy and spikes across base field
                        for pp = 1:size(lin_inds,1)
                            tempbins = nan(size(secdata,1),2);
                            tempbins(ismember(secdata(:,1),rate_components_perbasepx{pp}(:,1)),1) = rate_components_perbasepx{pp}(:,2); % Duration 
                            tempbins(ismember(secdata(:,1),rate_components_perbasepx{pp}(:,1)),2) = rate_components_perbasepx{pp}(:,3); % Spikes 
                            secdata(:,2) = nansum( [secdata(:,2) tempbins(:,1)] ,2); % Sum duration across session
                            secdata(:,3) = nansum( [secdata(:,3) tempbins(:,2)] ,2); % Sum spikes across session
                        end
                        secdata(:,4) = secdata(:,3)./secdata(:,2);

                        pseudodur = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
                        pseudodur(secdata(:,1)) = secdata(:,2);
                        pseudospk = zeros(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
                        pseudospk(secdata(:,1)) = secdata(:,3);
                        pseudomap = nan(size(data.(pairname_short).(msvar{2-oo+1}).basemapLrw));
                        pseudomap(secdata(:,1)) = secdata(:,4);

                        pseudosecdur(ff,:) = pseudodur;
                        pseudosecspk(ff,:) = pseudospk;
                        pseudosecmap_raw(ff,:) = pseudomap;
                        pseudobasebins{ff} = lin_inds;
                        pseudocondbins_out{ff,:} = nonoverlap;
                        pseudotertbins{ff,:} = tertlinbin;

                        ff = ff + 1;
                    end
                toc;
                pseudosecdataperfield{ii} = {pseudobasebins pseudosecdur pseudosecspk pseudosecmap_raw};
                secfieldnumbinperfield(ii) = secfieldnumbin;
                condfield_inbinsperfield{ii} = condfield_inbins;
                condfield_nonoverlapbinsperfield{ii} = pseudocondbins_out;
                cond_tertlinbinsperfield{ii} = pseudotertbins;
            end
            data.(pairname_short).(msvar{oo}).pseudosecdataperfield = pseudosecdataperfield;
            secfieldnumbinset{oo} = secfieldnumbinperfield;
            condfield_inbinsset{oo} = condfield_inbinsperfield;
            condfield_nonoverlapbinset{oo} = condfield_nonoverlapbinsperfield;
            condfield_tertlinbin{oo} = cond_tertlinbinsperfield;
        end

        % Smooth
        disp('Smoothing pseudopopulation of secondary maps ...');
        shuff = 1000;
        for oo = 1:size(msvar,2)
            pseudosecmaps_sm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            pseudosecSIC_sm = cell(data.(pairname_short).(msvar{oo}).sigfields,1);
            tertrate_orig = nan(shuff,data.(pairname_short).(msvar{oo}).sigfields);
            tertrate_pseudo = nan(shuff,data.(pairname_short).(msvar{oo}).sigfields);
            for ii = 1:size(data.(pairname_short).(msvar{oo}).seccomponents_perbasepx,1) % For each base field

                temp = data.(pairname_short).(msvar{oo}).pseudosecdataperfield{ii};
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
                        tic;
                        %%% Smoothing loop
                        disp('loop smooth');
                        secmap_sm = nan(size(seclindur'));
                        secdur_sm = nan(size(seclindur'));
                        for kk = 1:size(seclindur,1)
                            if mod(kk,100) == 0
                                disp(['smooth ' num2str(kk)]);
                            end
                            seclindursingle = seclindur(kk,:);
                            seclinspksingle = seclinspk(kk,:);
                            % Assign linear bin to grid bin - left to right, bottom to top
                            durG = lineartogrid(seclindursingle','view',binDepths);
                            spkG = lineartogrid(seclinspksingle','view',binDepths);
        %                             rateG = lineartogrid(seclinmap','view',binDepths);

                            % Pad sv map with 5 extra rows
                            n = 5;
                            padpillar = false;
                            [emptyfloorref_pad,~] = padsvmap(n,durG,gazeSections,padpillar);
                            padpillar = true;
                            [durGpad,retrievemap] = padsvmap(n,durG,gazeSections,padpillar);
                            [spkGpad,~] = padsvmap(n,spkG,gazeSections,padpillar);
        %                             [rateGpad,~] = padsvmap(n,rateG,gazeSections,padpillar);

                            % Adaptive smooth
                            maps_adsmGpad = cell(size(durGpad));
                            dur_adsmGpad = cell(size(durGpad));
                            for jj = 1:size(binDepths,1)
                                if jj == 1 || jj == 2
                                    maps_adsmGpad{jj} = spkGpad{jj}./durGpad{jj};
                                    dur_adsmGpad{jj} = durGpad{jj};
                                else
        %                             [maps_adsmGpad{jj},spk_adsmG,dur_adsmG] = adsmooth(durGpad{jj},spkGpad{jj},sv.Args.Alpha);
                                    [maps_adsmGpad{jj},spk_adsmG,dur_adsmGpad{jj}] = adsmooth(durGpad{jj},spkGpad{jj},1e2);
                                end
                            end

                            % Unpad smoothed map
                            maps_smG = unpadsvmap(maps_adsmGpad,retrievemap,durG);
                            dur_smG = unpadsvmap(dur_adsmGpad,retrievemap,durG);
                            % Convert grid map back to linear sv map
                            secmap_sm(:,kk) = gridtolinear(maps_smG,'view',binDepths);
                            secdur_sm(:,kk) = gridtolinear(dur_smG,'view',binDepths);
                            secdur_sm(isnan(secdur_sm(:,kk)),kk) = 0;

                        end
                        toc;
                    case 'headdirection'
                        n = 5;
                        % Smooth
        %                     rateG = seclinmap'; % raw
                        secmap_sm = smoothdir(seclinmap',n,cr.headdirectionbins);
        %                     maps_smG = secmap_sm; % no difference between linear and grid map for HD
                        secdur_sm = smoothdir(seclindur',n,cr.headdirectionbins);
                        toc;
                end

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

                % Get mean firing rate of bins outside of conditioned
                % fields of original map
                rate_orig = nan(shuff,1);
                rate_pseudo = nan(shuff,1);
                for kk = 1:shuff
                    px = condfield_tertlinbin{oo}{ii}{kk};
                    rate_orig(kk) = nanmean(data.(pairname_short).(msvar{oo}).secmaps_sm{ii}(px));
                    rate_pseudo(kk) = nanmean(secmap_sm(px,kk));
                end

                tertrate_orig(:,ii) = rate_orig;
                tertrate_pseudo(:,ii) = rate_pseudo;

        %                 
        %                 
        %                 temp = nan(size(data.(pairname_short).(msvar{oo}).seccond_fieldlinbin{ii},1),1);
        %                 for kk = 1:size(data.(pairname_short).(msvar{oo}).seccond_fieldlinbin{ii},1)
        %                     outpx = setdiff(data.(pairname_short).(msvar{oo}).secmaps_raw{ii}(:,1),data.(pairname_short).(msvar{oo}).seccond_fieldlinbin{ii}{kk});
        %                     temp(kk) = nanmean(data.(pairname_short).(msvar{oo}).secmaps_sm{ii}(outpx));
        %                 end
        %                 firingrate_out_orig = nanmean(temp);
        %                 
        %                 % Get mean firing rate of bins outside of conditioned
        %                 % fields of pseudo maps
        %                 nonoverlapset = condfield_nonoverlapbinset{oo}{ii};
        %                 for kk = 1:size(nonoverlapset,1)
        % %                     % find sec field that is in both pseudo map and original conditioned map
        % %                     px = randsample(nonoverlapset{kk},size(data.(pairname_short).(msvar{oo}).fieldlinbin{ii},1));
        %                     firingrate_out_pseudo(kk) = nanmean(secmapadsm(nonoverlapset{kk},kk));
        %                 end

            end

            % Store data
            data.(pairname_short).(msvar{oo}).pseudosecmaps_adsm = pseudosecmaps_sm;
            data.(pairname_short).(msvar{oo}).pseudosecSIC_adsm = pseudosecSIC_sm;
            data.(pairname_short).(msvar{oo}).tertrate_orig = tertrate_orig;
            data.(pairname_short).(msvar{oo}).tertrate_pseudo = tertrate_pseudo;
        end

        %         %% Get fields in pseudomaps
        %         
        %         for oo = 1:size(msvar,2)
        %             
        %             
        %             
        %             
        %         end
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
    case 'view'
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
    case 'headdirection'
        mapLdummy = 1:60;
        gridnum = 1;
        temp = flipud(reshape(mapLdummy, 60, 1)');
end
[y,x] = find(temp == px);
y = size(temp,1)-y+1;
