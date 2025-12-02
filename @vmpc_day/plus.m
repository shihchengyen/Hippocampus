function r = plus(p, q, varargin)
%@vmpc_day/plus Plus function for vmpc_day objects
%   R = PLUS(P, Q) combines two vmpc_day objects by vertically
%   concatenating their data fields for use with InspectGUI.

% Get classname
classname = mfilename('class');

% Boilerplate code to check if inputs are the correct class
if (~isa(p, classname))
    if (~isa(q, classname))
        % Both inputs are not vmpc_day objects, create empty
        r = feval(classname);
    else
        % Second input is a vmpc_day object, return that
        r = q;
    end
else
    if (~isa(q, classname))
        % p is a vmpc_day object but q is not, so just return p
        r = p;
    elseif (isempty(p))
        r = q;
    elseif (isempty(q))
        r = p;
    else
        % Both p and q are vmpc_day objects, so add them
        r = p;
   
        
        % --- 1. Combine Parent nptdata Properties ---
        r.nptdata = plus(p.nptdata, q.nptdata);

        % --- 2. Combine vmpc_day Data Properties ---
        if (p.data.numSets == 0)
             r = q;
             return;
        end
        r.data.numSets = p.data.numSets + q.data.numSets;
        
        % Raw data
        r.data.maps_raw = [p.data.maps_raw; q.data.maps_raw];
        r.data.dur_raw = [p.data.dur_raw; q.data.dur_raw];
        r.data.spk_raw = [p.data.spk_raw; q.data.spk_raw];
        r.data.maps_raw1 = [p.data.maps_raw1; q.data.maps_raw1];
        r.data.dur_raw1 = [p.data.dur_raw1; q.data.dur_raw1];
        r.data.spk_raw1 = [p.data.spk_raw1; q.data.spk_raw1];
        r.data.maps_raw2 = [p.data.maps_raw2; q.data.maps_raw2];
        r.data.dur_raw2 = [p.data.dur_raw2; q.data.dur_raw2];
        r.data.spk_raw2 = [p.data.spk_raw2; q.data.spk_raw2];
        
        % Smoothed maps (real)
        r.data.maps_adsm = [p.data.maps_adsm; q.data.maps_adsm];
        r.data.maps_bcsm = [p.data.maps_bcsm; q.data.maps_bcsm];
        r.data.maps_dksm = [p.data.maps_dksm; q.data.maps_dksm];
        r.data.maps_sm = [p.data.maps_sm; q.data.maps_sm];
        
        % Smoothed maps (shuffled)
        r.data.maps_adsmsh = [p.data.maps_adsmsh; q.data.maps_adsmsh];
        r.data.maps_bcsmsh = [p.data.maps_bcsmsh; q.data.maps_bcsmsh];
        r.data.maps_dksmsh = [p.data.maps_dksmsh; q.data.maps_dksmsh];
        r.data.maps_smsh = [p.data.maps_smsh; q.data.maps_smsh];

        % Smoothed duration & radii
        r.data.dur_adsm = [p.data.dur_adsm; q.data.dur_adsm];
        r.data.radii = [p.data.radii; q.data.radii];
        r.data.radiish = [p.data.radiish; q.data.radiish];

        % Metrics (real)
        r.data.crit_sm = [p.data.crit_sm; q.data.crit_sm];
        r.data.SIC_adsm = [p.data.SIC_adsm; q.data.SIC_adsm];
        r.data.SIC_bcsm = [p.data.SIC_bcsm; q.data.SIC_bcsm];
        r.data.SIC_dksm = [p.data.SIC_dksm; q.data.SIC_dksm];
        
        % Metrics (shuffled)
        r.data.critsh_sm = [p.data.critsh_sm; q.data.critsh_sm];
        r.data.SICsh_adsm = [p.data.SICsh_adsm; q.data.SICsh_adsm];
        r.data.SICsh_bcsm = [p.data.SICsh_bcsm; q.data.SICsh_bcsm];
        % --- THIS IS THE CORRECTED LINE ---
        r.data.SICsh_dksm = [p.data.SICsh_dksm; q.data.SICsh_dksm];
        
        % Other metrics
        r.data.sparsity = [p.data.sparsity; q.data.sparsity];
        r.data.sig2noise = [p.data.sig2noise; q.data.sig2noise];
        r.data.coherence = [p.data.coherence; q.data.coherence];
        r.data.coherence_sm = [p.data.coherence_sm; q.data.coherence_sm];
        r.data.intracorr = [p.data.intracorr; q.data.intracorr];
        
        % Quality flags
        r.data.discard = [p.data.discard; q.data.discard];
        r.data.rateok = [p.data.rateok; q.data.rateok];
        r.data.significant = [p.data.significant; q.data.significant];

        % Update counts
        r.data.filtspknum = [p.data.filtspknum; q.data.filtspknum];
        
        % Append cell arrays
        r.data.sessionNames = {p.data.sessionNames{:} q.data.sessionNames{:}};
        
        if (p.data.numSets == 1)
            % This is the FIRST addition (e.g., Cell 1 + Cell 2).
            % We take the first path from 'p' and the first path from 'q'.
            r.data.cellPaths = {p.data.cellPaths{1}, q.data.cellPaths{1}};
        else
            % 'p' is already a combined object (e.g., [Cell 1, Cell 2] + Cell 3).
            % We take all paths from 'p' and append the first path from 'q'.
            r.data.cellPaths = {p.data.cellPaths{:}, q.data.cellPaths{1}};
        end
        
    end
end