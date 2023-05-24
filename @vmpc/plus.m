function r = plus(p,q,varargin)
%@dirfiles/plus Overloaded plus function for dirfiles objects.
%   R = plus(P,Q) combines dirfiles objects P and Q and returns the
%   dirfiles object R.

% for vmpc

% get name of class
classname = mfilename('class');

% check if first input is the right kind of object
if(~isa(p,classname))
	% check if second input is the right kind of object
	if(~isa(q,classname))
		% both inputs are not the right kind of object so create empty
		% object and return it
		r = feval(classname);
	else
		% second input is the right kind of object so return that
		r = q;
	end
else
	if(~isa(q,classname))
		% p is the right kind of object but q is not so just return p
		r = p;
    elseif(isempty(p))
        % p is right object but is empty so return q, which should be
        % right object
        r = q;
    elseif(isempty(q))
        % p are q are both right objects but q is empty while p is not
        % so return p
        r = p;
	else
		% both p and q are the right kind of objects so add them 
		% together
		% assign p to r so that we can be sure we are returning the right
		% object
		r = p;
		% useful fields for most objects
		r.data.numSets = p.data.numSets + q.data.numSets;
        r.data.sessionTimeC = {p.data.sessionTimeC; q.data.sessionTimeC};
        r.data.stcfilt = {p.data.stcfilt; q.data.stcfilt};
        r.data.filtspknum = [p.data.filtspknum; q.data.filtspknum];
        r.data.discard = [p.data.discard; q.data.discard];

		% object specific fields
% 		r.data.detailed_fr = [p.data.detailed_fr; q.data.detailed_fr];
        r.data.maps_raw = [p.data.maps_raw; q.data.maps_raw];
        r.data.maps_raw1 = [p.data.maps_raw1; q.data.maps_raw1];
        r.data.maps_raw2 = [p.data.maps_raw2; q.data.maps_raw2];
%         r.data.maps_adsm = [p.data.maps_adsm; q.data.maps_adsm];
%         r.data.maps_adsm1 = [p.data.maps_adsm1; q.data.maps_adsm1];
%         r.data.maps_adsm2 = [p.data.maps_adsm2; q.data.maps_adsm2];
%         r.data.maps_adsmsh = [p.data.maps_adsmsh; q.data.maps_adsmsh];
%         r.data.maps_bcsm = [p.data.maps_bcsm; q.data.maps_bcsm];
%         r.data.maps_bcsm1 = [p.data.maps_bcsm1; q.data.maps_bcsm1];
%         r.data.maps_bcsm2 = [p.data.maps_bcsm2; q.data.maps_bcsm2];
%         r.data.maps_bcsmsh = [p.data.maps_bcsmsh; q.data.maps_bcsmsh];
%         r.data.maps_dksm = [p.data.maps_dksm; q.data.maps_dksm];
%         r.data.maps_dksm1 = [p.data.maps_dksm1; q.data.maps_dksm1];
%         r.data.maps_dksm2 = [p.data.maps_dksm2; q.data.maps_dksm2];
%         r.data.maps_dksmsh = [p.data.maps_dksmsh; q.data.maps_dksmsh];
        r.data.maps_sm = [p.data.maps_sm; q.data.maps_sm];
        r.data.maps_sm1 = [p.data.maps_sm1; q.data.maps_sm1];
        r.data.maps_sm2 = [p.data.maps_sm2; q.data.maps_sm2];

        r.data.dur_raw = [p.data.dur_raw; q.data.dur_raw];
        r.data.dur_raw1 = [p.data.dur_raw1; q.data.dur_raw1];
        r.data.dur_raw2 = [p.data.dur_raw2; q.data.dur_raw2];
        r.data.dur_adsm = [p.data.dur_adsm; q.data.dur_adsm];
        r.data.spk_raw = [p.data.spk_raw; q.data.spk_raw];
        r.data.spk_raw1 = [p.data.spk_raw1; q.data.spk_raw1];
        r.data.spk_raw2 = [p.data.spk_raw2; q.data.spk_raw2];
        
%         r.data.dur_adsm1 = [p.data.dur_adsm1; q.data.dur_adsm1];
%         r.data.dur_adsm2 = [p.data.dur_adsm2; q.data.dur_adsm2];
%         r.data.dur_adsm1 = [p.data.dur_adsm1; q.data.dur_adsm1];
%         r.data.dur_adsm2 = [p.data.dur_adsm2; q.data.dur_adsm2];
%         r.data.dur_adsmsh = [p.data.dur_adsmsh; q.data.dur_adsmsh];
%         r.data.SIC_adsm = [p.data.SIC_adsm; q.data.SIC_adsm];
%         r.data.SIC_adsm1 = [p.data.SIC_adsm1; q.data.SIC_adsm1];
%         r.data.SIC_adsm2 = [p.data.SIC_adsm2; q.data.SIC_adsm2];
%         r.data.SICsh_adsm = [p.data.SICsh_adsm; q.data.SICsh_adsm]; %
%         r.data.SIC_bcsm = [p.data.SIC_bcsm; q.data.SIC_bcsm];
%         r.data.SIC_bcsm1 = [p.data.SIC_bcsm1; q.data.SIC_bcsm1];
%         r.data.SIC_bcsm2 = [p.data.SIC_bcsm2; q.data.SIC_bcsm2];
%         r.data.SICsh_bcsm = [p.data.SICsh_bcsm; q.data.SICsh_bcsm];
%         r.data.SIC_dksm = [p.data.SIC_dksm; q.data.SIC_dksm];
%         r.data.SIC_dksm1 = [p.data.SIC_dksm1; q.data.SIC_dksm1];
%         r.data.SIC_dksm2 = [p.data.SIC_dksm2; q.data.SIC_dksm2];
%         r.data.SICsh_dksm = [p.data.SICsh_dksm; q.data.SICsh_dksm];
%         r.data.ISE_adsm = [p.data.ISE_adsm; q.data.ISE_adsm];
%         r.data.ISE_adsm1 = [p.data.ISE_adsm1; q.data.ISE_adsm1];
%         r.data.ISE_adsm2 = [p.data.ISE_adsm2; q.data.ISE_adsm2];
%         r.data.ISEsh_adsm = [p.data.ISEsh_adsm; q.data.ISEsh_adsm];
        r.data.radii = [p.data.radii; q.data.radii];
%         r.data.radii1 = [p.data.radii1; q.data.radii1];
%         r.data.radii2 = [p.data.radii2; q.data.radii2];
%         r.data.radiish = [p.data.radiish; q.data.radiish];
        r.data.crit_sm = [p.data.crit_sm; q.data.crit_sm];
        r.data.critsh_sm = [p.data.critsh_sm; q.data.critsh_sm];
        r.data.critthrcell = [p.data.critthrcell; q.data.critthrcell];
        r.data.crit_sm1 = [p.data.crit_sm1; q.data.crit_sm1];
        r.data.crit_sm2 = [p.data.crit_sm2; q.data.crit_sm2];

        r.data.origin = [p.data.origin; q.data.origin];
             
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
