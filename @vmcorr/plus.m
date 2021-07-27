function r = plus(p,q,varargin)
%@dirfiles/plus Overloaded plus function for dirfiles objects.
%   R = plus(P,Q) combines dirfiles objects P and Q and returns the
%   dirfiles object R.

% For vmcorr

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
        r.data.origin = [p.data.origin; q.data.origin];

		% object specific fields %% Missing some
        r.data.maps_raw_corrp = [p.data.maps_raw_corrp; q.data.maps_raw_corrp];
        r.data.maps_raw_corrpset = [p.data.maps_raw_corrpset; {q.data.maps_raw_corrpset}];
        r.data.maps_adsm_corrp = [p.data.maps_adsm_corrp; q.data.maps_adsm_corrp];
        r.data.maps_adsm_corrpset = [p.data.maps_adsm_corrpset; {q.data.maps_adsm_corrpset}];
        r.data.SIC_corrp = [p.data.SIC_corrp; q.data.SIC_corrp];
        r.data.SIC_corrpset = [p.data.SIC_corrpset; {q.data.SIC_corrpset}];
        r.data.ISE_corrp = [p.data.ISE_corrp; q.data.ISE_corrp];
        r.data.ISE_corrpset = [p.data.ISE_corrpset; {q.data.ISE_corrpset}];

        
        r.data.maps_raw_corrsv = [p.data.maps_raw_corrsv; q.data.maps_raw_corrsv];
        r.data.maps_raw_corrsvset = [p.data.maps_raw_corrsvset; {q.data.maps_raw_corrsvset}];
        r.data.maps_adsm_corrsv = [p.data.maps_adsm_corrsv; q.data.maps_adsm_corrsv];
        r.data.maps_adsm_corrsvset = [p.data.maps_adsm_corrsvset; {q.data.maps_adsm_corrsvset}];
        r.data.SIC_corrsv = [p.data.SIC_corrsv; q.data.SIC_corrsv];
        r.data.SIC_corrsvset = [p.data.SIC_corrsvset; {q.data.SIC_corrsvset}];
        r.data.ISE_corrsv = [p.data.ISE_corrsv; q.data.ISE_corrsv];
        r.data.ISE_corrsvset = [p.data.ISE_corrsvset; {q.data.ISE_corrsvset}];
        r.data.covmat = [p.data.covmat; {q.data.covmat}];
        r.data.covmat_norm = [p.data.covmat_norm; {q.data.covmat_norm}];
        r.data.norml1 = [p.data.norml1; q.data.norml1];
        r.data.norml2 = [p.data.norml2; q.data.norml2];
        
        r.data.convergewithsv = [p.data.convergewithsv; q.data.convergewithsv];
        r.data.convergewithp = [p.data.convergewithp; q.data.convergewithp];
        r.data.smoothpick = [p.data.smoothpick; q.data.smoothpick];
        r.data.llhpick = [p.data.llhpick; q.data.llhpick];
        r.data.NumIterLlh = [p.data.NumIterLlh; q.data.NumIterLlh];
        r.data.llhpicklabel = [p.data.llhpicklabel; {q.data.llhpicklabel}];
        r.data.llh = [p.data.llh; {q.data.llh}];
			
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
