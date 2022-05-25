function r = plus(p,q,varargin)
%@vmhd/plus Overloaded plus function for vmhd objects.
%   R = plus(P,Q) combines dirfiles objects P and Q and returns the
%   dirfiles object R.

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
        
        % Full session
        r.data.maps_raw = [p.data.maps_raw; q.data.maps_raw];
        r.data.dur_raw = [p.data.dur_raw; q.data.dur_raw];
        r.data.spk_raw = [p.data.spk_raw; q.data.spk_raw];
        r.data.maps_bcsm = [p.data.maps_bcsm; q.data.maps_bcsm];
        r.data.maps_bcsmsh = [p.data.maps_bcsmsh; q.data.maps_bcsmsh];
        r.data.RV_bcsm = [p.data.RV_bcsm; q.data.RV_bcsm];
        r.data.RVsh_bcsm = [p.data.RVsh_bcsm; q.data.RVsh_bcsm];
        
        % First half
        r.data.maps_raw1 = [p.data.maps_raw1; q.data.maps_raw1];
        r.data.dur_raw1 = [p.data.dur_raw1; q.data.dur_raw1];
        r.data.spk_raw1 = [p.data.spk_raw1; q.data.spk_raw1];
        r.data.maps_bcsm1 = [p.data.maps_bcsm1; q.data.maps_bcsm1];
        r.data.RV_bcsm1 = [p.data.RV_bcsm1;q.data.RV_bcsm1];
        
        % Second half
        r.data.maps_raw2 = [p.data.maps_raw2; q.data.maps_raw2];
        r.data.dur_raw2 = [p.data.dur_raw2; q.data.dur_raw2];
        r.data.spk_raw2 = [p.data.spk_raw2; q.data.spk_raw2];
        r.data.maps_bcsm2 = [p.data.maps_bcsm2; q.data.maps_bcsm2];
        r.data.RV_bcsm2 = [p.data.RV_bcsm2;q.data.RV_bcsm2];

		% object specific fields
% 		r.data.dlist = [p.data.dlist; q.data.dlist];
% 		r.data.setIndex = [p.data.setIndex; (p.data.setIndex(end) ...
% 			+ q.data.setIndex(2:end))];

        r.data.origin = [p.data.origin; q.data.origin];
			
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
