function r = plus(p,q,varargin)
%@dirfiles/plus Overloaded plus function for dirfiles objects.
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
		
		% object specific fields
		r.data.maps_raw = [p.data.maps_raw q.data.maps_raw];
        r.data.maps_boxsmooth = [p.data.maps_boxsmooth q.data.maps_boxsmooth];
        r.data.maps_adsmooth = [r.data.maps_adsmooth q.data.maps_adsmooth];
		r.data.mapsSEM_raw = [p.data.mapsSEM_raw q.data.mapsSEM_raw];
		r.data.SIC = [p.data.SIC; q.data.SIC];
        r.data.SICsh = [p.data.SICsh q.data.SICsh];
        r.data.maps_raw1sthalf = [p.data.maps_raw1sthalf q.data.maps_raw1sthalf];
        r.data.maps_boxsmooth1sthalf = [p.data.maps_boxsmooth1sthalf q.data.maps_boxsmooth1sthalf];
        r.data.maps_adsmooth1sthalf = [p.data.maps_adsmooth1sthalf q.data.maps_adsmooth1sthalf];
        r.data.mapsSEM_raw1sthalf = [p.data.mapsSEM_raw1sthalf q.data.mapsSEM_raw1sthalf];
        r.data.SIC1sthalf = [p.data.SIC1sthalf; q.data.SIC1sthalf];
        r.data.maps_raw2ndhalf = [p.data.maps_raw2ndhalf q.data.maps_raw2ndhalf];
        r.data.maps_boxsmooth2ndhalf = [p.data.maps_boxsmooth2ndhalf q.data.maps_boxsmooth2ndhalf];
        r.data.maps_adsmooth2ndhalf = [p.data.maps_adsmooth2ndhalf q.data.maps_adsmooth2ndhalf];
        r.data.mapsSEM_raw2ndhalf = [p.data.mapsSEM_raw2ndhalf q.data.mapsSEM_raw2ndhalf];
        r.data.SIC2ndhalf = [p.data.SIC2ndhalf; q.data.SIC2ndhalf];
			
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
