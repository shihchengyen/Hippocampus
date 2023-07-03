function r = plus(p,q,varargin)
%@dirfiles/plus Overloaded plus function for dirfiles objects.
%   R = plus(P,Q) combines dirfiles objects P and Q and returns the
%   dirfiles object R.

% For vmms

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
		r.data.dlist = [p.data.dlist; q.data.dlist];
		r.data.setIndex = [p.data.setIndex; (p.data.setIndex(end) ...
			+ q.data.setIndex(2:end))];
        
        % object specific fields
        % r.data.filtspkcount = [p.data.filtspkcount; q.data.filtspkcount];
%         r.data.mixsel = [p.data.mixsel; q.data.mixsel];
        % r.data.placesel = [p.data.placesel; q.data.placesel];
        % r.data.viewsel = [p.data.viewsel; q.data.viewsel];
        % r.data.headdirectionsel = [p.data.headdirectionsel; q.data.headdirectionsel];
        r.data.discard = [p.data.discard; q.data.discard];
        r.data.pv = [p.data.pv; q.data.pv];
        r.data.ph = [p.data.ph; q.data.ph];
        r.data.hv = [p.data.hv; q.data.hv];
%         r.data.place = [p.data.place; q.data.place];
%         r.data.view = [p.data.view; q.data.view];
%         r.data.headdirection = [p.data.headdirection; q.data.headdirection];
        r.data.Args = p.data.Args;
        r.data.Args.UseCorr = [p.data.Args.UseCorr; q.data.Args.UseCorr]; % overwrite
        r.data.Args.FieldThr = [p.data.Args.FieldThr; q.data.Args.FieldThr]; % overwrite
        r.data.Args.FieldSplitThr = [p.data.Args.FieldSplitThr; q.data.Args.FieldSplitThr]; % overwrite
        r.data.Args.NumShuffles = [p.data.Args.NumShuffles; q.data.Args.NumShuffles]; % overwrite
			
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
