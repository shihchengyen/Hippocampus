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
		r.data.meanFRs = [p.data.meanFRs q.data.meanFRs];
		r.data.semFRs = [p.data.semFRs q.data.semFRs];
		r.data.SIC = [p.data.SIC; q.data.SIC];
        r.data.SICsh = [p.data.SICsh q.data.SICsh];
        r.data.half1stmeanFRs = [p.data.half1stmeanFRs q.data.half1stmeanFRs];
        r.data.half1stsemFRs = [p.data.half1stsemFRs q.data.half1stsemFRs];
        r.data.half1stSIC = [p.data.half1stSIC; q.data.half1stSIC];
        r.data.half2ndmeanFRs = [p.data.half2ndmeanFRs q.data.half2ndmeanFRs];
        r.data.half2ndsemFRs = [p.data.half2ndsemFRs q.data.half2ndsemFRs];
        r.data.half2ndSIC = [p.data.half2ndSIC; q.data.half2ndSIC];
			
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
