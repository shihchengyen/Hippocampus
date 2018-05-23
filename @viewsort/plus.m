function r = plus(p,q,varargin)
%@viewsort/plus Overloaded plus function for viewsort objects.
%   R = plus(P,Q) combines viewsort objects P and Q and returns the
%   viewsort object R.

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
		r.data.spikeForms = [p.data.spikeForms; q.data.spikeForms];
        r.data.Noise = [p.data.Noise;q.data.Noise];
        r.data.numChannels = p.data.numChannels + q.data.numChannels;
		% object specific fields
		r.data.ChannelIndex = [p.data.ChannelIndex; (p.data.ChannelIndex(end) ...
			+ q.data.ChannelIndex(2:end))];


		% check if p and q are from the same array
		if(~strcmp(p.data.arrstr(end,:),q.data.arrstr))
            % concatenate p and q
			r.data.arrstr = strvcat(p.data.arrstr, q.data.arrstr);
			r.data.ArrayIndex = [p.data.ArrayIndex; (p.data.ArrayIndex(end) ...
				+ q.data.ArrayIndex(2:end))];
		else
			% just update array index
			r.data.ArrayIndex(2) = p.data.ArrayIndex(2)+ q.data.ArrayIndex(2);
		end

		% check if p and q are from the same session
		if(~strcmp(p.data.sesstr(end,:),q.data.sesstr))
            % concatenate p and q
			r.data.sesstr = strvcat(p.data.sesstr, q.data.sesstr);
			r.data.SessionIndex = [p.data.SessionIndex; (p.data.SessionIndex(end) ...
				+ q.data.SessionIndex(2:end))];
		else
			% just update array index
			r.data.SessionIndex(2) = p.data.SessionIndex(2)+ q.data.SessionIndex(2);
		end

		% check if p and q are from the same day
		if(~strcmp(p.data.daystr(end,:),q.data.daystr))
            % concatenate p and q
			r.data.daystr = strvcat(p.data.daystr, q.data.daystr);
			r.data.DayIndex = [p.data.DayIndex; (p.data.DayIndex(end) ...
				+ q.data.DayIndex(2:end))];
		else
			% just update array index
			r.data.DayIndex(2) = p.data.DayIndex(2)+ q.data.DayIndex(2);
		end
					
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
