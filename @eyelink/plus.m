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
        r.data.noOfTrials = concat(p.data.noOfTrials, q.data.noOfTrials, 'Columnwise');
        r.data.noOfSessions = p.data.noOfSessions + q.data.noOfSessions;

		% object specific fields
		r.data.trial_timestamps = concat (p.data.trial_timestamps, q.data.trial_timestamps, 'Columnwise');
        r.data.sacc_event = concat (p.data.sacc_event, q.data.sacc_event, 'Columnwise');
        r.data.fix_event = concat (p.data.fix_event, q.data.fix_event, 'Columnwise');
        %Timestamps takes the max time that the exp were run for  
        if (size(p.data.timestamps,1) > size (q.data.timestamps))
            r.data.timestamps = p.data.timestamps;
        else 
            r.data.timestamps = q.data.timestamps;
        end 
        
        %contains all the eye positions for all the sessions 
        r.data.eye_pos = concat(p.data.eye_pos, q.data.eye_pos);  
        %OOF
        r.data.timeouts = concat(p.data.timeouts, q.data.timeouts, 'Rowise');
        
		% add nptdata objects as well
		r.nptdata = plus(p.nptdata,q.nptdata);
	end
end
