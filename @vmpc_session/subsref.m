function [b,res] = subsref(obj,index)
%@dirfiles/subsref Index function for DIRFILES object.

res = 1;
unknown = 0;

indlength = length(index);

switch index(1).type
case '.'
	switch index(1).subs
	case 'data'
		% if this is the only index, just return the entire thing
		if(indlength==1)
			b = obj.data;
		else
			% there are more than 1 index so pass on the subsref call
			b = subsref(obj.data,index(2:end));
		end
	otherwise
		unknown = 1;
	end
otherwise
	unknown = 1;
end

if(unknown)
	[b,res] = subsref(obj.nptdata,index);
end

if(res == 0)
	if(isempty(inputname(1)))
		% means some other function is calling this function, so
		% just return error instead of printing error
		res = 0;
		b = 0;
	else
		error('Invalid field name');
	end
end
