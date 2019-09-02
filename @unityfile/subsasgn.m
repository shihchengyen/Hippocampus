function [obj,res] = subsasgn(obj,index,value)
%@dirfiles/subsasgn Assignment function for DIRFILES object

res = 1;
unknown = 0;

indlength = length(index);

switch index(1).type
case '.'
	switch index(1).subs
	case 'data'
		% if this is the only index, just replace the entire thing
		if(indlength==1)
			obj.data = value;
		else
			% there are more than 1 index so pass on the subsasgn call
			obj.data = subsasgn(obj.data,index(2:end),value);
		end
	otherwise
		unknown = 1;
	end
otherwise
	unknown = 1;
end

if(unknown)
	[obj.nptdata,res] = subsasgn(obj.nptdata,index,value);
end

if(res == 0)
	if(isempty(inputname(1)))
		% means some other function is calling this function, so
		% just return error instead of printing error
		res = 0;
	else
		error('Invalid field name');
	end
end
