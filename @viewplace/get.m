function [r,varargout] = get(obj,varargin)
%viewsort/get Get function for viewsort objects
%viewsort/GET Returns object properties
%   VALUE = GET(OBJ,PROP_NAME) returns an object 
%   property.
%   In viewsort, PROP_NAME can be one of the following:
%      'ObjectLevel'
%	 'AnalysisLevel'
%
%   Dependencies: 

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'GroupPlotProperties',0, ...
				'Number',0, 'Channel',0, 'Array',0, 'Session',0, 'Day',0);
Args.flags ={'ObjectLevel','AnalysisLevel','GroupPlotProperties','Number', ...
				'Channel', 'Array', 'Session', 'Day'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.ObjectLevel)
	% specifies that the object should be created in the session directory
	r = 'Channel';
elseif(Args.AnalysisLevel)
	% specifies that the AnalysisLevel of the object is 'AllIntragroup'
	r = 'Single';
elseif(Args.GroupPlotProperties)
	r = 'Vertical';
elseif(Args.Number)
	if(Args.Array)
		r = size(obj.data.arrstr,1);
	elseif(Args.Session)
		r = size(obj.data.sesstr,1);
	elseif(Args.Day)
		r = size(obj.data.daystr,1);
	elseif(Args.Channel)
		r = obj.data.numChannels;
	else
		r = size(obj.data.locSI,1);
	end
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end