function [r,varargout] = get(obj,varargin)
%dirfiles/get Get function for dirfiles objects
%dirfiles/GET Returns object properties
%   VALUE = GET(OBJ,PROP_NAME) returns an object 
%   property.
%   In dirfiles, PROP_NAME can be one of the following:
%      'ObjectLevel'
%	 'AnalysisLevel'
%
%   Dependencies: 

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'GroupPlotProperties',0, ...
				'Number',0, 'Trial',0,'NumProcessTrials',0,'ProcessTrialsDiff',-100);
Args.flags ={'ObjectLevel','AnalysisLevel','GroupPlotProperties','Number', ...
    'Trial','NumProcessTrials'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.ObjectLevel)
	% specifies that the object should be created in the session directory
	r = 'Session';
elseif(Args.AnalysisLevel)
	% specifies that the AnalysisLevel of the object is 'AllIntragroup'
	r = 'Single';
elseif(Args.GroupPlotProperties)
	r = 'Vertical';
elseif(Args.Number && Args.Trial)
	r = obj.data.setIndex(end);
elseif(Args.NumProcessTrials)
    % get processTrial field
    pT = obj.data.processTrials;
    % take the difference between consecutive points to find the transition
    % between sessions
    di = diff(pT);
    % find the indices where the difference is more negative than -100
    dim = find(di<Args.ProcessTrialsDiff);
    % find the difference between those indices to get the number of trials
    % per session
    dimd = diff(dim);
    % this will be missing the number of trials from the first and last 
    % sessions so we will need to add and compute those separately
    r = [dim(1); dimd; size(pT,1)-dim(end)];
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end