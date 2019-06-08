function [r,varargout] = get(obj,varargin)
%eyelink/get Get function for eyelink objects
%eyelink/GET Returns object properties
%   VALUE = GET(OBJ,PROP_NAME) returns an object 
%   property.
%   In eyelink, PROP_NAME can be one of the following:
%      'ObjectLevel'
%      'AnalysisLevel'
%
%   Dependencies: 

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'Number', 0, 'Trial', 0);
Args.flags ={'ObjectLevel','AnalysisLevel', 'Number', 'Trial'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.ObjectLevel)
	% specifies that the object should be created in the session directory
	r = 'Day';
elseif(Args.AnalysisLevel)
	% specifies that the AnalysisLevel of the object is 'AllIntragroup'
	r = 'Single';
elseif (Args.Number && Args.Trial)
    %get the number of trials that are in the object
    %add every third column of the matrix 
    
   n = obj.data.noOfTrials;
   disp (n);
   r = sum (n,2);
    
else 
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end