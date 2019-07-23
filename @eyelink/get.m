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

Args = struct('ObjectLevel',0, 'AnalysisLevel',0, 'Number', 0, 'Trial', 0, 'XY', 0, 'Calibration', 0);
Args.flags ={'ObjectLevel','AnalysisLevel', 'Number', 'Trial', 'XY', 'Calibration'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.ObjectLevel)
	% specifies that the object should be created in the session directory
	r = 'Session';
elseif(Args.AnalysisLevel)
	% specifies that the AnalysisLevel of the object is 'AllIntragroup'
	r = 'Single';
elseif (Args.Number && (Args.Trial || Args.XY))
    %get the number of trials that are in the object
    %add every third column of the matrix 
    
   n = obj.data.noOfTrials;
   disp (n);
   r = sum (n,2);
elseif (Args.Number && Args.Calibration) 
    n = size(obj.data.trial_timestamps,1); 
else 
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end