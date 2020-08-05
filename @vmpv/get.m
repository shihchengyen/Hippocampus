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

Args = struct('ObjectLevel',0, 'AnalysisLevel',0,'SpeedLimit',0);
Args.flags ={'ObjectLevel','AnalysisLevel'};
Args = getOptArgs(varargin,Args);

% set variables to default
r = [];

if(Args.ObjectLevel)
	% specifies that the object should be created in the session directory
	r = levelConvert('levelNo',1);
elseif(Args.AnalysisLevel)
	% specifies that the AnalysisLevel of the object is 'AllIntragroup'
	r = 'Single';
    
elseif(Args.SpeedLimit)
    
    % get row by row filter for which sessionTimeC rows are valid
    speeding_checker = [obj.data.unityTime [0; obj.data.unityData(:,7)./obj.data.unityData(:,2)]];
    speeding_checker(:,2) = speeding_checker(:,2) > Args.SpeedLimit;
    speeding_checker(:,3) = [0; diff(speeding_checker(:,2))];
    if speeding_checker(1,2) == 0
        speeding_checker(1,3) = -1;
    else
        speeding_checker(1,3) = 1;
    end
    stop_intervals = nan(max([sum(speeding_checker(:,3)==1) sum(speeding_checker(:,3)==-1)]), 2);
    % stores intervals that have low speed
    stop_intervals(1:sum(speeding_checker(:,3)==-1),1) = speeding_checker(find(speeding_checker(:,3)==-1),1);
    stop_intervals(1:sum(speeding_checker(:,3)==1),2) = speeding_checker(find(speeding_checker(:,3)==1),1);
    % remove intervals with same time stamp
    stop_intervals(find((stop_intervals(:,2)-stop_intervals(:,1))==0),:) = [];
    stop_intervals(:,3) = [stop_intervals(2:end,1) ;NaN];
    nextrowsame = find(stop_intervals(:,2) == stop_intervals(:,3));
    stop_intervals(nextrowsame,2) = stop_intervals(nextrowsame+1,2);
    stop_intervals(nextrowsame+1,:) = [];
    stop_intervals(:,3) = [];
    
    stc = obj.data.sessionTimeC;
    
    delta_t = find(diff(stc(:,1))>0);
    delta_t = [[1; delta_t+1] [delta_t; size(stc,1)]];
    delta_t = [delta_t stc(delta_t(:,1),1)];

    histbins = stop_intervals';
    histbins = histbins(:);
    % Remove intervals outside of end of sessionTimeC
    if histbins(find(isnan(histbins))-1,1) > max([stc(end,1) delta_t(end,3)])
        histbins(find(isnan(histbins))-1:find(isnan(histbins))) = [];
    else
        histbins(isnan(histbins)) = max([stc(end,1) delta_t(end,3)]);
    end

    [~,~,binned] = histcounts(delta_t(:,3), histbins);

    delta_t(:,4) = rem(binned,2); % ones indicate within interval

    r = ~repelem(delta_t(:,4), delta_t(:,2)-delta_t(:,1)+1); 
    
else
	% if we don't recognize and of the options, pass the call to parent
	% in case it is to get number of events, which has to go all the way
	% nptdata/get
	r = get(obj.nptdata,varargin{:});
end