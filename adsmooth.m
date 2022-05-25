function [smoothedRate,smoothedSpk,smoothedDur]=adsmooth(dur,spk,alpha)
% Adaptive smoothing of rate maps.
% *** Modified to work for n shuffles. dur = b x b x n for each grid
%
%       [smoothedRate,smoothedSpk,smoothedPos]=rates_adaptivesmooth(posMap,spkMap,alpha)
%
% Each bin is smoothed using a flat, circular kernal. The circle radius 
% is set for each bin, indivdually, such that 
%
%   radius => alpha ./ ( sqrt(nSpike) .* nDwell )
%
% where nSpike and nDwell are the number of spikes, and the amount of dwell time (in s) within the kernel.
%
% smoothedRate, smoothedSpk, smoothedPos are the smoothed maps (spike and pos maps are smoothed 
% with the same kernal group as for the rate map.

% % Check for empty spk maps %
% for ii = 1:size(spk,3)
%     if sum(sum(spk(:,:,ii)))==0
%         smoothedDur = zeros(size(dur));
%         smoothedSpk = zeros(size(spk));
%         smoothedRate = zeros(size(spk));
%         smoothedDur(:,:,ii)=dur(:,:,ii);    smoothedDur(dur(:,:,ii)==0,ii)=nan;
%         smoothedSpk(:,:,ii)=spk(:,:,ii);    smoothedSpk(dur(:,:,ii)==0,ii)=nan;
%         smoothedRate(:,:,ii)=spk(:,:,ii);   smoothedRate(dur(:,:,ii)==0,ii)=nan;
%         return
%     end
% end
% Check for empty spk maps % Only for no shuffles
if size(spk,3)==1 && sum(sum(spk))==0
    smoothedDur=dur;    smoothedDur(dur==0)=nan;
    smoothedSpk=spk;    smoothedSpk(dur==0)=nan;
    smoothedRate=spk;   smoothedRate(dur==0)=nan;
    return
end
% Pre-assign output %
smoothedDur=zeros(size(dur));
smoothedSpk=zeros(size(dur));
% Visited env template: use this to get numbers of visited bins in filter at edge of environemnt %
vis=zeros(size(dur));
vis(dur>0)=1;
% Pre-assign map which records which bins have passed %
smoothedCheck=false(size(dur));
smoothedCheck(dur==0)=true; % Disregard unvisited - mark as already done.
% Pre-assign list of radii used (this is for reporting purposes, not used for making maps) %
radiiUsedList=nan(1,sum(sum(sum(dur>0)))); %%% FIX
% radiiUsedList=nan(1,sum(sum(dur>0)));
radiiUsedCount=1;

%%% Run increasing radius iterations %%%
r=1; % Circle radius
boundary=0; % IMFILTER boundary condition
while any(any(any(~smoothedCheck)))
    % Check radius isn't getting too big (if >map/2, stop running) %
    if r>max(size(dur(:,:,1)))/2
%     if r>max(size(dur))/2
%     if r>20
        smoothedSpk(~smoothedCheck)=nan;
        smoothedDur(~smoothedCheck)=nan;
        break
    end
    % Construct filter kernel ...
    % Place: Flat disk, where r>=distance to bin centre %
    f=fspecial('disk',r); 
    f(f>=(max(max(f))/3))=1;
    f(f~=1)=0;   
    % Filter maps (get N spikes and pos sum within kernel) %
    fSpk=imfilter(spk,f,boundary);
    fDur=imfilter(dur,f,boundary);
    fVis=imfilter(vis,f,boundary);
    % Which bins pass criteria at this radius? %
    warning('off', 'MATLAB:divideByZero');
    binsPassed=alpha./(sqrt(fSpk).*fDur) <= r;
    warning('on', 'MATLAB:divideByZero');
    binsPassed=binsPassed & ~smoothedCheck; % Only get the bins that have passed in this iteration.
    % Add these to list of radii used %
    nBins=sum(sum(sum(binsPassed)));
%     nBins=sum(binsPassed(:));
    radiiUsedList(radiiUsedCount:radiiUsedCount+nBins-1)=r;
    radiiUsedCount=radiiUsedCount+nBins;
    % Assign values to smoothed maps %
    smoothedSpk(binsPassed)=fSpk(binsPassed)./fVis(binsPassed);
    smoothedDur(binsPassed)=fDur(binsPassed)./fVis(binsPassed);
    % Record which bins were smoothed this iteration %
    smoothedCheck(binsPassed)=true;
    % Increase circle radius (half-bin steps) %
    r=r+0.5; % Increase radius in 0.5 bin steps.
end

% Assign Output %
warning('off', 'MATLAB:divideByZero');
smoothedRate=smoothedSpk./smoothedDur;
warning('on', 'MATLAB:divideByZero');
smoothedRate(dur==0)=nan;
smoothedDur(dur==0)=nan;
smoothedSpk(dur==0)=nan;

% Report radii sizes %
if 0
    hAllFigs = get(0, 'children');
    hFig = findobj(hAllFigs, 'flat', 'tag', 'adaptiveSmoothPlotWindow');
    if isempty(hFig);
        hFig=figure;
        set(hFig,'tag','adaptiveSmoothPlotWindow');
    else
        figure(hFig);
    end
    hist(radiiUsedList,1:10);
    uiwait(hFig,1.5);
end