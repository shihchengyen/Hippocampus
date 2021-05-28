% get each cell's directory with the spiketrain.mat
% calculate each vmsv
    
function [obj, varargout] = test_ise(savefig,varargin)
% Batch 
    % Load cell list
    cwd = '/Users/yuhsuan/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis';
%old: cwd = '/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current
%     Analysis'
    fid = fopen([cwd '/cell_list.txt'],'rt');
%     fid = fopen([cwd '/cell_list_1pxFiltAll.txt'],'rt');
    cellList = textscan(fid,'%s','Delimiter','\n');
    cellList = cellList{1};

% Load vmpv object for each session
for ss = 1:size(cellList,1)

% % % %     generate new vmsv
    cd(strrep(cell2str(cellList(ss))),'Volumes','Users/yuhsuan');%directory specific to my workstation
    sv=vmsv('auto');
%     get each cell's ise
    ise=[ise(:) sv.data.ISE];
%     get each cell's ise_threshold 95% percentile
    ise_thr =[ise_thr(:) prctile([sv.data.ISE; sv.data.ISEsh(:)],95)];
end
% average the each cell's ise threshold to get overall threshold
    ise_thr_avg=mean(ise_thr(:));
% find the difference between the ise and overall threshold
%     each cell's ise-overall threshold
       difference=ise-ise_thr_avg;
%     each cell's ise/overall threshold
        ratio=ise/ise_thr_avg;
 save('/Desktop/test_ise.mat');
end



