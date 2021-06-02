% get each cell's directory with the spiketrain.mat
% calculate each vmsv

% % % % % % % need to deal with missing file

function [obj, varargout] = test_ise(savefig,varargin)
% Batch 
    % Load cell list
    cwd = '/Users/yuhsuan/Desktop';
%     Analysis'
    fid = fopen([cwd '/cell_list copy 2.txt'],'rt');
    cellList = textscan(fid,'%s','Delimiter','\n');
    cellList = cellList{1};

% Make sure no empty cells
notempty = ~cellfun(@isempty,cellList);
cellList = cellList(notempty,:);
% Generate unique identifier for each cell
s = regexp(cellList{1},'session');
identifiers = zeros(size(cellList,1),5);
cellid = cell(size(cellList,1),1);
missing = [];
    
% Load vmpv object for each session
ise=[];
ise_thr =[];
for ss = 1:size(cellList,1)

% % % %     generate new vmsv
    cd(strrep(cellList{ss},'/spiketrain.mat',''));
%     cd ..\ %up one to cell01
    sv=vmsv('auto');
%     get each cell's ise
    ise=[ise; sv.data.ISE];
%     get each cell's ise_threshold 95% percentile
    ise_thr =[ise_thr; prctile([sv.data.ISE; sv.data.ISEsh(:)],95)];
end
% average the each cell's ise threshold to get overall threshold
    ise_thr_avg=mean(ise_thr(:));
% find the difference between the ise and overall threshold
%     each cell's ise-overall threshold
       difference=ise-ise_thr_avg;
%     each cell's ise/overall threshold
        ratio=ise/ise_thr_avg;
 save('Users/yuhsuan/Desktop/test_list2ise.mat');
end



