% get each cell's directory with the spiketrain.mat
% calculate each vmsv

% % % % % % % need to deal with missing file

function [obj, varargout] = test_ise(savefig,varargin)
% Batch 
    % Load cell list
    cwd = '/Users/yuhsuan/Desktop';
%     Analysis'
    fid = fopen([cwd '/cell_list copy 3.txt'],'rt');
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
for ii = 1:size(cellList,1)
    if exist(cellList{ii},'dir') == 2 
        %also want to check vmpv.ma file
        s=cellList{ii};
        array_p= strfind(s,'array');%position of 'array'
        filename= [s(1:array_p-2),'vmpv.mat'];%vmpv filename
        if exist(filename,'dir') == 2 
            % Collect date, session, array, channel, cell
            identifiers(ii,:) = [str2double(cellList{ii}(s-9:s-2)) str2double(cellList{ii}(s+7:s+8)) ...
                str2double(cellList{ii}(s+15:s+16)) str2double(cellList{ii}(s+25:s+27)) str2double(cellList{ii}(s+33:s+34))];
            % Cell identifier
            cellid{ii} = horzcat(num2str(identifiers(ii,4)),'-',num2str(identifiers(ii,5)));
        else
        missing = [missing ii];
        end
    end
end
% Remove missing cells
identifiers(missing,:) = [];
cellid(missing) = [];
cellList(missing) = [];
setsessions = unique(identifiers(:,1));
% setcells = unique(cellid);

% %use identifier to make file name tag
% for n=1:size(identifiers,1)
%     trail_name=[trail_name ;[string(identifiers(n,1)), '-session',string(identifiers(n,2)),'array',string(identifiers(n,3)),'channel',string(identifiers(n,4)),'cell',string(identifiers(n,5))];
% end

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
%     ise_result=[identifiers ise ise_thr difference ratio];
%  save('/Users/yuhsuan/Desktop/test_list2ise_ted.mat');
 save('/Users/yuhsuan/Desktop/test_list2ise_ted.mat','identifiers','ise','ise_thr','ise_thr_avg','difference','ratio');
end



