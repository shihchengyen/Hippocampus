% get each cell's directory with the spiketrain.mat
% calculate each vmsv

% % % % % % % need to deal with missing file

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
% % change directory
% for ss = 1:size(cellList,1)
% cellList(ss)=str2cell(strrep(cell2str(cellList(ss)),'Volumes','Users/yuhsuan'));
% end    

% % Make sure no empty cells
% notempty = ~cellfun(@isempty,cellList);
% cellList = cellList(notempty,:);
% % Generate unique identifier for each cell
% s = regexp(cellList{1},'session');
% identifiers = zeros(size(cellList,1),5);
% cellid = cell(size(cellList,1),1);
% missing = [];
% for ii = 1:size(cellList,1)
%     if exist(cellList{ii},'dir') == 7
%         % Collect date, session, array, channel, cell
%         identifiers(ii,:) = [str2double(cellList{ii}(s-9:s-2)) str2double(cellList{ii}(s+7:s+8)) ...
%             str2double(cellList{ii}(s+15:s+16)) str2double(cellList{ii}(s+25:s+27)) str2double(cellList{ii}(s+33:s+34))];
%         % Cell identifier
%         cellid{ii} = horzcat(num2str(identifiers(ii,4)),'-',num2str(identifiers(ii,5)));
%     else
%         missing = [missing ii];
%     end
% end
% % Remove missing cells
% identifiers(missing,:) = [];
% cellid(missing) = [];
% cellList(missing) = [];
% setsessions = unique(identifiers(:,1));
% % setcells = unique(cellid);
%     
    
% % % % %     generate new vmsv
for ss = 1:size(cellList,1)
    filename=strrep(cell2str(cellList(ss)),'Volumes','Users/yuhsuan');
    if isfile(filename)
        if isfile('spiketrain.mat')
            cd(filename);%directory specific to my workstation
            sv=vmsv('auto');
        %     get each cell's ise
            ise=[ise(:) sv.data.ISE];
        %     get each cell's ise_threshold 95% percentile
            ise_thr =[ise_thr(:) prctile([sv.data.ISE; sv.data.ISEsh(:)],95)];
        end
    end
end
% % % % % % % % old
% % Load vmpv object for each session
% for ss = 1:size(setsessions,1)
% 
%     cd(['/Users/yuhsuan/Hippocampus/Data/picasso-misc/' num2str(setsessions(ss)) '/session01']);
%     pv = load([num2str(pix) 'vmpv.mat']);
%     pv = pv.pv;
%     setcells = find(identifiers(:,1) == setsessions(ss));
%        
%     % Process placebyspatialview cell by cell
%     for cc = 1:size(setcells,1)
%         
%         cd(cellList{setcells(cc)});
%         disp(cellList{setcells(cc)});
%         sv=vmsv('auto');
% %       get each cell's ise
%         ise=[ise(:) sv.data.ISE];
% %       get each cell's ise_threshold 95% percentile
%         ise_thr =[ise_thr(:) prctile([sv.data.ISE; sv.data.ISEsh(:)],95)];
%    
%         close all;%from old file
%         
%     end
%     
% end




% % % % % % % % % 
% average the each cell's ise threshold to get overall threshold
    ise_thr_avg=mean(ise_thr(:),'All');
% find the difference between the ise and overall threshold
%     each cell's ise-overall threshold
       difference=ise-ise_thr_avg;
%     each cell's ise/overall threshold
        ratio=ise/ise_thr_avg;
 save('/Desktop/test_ise.mat');
end



