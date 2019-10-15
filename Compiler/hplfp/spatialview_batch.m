function [] = spatialview_batch()

% cwd = pwd;
cwd = '/Volumes/User/huimin/Desktop/condor_shuffle/';
% load celllist from text file
% celllist = textread('/Volumes/User/huimin/Desktop/condor_shuffle/cell_list copy.txt','%c');
fid = fopen([cwd '/cell_list sv.txt'],'rt');
cellList = textscan(fid,'%s','Delimiter','\n');
cellList = cellList{1};
% Make sure no empty cells
notempty = ~cellfun(@isempty,cellList);
cellList = cellList(notempty,:);

for ii = 1:size(cellList,1)
    
    if exist(cellList{ii},'dir') == 7
        cd(cellList{ii});
    else
        continue;
    end
%     if exist('spatialview.mat','file') == 2
%         delete('spatialview.mat');
%     end
    spatialview('auto','RedoLevels',1,'SaveLevels',1,'NumShuffles',10,'GridSteps',40,'UseAllTrials',1,'FiltLowOcc',0);
    
end

cd(cwd);
