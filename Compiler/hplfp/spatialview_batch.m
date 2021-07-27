function [] = spatialview_batch()

% cwd = pwd;
cwd = '/Volumes/User/huimin/Desktop/condor_shuffle/';
% load celllist from text file
% celllist = textread('/Volumes/User/huimin/Desktop/condor_shuffle/cell_list copy.txt','%c');
fid = fopen([cwd '/cell_list.txt'],'rt');
cellList = textscan(fid,'%s','Delimiter','\n');
cellList = cellList{1};
% Make sure no empty cells
notempty = ~cellfun(@isempty,cellList);
cellList = cellList(notempty,:);

% % Place
% for ii = 1:size(cellList,1)
%     
%     if exist(cellList{ii},'dir') == 7
%         cd(cellList{ii});
%         
%     else
%         continue;
%     end
% %     if exist('spatialview.mat','file') == 2
% %         delete('spatialview.mat');
% %     end
%     disp([cellList{ii} ': vmpc']);
%     vmpc('auto','RedoLevels',1,'SaveLevels',1);
%     
% end

% Spatial view
% day = strsplit(cellList{1},'picasso-misc/');
% day = day{2}(1:8);
for ii = 1:size(cellList,1)
    
    if exist(cellList{ii},'dir') == 7
        cd(cellList{ii});
        
    else
        continue;
    end
%     if exist('spatialview.mat','file') == 2
%         delete('spatialview.mat');
%     end
    
%     day2 = strsplit(cellList{1},'picasso-misc/');
%     day2 = day2{2}(1:8);

    disp([cellList{ii} ': vmsv']);
    vmsv('auto','RedoLevels',1,'SaveLevels',1);
    
end

cd(cwd);
