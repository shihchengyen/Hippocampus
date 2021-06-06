% Load cell list
    cwd = '/Users/yuhsuan/Desktop';
%     Analysis'
    fid = fopen([cwd '/cell_list copy.txt'],'rt');
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
    if exist(cellList{ii},'file') == 2 
        %also want to check vmpv.ma file
        st=cellList{ii};
        array_p= strfind(st,'array');%position of 'array'
        filename= [st(1:array_p-1),'vmpv.mat'];%vmpv filename
        if exist(filename,'file') == 2 
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
% remove already genreated vmsv.m


% Remove missing cells
identifiers(missing,:) = [];
cellid(missing) = [];
cellList(missing) = [];
% writecell(cellList,'/Users/yuhsuan/Documents/MATLAB/cell_list.txt');
writecell(identifiers(:,1),'/Users/yuhsuan/Documents/MATLAB/cell_list.txt');