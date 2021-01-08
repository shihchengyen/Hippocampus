function combine_spiketrain(cell_dirs_file)
%
%   cell_dirs_file (text file) contains directories of cells to perform the combine
%   for.
%   example line: '/Volumes/hippocampus/data/picasso-misc/20181101/session01/array01/ch019/cell01'
%   The session number does not matter, and only a single session for each cell
%   needs to be specified.
%

fid = fopen(cell_dirs_file,'rt');
cellList = textscan(fid,'%s','Delimiter','\n');
cellList = cellList{1};
fclose(fid);

done_cells = {};

for c = 1:size(cellList,1)
    disp(['Reading ' cell_dirs_file ': ' cellList{c}])
    done_cells{end+1} = combine_spiketrain_single(cellList{c}, done_cells);
end

end

function output_dir = combine_spiketrain_single(cell_dir, done_cells_input)

cd(cell_dir)
cd ../../../..
day_dir = pwd;
sessions = dir('session0*');
cd session01

rp = rplparallel('auto'); % load only the combined rplparallel
if rp.data.numSets ~= size(sessions,1)
    error('numSets in rplparallel does not match number of session folders found.')
elseif size(sessions,1) == 1
    error('There is only 1 session folder.')
end

a = regexp(cell_dir,'array');
cell_subdir = cell_dir(a:end);

cd(cell_subdir)
output_dir = pwd;
if any(strcmp(done_cells_input, output_dir))
    disp(['Duplicate entry found, skipping line: ' cell_dir])
    return
end

spiketrain = load('spiketrain.mat');
full_ts = spiketrain.timestamps;
session_start_markers = find(rp.data.Args.Data.markers == 84);

for s = 2:size(sessions,1)
    cd([day_dir '/session0' num2str(s) '/' cell_subdir])
    spiketrain = load('spiketrain.mat');
    front_shift = 1000*rp.data.Args.Data.timeStamps(session_start_markers(s) - 1); % corresponds to the timestamp of the last marker of the previous session
    full_ts = [full_ts (spiketrain.timestamps + front_shift - 1000*rp.data.session_start_all(s))]; % similar to how the timestamps of rplparallel are combined
end

cd([day_dir '/session01/' cell_subdir])
timestamps = full_ts;
copyfile spiketrain.mat spiketrain_original.mat % backup original
save('spiketrain.mat','timestamps')


end

