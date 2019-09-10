function mountain_batch(target)

    if ~exist('target','var')
        disp('error function, specify arguments - all, hippo, <numbers>');
        return;
    end
    
%     if ~exist('mountains', 'dir')
    mkdir('mountains');
%     end

    full_cell = comb_channels_ms;

    for i = 1:size(full_cell, 1)

        channel_no = extractAfter(full_cell{i,1}, 'channel');
        channel_no = string(str2double(channel_no));
        
        if string(target) == 'all'
            mountain_channel(full_cell, i);
        elseif string(target) == 'hippo'
            for j = 3:2+full_cell{i,2}
                session_extract = extractAfter(full_cell{i,j}, 'session');
                if startsWith(session_extract, '0')
                    mountain_channel(full_cell, i);
                    break;
                end
            end
        else
            splits = strsplit(string(target), ' ');
            for j = 1:length(splits(1,:))
                if splits{1,j} == channel_no
                    mountain_channel(full_cell, i);
                    break;
                end
            end
        end
    end

end

function mountain_channel(full_cell, index)

    disp('processing channel');
    disp(full_cell{index,1});
    
    origin = pwd;
    cd('mountains');
    
    if ~exist(full_cell{index,1}, 'dir')
        mkdir(full_cell{index,1});
    else
        disp('overwriting existing folder in mountains');
        cd(full_cell{index,1});
        unix('rm -r *');
        cd('..');
    end
    
    cd(full_cell{index,1});
    start_indices = cell(2,full_cell{index,2});
    current = 1;
    index1 = 1;
    mc_path = pwd;
    for i = 3:2+full_cell{index,2}
        cd(full_cell{index,i});
        cd(full_cell{index,1});
        disp('debugging clones!');
        disp(pwd);
        start_indices{1,index1} = current;
        cut = strsplit(full_cell{index,i}, '/');
        start_indices{2,index1} = cut{1,length(cut)-1};
        if i == 3
            data = rplhighpass('auto');
            data = data.data.analogData';
        else
            temp = rplhighpass('auto');
            temp = temp.data.analogData';
            data = [data temp];
        end
        index1 = index1 + 1;
        current = length(data) + 1;
    end
    
    cd(mc_path);
    save('start_times.mat','start_indices');
    mkdir('dataset');
    cd('dataset');
    disp('debugging clones');
    disp(max(data));
    disp(min(data));
    
    writemda(data, 'raw_data.mda', 'float32');
    unix('cp $GITHUB_MATLAB/Hippocampus/mountains/geom.csv .');
    cd('..');
    unix('cp $GITHUB_MATLAB/Hippocampus/mountains/sort.sh .');
    unix('source ~/.bash_profile; conda activate mountainlab; sh sort.sh')
    
    disp('finished for this channel');
    
    cd(origin);

end


function [channels_identified] = comb_channels_ms()

    origin = pwd;
    full_list = dir();
    dirFlags = [full_list.isdir];

    top_folders = full_list(dirFlags);
    top_folders = top_folders(3:length(top_folders));
    
    
    session_names = cell(1,length(top_folders));
    count = 1;
    for i = 1:length(session_names)
        if strncmpi(top_folders(i).name, 'session0', 8) || strncmpi(top_folders(i).name, 'sessioneye', 10)
            session_names{1,count} = top_folders(i).name;
            count = count + 1;
        end
    end
    
    session_names = session_names(1,1:count-1);
        
    day_path = pwd;
    channels_identified = cell(150,count+1);
    identified_count = 0;
    
    for i = 1:length(session_names)
        
        cd(strcat('./',session_names{1,i}));
        
        sub_list = dir();
        dirFlags = [sub_list.isdir];
        sub_folders = sub_list(dirFlags);
        sub_folders = sub_folders(3:length(sub_folders));
        
        array_folders = cell(1,length(sub_folders));
        count = 1;
        
        for j = 1:length(sub_folders)
            if strncmpi(sub_folders(j).name, 'array', 5) == 1
                array_folders{1,count} = sub_folders(j).name;
                count = count + 1;
            end
        end
        array_folders = array_folders(1,1:count-1);
        
        for a = 1:length(array_folders(1,:))
            cd(strcat('./',array_folders{1,a}));
            
            channel_list = dir();
            dirFlags = [channel_list.isdir];
            channel_list = channel_list(dirFlags);
            channel_list = channel_list(3:length(channel_list));

            for k = 1:length(channel_list)
%                 fprintf('%s %s %s\n',session_names{1,i}, array_folders{1,a}, channel_list(k).name);
                if strncmpi(channel_list(k).name, 'channel', 7) == 1
                    found = 0;
                    if identified_count == 0
                        channels_identified{1,1} = channel_list(k).name;
                        channels_identified{1,2} = 1;
                        channels_identified{1,3} = channel_list(k).folder;
                        identified_count = identified_count + 1;
                        found = 1;
                    else
                        for j = 1:identified_count
                            if channels_identified{j,1} == channel_list(k).name
                                channels_identified{j,channels_identified{j,2}+3} = channel_list(k).folder;
                                channels_identified{j,2} = channels_identified{j,2} + 1;
                                found = 1;
                            end
                        end
                    end
                    if found == 0
                        channels_identified{identified_count+1,1} = channel_list(k).name;
                        channels_identified{identified_count+1,2} = 1;
                        channels_identified{identified_count+1,3} = channel_list(k).folder;
                        identified_count = identified_count + 1;
                    end
                end
            end
            
            cd('..');
        end
            
        cd(day_path);
        
    end
    
    channels_identified = channels_identified(1:identified_count,:);
    
    cd(origin);
    
end