function [output] = placeandview_sic_transfer_in(current_path)
    
    cd(current_path);
    disp(current_path);
    split_path = strsplit(current_path, '/');
    
    partial_path = strjoin(split_path(end-4:end),'/'); % 20181102/session01/array01/channel019/cell01
    command = 'scp -P 8398 hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
    command = strcat(command, partial_path);
    command = strcat(command, '/vmsv.mat .');
    unix(command);
    
    partial_path = strjoin(split_path(end-4:end),'/'); % 20181102/session01/array01/channel019/cell01
    command = 'scp -P 8398 hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
    command = strcat(command, partial_path);
    command = strcat(command, '/spiketrain.mat .');
    unix(command);
    
    partial_path = strjoin(split_path(end-4:end-3),'/'); % 20181102/session01
    command = 'scp -P 8398 hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
    command = strcat(command, partial_path);
    command = strcat(command, '/umaze.mat .');
    unix(command);
    
    partial_path = strjoin(split_path(end-4:end-3),'/'); % 20181102/session01
    command = 'scp -P 8398 hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
    command = strcat(command, partial_path);
    command = strcat(command, '/rplparallel.mat .');
    unix(command);
    
    command = 'scp -P 8398 hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
    command = strcat(command, 'pillars.mat .');
    unix(command);

%     placeandview_shuffle(1000,100);
%     
%     partial_path = strjoin(split_path(end-4:end),'/'); % 20181102/session01/array01/channel019/cell01
%     command = 'scp -P 8398 pnv_shuffle.mat hippocampus@cortex.nus.edu.sg:/volume1/Hippocampus/Data/picasso-misc/';
%     command = strcat(command, partial_path);
%     unix(command);
%     unix('rm *.mat');
    
end



