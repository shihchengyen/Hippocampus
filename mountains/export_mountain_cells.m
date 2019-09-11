function export_mountain_cells(pruned)

    if ~exist('pruned','var')
        original = readmda('prespiketrain.mda');
    else
        original = pruned;
    end
        
    cutting = load('start_times.mat');
    cutting = cutting.start_indices;

    ch = strcat('channel', extractAfter(pwd, 'channel'));
    total_ind = length(original(2,:));
    
    marker = 1;
    count = 1;
    for i = 1:length(original(2,:))
        if count == size(cutting, 2)
            exporting_arr = original(2:3,marker:total_ind);
            location = pwd;
            cd('../../');
            cd(cutting{2,count});
            command = ['for j in `find . -name "', ch, '"`; do echo $j; done'];
            command = command(1) + command(2) + command(3);
            [~, found] = unix(command{:});     
            cd(found(1:10));
            split_into_cells_intra_session(ch, exporting_arr, cutting{1,count});
            cd(location);            
            break;
        end
        if original(2,i) >= cutting{1,count+1}
            exporting_arr = original(2:3,marker:i-1);
            location = pwd;
            cd('../../');
            cd(cutting{2,count});
            command = ['for j in `find . -name "', ch, '"`; do echo $j; done'];
            command = command(1) + command(2) + command(3);
            [~, found] = unix(command{:});
            cd(found(1:10));
            split_into_cells_intra_session(ch, exporting_arr, cutting{1,count});
            marker = i;
            count = count + 1;
            cd(location);
        end
    end
end

function split_into_cells_intra_session(channel, two_layer_chunk, start_ind)

    disp(pwd);
    cd(channel{:});
    disp(pwd);
    disp(size(two_layer_chunk));
    disp(start_ind);
    time_chunk = two_layer_chunk(1,:);
    time_chunk = time_chunk - start_ind + 1; % realigned to start of session
    time_chunk = time_chunk / 30000; % now s
    time_chunk = time_chunk * 1000; % now ms

    unique_cells_here = unique(two_layer_chunk(2,:));
    disp(unique_cells_here);
    assignment_layer = two_layer_chunk(2,:);
    unix('rm -r cell*');
    
    for i = 1:length(unique_cells_here)
        if i < 10
            cell_name = ['cell0', num2str(i)];
        else
            cell_name = ['cell', num2str(i)];
        end
        mkdir(cell_name);
        cd(cell_name);
        disp(pwd);
        strain.timestamps = time_chunk(assignment_layer==unique_cells_here(i));
        strain.components = unique_cells_here(i);
        strain.reference = 'from firings.curated2.mda';
        disp(length(strain.timestamps));
        save('spiketrain.mat', '-struct', 'strain');
        cd('..');    
    end
end


function A=readmda(fname)

    %READMDA - read the contents of a .mda file. MDA stands for
    %multi-dimensional array.
    %
    % See http://magland.github.io//articles/mda-format/
    %
    % Syntax: A=readmda(fname)
    %
    % Inputs:
    %    fname - path to the .mda file
    %
    % Outputs:
    %    A - the multi-dimensional array
    %
    % Other m-files required: none
    %
    % See also: writemda

    % Author: Jeremy Magland
    % Jan 2015; Last revision: 15-Feb-2106

    if (strcmp(fname(end-4:end),'.csv')==1)
        A=textread(fname,'','delimiter',',');
        return;
    end

    F=fopen(fname,'rb');

    try
    code=fread(F,1,'int32');
    catch
        error('Problem reading file: %s',fname);
    end
    if (code>0) 
        num_dims=code;
        code=-1;
    else
        fread(F,1,'int32');
        num_dims=fread(F,1,'int32');    
    end;

    dim_type_str='int32';
    if (num_dims<0)
        num_dims=-num_dims;
        dim_type_str='int64';
    end;

    S=zeros(1,num_dims);
    for j=1:num_dims
        S(j)=fread(F,1,dim_type_str);
    end;
    N=prod(S);

    if num_dims == 1,
      A = zeros(1,S);
    else
      A=zeros(S);
    end

    if (code==-1)
        M=zeros(1,N*2);
        M(:)=fread(F,N*2,'float');
        A(:)=M(1:2:prod(S)*2)+i*M(2:2:prod(S)*2);
    elseif (code==-2)
        A(:)=fread(F,N,'uchar');
    elseif (code==-3)
        A(:)=fread(F,N,'float');
    elseif (code==-4)
        A(:)=fread(F,N,'int16');
    elseif (code==-5)
        A(:)=fread(F,N,'int32');
    elseif (code==-6)
        A(:)=fread(F,N,'uint16');
    elseif (code==-7)
        A(:)=fread(F,N,'double');
    elseif (code==-8)
        A(:)=fread(F,N,'uint32');
    else
        error('Unsupported data type code: %d',code);
    end;

    fclose(F);
    
end







