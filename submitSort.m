function rval = submitSort(varargin)
% submitSort Submits channels for spike sorting
%   This function checks a spreadsheet to submit channels for spike sorting. This function should be called from
%   the session directory, or can be used in the following manner to submit several days of data for spike sorting:
%       ProcessLevel(nptdata,'Levels','Days','Exclude',{'analog'},'nptLevelCmd',{'Session','submitSort'})
%   This function uses the following optional arguments:
%       CellLog	- 'CellActivityLog.xlsx'
%       ChannelsPerArray - 32
%       SortCmd - 'source ~/.bash_profile; cp ~/Dropbox/Work/Matlab/hmmsort/hmmsort5.dag .; condor_submit_dag -maxpre 10
%                  hmmsort5.dag'

Args = struct('CellLog','CellActivityLog.xlsx','ChannelsPerArray',32, ...
    'ChannelRows',3:112, 'HPC', 0, 'UseHPC', 0, 'HPCUser', 'eleys', 'HPCMachine','atlas8', ...
    'HPCDir', '~/hpctmp/Data/', 'HPCCommand','~/hmmsort/hmmsort_pbs.py ~/hmmsort', ...
    'SkipMarker', 0, 'SkipSortName', 'skipsort.txt', 'SelectiveSort',0,...
    'SortCmd','source ~/.bash_profile; cp ~/Dropbox/Work/Matlab/hmmsort/hmmsort5.dag .; condor_submit_dag -maxpre 10 hmmsort5.dag');
Args.flags = {'HPC','SkipMarker','UseHPC','SelectiveSort'};

[Args,modvarargin] = getOptArgs(varargin,Args);

% set default return value
rval = 0;

% check if we want to sort by checking the spreadsheet
if(Args.SelectiveSort)
    % get day of current directory
    [p,daydirstr,e] = fileparts(getDataOrder('day'));
    
    % change to the days directory to read the spreadsheet
    [p,cwd] = getDataOrder('days','CDNow');
    
    % read spreadsheet indicating channels with possible single units
    num = xlsread(Args.CellLog);
    
    % return to original directory
    cd(cwd)
    % get index of current directory
    dayidx = find(num(1,:)==str2num(daydirstr));
    
    if(~isempty(dayidx))
        % find rows with single units
        ai1 =find(num(Args.ChannelRows,dayidx)==1);
        % if there are no single units, sort everything aside from the channels that are broken (marked by 0)
        % if(isempty(ai))
        ain = find(isnan(num(Args.ChannelRows,dayidx)));
        ai = sort([ai1; ain]);
        % end
        % get the channel numbers
        ch_nums = num(Args.ChannelRows(ai),1);
        
        
        if(Args.HPC)
            % assume we are in session directory
            [p1,sesstr] = nptFileParts(pwd);
            % get day string
            % [p2,daydirstr] = nptFileParts(p1);
            % create directory on HPC
            % daydirstr = [Args.HPCDir daystr];
            sesdirstr = [daydirstr '/' sesstr];
            usermachinestr = [Args.HPCUser '@' Args.HPCMachine];
            array_nums = floor((ch_nums-1)/Args.ChannelsPerArray)+1;
            if(Args.SkipMarker)
                % create array directories
                uarnum = unique(array_nums);
                % add array prefix
                arstring = sprintf('array%02d ', uarnum);
                % send ssh command to create array directories
                if ~Args.UseHPC
                    syscmd = ['ssh ' usermachinestr ' "cd ' Args.HPCDir '; mkdir ' daydirstr ' ' sesdirstr '; cd ' sesdirstr '; mkdir ' arstring '"'];
                    system(syscmd);
                end
                
                % find channels marked with 0's indicating that they are not working
                ai0 =find(num(Args.ChannelRows,dayidx)==0);
                skip_ch_nums = num(Args.ChannelRows(ai0),1);
                skip_array_nums = floor((skip_ch_nums-1)/Args.ChannelsPerArray)+1;
                for chidx = 1:size(skip_ch_nums,1)
                    % convert channel number to array number
                    array_dir = sprintf('array%02d', skip_array_nums(chidx));
                    chan_dir = sprintf('channel%03d', skip_ch_nums(chidx));
                    display(['Skip sorting ' array_dir filesep chan_dir])
                    cd(array_dir);
                    cd(chan_dir);
                    system(['touch ' Args.SkipSortName]);
                    cd(cwd)
                end  % for chidx = 1:size(skip_ch_nums,1)
            else  % if(Args.SkipMarker)
                system(['ssh ' usermachinestr ' "mkdir ' daydirstr ' ' sesdirstr '"']);
                fid = fopen('tarlist.txt','w');
                fprintf(fid,'array%02d/channel%03d/rplhighpass.mat\n',([array_nums ch_nums])');
                fclose(fid);
                system('cat tarlist.txt | xargs tar --options gzip:compression-level=9 -cvzf sort.tar');
                system(['scp sort.tar ' usermachinestr ':' sesdirstr]);
                system(['ssh ' usermachinestr ' "cd ' sesdirstr '; tar -xvzf sort.tar; rm sort.tar; ' Args.HPCCommand '"']);
                system('rm tarlist.txt sort.tar');
            end  % if(Args.SkipMarker)
        else  % if(Args.HPC)
            for chidx = 1:size(ch_nums,1)
                chan_num = ch_nums(chidx);
                % convert channel number to array number
                array_num = floor((chan_num-1)/Args.ChannelsPerArray)+1;
                array_dir = sprintf('array%02d',array_num);
                chan_dir = sprintf('channel%03d',chan_num);
                display([array_dir filesep chan_dir])
                cd(array_dir);
                cd(chan_dir);
                system(Args.SortCmd);
                cd(cwd)
            end  % for chidx = 1:size(ch_nums,1)
        end  % if(Args.HPC)
    else  % if(~isempty(dayidx))
        display('Day not found!')
        % set return value to -1 to make sure rsCreateObject does not launch hplfp scripts
        rval = -1;
    end  % if(~isempty(dayidx))
else  % if(Args.SelectiveSort)
    % ignore spreadsheet and sort all channels
    if(Args.HPC) && ~Args.UseHPC
        % get array directories string
        arstring = sprintf('array%02d ', 1:4);
        % assume we are in session directory
        [p1,sesstr] = nptFileParts(pwd);
        % get day string
        [p2,daydirstr] = nptFileParts(p1);
        % get day+session string
        sesdirstr = [daydirstr '/' sesstr];
        % get username and machine
        usermachinestr = [Args.HPCUser '@' Args.HPCMachine];
        % create day, session, and array directories on the HPC
        syscmd = ['ssh ' usermachinestr ' "cd ' Args.HPCDir '; mkdir ' daydirstr ' ' sesdirstr '; cd ' sesdirstr '; mkdir ' arstring '"'];
        system(syscmd);
        % find all rplhighpass files and create a tar archive
        system('find . -name "rplhighpass.mat" | xargs tar --options gzip:compression-level=9 -cvzf sort.tar');
        system(['scp sort.tar ' usermachinestr ':' Args.HPCDir sesdirstr]);
        system(['ssh ' usermachinestr ' "cd ' Args.HPCDir '; cd ' sesdirstr '; tar -xvzf sort.tar; rm sort.tar; ' Args.HPCCommand '"']);
        system('rm sort.tar');
    else  % if(Args.HPC) && ~Args.UseHPC
        % find directories named channel???, and run the SortCmd inside each directory
        system(['cwd=`pwd`; for i in `find . -name "channel???"`; do echo $i; cd $i; ' Args.SortCmd '; cd $cwd; done'])
    end  % if(Args.HPC)
end   % if(Args.SelectiveSort)
