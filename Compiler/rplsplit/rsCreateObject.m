function rsCreateObject()

if ~isdeployed
    addpath('~/matlab/DPV')
    addpath('~/matlab/newNpt')
    addpath('~/matlab/Hippocampus')
    addpath('~/matlab/neuroshare')
    addpath('~/hmmsort')
end

% load arguments
% we are doing this so this program can run in compiled form
load('rsData');
% delete rsData.mat
varargin = modvarargin;

if(~Args.SkipSplit)
    % these are object specific fields
    rawfname = dfile(1).name;
    fprintf('Splitting %s\n',rawfname);
    % open the file, and read the information
    [ns_status, hFile] = ns_OpenFile(dfile(1).name);
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    % get number of EntityCount
    nec = nsFileInfo.EntityCount;
    % go through and create the appropriate subdirectory for entity
    for ni = 1:nec
        % get info on the type of entity
        [ns_status, nsEI] = ns_GetEntityInfo(hFile, ni);
        numSamples = nsEI.ItemCount;
        switch(nsEI.EntityType)
            case 1
                % parallel input
                if(~Args.SkipParallel)
                    % parallel input
                    % Get events and time stamps
                    tData.markers = NaN(1, numSamples);
                    tData.timeStamps = NaN(1, numSamples);
                    for i = 1:numSamples
                        [~, tData.timeStamps(i), tData.markers(i)] = ns_GetEventData(hFile, ni, i);
                    end
                    % create and save obj
                    rplparallel('auto','Data',tData,'save',varargin{:});
                    clear tData
                end
            case 2
                b_lfp = 0;
                b_raw = 0;
                % get label
                eLabel = nsEI.EntityLabel;
                % check if it is analog data
                if(~isempty(strfind(eLabel,'analog')))
                    if(~Args.SkipAnalog)
                        chan_num = sscanf(eLabel,'analog %d');
                        % read data
                        [ns_RESULT, tData.analogInfo] = ns_GetAnalogInfo(hFile, ni);
                        [ns_RESULT, ~, tData.analogData] = ns_GetAnalogData(hFile, ni, 1, numSamples);
                        cwd = pwd;
                        nptMkDir('analog');
                        cd('analog');
                        chan_dir = sprintf('channel%02d',chan_num);
                        display(['analog' filesep chan_dir])
                        nptMkDir(chan_dir);
                        cd(chan_dir);
                        tData.analogInfo.NumberSamples = numSamples;
                        rplraw('auto','Data',tData,'save',varargin{:});
                        
                        if ~Args.UseHPC
                            submitSort('HPC','SkipMarker') % do scp to transfer the files to HPC
                        end
                        
                        system('transfersession.sh');
                        fid = fopen('transferred.txt','w');
                        fclose(fid);
                        % submitJob(Args); % submit job onto PBS queue
                        
                        cd(cwd)
                        clear tData
                    end
                    % check if it is raw data
                elseif(~isempty(strfind(eLabel,'raw')))
                    chan_num = sscanf(eLabel,'raw %d');
                    b_raw = 1;
                    b_lfp = 0;
                elseif(~isempty(strfind(eLabel,'lfp')))
                    chan_num = sscanf(eLabel,'lfp %d');
                    b_lfp = 1;
                    b_raw = 0;
                end
                % check if we should process this channel
                chanArgs = isempty(Args.Channels);
                if( chanArgs | (~chanArgs && ~isempty(find(Args.Channels==chan_num)) ) )
                    if( (b_raw * ~Args.SkipRaw) | (b_lfp * ~Args.SkipLFP) )
                        % entity is raw data, so create a channel directory
                        % read data
                        [ns_RESULT, tData.analogInfo] = ns_GetAnalogInfo(hFile, ni);
                        [ns_RESULT, ~, tData.analogData] = ns_GetAnalogData(hFile, ni, 1, numSamples);
                        % check how channels are arranged with respect to arrays
                        array_num = floor((chan_num-1)/Args.ChannelsPerArray)+1;
                        array_dir = sprintf('array%02d',array_num);
                        cwd = pwd;
                        nptMkDir(array_dir);
                        cd(array_dir);
                        %
                        % restart channel numbers from 001 for each array
                        % chan_dir = sprintf('channel%03d',rem(chan_num-1,Args.ChannelsPerArray)+1);
                        %
                        % change to keep original channel numbers so we don't
                        % have to switch between channel numbers
                        chan_dir = sprintf('channel%03d',chan_num);
                        display([array_dir filesep chan_dir])
                        nptMkDir(chan_dir);
                        cd(chan_dir);
                        tData.analogInfo.NumberSamples = numSamples;
                        % create and save obj
                        if(b_raw)
                            rplraw('auto','Data',tData,'save',varargin{:});
                        elseif(b_lfp)
                            rpllfp('auto','Data',tData,'save',varargin{:});
                        end
                        
                        if ~Args.UseHPC
                            submitSort('HPC','SkipMarker') % do scp to transfer the files to HPC
                        end
                        
                        submitJob(Args); % submit job onto PBS queue
                        
                        cd(cwd);
                        clear tData
                    end % if( (b_raw * ~Args.SkipRaw) | (b_lfp * ~Args.SkipLFP) )
                end % if( chanArgs && (~chanArgs && ~isempty(find(Args.Channels==chan_num)) )
        end
    end
    ns_status = ns_CloseFile(hFile);
end  % if(~Args.SkipSplit)


display('Done!')

end

function [] = submitJob(Args)
if ~Args.UseHPC % swap between the HPC and HTCondor
    cmdPath = ['condor_submit ',fullfile('~','cbin')];
    cmdScript = '';
else
    cmdPath = ['qsub ',fullfile('$GITHUB_MATLAB','Hippocampus','Compiler','hplfp')];
    cmdScript = 'HPC';
end
disp('Launching eyehplfp scripts...')
cmdSubmit = [cmdPath, filesep, 'eyehplfp', cmdScript, '_submit_file.txt'];
system(['source ',fullfile('~','.bash_profile'),'; source /etc/profile.d/rec_modules.sh; module load pbs; ', cmdSubmit]);
end
