function [obj, varargout] = rplsplit(varargin)
%@rplsplit Constructor function for rplsplit class
%   OBJ = rplsplit(varargin) extracts LFPs from a RIPPLE recording
%
%   OBJ = rplsplit('auto') attempts to create a rplsplit object by ...
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Instructions on rplsplit %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%example [as, Args] = rplsplit('save','redo')
%
%dependencies: 

Args = struct('RedoLevels',0, 'SaveLevels',0, 'Auto',0, 'ArgsOnly',0, ...
	'SkipRaw',0, 'SkipLFP',0, 'SkipParallel',0, 'SkipAnalog',0, 'SkipSort', 0, ...
	'Channels',[], 'ChannelsPerArray',32, 'HPCCmd','', 'SkipSplit',0, ...
    'HPCInputFilename','rsData.mat','UseHPC',0, 'SelectiveSort', 0,'skipCheckingRplsplit',0, ...
    'SkipOSort',0);
Args.flags = {'Auto','ArgsOnly','SkipRaw','SkipLFP','SkipParallel', 'SkipSort', ...
	'SkipAnalog','SkipSplit','UseHPC','SelectiveSort','skipCheckingRplsplit', 'SkipOSort'};
% The arguments which can be neglected during arguments checking
Args.DataCheckArgs = {'Channels'};                            

[Args,modvarargin] = getOptArgs(varargin,Args, ...
	'subtract',{'RedoLevels','SaveLevels'}, ...
	'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
	'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'rplsplit';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'rs';

% To decide the method to create or load the object
if Args.skipCheckingRplsplit && exist(fullfile(pwd,Args.matname),'file')
    delete(Args.matname)
end
command = checkObjCreate('ArgsC',Args,'narginC',nargin,'firstVarargin',varargin,'saverplsplit',0);

if(strcmp(command,'createEmptyObjArgs'))
    varargout{1} = {'Args',Args};
    obj = createEmptyObject(Args);
elseif(strcmp(command,'createEmptyObj'))
    obj = createEmptyObject(Args);
elseif(strcmp(command,'passedObj'))
    obj = varargin{1};
elseif(strcmp(command,'loadObj'))
    l = load(Args.matname);
    obj = eval(['l.' Args.matvarname]);
elseif(strcmp(command,'createObj'))
    % IMPORTANT NOTICE!!! 
    % If there is additional requirements for creating the object, add
    % whatever needed here
	dfile = dir('*.ns5');
	dnum = size(dfile,1);

	% check if the right conditions were met to create object
	if(dnum>0)
		% this is a valid object
		if(dnum>1)
			% there is more than 1 ns5 file in the current directory
			% print a warning that only the 1st file will be used to create the object
			warning('More than 1 ns5 file found')
			warning(['Creating object only from ' dfile(1).name])
		end
	
		% save Args and modvarargin so that compiled program can retrieve them
		% look for ripple file
		save(Args.HPCInputFilename,'Args','modvarargin','dfile');
		if(~isempty(Args.HPCCmd))
			% use HTCondor to run rplsplitcreateObject
			[s,w] = system(Args.HPCCmd);
			% make sure shell script ran without problems
			if(s~=0)
				error([Args.classname ': Error splitting files!'])
            else
                % display output
                fprintf('%s\n',w);

                % create plot unity maze job
				cwd = pwd;
                if(~isempty(strfind(cwd,'session01')))
                    fid = fopen('unitymaze_job0000.pbs','w');
                    fprintf(fid,...
                        ['#!/bin/bash\n'...
                        '#PBS -q serial\n'...
                        '#PBS -W depend=afterok:',w,'\n'...
                        'cd "$PBS_O_WORKDIR"\n'...
                        'cwd=$PWD\n'...
                        'matlab2016b2 -nosplash -r "figure; um=plot(unitymaze(''auto''),1); saveas(gcf,''um.png''); close; exit"']);
                    fclose(fid);
                    
                    % submit transfer job
                    [s,w] = system('source ~/.bash_profile; source /etc/profile.d/rec_modules.sh; module load pbs; qsub unitymaze_job0000.pbs');
                end  % if(~isempty(strfind(cwd,'session01')))
                
                % transfer raw data, rplparallel.mat, um.fig, um.png
                fid = fopen('transferRaw_job0000.pbs','w');
                fprintf(fid,...
                    ['#!/bin/bash\n'...
                    '#PBS -q serial\n'...
                    '#PBS -W depend=afterok:',w,'\n'...
                    'cd "$PBS_O_WORKDIR"\n'...
                    'cwd=$PWD\n'...
                    'b="2018"\n'...
                    'index="${cwd%%$b*}"\n'...
                    'channelStr=${cwd:${#index}}\n'...
                    'targetDir=''/volume1/Hippocampus/Data/picasso/''\n'...
                    'targetDir+=$channelStr\n'...
                    'ssh -p 8398 hippocampus@cortex.nus.edu.sg mkdir -p $targetDir\n'...
                    'scp -P 8398 -rv {RawData*,18*,unity*,rplparallel.mat,um*} hippocampus@cortex.nus.edu.sg:$targetDir &&\n'...
                    'rm -rv {RawData*,18*,unity*,rplparallel.mat,um*} && \n'...
                    'touch transferredRaw.txt']);
                fclose(fid);
                
                % submit transfer job
                system('source ~/.bash_profile; source /etc/profile.d/rec_modules.sh; module load pbs; qsub transferRaw_job0000.pbs');

			end
		else
			rsCreateObject();
		end
	
		% create nptdata so we can inherit from it    
		data.Args = Args;
		data.rawfname = dfile(1).name;
		data.numSets = 1;
		n = nptdata(data.numSets,0,pwd);
		d.data = data;
		obj = class(d,Args.classname,n);
		saveObject(obj,'ArgsC',Args);	
		
	else  % if(dnum>0)
		% create empty object
		obj = createEmptyObject(Args);
	end
end

function obj = createEmptyObject(Args)

% create nptdata so we can inherit from it
data.Args = Args;
data.rawfname = '';
data.numSets = 0;
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
