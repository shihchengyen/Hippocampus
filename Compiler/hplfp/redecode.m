function redecode(varargin)
	% check to see if we should sort this channel
	if(isempty(dir('skipsort.txt')))
		if(size(h5read('hmmsort/spike_templates.hdf5','/after_combine/p'),1)>size(h5read('hmmsort/spike_templates.hdf5','/after_noise/p'),1))
			display('Launching spike sorting ...')
			% get channel string
			[p1, chstr] = nptFileParts(pwd);
			% get array string
			[p2, arrstr] = nptFileParts(p1);
			% get session string
			[p3, sesstr] = nptFileParts(p2);
			% get day string
			[p4, daystr] = nptFileParts(p3);
	
			% make channel direcory on HPC, copy to HPC, cd to channel directory, and then run hmmsort
			display('Creating channel directory ...')
			syscmd = ['ssh eleys@atlas7 "cd ~/hpctmp/Data; mkdir ' daystr '; cd ' daystr '; mkdir ' sesstr '; cd ' sesstr '; mkdir ' arrstr '; cd ' arrstr '; mkdir ' chstr '"'];
			display(syscmd)
			system(syscmd);
			display('Transferring rplhighpass file ...')
			syscmd = ['scp -r rplhighpass.mat hmmsort eleys@atlas7:~/hpctmp/Data/' daystr '/' sesstr '/' arrstr '/' chstr];
			display(syscmd)
			rval=1;
			while(rval~=0)
				rval=system(syscmd);
			end
			display('Running spike sorting ...')
			syscmd = ['ssh eleys@atlas7 "cd ~/hpctmp/Data/' daystr '/' sesstr '/' arrstr '/' chstr '; cp ~/hmmsort/decode_job0000.pbs .; qsub ./decode_job0000.pbs"'];
			display(syscmd)
			system(syscmd);
		end
	end