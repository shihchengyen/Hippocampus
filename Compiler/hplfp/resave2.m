function resave2

rawname = 'rplraw.mat';
hpname = 'rplhighpass.mat';
lfpname = 'rpllfp.mat';
rawname2 = 'rplraw2.mat';
hpname2 = 'rplhighpass2.mat';
lfpname2 = 'rpllfp2.mat';

% check if rplraw exists
if(~isempty(dir(rawname)))
	l = load(rawname);
	oname = fieldnames(l);
	rw = eval(['l.' oname{:}]);
    % rw = rplraw('auto');
    % copy all fields except analogTime
    % rw.data.analogInfo = rw2.data.analogInfo;
    % rw.data.analogData = rw2.data.analogData;
    % rw.data.numSets = rw2.data.numSets;
    % rw.data.Args = rw2.data.Args;
    bmod = 0;
    if(isfield(rw.data,'analogTime'))
        display('rplraw: Removing analogTime')
	    rw.data = rmfield(rw.data,'analogTime');
	    bmod = 1;
	end
	if(isa(rw.data.analogData,'double'))
        display('rplraw: Converting to single')
		rw.data.analogData = single(rw.data.analogData);
	    bmod = 1;
	end
	if(bmod==1)
        display('rplraw: Saving modified object')
	    save(rawname2,'rw');
	end
	clear l oname rw
end

% check if rplhighpass exists
if(~isempty(dir(hpname)))
	l = load(hpname);
	oname = fieldnames(l);
	rh = eval(['l.' oname{:}]);
    % rh = rplhighpass('auto');
    % copy all fields except analogTime
    % rh.data.analogInfo = rh2.data.analogInfo;
    % rh.data.analogData = rh2.data.analogData;
    % rh.data.numSets = rh2.data.numSets;
    % rh.data.Args = rh2.data.Args;
    bmod = 0;
	if(isfield(rh.data,'analogTime'))
        display('rplhighpass: Removing analogTime')
		rh.data = rmfield(rh.data,'analogTime');
		bmod = 1;
	end
	if(isa(rh.data.analogData,'double'))
        display('rplhighpass: Converting to single')
		rh.data.analogData = single(rh.data.analogData);
	    bmod = 1;
	end
	if(bmod==1)	
        display('rplhighpass: Saving modified object')
		save(hpname2,'rh');
	end
	clear l oname rh
end

% check if rpllfp exists
if(~isempty(dir(lfpname)))
	l = load(lfpname);
	oname = fieldnames(l);
	df = eval(['l.' oname{:}]);
    % df = rpllfp('auto');
    % copy all fields except analogTime
    % df.data.analogInfo = df2.data.analogInfo;
    % df.data.analogData = df2.data.analogData;
    % df.data.numSets = df2.data.numSets;
    % df.data.Args = df2.data.Args;
    bmod = 0;
	if(isfield(df.data,'analogTime'))
        display('rpllfp: Removing analogTime')
		df.data = rmfield(df.data,'analogTime');
		bmod = 1;
	end
	if(isa(df.data.analogData,'double'))
        display('rpllfp: Converting to single')
		df.data.analogData = single(df.data.analogData);
	    bmod = 1;
	end
	if(bmod==1)	
        display('rpllfp: Saving modified object')
		save(lfpname2,'df');
	end
	clear l oname df
end
