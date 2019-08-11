function resave3
%ProcessLevel(nptdata,'Levels','Day','nptLevelCmd',{'Channel','resave2'})
rawname = 'rplraw.mat';
hpname = 'rplhighpass.mat';
lfpname = 'rpllfp.mat';
rawname2 = 'rplraw2.mat';
hpname2 = 'rplhighpass2.mat';
lfpname2 = 'rpllfp2.mat';

% check if rplraw exists
if(~isempty(dir(rawname)))
    try
        l = load(rawname);
        oname = fieldnames(l);
        rwOriginal = eval(['l.' oname{:}]);
        rw = eval(['l.' oname{:}]);
    catch
        return;
    end
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
        
        %Compare original data to saved modified data
        l = load(rawname2);
        oname = fieldnames(l);
        rw2 = eval(['l.' oname{:}]);
        if (isequal(rwOriginal.data.analogInfo, rw2.data.analogInfo)) 
            display('rplraw2: Data match')
        else
            display('rplraw2: Data mismatch')
            return;
        end
	end
	clear l oname rw
end

% check if rplhighpass exists
if(~isempty(dir(hpname)))
    try
        l = load(hpname);
        oname = fieldnames(l);
        rhOriginal = eval(['l.' oname{:}]);
        rh = eval(['l.' oname{:}]);
    catch
        return;
    end
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
        
        l = load(hpname2);
        oname = fieldnames(l);
        rh2 = eval(['l.' oname{:}]);
        if (isequal(rhOriginal.data.analogInfo, rh2.data.analogInfo)) 
            display('rplhighpass2: Data match')
        else
            display('rplhighpass2: Data mismatch')
            return;
        end
	end
	clear l oname rh
end

% check if rpllfp exists
if(~isempty(dir(lfpname)))
    try
        l = load(lfpname);
        oname = fieldnames(l);
        dfOriginal = eval(['l.' oname{:}]);
        df = eval(['l.' oname{:}]);
    catch
        return;
    end
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
        
        l = load(lfpname2);
        oname = fieldnames(l);
        df2 = eval(['l.' oname{:}]);
        if (isequal(dfOriginal.data.analogInfo, df2.data.analogInfo)) 
            display('rpllfp2: Data match')
        else
            display('rpllfp2: Data mismatch')
            return;
        end
	end
	clear l oname df
end

saveFile = fopen('size.txt','w');

if(~isempty(dir(rawname2)))
    filename1=dir(rawname2);
    filesize1=filename1.bytes;
    fprintf(saveFile,'rplraw2.mat: %d\n',filesize1);
end

if(~isempty(dir(hpname2)))
    filename2=dir(hpname2);
    filesize2=filename2.bytes;
    fprintf(saveFile,'rplhighpass2.mat: %d\n',filesize2);
end

if(~isempty(dir(lfpname2)))
    filename3=dir(lfpname2);
    filesize3=filename3.bytes;
    fprintf(saveFile,'rpllfp2.mat: %d\n',filesize3);
end

fclose(saveFile);
