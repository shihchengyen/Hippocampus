function resave

rawname = 'rplraw.mat';
hpname = 'rplhighpass.mat';
lfpname = 'rpllfp.mat';

% check if rplraw exists
if(~isempty(dir(rawname)))
    rw = rplraw('auto');
    % convert to single precision float
    rw.data.analogData = single(rw.data.analogData);
    rw.data.analogTime = single(rw.data.analogTime);
    % resave in single precision float
    save(rawname,'rw');
end

% check if rplhighpass exists
if(~isempty(dir(hpname)))
    rh = rplhighpass('auto');
    % convert to single precision float
    rh.data.analogData = single(rh.data.analogData);
    rh.data.analogTime = single(rh.data.analogTime);
    % resave in single precision float
    save(hpname,'rh');
end

% check if rpllfp exists
if(~isempty(dir(lfpname)))
    df = rpllfp('auto');
    % convert to single precision float
    df.data.analogData = single(df.data.analogData);
    df.data.analogTime = single(df.data.analogTime);
    % resave in single precision float
    save(lfpname,'df');
end
