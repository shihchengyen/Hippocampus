function IsDataCompressed

%ProcessLevel(nptdata,'Levels','Days','nptLevelCmd',{'Channel','IsDataCompressed'});
rawname = 'rplraw.mat';
hpname = 'rplhighpass.mat';
lfpname = 'rpllfp.mat';

saveFile = fopen('~/Documents/dataCheckList.txt','a+');
printDir = true;

% check if rplraw exists
if(~isempty(dir(rawname)))
    try 
        l = load(rawname);
        oname = fieldnames(l);
        rw = eval(['l.' oname{:}]);
    catch
       if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
       end
        fprintf(saveFile,'rplraw: error using load\n');
        return;
    end
    
    if(isfield(rw.data,'analogTime'))
       % analogTime is not removed
       if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
       end
       fprintf(saveFile,'rplraw: AnalogTime not deleted\n');
    end
	if(~isa(rw.data.analogData,'single'))
        % data is not single precision
        if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
        end
        fprintf(saveFile,'rplraw: Precision not single\n');
    end
	clear l oname rw
end

% check if rplhighpass exists
if(~isempty(dir(hpname)))
    try 
        l = load(hpname);
        oname = fieldnames(l);
        rh = eval(['l.' oname{:}]);
    catch
        if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
        end
        fprintf(saveFile,'rplhighpass: error using load\n');
        return;
    end
    
    if(isfield(rh.data,'analogTime'))
       % analogTime is not removed
       if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
       end
       fprintf(saveFile,'rplhighpass: AnalogTime not deleted\n');
    end
	if(~isa(rh.data.analogData,'single'))
        % data is not single precision
       if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
       end
       fprintf(saveFile,'rplhighpass: Precision not single\n');
    end
    
	clear l oname rh
end

% check if rpllfp exists
if(~isempty(dir(lfpname)))
    try 
        l = load(lfpname);
        oname = fieldnames(l);
        df = eval(['l.' oname{:}]);
    catch
       if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
       end
       fprintf(saveFile,'rpllfp: error using load\n');
        return;
    end
    
    if(isfield(df.data,'analogTime'))
       % analogTime is not removed
       if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
       end
       fprintf(saveFile,'rpllfp: AnalogTime not deleted\n');
    end
	if(~isa(df.data.analogData,'single'))
        % data is not single precision
       if (printDir) 
           fprintf(saveFile,'%s\n',pwd);
           printDir = false;
       end
       fprintf(saveFile,'rpllfp: Precision not single\n');
    end
    
	clear l oname df
end

fclose(saveFile);


