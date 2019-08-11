function rerpll
% Recreates rplparallel object by calling arrangeMarkers to re-arrange the markers
% This function loads an existing rplparallel object, renames the old object to
% rplparallel0.mat, adds a SampleRate field, and then instantiates a new rplparallel
% object that re-arranges the markers before saving the new rplparallel object.
%
% The function can be used in the following way:
%    ProcessLevel(nptdata,'Include',{'20171123','20171124','20171127','20171130'}, ...
%        'Levels','Days','nptLevelCmd',{'Session','rerpll'})

rp = rplparallel('auto');
!mv rplparallel.mat rplparallel0.mat
rp.data.SampleRate = 30000;
df = rplparallel('Data',rp.data);
save rplparallel.mat df
