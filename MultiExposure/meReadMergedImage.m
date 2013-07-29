function [m mergeFullPath]=meReadMergedImage(mi,doNorm)
% function [m mergeFullPath]=ReadMergedImage(mi)
% Try to load the merged image m, and normalize it to total dose (if doNorm = 1;
% default is no normalization).
% If the image file isn't found, put up a file selector, and return
% both the path where the image is found.
% m=0 if nothing is found.
if nargin<2
    doNorm=0;
end;
mergeFullPath='';
m=0;
iname=[mi.basePath mi.procPath mi.baseFilename 'm.mrc'];
inamec=[mi.basePath mi.procPath mi.baseFilename 'mc.mrc'];
if FileExists(iname)
    m=ReadEMFile(iname);
elseif FileExists(inamec)
    m=ReadEMFile(inamec);
else
    %%
    [iname mergeFullPath]=uigetfile('*m.mrc','Find the merged image');
    if numel(iname)>1
        m=ReadEMFile([mergeFullPath iname]);
    end;
end;
if doNorm
m=m/mi.doses(1);
end;