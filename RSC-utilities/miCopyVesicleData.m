% miCopyVesicleData
% Copy the vesicleModel and vesicle fields of one set of mi files to
% another.  This is useful in cases where images have been re-merged e.g.
% to change the downsampling factor of the merged images.  In this case the
% vesicle parameters can remain the same and don't need re-fitting.

[fname, infoPath]=uigetfile('*mi.mat','Select source mi files','multiselect','on');
if isnumeric(fname)  % Cancel
    return
end;
if ~iscell(fname)
    fname={fname};
end;
cd(infoPath);

% Put up a second file selector
[fname2, infoPath2]=uigetfile('*mi.mat','Select target mi files','multiselect','on');
if isnumeric(fname2)  % Cancel
    return
end;
%%

for i=1:numel(fname);
    miName=fname{i};
    mi2Name=[infoPath2 fname{i}];
    if exist(miName,'file') && exist(mi2Name,'file')
        mi1=load(miName);
        mi1=mi1.mi;
        mi2=load(mi2Name);
        mi2=mi2.mi;
        
        mi2.vesicleModel=mi1.vesicleModel;
        mi2.vesicle=mi1.vesicle;
        
        mi=mi2;
        save(mi2Name,'mi');
        disp(['Updated: ' mi2Name]);
    end;
 end;
        