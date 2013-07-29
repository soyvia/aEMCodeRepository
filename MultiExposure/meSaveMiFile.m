function fullName=meSaveMiFile(mi)
% function fullName=meSaveMiFile(mi)
fullName=[mi.basePath mi.infoPath mi.baseFilename 'mi.mat'];
% disp(['saving:  ' fname]);
save(fullName,'mi');
