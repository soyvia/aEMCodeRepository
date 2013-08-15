function sp=spCreateStackParams(stack,pixA)
% Create the sp structure to go with a particle-image stack
[n, ny, nim]=size(stack);
sp=struct;
sp.pixA=pixA;
sp.boxSize=n;
sp.trans=single(zeros(nim,2));
sp.rot=single(zeros(nim,1));
sp.class=single(zeros(nim,1));
sp.thetaPhi=single(zeros(nim,2));
sp.cc=single(zeros(nim,1));
sp.amp=single(zeros(nim,1));
sp.active=true(nim,1);
sp.flip=false(nim,1);

