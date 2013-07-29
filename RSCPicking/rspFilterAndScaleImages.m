function imgs=rspFilterAndScaleImages(mi,dis,m0,mVes)
% Load the first two images of the imgs stack.

fc=mi.pixA*dis.ds./dis.filter; % Convert from A to A^-1
fc(dis.filter==0)=0;           % 0 A -> 0 A^-1
imgs=single(zeros(dis.ndis,dis.ndis,2));
mf=GaussFilt(GaussHP(m0,fc(1)),fc(2));
mVesf=GaussFilt(GaussHP(mVes,fc(1)),fc(2));
imgs(:,:,1)=imscale(mf,256,dis.contrast);
imgs(:,:,2)=imscale(mf-mVesf,256,dis.contrast);
