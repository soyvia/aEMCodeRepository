function eigenSet=rsMakeEigenrefs(templates,nterms)
% function eigenSet=rsMakeEigenrefs(templates,nterms)
% Given a set of templates (nt x nt x nim) where nim is the product of the
% three higher dimensions nHemi, nGamma and nAngs, create an eigenimage
% expansion using nterms.  A typical value for nterms is 20.
% The result is the structure eigenSet with the fields...
% eigenSet.imgs: Eigenimages, nt x nt x nterms, each with unity power.
% eigenSet.vList: Expansion of each template in factor space, nterms x nim.
% eigenSet.vListNorm: the same, but each vector (column) is normalized.
% eigenSet.ampList:  rms amplitude of each vList vector, nim x 1.
%   So that vListNorm = vList / ampList.
% eigenSet.termVar:  Shows the power as a function of the number of terms
%   by computing the cumulative sum of the squared amplitudes, normalized
%   by the original template variances.  nterms x nim

[nt nt1 nHemi nGamma nAngs]=size(templates);  % nHemi should be 2.
nn=nt^2;  % number of pixels
nim=nHemi*nGamma*nAngs;

templateVector=double(reshape(templates,nt^2,nim));
tplVar=sum(templateVector.^2)';

disp('starting svd..');
tic
% [u s v]=svd(p3,0); nn=nt^2;  % full svd
[u s v]=svds(templateVector,nterms);  % takes 5x less time for 80 than for all 4096.
disp('...done.');
toc
% u is the np x nn matrix of eigenimages (each np pixels image is a column)
% s is the nn x nn diagonal matrix of singular values
% v is the ntot x nn matrix of coefficients

% scale the projection vector
vp=v*s';  % Get the absolute scaling for each projection. u*vp gives the templates back.
vpvar=sum(vp'.^2)';  % Get the sum of each row squared (=var of each projection).
projamps=sqrt(vpvar);
vpn=vp./repmat(projamps,1,nterms);  % Normalized projection coordinates.
% the sum of squares of each row of vpn is 1

% Returned quantities
% Eigenimages, nt x nt x nterms, each with unity power.
eigenSet.imgs=single(reshape(single(u),nt,nt,nterms));
% Expansion of each template in factor space, nterms x nim
eigenSet.vList=reshape(vp',nterms,nHemi,nGamma,nAngs);
% Expansion of each template, but each vector (column) is normalized.
eigenSet.vListNorm=reshape(vpn',nterms,nHemi,nGamma,nAngs); % nterms x nim
% rms amplitude of each vList vector.  vListNorm = vList / ampList.
if numel(projamps)>1
    eigenSet.ampList=reshape(projamps,nHemi,nGamma,nAngs);  % nim elements
else
    eigenSet.ampList=projamps*ones(nHemi,nGamma,nAngs);
end;
% Show the power as a function of the number of terms by computing the
% cumulative sum of the squared amplitudes, normalized by the original
% template variances.
eigenSet.termVar=(cumsum(vp'.^2)'./repmat(tplVar,1,nterms))';  % nterms x nim
