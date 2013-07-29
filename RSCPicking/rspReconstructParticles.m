function img=rspReconstructParticles(dis,mi,picks,ptrs,rscc)
% For SimpleRSPicker, create an image of autopicked particles, based on the
% best-matching template and the particle amplitudes.

img=single(zeros(dis.ndis));
if ~isfield(rscc,'eigenImgs') || ~isfield(rscc,'vList') % data not there
    return
end;
apFlags=16:47;  % Flag for picked particles
apPtrs=[2 3 4 6];    % pointers for picked particles
n=size(rscc.mxCC,1);
ds1=mi.imageSize(1)/n;
[ne, ne1, nterms]=size(rscc.eigenImgs);
vsz=numel(rscc.vList)/nterms;
eigs=reshape(rscc.eigenImgs,ne*ne,nterms);
vlst=reshape(rscc.vList,nterms,vsz);

for ptrIndex=apPtrs
    for ip=1:ptrs(ptrIndex)  % loop over all auto-picked particles
        c=squeeze(picks(ptrIndex,ip,:))';
        pos=c(1:2)/ds1+1;
        flag=c(3);
        vesInd=c(4);
        amp=-c(5);  % invert the contrast
        templ=c(6);
        if any(flag==apFlags) && templ>0 && vesInd>0
            s=mi.vesicle.s(vesInd);
            modelParticle=reshape(eigs*vlst(:,templ),ne,ne);
            img=img+amp*s*ExtractImage(modelParticle,round(pos),n,1);
        end;
    end;
end;
img=imscale(img,256);
% dis.mode=3;
% imgs(:,:,3)=img;
%             rspUpdateDisplay(mi,dis,imgs,masks,picks,ptrs);

