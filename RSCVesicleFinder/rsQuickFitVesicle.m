function [miNew, diffIm, vesFit]=rsQuickFitVesicle(img,mask,mi,mi0,vindex,effCTF,hfVar,mode,displayOn)
% function miNew=rsQuickFitVesicle(img,mi,mi0,vindex,effCTF,displayOn)
% Given a subtracted image, restore one vesicle using the model from mi0
% and the vesicle parameters from mi. With one vesicle restored, fit by
% least-squares a vesicle model with fractional shifts. 
% effCTF has zero frequency at the origin.
% A small image of size ndis=size(effCTF) is extracted from img and this is the portion
% that is fitted.  The default for displayOn = 1.
% The starting vesicle info is taken from entry vindex in the mi.vesicle
% arrays, and the refined information is updated into the mi copy, miNew.
% mode=1 fit amplitude only
% mode=2 fit position and amplitude (default)
% mode=3 fit radius too.

ndis=size(effCTF,1);
if nargin<9
    displayOn=1;
end;
if nargin<8
    mode=2;  % don't fit radius
end;

mask=single(mask);

if mode>2
    niters=20;  % iterate fitting of the vesicle radius
else
    niters=2;
end;
sigmaT=2;  % SD of prior for shifts
n=size(img,1);
ds=mi.imageSize(1)/n;

vx=(mi.vesicle.x(vindex))/ds+1;  % assume zero-based coordinates
vy=(mi.vesicle.y(vindex))/ds+1;
vr=mi.vesicle.r(vindex)/ds;
vs=mi.vesicle.s(vindex);

% Get the membrane cross-section density
vd1=mi.vesicleModel;
vd=meDownsampleVesicleModel(vd1,ds)*mi.pixA*ds;

vdOrig=mi0.vesicleModel;
vd0=meDownsampleVesicleModel(vdOrig,ds)*mi0.pixA*ds;

approxLoc=round([vx vy]);
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get a centered, cropped portion of the image.
% imgc0=Crop(circshift(img,-approxShift),ndis);
imgc0=ExtractImage(img,approxLoc,ndis);
mask0=ExtractImage(mask,approxLoc,ndis);
% add back the original vesicle
v=-vs*VesicleFromModel(ndis,vr,vd0,[vx vy]-approxLoc+ndis/2+1);
vfilt=real(ifftn(fftn(v).*effCTF));
% We assume that img has been filtered equivalently to effCTF.  That is, if
% prewhitening has been done to img, it should also be included in effCTF.
imgc=(imgc0+vfilt).*mask0;

if displayOn
    subplot(2,2,1);
    imacs(imgc);
    title(vindex(1));
    
    subplot(2,2,2);
    imacs(imgc0);
    drawnow;
end;

logP=-Radius(ndis).^2./(2*sigmaT^2);  % unscaled log prior
ctr=floor(ndis/2+1);

% initial value for location
const=ones(ndis,ndis);
F(:,2)=const(:);
P=vr;
%
sumT=[ctr ctr];
shifts=zeros(niters+1,4);
if mode>2
    vr=Simplex('init',P);
end;
for iRad=1:niters
    v=-VesicleFromModel(ndis,vr,vd,sumT);
    fvfilt=fftn(v).*effCTF;
    vfilt=real(ifftn(fvfilt)).*mask0;  % f
    %     eRef=vfilt(:)'*vfilt(:);  % approx. power in the reference.
    % Compute the CCF
    cc=fftshift(real(ifftn(fftn(imgc).*conj(fvfilt))));
    logPScaled=logP*hfVar/vs;  % scale up to match ccf
    [mxc xi yi]=max2di(cc+logPScaled);
    sumT=sumT+[xi yi]-ctr;
    %         v=-VesicleFromModel(ndis,vr,vd,sumT);
    %         fvfilt=fftn(v).*effCTF;
    %         vfilt=real(ifftn(fvfilt));
    F(:,1)=vfilt(:);
    warning('off','MATLAB:singularMatrix');
    a=LinLeastSquares(F,imgc(:));
    warning('on','MATLAB:singularMatrix');
    vs=a(1);
    vfit=F*a;
%     shifts(iRad,:)=[vs*1000 sumT vr];

    if mode>2 && iRad<niters
        diff=imgc(:)-vfit;
        err=diff'*diff;
        vr=Simplex(err);
    end;
end;
if mode>2  % Get the radius centroid, then do the linear fit
    vr=Simplex('centroid');
    v=-VesicleFromModel(ndis,vr,vd,sumT);
    vfilt=real(ifftn(fftn(v).*effCTF)).*mask0;
    cc=fftshift(real(ifftn(fftn(imgc).*conj(fftn(vfilt)))));
    logPScaled=logP*hfVar/vs;  % scale up to match ccf
    [mxc xi yi]=max2di(cc+logPScaled);
    sumT=sumT+[xi yi]-ctr;
    F(:,1)=vfilt(:);
    warning('off','MATLAB:singularMatrix');
    a=LinLeastSquares(F,imgc(:));
    warning('on','MATLAB:singularMatrix');
    vfit=F*a;
%     shifts(niters+1,:)=[vs*1000 sumT vr];
end;


diff=imgc(:)-vfit;
vesFit=reshape(vfit,ndis,ndis);
diffIm=reshape(diff,ndis,ndis);
%     subplot(2,2,3);
%     imacs(vesfit);

if displayOn
    %     disp(xShift);
    subplot(2,2,4);
    imacs(diffIm);
    subplot(2,2,3);
    %     imacs(vesFit);
    plot([sect(cc) sect(logPScaled)]);  % plot the CC and the prior
    drawnow;
end;

% disp(shifts)


%% Construct the model image
% v=-VesicleFromModel(ndis,vr,vd,sumT);
% vfilt=real(ifftn(fftn(v).*effCTF));
% F=[vfilt(:) const(:)];
% a=LinLeastSquares(F,imgc(:));
% vfit=F*a;  % Ignore the constant term
% diff=imgc(:)-vfit;
% vesFit=reshape(vfit,ndis,ndis);
% diffIm=reshape(diff,ndis,ndis);  % residual image

% update the info structure
miNew=mi;
sh=approxLoc-ndis/2-1;
miNew.vesicle.r(vindex)=vr*ds;  % radius isn't changed
miNew.vesicle.x(vindex)=(sumT(1)+sh(1)-1)*ds;  % zero-based position in image
miNew.vesicle.y(vindex)=(sumT(2)+sh(2)-1)*ds;
miNew.vesicle.s(vindex)=a(1);



% disp([P(1) P(2) P(3) a']);  % show the main fit

end
