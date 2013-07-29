function mi1=rsRefineVesicleFitsSub(mi,m,fittingMode)
% modes are:  1 - fit overall amplitude
%             2 - fit amplitude and shift only, don't change radius
%             3 - fit amplitude, shift and radius
if nargin<3
    fittingMode=3;
end;
maxVesiclesToFit=inf;
useOkField=1;      % refine every vesicle for which ok is true.
doDownsampling=1;  % Downsample for speed
disA=800;          % size of displayed/fitted window in angstroms
sRange=2;          % Outlier amplitudes must lie within this factor of the
% median of vesicle.s, otherwise set equal to median.
% Set to 1 to force median for all vesicles.
displayPeriod=5;
%       faster than RadialPowerSpectrum:
n=size(m);
annulus=fuzzymask(n,2,0.225*n,.05*n)-fuzzymask(n,2,0.15*n,.05*n);
spc=annulus.*fftshift(abs(fftn(m)).^2)/(n(1)*n(2));
hfVar0=sum(spc(:))/sum(annulus(:));

% Handle cases where the ok field is old style
if ~isfield(mi.vesicle,'ok') || numel(mi.vesicle.ok)<numel(mi.vesicle.x) % no ok field at all
    mi.vesicle.ok=true(numel(mi.vesicle.x),1);
end;
[nv, ne]=size(mi.vesicle.ok);
if ne<4
    disp(['Expanding the vesicle.ok array from ' num2str([nv ne])]);
    for i=ne+1:4
    mi.vesicle.ok(:,i)=mi.vesicle.ok(:,ne);
    end;
end;
mi1=mi;

% Get image and pixel sizes
n=size(m,1);
ds0=mi.imageSize(1)/n;  % downsampling factor of m
pixA0=mi.pixA*ds0;    % pixel size of m
if doDownsampling
    % downsample the merged image to about 10A per pixel, yielding the image ms
    targetPixA=12;  % maximum pixel size
    ns=NextNiceNumber(n*pixA0/targetPixA,5,4);  % multiple of 4, max factor 5.
    if ns<n
        disp(['Downsampling to ' num2str(ns) ' pixels.']);
        ms=Downsample(m,ns);
    else
        ns=n;
        ms=m;
    end;
    ds=ds0*n/ns;  % downsampling factor of ms relative to original images.
    pixA=ds*mi.pixA;  % pixA in the image ms.
else  % use the original merged image scale
    ds=ds0;
    pixA=pixA0;
    ns=n;
    ms=m;
end;
hfVar=hfVar0*(ds0/ds)^2;  % hf spectral density after downsampling
ndis=NextNiceNumber(disA/pixA);  % size of display/fitting image


%%  Get the original subtraction, and modify the amplitudes if necessary.\
% subtract every vesicle that was originally fit
mi1.vesicle.s(isnan(mi.vesicle.s))=0;
sMedian=median(mi1.vesicle.s);
%   Force ridiculous amplitude values to the median.
q=mi1.vesicle.s<sMedian/sRange | mi1.vesicle.s > sMedian*sRange;
mi1.vesicle.s(q)=sMedian;
%   Compute the basic (old) subtraction using mi1 with new model
vs=meMakeModelVesicles(mi1,ns,find(mi.vesicle.ok(:,1)));
%   Use linear least-squares to fine-tune the scaling
msAmpScale=(vs(:)'*ms(:))/(vs(:)'*vs(:))
if msAmpScale>1e-3 % don't allow ridiculous values
    vs=vs*msAmpScale;
    mi1.vesicle.s=mi1.vesicle.s*msAmpScale;
end;
msub=ms-vs;  % This is the subtracted image we'll use

%         Whiten and mask the subtracted image
mmask=meGetMask(mi1,ns);
pwH=meGetNoiseWhiteningFilter(mi1,ns);
msubf=mmask.*real(ifftn(fftn(msub).*ifftshift(pwH)));
figure(1);
SetGrayscale;
subplot(2,3,1);
imacs(msubf);
title(['Old subtraction, scaling ' num2str(msAmpScale)]);

%%
switch fittingMode
    %             case 1 hasn't been tested yet-------------
    case 1  % normalize the amplitude of the whole image only
        %                 vs1=meMakeModelVesicles(mi1,ns).*mmask;  % new subtraction
        %                 ms1AmpScale=(vs1(:)'*ms(:))/(vs1(:)'*vs1(:))
        %                 if ms1AmpScale>1e-3 % don't allow ridiculous values
        %                     vs1=vs1*ms1AmpScale;
        %                 end;
        %                 msub1=ms-vs1;
        subplot(2,3,2);
        imacs(msub);
        title('New subtraction');
        %                 mi1.vesicle.s=mi1.vesicle.s*ms1AmpScale;
        q=isnan(mi.vesicle.s);
        mi1.vesicle.s(q)=0;
        mi1.vesicle.ok(:,3)=mi1.vesicle.ok(:,1) & ~q;
        
    case {2 3}  % tune up the x, y and s for each vesicle.
        %
        %                 if ~isfield(mi.vesicle,'ok')
        %                     mi1.vesicle.ok=ones(size(mi1.vesicle.x,3));
        %                 end;
        %%  Actual fitting is done here
        vfit=zeros(ns,ns);
        figure(1)
        disp('Fine fitting translation and amplitude.');
        drawnow;
        nVesicles=sum(mi1.vesicle.ok(:,1)>0);
        nVesicles=min(nVesicles,maxVesiclesToFit);
        nVesicles
        effCTF=ifftshift(meGetEffectiveCTF(mi1,ndis,ds));
        pwFilter=ifftshift(meGetNoiseWhiteningFilter(mi1,ndis,ds));
        figure(1);
        mi1.vesicle.ok(:,3)=true;  % we'll mark unfitted vesicles here.
        for ind=1:nVesicles
            ok=mi1.vesicle.ok(ind,1)>0;  % The vesicle exists
            if ~useOkField || ok
                doDisplay=mod(ind,displayPeriod)==0;
                [mi1, diffIm, vesFit]=rsQuickFitVesicle(msubf,mmask,mi1,mi1,...
                    ind,effCTF.*pwFilter,hfVar,fittingMode,doDisplay);
                
                %                         subplot(2,2,2);
                %                         title(num2str([ind mi1.vesicle.s(ind) mi.vesicle.ok(ind,2)]));
                %                         pause;
                %                         title(fname{fileIndex},'interpreter','none');
            else
                mi1.vesicle.s(ind)=NaN;
            end;
        end;
        q=isnan(mi1.vesicle.s);  % blank the unfittable vesicles.
        mi1.vesicle.s(q)=0;
        mi1.vesicle.ok(:,3)=~q;  % unrefinable vesicles are marked 0
        %                 medianS=median(mi1.vesicle.s);
        %                 good=(mi1.vesicle.s>medianS/(1+sAcceptSd)) & (mi1.vesicle.s<(1+sAcceptSd)*medianS)...
        %                     & (mi1.vesicle.r < maxR) & (mi1.vesicle.r > minR);
        %              oldOk=mi.vesicle.ok;
        %                 mi1.vesicle.ok=min(mi1.vesicle.ok,1+good);
        numberRefined=sum(mi1.vesicle.ok(:,3))
        numberGood=sum(all(mi1.vesicle.ok(:,1:3),2))  %    exists, in range, refined
        mi1.vesicle.refined=1;
        
end; % switch
