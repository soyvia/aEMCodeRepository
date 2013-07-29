% meInverseFilterAuto
% Given a good subtraction of vesicles, we search for regions of the
% micrograph with low variance, and mark them as background boxes.
% Then we compute an inverse filter.

writeOut=0;  % write the pre-whitened image.
f0=.002; % A^1  Lower frequency cutoff
f1=.06; % A^1  Upper
nb=64;  % box size
nBoxes=20;

% Select one or mi files
[fname, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);

for fileIndex=1:numel(fname)
    disp(['Reading ' fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);
    iname=[mi.procPath mi.baseFilename 'm.mrc'];
    inamec=[mi.procPath mi.baseFilename 'mc.mrc'];
    if FileExists(iname)
        disp(['        ' iname]);
        m=ReadEMFile(iname);
    elseif FileExists(inamec)
        disp(['        ' iname]);
        m=ReadEMFile(inamec);
        iname=inamec;
    else
        [iname pa]=uigetfile('*m.mrc','Find the merged image');
        cd(pa);
        m=ReadEMFile(iname);
    end;
    
    %%
    
    % Get image and pixel sizes
    n=size(m);
    ds=mi.imageSize(1)/n(1);  % downsampling factor of m
    pixA=mi.pixA*ds;    % pixel size of m
    
    % get the vesicle subtraction
    disp('Subtracting vesicles');
    vfile=[mi.tempPath mi.baseFilename 'v.mrc'];
    ok=0;
%     if isfield(mi,'tempPath') && exist(vfile,'file')
%         vs=ReadMRC(vfile);
%         ok=all(size(vs)==n);
%     end;
%     if ~ok
        vs=meMakeModelVesicles(mi,n);
%     end;
%   Do a least-squares amp correction
    vesicleAmpCorrection=(vs(:)'*m(:))/(vs(:)'*vs(:))
    vs=vs*vesicleAmpCorrection;
    msub=m-vs;
    
    figure(1);
    SetGrayscale;
    subplot(2,2,1);
    imacs(GaussFilt(msub,.1));
    title(iname,'interpreter','none');
    drawnow;
    
    % Get the effective CTF from the merging.
    H=ifftshift(meGetEffectiveCTF(mi,n,ds));
    
    %% Get the variance map
    disp('Computing the variance map');
    %     1.  Filter the image
    df=1./(n*pixA);  % frequency step
    filt=fuzzymask(n,2,f1./df,f1/(4*df(1)))-fuzzymask(n,2,f0./df,f0/(4*df(1)));
    mfilt=real(ifftn(fftn(msub).*ifftshift(filt)));
    %     2.  Get the convolution box
    box=Crop(SquareWindow(nb+2,2),n);
    boxo=ifftshift(box);
    localVar=real(ifftn(fftn(mfilt.^2).*fftn(boxo)));
    mask=meGetMask(mi,n);
    mask=Crop(Crop(mask,n-ceil(max(64/ds+nb/2,n/6))),n);  % Force a band of zeros at edges
    mxVar=max(localVar(:));
    lVar=localVar.*mask+mxVar*(1-mask);
    subplot(2,2,1);
    imacs(lVar)
    colormap jet
    drawnow;
    %% Get some boxes
    lv2=lVar;
    blankBox=ifftshift(Crop(ones(nb*2),n));
    boxX=zeros(1,nBoxes);
    boxY=zeros(1,nBoxes);
    imgs=zeros(nb,nb,nBoxes);
    for i=1:nBoxes
        [val ix iy]=max2d(-lv2);
        lv2=max(lv2,circshift(blankBox,[ix iy])*mxVar);  % blank the region
        imgs(:,:,i)=ExtractImage(msub,[ix iy],nb);
        %         disp([-val ix iy]);
        boxX(i)=(ix-1)*ds;
        boxY(i)=(iy-1)*ds;
    end;
    SetGrayscale;
    subplot(2,2,1);
    ShowImageAndBoxes(GaussFilt(msub,.1),[boxX' boxY']/ds,nb,8,[1 1 0]);
    drawnow;
       
    sp=mean(RadialPowerSpectrum(imgs),2);

    f=(0:nb/2-1)'/(nb*pixA);
    c=sectr(meGetEffectiveCTF(mi,nb,ds));
    
    subplot(2,2,2);
    plot(f,[sp/40 c]);
    %%
    % Fit the power spectrum
    
    noiseModelFcn='NoiseModel1';
        
    niters=1000;
    
    df=1/(pixA*nb);
    freqs=(0:df:(nb/2-1)*df)';  % frequencies in spectrum of a single box
    freqsx=(0:df/10:(nb/2-1)*df)';  % 10x oversampled frequencies
    
    
    % Get the effective CTFs
    [c effctf]=meComputeMergeCoeffs(freqs, mi.ctf, mi.doses);
    [c effctfx]=meComputeMergeCoeffs(freqsx, mi.ctf, mi.doses);
    
    % fit the noise model

    nPSets=2;  % try this many sets of parameters
    % Initialize parameters:
    q=sp(nb/4);
    
    %        af1 af2 ag    sigma  bf  s0 f1exp f2exp
    p0=[     1*q 1*q 10*q .007  100  q   1.5   1.5];
    ac=[     1   1    1    1    0   1    1     1 ];

    p0(2,:)=[1*q 1*q  0*q .007  100  q   1.5   1.5];
    ac(2,:)=[1   1    0    0    1   1    1     1 ];

    ps=p0;   % Store the final parameters
    errs=zeros(size(p0,1),1); % store the final errors.
    
    for jp=1:size(p0,1);  % loop over sets of parameters
        p=Simplex('init',p0(jp,:),ac(jp,:));
        for i=1:niters
            [spec, shot]=eval([noiseModelFcn '(freqs,p)']);
            model=spec.*effctf.^2+shot;
            d=(model-sp);
            err=d'*d;
            p=Simplex(err,p);
            if mod(i,50)==0
                plot(freqs,sp,'k.',freqs,model,'b-');
                title(i);
                xlabel('Spatial frequency, A^{-1}');
                drawnow;
            end;
        end;
        errs(jp)=err;
        ps(jp,:)=Simplex('centroid');
    end;
%     Pick the one which converged better
    [minErr jBest]=min(errs);
    p=ps(jBest,:);
    
    
    mi=meStoreNoiseModel(p,noiseModelFcn,mi);
    save([rootPath infoPath fname{fileIndex}],'mi');
    disp(['Updated ' fname{fileIndex}]);

    Ti=meGetNoiseWhiteningFilter(mi,n);
    f2d=RadiusNorm(n);
%     % Compute the effective ctf and inverse filter
%     f2d=RadiusNorm(n)/pixA;
%     [coeffs effct]=meComputeMergeCoeffs(f2d,mi.ctf,mi.doses);
%     [spec shot]=meEvalNoiseModel(f2d,mi);
%     %             [spec shot]=eval([noiseModelFcn '(f2d,p)']);
%     T=(spec.*effct.^2+shot)./shot;
%     Ti=1./sqrt(T);
    subplot(2,2,3);                         %%
    plot(sectr(f2d),sectr(Ti));             %%
%     axis([0 inf 0 1.02*max(sectr(Ti))]);    %%
    title('Prewhitening filter');           %%
    xlabel('Normalized spatial frequency');    %%
    mf=real(ifftn(fftn(m).*fftshift(Ti)));
    subplot(2,2,4);                         %%
    imacs(BinImage(mf,4));
    axis off;
    disp(mi.baseFilename);
    title(mi.baseFilename,'interpreter','none');
    drawnow;
    if writeOut
        outName=[mi.baseFilename 'mw.mrc'];
        disp(['Writing ' outName]);
        WriteMRC(mf,pixA,[mi.procPath outName]);
    end;
    disp(' ');
end;