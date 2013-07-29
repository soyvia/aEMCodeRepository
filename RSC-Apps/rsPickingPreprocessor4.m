% rsPickingPreprocessor

membraneOffsetA = -70;  % Membrane center is this distance from particle center
localVarRadius=100;  % angstroms
maskPaddingA=20;     % extra space around outer radius, should be greater than maxBob in picker.

outputImageSize=1024;  % size of output images.
nterms=32;


showTemplates=0;
simulateImage=0;

mapName='/Volumes/TetraData/Structures/AMPAR/3KG2map58.mrc';
% mapName='/Volumes/TetraData/Structures/AMPAR/3KG2mapsub5.8A.mrc';

% Have the user select some mi files
[fname, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
if isnumeric(pa) % File selection cancelled
    return
end;
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(rootPath);

oldPixA=0;

for fileIndex=1:numel(fname); % Operate on a single micrograph
    miName=[infoPath fname{fileIndex}];
    disp(['Reading ' miName]);
    load(miName);
    mi.basePath=rootPath;
    if numel(mi.vesicle.x)<1
        disp('No vesicles found.  Terminating the program.');
        return
    end;
    %     Pick up the merged image
    [origImg, mergeFullPath]=meReadMergedImage(mi);
    %%
    %     Possibly downsample it
    n0=size(origImg);
    dsm=n0(1)/outputImageSize;
    if dsm~=1  % further downsampling
        m0=Downsample(origImg,outputImageSize);
    else
        m0=origImg;
    end;
     n=size(m0);
    ds=mi.imageSize(1)/n(1);  % downsampling relative to original image
    pixA=mi.pixA*ds;
    
    disp('Making model vesicles');
    ves=meMakeModelVesicles(mi,n,find(mi.vesicle.ok(:,3))); % everything that could be refined.
    m1=m0-ves;
    
    useNoiseWhitening=(numel(mi.noiseModelPars)>0);
    if  useNoiseWhitening
        disp('Noise whitening');
        H=meGetNoiseWhiteningFilter(mi,size(m0));
        m2=real(ifftn(fftn(m1).*ifftshift(H)));
    else
        disp('No specimen-noise whitening');
        m2=m1;
    end;
    %     m=m1-ves;  % vesicle subtraction
    %%
        figure(1); clf;
        SetGrayscale;
    imac(imscale(GaussFilt(m2,.2),256,.001));
    drawnow;
    %%
    if pixA~=oldPixA  % We haven't already made templates of the correct size
        oldPixA=pixA;
        % Load the 3D map
        disp('Loading the 3D map');
        [origMap, s]=ReadMRC(mapName);
        mpixA=s.pixA;
        nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
        nt=ceil(nt1/8)*8;
        [map, finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
        map=map*pixA;  % approx amplitude correction (V-A scaling)
        
        % magnifications=[mpixA/pixA finalmag]
        
        %% Create the list of angles for the templates
        membraneOffset=membraneOffsetA/pixA;
        nAlpha=32; % about 10 degrees
        nBeta=12;  % even is best.  Number on hemisphere.
        nGamma=8;
        symmetry=2;
        gammaStep=360/(symmetry*nGamma);
        
%         hemiAngles run from alpha=[0..360) and beta=[0..90)
        [hemiAngles, angleInds]=rsListHemisphereAngles(nAlpha, nBeta);
        nHemiAngles=size(hemiAngles,1);
        nHemi=2;  % both hemispheres
        
        angleList=zeros(nHemi,nGamma,nHemiAngles,3);
        for j=1:nGamma;
            gamma=(j-1)*gammaStep;
            for k=1:nHemiAngles
                angleList(1,j,k,:)=[hemiAngles(k,:) gamma];
                angleList(2,j,k,:)=[hemiAngles(k,1) 180-hemiAngles(k,2) gamma];
%                 angleList(2,j,k,:)=[[0 180]-hemiAngles(k,:) gamma];
            end;
        end;
        nAngles=numel(angleList)/3;
        
        % angle list is of size
        % (nHemi x nGamma x nHemiAngles, 3)  where nHemi=2 is the number of hemispheres.
        % Note that the beta angles alternate such that, if betastep=1, they are
        % are (0, 180) for each gamma, then (1, 179) for each gamma, up to 89, 91.
        % that is, the same projected position is described twice.
        
        %% Make the templates
        
        disp(['Making ' num2str(nAngles) ' templates']);
        
        tic
        % allTemplates=rsMakeTemplatesQuick(angleList,map);
        allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),map);
        toc
        allTemplates=reshape(allTemplates,nt,nt,nHemi,nGamma,nHemiAngles);
        
    end;
    
    %% Filter the templates according to the CTF
    [nt, nt, nHemi, nGamma, nHemiAngles]=size(allTemplates);
    nAngles=nHemi*nGamma*nHemiAngles;
    ne=NextNiceNumber(nt*1.2);  % increase the size to allow CTF rings
    ctf=meGetEffectiveCTF(mi,ne,ds);  % put in dqe, pw filter.
    if useNoiseWhitening
        h=meGetNoiseWhiteningFilter(mi,ne,ds);
        ctf=ctf.*h;
    end;
    %     % evaluate a generic inverse filter
    %     f=RadiusNorm(ne)/(mi.pixA*ds);  % frequencies for evaluating the CTF
    %     hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
    %     hinv(f==0)=0;

    % Pad the templates to avoid ctf artifacts
    xTemplates=Crop(reshape(allTemplates,nt,nt,nAngles),ne,1);
    nim=size(xTemplates,3);
%     operate with the CTF and mask
    H=ifftshift(ctf);
    msk=fuzzymask(ne,2,0.45*ne,.1*ne);
    for i=1:nim
        xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*H)).*msk;
    end;
    xTemplates=reshape(xTemplates,ne,ne,nHemi,nGamma,nHemiAngles);
    ntstr=num2str(nt);
    nestr=num2str(ne);
    disp(['Templates expanded from ' ntstr 'x' ntstr ' to ' nestr 'x' nestr]);
    
    %% Make the eigenreferences
    eigenSet=rsMakeEigenrefs(xTemplates,nterms);
    figure(2);
    SetGrayscale;
    plot(1-eigenSet.termVar);
    % figure(1);
    % ImagicDisplay(eigenSet.imgs,2);
    
    
    if showTemplates
        %%  % Show the templates and the reconstructions
        timgs=reshape(eigenSet.imgs,ne*ne,nterms);
        nG=1;
        nH=2;
        nAngs=nHemiAngles;
        rImg=single(zeros(ne,ne,2*nH,nG,nAngs));
        for k=1:nHemiAngles
            for j=1:nG
                for i=1:nH
                    %             rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne).*squeeze(eigenSet.ampList(i,j,k));
                    rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne);
                    rImg(:,:,2*i,j,k)=xTemplates(:,:,i,j,k);
                end;
            end;
        end;
        figure(1);
        ImagicDisplay2(rImg,2);
        
        %%  % Make a complete set of reconstructions for comparison
        %   and make the average power spectrum
        nTotal=nAngs*nGamma*nHemi;
        timgs=reshape(eigenSet.imgs,ne*ne,nterms);
        rImg=single(zeros(ne,ne,nTotal));
        vL=reshape(eigenSet.vList,nterms,nTotal);
        sp=zeros(ne,ne);
        sp0=zeros(ne,ne);
        for i=1:nTotal
            %     img=reshape(timgs*vL(:,i),ne,ne)*eigenSet.ampList(i);
            img=reshape(timgs*vL(:,i),ne,ne);
            sp0=sp0+abs(fftn(xTemplates(:,:,i))).^2;
            sp=sp+abs(fftn(img)).^2;
            rImg(:,:,i)=img;
        end;
        spr=Radial(fftshift(sp/nTotal));
        spr0=Radial(fftshift(sp0/nTotal));
        figure(2);
        semilogy([spr0 spr]);
        figure(3);
        plot([cumsum(spr0) cumsum(spr)]);
        ImagicDisplay2(rImg);
    end;
    % mic=rsSortVesicles(mi);  % make a copy with the vesicles sorted by position.
    mic=mi;
    
    %%  Evaluate the cc for each vesicle
    disp('Evaluating cross-correlations');
    disp([' using ' num2str(nterms) ' terms']);
    figure(1); clf; SetGrayscale;
    imac(imscale(GaussHP(GaussFilt(origImg,.1),.005),256,.003));    
    
    %% Pick vesicles to work on
%   definition: vesicles.ok(i,:) = [aVesicle inRange refined -- ]
    okVesicles=all(mic.vesicle.ok(:,2:3),2); % in-range and refined
    medianS=median(mic.vesicle.s);
    mic.vesicle.s(okVesicles)=medianS;  % force amplitudes to the median
    
    allVesicles=mic.vesicle.ok(:,1);  % every vesicle that was found
    goodVesInds=find(okVesicles);
    badVesInds=find(allVesicles & ~okVesicles);  % exist, but not ok.

    gMaxVals=single(zeros(n));
    gMaxNCC=single(zeros(n));
    mxTemplInds=uint16(zeros(n));
    mxVesInds=uint16(zeros(n));
    mxVars=single(zeros(n));
    mxRsos=single(zeros(n));
    mxDist=single(ones(n))*max(n);
    figure(3);
    SetGrayscale;
    
    msklRadius=round(localVarRadius/pixA);
    nl=2*ceil(1.2*msklRadius);
    mskl=fuzzymask(nl,2,msklRadius,0.12*msklRadius);  % local mask
    npts=sum(mskl(:));
    
    nves=numel(goodVesInds);
    disp(['Total vesicles: ' num2str(nves)]);
    maskPadding=ceil(maskPaddingA/pixA);
    figure(3);
    
    for j=1:nves
        i=goodVesInds(j);
        maskRadii=mic.vesicle.r(i)/ds+membraneOffset*[-1 1 0]+maskPadding;
        maskRadii(3)=maskRadii(3)+abs(membraneOffset)+msklRadius+maskPadding;
        
        % Get the single-vesicle correlation function
        [mxVals, mxInds, mxNCC, mxValsV, mxRso, localVar]=...
            rsVesicleCorrelation4(-m2,mic,i,membraneOffset,...
                                  maskRadii,angleInds,eigenSet);
        localDist=RadiusNorm(size(mxVals));
        ctr=round([mic.vesicle.x(i) mic.vesicle.y(i)]/ds+1);  % shift to 1-based coords
        nv=size(mxVals,1);
        
        % Do a local averaging of the squared CC in the whole vesicle.
        var=mxVals.^2;
        h=ifftshift(Crop(mskl,nv));  % average over the local mask
        filtVar=real(ifftn(fftn(var).*fftn(h))).*fuzzymask(nv,2,max(maskRadii(1:2)),.5);
        
        %     Pad to the full-sized image
        xVals=ExtractImage(mxVals,ctr,n,1);
        xTemplInds=ExtractImage(mxInds,ctr,n,1);
        xVar=ExtractImage(filtVar,ctr,n,1);
        %         xVar=ExtractImage(localVar,ctr,n,1);
        xNCC=ExtractImage(mxNCC,ctr,n,1);
        xRso=ExtractImage(mxRso, ctr,n,1);
        xDist=Radius(n,ctr);  % Get the distance map too.
        %       Incorporate into the composite image
        q=xVar>mxVars;
        mxVars(q)=xVar(q);
        
        q=xDist<=mxDist;  % all points closer to this vesicle than any other
        mxDist(q)=xDist(q);
        %           q=xVals>gMaxVals;  % find maximum cc values
        gMaxVals(q)=xVals(q); % update the values where these ones are higher.
        %         q=xNCC>gMaxNCC;  % find maximum cc values
        gMaxNCC(q)=xNCC(q); % update the values where these ones are higher.
        mxVesInds(q)=i;  % the related vesicle indices
        mxTemplInds(q)=xTemplInds(q);  % the template indices
        mxRsos(q)=xRso(q); % the right-side-out flags.
        q1=xVar>mxVars;  % find maximum variance values
        mxVars(q1)=xVar(q1);
        
        if mod(i,10)==0
            imacs(gMaxVals);  % Show the CC function
            %         imacs(mxVars);  % Show the local variances
            title([num2str(i) ' / ' num2str(nves)]);
            drawnow;
        end;
        
    end;  % loop over vesicles
    %%  Zero the cc of bad vesicles
    badVesMask=single(zeros(n));
    nbad=numel(badVesInds);
    if nbad>0
        for i=1:nbad
            badVesMask=(badVesMask | mxVesInds==badVesInds(i));
        end;
        gMaxVals(badVesMask)=0;
        gMaxNCC(badVesMask)=0;
    end;
    %%
    
    figure(2);
    imacs(sqrt(mxVars));
    title('mxVars');
    figure(3);
    imacs(min(1,gMaxVals));
    title('gMaxVals');
%     figure(4); SetGrayscale;
%     imacs(gMaxNCC);
%     title('gMaxNCC');
    
    mxCC=gMaxVals;
    eigenImgs=eigenSet.imgs;
    vList=eigenSet.vList;
    
    partRadius=18;
    save([mi.procPath mi.baseFilename 'rscc.mat'],'mxCC','mxVars','mxVesInds',...
        'mxTemplInds','mxRsos','partRadius','ds',...
        'gMaxNCC','gMaxVals','badVesMask','eigenImgs','vList','angleList');
    disp(['written: ' mi.procPath mi.baseFilename 'rscc.mat']);
    disp(' ');
end;
