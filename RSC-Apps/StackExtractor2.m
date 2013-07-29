% StackExtractor2
% Version 2 uses vesicle shifts in subtraction and particle extraction.
% For each info file, read the micrographs and weight them according to the
% merge coefficients. For each particle, extract its region from each
% micrograph, shifted according to the mi.vesicle.shiftX and shiftY values.
% and sum the extracted images.
% Rotate the merged particle to
% the standard position, and write *st.mrc and *stu.mrc files (stack and
% unsubtracted stack) into the Stack/ directory.  Also store *stall.mrc the
% merged stack from all images, and *tsi.mat which contains the si structure.
% Contrast is reversed so protein is white.
% The weights variable determines which exposures are used.  If the weights
% are not all 1s, e.g. [0 1 0] for the second exposure only, the saved files
% have the form *w010st.mrc and *w010tsi.mat.

boxSize=48;  % Size of boxes to be extracted from merged images.
ds=4;        % downsampling of boxed particles from original micrograph
shiftMax=20; % maximum allowed vesicle shift (orig. pixels) between exposures
cpe0=16;     % default value, if not already in mi file
iCamera=1;

dds=2;       % further downsampling for micrograph display
% weights=[0 1 0]  % which exposures are used.
weights=[1 1 1]  % which exposures are used.
nZeros=1;    % frequency limit for merging 2nd and 3rd exposures
vindex=0;
padSize=NextNiceNumber(boxSize*1.5);
% Make the upsampled pad mask, and the particle fourier mask
padMask=fuzzymask(padSize,2,padSize*.48,padSize*.04);
fmasko=ifftshift(fuzzymask(padSize,2,.45*padSize,.05*padSize));

[fname, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;

cd(rootPath);
if ~exist('Stack','dir');
    mkdir('Stack');
end;

modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD
%%
totalNParts=0;
np=1e6;  % preliminary stack size
si=struct;
si.miIndex=     uint16(zeros(np,1));
si.miParticle=  uint16(zeros(np,1));
si.alpha0=      single(zeros(np,1));
si.yClick=      single(zeros(np,1));  % in units of si.pixA
si.rVesicle=    single(zeros(np,1));
si.sVesicle=    single(zeros(np,1));
si.mi=cell(0);
nfiles=numel(fname);
ctfs=single(zeros(boxSize,boxSize,nfiles));
pixA0=0;  % unassigned value


for fileIndex=1:numel(fname)
    disp(['Reading ' fname{fileIndex}]);
    load([infoPath fname{fileIndex}]);  % Load the mi file
    % Make new entries into the mi file
    mi.boxSize=boxSize;
    if ~isfield(mi,'stackPath')
        mi.stackPath='';
    end;
    mi.stackPath='Stack/';
    if ~isfield(mi,'cpe')
        mi.cpe=cpe0;
    else
        cpe=mi.cpe;
    end;
    origWeights=mi.weights;
    mi.weights=weights;  % Modify the weights.
    
    if isfield(mi.particle,'picks')
        nPicks=size(mi.particle.picks,1);
    else
        nPicks=0;
    end;
    n0=mi.imageSize/ds;
    %     Get the final pixel size
    pixA=mi.pixA*ds;
    if pixA0==0
        pixA0=pixA;
    end;
    if pixA0~=pixA
        warning(['Change in pixA values: ' num2str([pixA0 pixA]) '  ' fname{fileIndex}]);
    end;
    
    if nPicks>0 % there are particles
        %     Noise-whitening filter
        if numel(mi.noiseModelCode)>0
            disp('Noise-whitening filter');
            Ti=meGetNoiseWhiteningFilter(mi,n0,ds);  % for micrograph
            Tic=meGetNoiseWhiteningFilter(mi,boxSize,ds);  % for particle ctf
            %             msubFilt=real(ifftn(fftn(msub).*ifftshift(Ti)));
            %             mFilt=real(ifftn(fftn(m).*ifftshift(Ti)));
        else
            disp('No noise-whitening filter.')
            Ti=1;
            Tic=1;
        end;
        disp('Compute image set');
        mSet=meMakeMergeImageSet(mi,cpe,ds,Ti,nZeros); % downsampling by ds
        mMerge=sum(mSet,3);  % merged image for checking
        nim=size(mSet,3);
        imSet=find(mi.weights>0);  % The images we'll actually use.
        figure(1);
        SetGrayscale;
        subplot(2,1,1);
        imacs(BinImage(mMerge,dds));
        title(mi.baseFilename,'interpreter','none');
        drawnow;
        
        disp('Compute vesicle models');
        doCTF=1;
        vSet=meMakeModelVesicleSet(mi,mMerge,vindex,doCTF,Ti,nZeros);
        mvSet=mSet-vSet;  % subtracted images for extraction.
        mvMerge=sum(mvSet,3);
        subplot(2,1,2);
        imacs(BinImage(mvMerge,dds));
        drawnow;
        % mtf=TestMarkParticles(mi,mMerge);       
        % ctf, including the prewhitening filter effects, for individual
        % particles
        freqs=RadiusNorm(boxSize)/pixA;
        [pcoeffs, pctf]=meComputeMergeCoeffs(freqs,mi.ctf,mi.doses);
        ctfs(:,:,fileIndex)=pctf.*Tic;
        %%
        nPickerEntries=nPicks
        if totalNParts==0  % the first time, allocate part of the total stack.
            totalStack=single(zeros(boxSize,boxSize,nPicks));
        end;
        nv=numel(mi.vesicle.x);
        vCoords=[mi.vesicle.x(:) mi.vesicle.y(:)];
        vRadii=mi.vesicle.r(:);
        stackSub=single(zeros(boxSize,boxSize,nPicks));
        stackImg=single(zeros(boxSize,boxSize,nPicks));
        sumImg=zeros(boxSize,boxSize);
        sumSub=zeros(boxSize,boxSize);
        if isfield(mi.vesicle,'shiftX') && numel(mi.vesicle.shiftX)>0
            shiftX=mi.vesicle.shiftX;
            shiftY=mi.vesicle.shiftY;
        else
            disp('No shift values have been set.');
            shiftX=zeros(nv,nim);
            shiftY=zeros(nv,nim);
        end;
        bad=(abs(shiftX)>shiftMax) | (abs(shiftY)>shiftMax);
        shiftX(bad)=0;
        shiftY(bad)=0;
        disp('Extracting');
        figure(2);
        SetGrayscale;
        j=0;
        for i=1:nPicks  % scan the particle coordinates
            type=mi.particle.picks(i,3);
            if (type>=16 && type<=32)  % valid particle flags
                j=j+1;  % counter
                pCoords=mi.particle.picks(i,1:2);
                vInd=mi.particle.picks(i,4);
                if vInd<1  % no vesicle marked.  Find the nearest vesicle
                    dists=sqrt(sum((vCoords-repmat(pCoords,nv,1)).^2,2));
                    [minDist, vInd]=min(dists);
                    if minDist > vRadii(vInd) % outside of membrane, check distance from mbn.
                        distm=abs(dists-vRadii);
                        [minDist, vInd]=min(distm);
                    end;
                    mi.particle.picks(i,4)=vInd;  % insert our estimated vesicle index.
                end;
                
                xves=zeros(padSize,padSize,nim);
                ximg=zeros(padSize,padSize,nim);
                %         Up to this point, all coordinates are in original pixel size.
                for k=imSet  % loop over active images
                    intCoords=round(pCoords/ds);  % downsampled coordinates
                    fraCoords=pCoords/ds-intCoords;
                    xm=ExtractImage(mSet(:,:,k),intCoords+1,padSize);
                    xv=ExtractImage(vSet(:,:,k),intCoords+1,padSize);
                    shifts=[shiftX(vInd,k) shiftY(vInd,k)]/ds+fraCoords;
                    ximg(:,:,k)=real(ifftn(fftn(xm).*fmasko.*FourierShift(padSize,-shifts)));
                    xves(:,:,k)=real(ifftn(fftn(xv).*fmasko.*FourierShift(padSize,-shifts)));
                end;
                %%
                %             merge the padded particle images, reverse the contrast.
                ximg2=-sum(ximg,3).*padMask;
                xsub2=-sum(ximg-xves,3).*padMask;
                %             Compute the angle relative to the normal
                pVec=pCoords-vCoords(vInd,:);
                alpha=atan2(pVec(2),pVec(1))-pi/2;
                %        alpha=0;
                %               rotate them
                pSub=Crop(grotate(xsub2,-alpha),boxSize);
                stackSub(:,:,j)=pSub;
                sumSub=sumSub+pSub;
                
                pImg=Crop(grotate(ximg2,-alpha),boxSize);
                stackImg(:,:,j)=pImg;
                sumImg=sumImg+pImg;
                
                subplot(2,2,1);
                imacs(pSub);
                title(mi.baseFilename,'interpreter','none');
                subplot(2,2,2);
                imacs(pImg);
                subplot(2,2,3)
                imacs(sumSub);  % unrotated image
                title(j);
                subplot(2,2,4);
                imacs(sumImg);
                drawnow;
                
                jt=j+totalNParts;  % pointer to overall stack index
                si.miIndex(jt)=     fileIndex;
                si.miParticle(jt)=  i;
                si.alpha0(jt)=      alpha*180/pi;
                si.yClick(jt)=      hypot(pVec(1),pVec(2))/ds;  % in units of si.pixA
                si.rVesicle(jt)=    mi.vesicle.r(vInd)/ds;
                si.sVesicle(jt)=    mi.vesicle.s(vInd);
            end; % if type
        end; % for i
        nparts=j;
        stackSub=stackSub(:,:,1:nparts);  % truncate the stack
        stackImg=stackImg(:,:,1:nparts);
        
        totalStack(:,:,totalNParts+1:totalNParts+nparts)=stackSub;
        totalNParts=totalNParts+nparts
        si.mi{fileIndex}=mi;  % store a copy of the micrograph info
        mi.weights=origWeights;   % Replace the original weights.
        save([infoPath fname{fileIndex}],'mi');
        disp([infoPath fname{fileIndex} ' updated']);
        
        
        outname=[mi.stackPath mi.baseFilename 'st.mrc'];
        WriteMRC(stackSub,pixA,outname);
        outname=[mi.stackPath mi.baseFilename 'stu.mrc'];
        WriteMRC(stackImg,pixA,outname);
    else
        disp(' -no particle picks found');
    end;  % if np
    
end; % for fileIndex
%%
if totalNParts>0
    
    si.pixA=pixA;
    np=totalNParts;
    % truncate the si arrays.
    si.miIndex=     si.miIndex(1:np,1);
    si.miParticle=  si.miParticle(1:np,1);
    si.alpha0=      si.alpha0(1:np,1);
    si.yClick=      si.yClick(1:np,1);
    si.rVesicle=    si.rVesicle(1:np,1);
    si.sVesicle=    si.sVesicle(1:np,1);
    si.ctfs = ctfs;
    %
    if all(weights>0)
        weightString='';
    else
        weightString=['w' sprintf('%1d',weights>0)];
    end;
    outname=[mi.stackPath mi.baseFilename weightString 'stall.mrc'];
    WriteMRC(totalStack,pixA,outname);
    disp(['Wrote the stack ' outname]);
    
    % save the stackInfo structure
    outname=[mi.stackPath mi.baseFilename weightString 'tsi.mat'];
    save(outname,'si');
    disp(['Wrote the total stack info ' outname]);
else
    disp('--no particles, nothing written.');
    disp(' ');
end;