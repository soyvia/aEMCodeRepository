function [mergeMats corrCTPars Ps]=meAlignMultiExposures(m,pixA,rawCTFs,doses,ds, Ps)
% function [mergeMats corrCTPars]=meAlignMultiExposures(m,pixA,CTFitPars,doses,ds)
% Given the stack of images m (n x n x nim), align them all, returning the
% transformation matrices mergeMats (3 x 3 x nim).  Alignments are carried
% out sequentially (m(:,:,2) aligned to m(:,:,1), m(:,:,3) aligned to
% m(:,:,2) etc.) and then the products of the transformations are taken.
% In the end, if you compute
% ma(:,:,j)=AffineTransform(m(:,:,j), mergeMats(:,:,j)) for each j, the
% resulting ma images will all be aligned.
% pixA is the angstroms/pixel.  CTFitPars is an array (nim) of the CTF
% parameter structures. ds is the optional downsample factor for aligning;
% default is 4.
% A fixed noise correction is performed.  A cross-correlation is computed
% between the first and last image.  It is assumed that the defocus
% difference between these will be large enough that the actual CC peak
% will be separate from the zero-lag peak.  We take the zero-lag peak to be
% spurious (from common gain references) and we null that peak out by
% adding appropriate anti-correlated noise.

% Variable naming conventions:
%   d means downsampled
%   z means the zero-frequency point is in the center.
antiNoiseFactor=1;  % should be 1, I don't know why not..


if nargin<5
    ds=4;  % downsampling factor for alignment, relative to original micrograph
end;

if nargin>5 && any(Ps(:))
    refineOnly=1;
else
    refineOnly=0;
    Ps=[];
    disp('Aligning exposures for merging');
end;

doFixedNoiseCorrection=1;
displayFNC=0;  % Don't display it.
windowFractionalWidth=16;  % used to be 128, but this seems better.


%%
[nx ny nim]=size(m);

n0=[nx ny];
nxd=nx/ds;
nyd=ny/ds;
nd=[nxd nyd];

pixAd=pixA*ds;  % Angstroms per pixel

win=single(SquareWindow(n0,nx/windowFractionalWidth));
winFraction=sum(win(:))/numel(win);
% fy=fftn(y.*win);  % FT of compensation noise


nit=[40 30 100];  % iterations of the Simplex fitting.
nit=[60 40 200];  % iterations of the Simplex fitting.
nit0=100;  % initial fit

% allocate arrays

Tmats=zeros(3,3,nim);     % transform matrices
Tmats(:,:,1)=eye(3);
md=single(zeros(nxd,nyd,nim));  % downsampled images
ctz=single(zeros(nxd,nyd,nim)); % corresponding CTFs
fd=complex(single(zeros(nxd,nyd,nim)));  % FTs
fd0=fd;
dsmaskz=single(zeros(nx,ny,2));  % downsampling masks
dsmaskz(:,:,1)=single(fuzzymask(n0,2,nd*.45,nd*.05));
dsmaskz(:,:,2)=single(fuzzymask(n0,2,nd*.225,nd*.05));

% Do fixed-noise correction on the two highest-defocus images.
if nim>1 && doFixedNoiseCorrection
    [fy cc]=meMakeCompNoiseAuto(m(:,:,nim-1),m(:,:,nim));
    dose0=sqrt(prod(doses(nim-1:nim)));
else
    fy=0;
    dose0=1;
end;

% Process and downsample the images.
figure(1); clf;
SetGrayscale;
ndis=64;
%%
for i=1:nim
    f0=fftn(m(:,:,i).*win);
    f1=f0+(-1)^i*doses(i)/dose0*fy*winFraction*antiNoiseFactor;  % put in fixed-noise correction
    %     make the FT of downsampled images with anti-noise
    fd(:,:,i) =ifftshift(Crop(fftshift(f1).*dsmaskz(:,:,1+(i>1)),nd));
    
    % Optional: display the residual correlations to check fixed-noise.
    if doFixedNoiseCorrection && displayFNC
        %         same as fd, but with no anti-noise
        fd0(:,:,i)=ifftshift(Crop(fftshift(f0).*dsmaskz(:,:,1+(i>1)),nd));
        if i>1
            ccs0=Crop(fftshift(real(ifftn(fd0(:,:,i-1).*conj(fd0(:,:,i))))),ndis);
            ccs=Crop(fftshift(real(ifftn(fd(:,:,i-1).*conj(fd(:,:,i))))),ndis);
            subplot(3,nim,i+nim);  % second row of figure, column i
            imacs(ccs0);
            title('CC not corrected');
            drawnow;
            subplot(3,nim,i+2*nim); % 3rd row of figure
            imacs(ccs);
            title('CC corrected');
            drawnow;
        end;
    end;
    % Construct the downsampled images
    md(:,:,i)=real(ifftn(fd(:,:,i)));
    if displayFNC
        subplot(3,nim,i);  % top row of figure
        imacs(md(:,:,i));  % original image
        axis off
        title(['Image ' num2str(i)]);
        drawnow;
    end;
    ctz(:,:,i)=single(CTF(nd,pixAd,rawCTFs(i)));
end;
fdz=RadiusNorm(nd);
%%
for i=2:nim  % find Tmat such that Tmat(image i) matches (image i-1)
    figure(2);
    clf;
    SetGrayscale;
    
    % Matching filter.  Empirically a sqrt(c1*c2) filter looks good.
    %     filt=ifftshift( sign(ctz(:,:,i-1)).*sign(ctz(:,:,i))...
    %         .*sqrt( abs(ctz(:,:,i-1).*ctz(:,:,i)).*fdz ) );  % added sqrt(fdz)
    
    filt=ifftshift( sign(ctz(:,:,i-1)).*sign(ctz(:,:,i))...
        .*( abs(ctz(:,:,i-1).*ctz(:,:,i)).*fdz ) );  % added fdz factor
    
    %     disp(' aligning');
    %     disp(' displayed parameters are: theta (radians), mag, stx, sty')
    deltadef=rawCTFs(i).defocus-rawCTFs(i-1).defocus;
    theta=-2e-4*deltadef;  % guess at rotation
    mag=1+1.6e-3*deltadef;  % guess at magnification
    if ~refineOnly  % Make a de novo fit.

        %     First, do two initial alignments at very low resolution, and
        %     a rought refinement of each
        ds1=4;  % another 4x downsampling
        filt1=Downsample(filt,nd/ds1);
        md1=Downsample(md,nd/ds1,1);
        P=[theta mag 0 0 0 0];  % keep the skews at zero.
        Psteps0=[.01 .003*deltadef+.01 0 0 0 0];
        Psteps1=[.01 .01 0 0 0 0];
        
        P01=meMergeAligner2(md1(:,:,i-1),md1(:,:,i),filt1,P,Psteps0,nit0,1);
        [P01 T cc01]=meMergeAligner2(md(:,:,i-1),md(:,:,i),filt,P01,Psteps1,nit(1),1);
        
        Psteps0=[.04 .005*deltadef+.04 0 0 0 0];
        P02=meMergeAligner2(md1(:,:,i-1),md1(:,:,i),filt1,P,Psteps0,nit0,1);
        [P02 T cc02]=meMergeAligner2(md(:,:,i-1),md(:,:,i),filt,P02,Psteps1,nit(1),1);
        
        mx1=max(cc01(:))/std(cc01(:));
        mx2=max(cc02(:))/std(cc02(:));
        if mx1 > mx2
            P=P01;
        else
            P=P02;
        end;
        disp(['Initial alignment SNR = ' num2str(max(mx1,mx2))]);

        %     %% Optimize the alignment
        % %         Psteps1=[.01 .01 0 0 0 0];
        % %         % We vary the transformation of image i to match image i-1.
        % %         P=meMergeAligner2(md(:,:,i-1),md(:,:,i),filt,P,Psteps1,nit(1),1);
        % %         %     P=MergeAligner2(md(:,:,i-1),md(:,:,i),filt,P,Psteps1,10,1);
        % %
        % Next, vary stretch only
        P(5:6)=0;  % Get rid of translation
        Psteps2=[0 0 10/nx 10/ny 0 0];
        P=meMergeAligner2(md(:,:,i-1),md(:,:,i),filt,P,Psteps2,nit(2),1);
        
    else
        P=Ps(i,:);
        Psteps1=[.002 .002   0    0    0 0];
        Psteps2=[  0   0  .0002 .0002  0 0];
    end;  % ~refineOnly
    
    % Finally, vary all parameters
    P(5:6)=0;
    Psteps3=(Psteps1+Psteps2);
    [P T]=meMergeAligner2(md(:,:,i-1),md(:,:,i),filt,P,Psteps3,nit(3),1);
    
    % Store the transformation matrix.
    Tmats(:,:,i)=T;
    Ps(i,:)=P;
end;

% Form the composite matrices.  to match image 3 to image 1, we operate
% with T1*T2*T3, while image 1 is just operated on by T1.
mergeMats=Tmats;  % Allocate it, and copy the first matrix.
for i=2:nim % form the composite matrices
    mergeMats(:,:,i)=mergeMats(:,:,i-1)*Tmats(:,:,i);
end;

% Update the CTF parameters
% in view of the magnification changes.
corrCTPars=rawCTFs;  % Most fields are unchanged.
for i=2:nim
    q=mergeMats(1:2,1:2,i);   % Get the rotate/mag part
    scalesq=prod(sqrt(sum(q.^2)));  % magnification squared
    corrCTPars(i).defocus=rawCTFs(i).defocus*scalesq;
    corrCTPars(i).deltadef=rawCTFs(i).deltadef*scalesq;
end;

