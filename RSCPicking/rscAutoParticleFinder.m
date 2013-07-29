function [coords, ovMask, endCC]=rscAutoParticleFinder(mi,rscc,pars)
% %
% % We search for peaks in rscc.mxCC and check rscc.mxVars for excessive
% % variance.
% % rscc.MxVesInds gives the vesicle number for each particle; 
% % rscc.mxTempInds is used to get the template number.
% % each row of coords contains [x y flag vesIndex amp templ 0 0]
% % where flag = 32 (autopicked particle)
% rsoOffset (=pars(4) is how far (in angstroms) the particle center may lie
% inside the membrane center and still be counted as a bona fide right side
% out particle.   rsoOffset=0
% indicates that we take both rso and 


minAmp=pars(1);
maxAmp=pars(2);
maxVar=pars(3);
rsoOffset=pars(4)/mi.pixA;  % offset of particle center from tip (for RSO selection)
partRadius=pars(5)/mi.pixA; % blanking radius around particle, in A
overlapRadius=pars(6); % "forbidden zone" for overlap around each vesicle
maxBob=pars(7)/mi.pixA; % maximum distance of particle center outside the vesicle.
border=pars(8)/mi.pixA;   % border in A

n=size(rscc.mxCC,1);
ovMask=single(zeros(n,n));
ds1=mi.imageSize(1)/n;
n0=mi.imageSize(:)';

nb=2*ceil(partRadius/ds1)+4;  % size of the blanking mask box
blankMask=1-fuzzymask(nb,2,partRadius/ds1,1);

coords=single(zeros(1000,8));
k=0;
nFound=0;

if max(rscc.mxVesInds(:))>numel(mi.vesicle.x)  % out of bounds vesicle numbers
    disp('Autopicking aborted: mismatch between vesicle picker and preprocessing');
    return
end;

% Draw a disc of r+overlapRadius around each vesicle center.  Blank the CC
% wherever these discs overlap.
[mxCC2, ovMask]=rscBlankOverlaps(mi,rscc.mxCC,overlapRadius);

% Search for cc peaks
[amp, ix, iy]=max2d(mxCC2);
while amp>minAmp
    iv=rscc.mxVesInds(ix,iy);
    ccVar=rscc.mxVars(ix,iy);
    if (amp < maxAmp) && (ccVar < maxVar)
        [ampi, xi, yi]=max2di(mxCC2);  % get the interpolated values
        partCtr=([xi yi]-1)*ds1;
        if iv>0 % There is a matching vesicle.
            vesCtr=[mi.vesicle.x(iv) mi.vesicle.y(iv)];
            vesR=mi.vesicle.r(iv);
            rPart=sqrt(sum((partCtr-vesCtr).^2)); % dist particle to ves ctr
%             Geometry check:
%               -particle center is beyond r-rsoOffset, or else rsoOffset=0
%               -particle center is within r+maxBob
%               -particle center is beyond the border
            if (rsoOffset==0 || rPart > vesR-rsoOffset)...
                    && (rPart <= vesR+maxBob)...
                    && all(abs(partCtr-n0/2) < (n0/2-border))
            nFound=nFound+1;
            templ=single(rscc.mxTemplInds(ix,iy));
            coords(nFound,:)=[ partCtr 32 single(iv) amp templ 0 0];
            end;
        end;
    end;
    % Blank the vicinity of the found peak
    mxCC2=Mask(mxCC2,[ix iy],blankMask);
    [amp, ix, iy]=max2d(mxCC2);
end;
coords=coords(1:nFound,:);
endCC=mxCC2;
