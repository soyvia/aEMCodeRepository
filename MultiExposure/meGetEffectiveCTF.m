function c=meGetEffectiveCTF(mi,n,ds,noflip,nzeros)
% function c=meGetEffectiveCTF(mi,n,ds,noflip)
% Given the micrograph info structure mi, the size n of the desired 2D ctf
% and the downsampling ratio ds, return the effective ctf including the
% effects of merging and CCD (after doing CCD pre-whitening).
% By default, ds=mi.imageSize(1)/n.  You shouldn't use the default if n
% represents a cropped image, e.g. a particle image. The default is used if
% ds is given as 0.  The returned c has zero frequency in the center.  By
% default c is non-negative.  If noflip=1, c alternates sign, following the
% polarity of first image's ctf. (Merged images, or single images processed
% by the merging routines, are always phase-flipped.)
% fs 3 Sep 11

if numel(n)<2
    n=[1 1]*n;
end;
if nargin<3 || ds==0
    ds=mi.imageSize(1)/n(1);
end;
if nargin<4
    noflip=0;
end;
if nargin<5
    nzeros=1;  % default is to include only the 1st zero in later images
end;
if ~isfield(mi,'weights') || numel(mi.weights)<numel(mi.doses)
    mi.weights=single(mi.doses>0);
end;
effPixA=mi.pixA*ds;
freqs=RadiusNorm(n)/effPixA;  % frequencies for evaluating the CTF
% Change the weights to include the ampFactors
weights=mi.weights;
if isfield(mi.ctf,'ampFactor')
    for i=1:numel(weights)
        weights(i)=weights(i)*mi.ctf(i).ampFactor;
    end;
end;
[coeffs, effctf]=meComputeMergeCoeffs(freqs, mi.ctf, mi.doses,nzeros,weights);
if noflip
    effctf=effctf.*sign(coeffs(:,:,1));
end;
c=effctf.*CCDEffCTF(n,ds);
