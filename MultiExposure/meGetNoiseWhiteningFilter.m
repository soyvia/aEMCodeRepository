function H=meGetNoiseWhiteningFilter(mi,n,ds)
% function H=meGetNoiseWhiteningFilter(mi,n,ds)
% Given the micrograph info structure mi, returns a filter frequency
% response H of size n, f=0 in the center, that is the inverse of the
% mi.noiseModel as filtered by the effective CTF.
% The result is for an image downsampled by ds (default ds=mi.imageSize/n).
% No need to give ds if the image wasn't cropped.
%    H=sqrt(shot/(excessSpectrum * ctf^2 + shotSpectrum)).
if nargin<3
    ds=mi.imageSize(1)/n(1);
end;
if ~isfield(mi,'weights')
    mi.weights=mi.doses>0;
end;
if numel(mi.noiseModelPars)<1 || any(isnan(mi.noiseModelPars))
    H=1;
%     warning('No noise model found in mi file.');
    return
end;
pixA=mi.pixA*ds;  % Downsampled pixel size
f2d=RadiusNorm(n)/pixA;  % frequency in A^-1
weights=mi.weights;
for i=1:numel(mi.ctf)
    weights(i)=mi.weights(i)*mi.ctf(i).ampFactor;
end;
% Get the effective CTF and noise spectrum
[coeffs effct]=meComputeMergeCoeffs(f2d,mi.ctf,mi.doses,1,weights);
[spec shot]=meEvalNoiseModel(f2d,mi);
spec=max(spec,0);  % This excess noise should be positive
H=sqrt(mean(shot(:))./(spec.*effct.^2+shot)); % =1 at high frequencies
