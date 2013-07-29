function [angleList listIndices]=rsListHemisphereAngles(nAlpha, nBeta)
% Create an array t x 2 of angle values which are ~uniform on a hemisphere,
% 
% If nBeta is odd, we include a point at 90 degees; if even, not.
% e.g. nBeta = 2 yields 0 and 60 degrees.

corr=(mod(nBeta,2)+1)/2; % .5 even, 1 odd
betaStep=90/(nBeta-corr); % odd nBeta gives a final value at 90 degrees;
% even nBeta gives a final value of 90*(nBeta-1)/(nBeta-.5)), so that
% thetas will be evenly spaced when the lower hemisphere is included.
k=0;
angleList=zeros(nAlpha*nBeta,2);  % allocate an array bigger than we need.
listIndices=int32(zeros(nAlpha,nBeta));
for i=1:nBeta
  beta=(i-1)*betaStep;
  na=min(round(sind(beta)*nAlpha+.5),nAlpha);
  alphaStep=360/na;
%   [beta na]
  for j=1:na
      alpha=(j-1)*alphaStep;
      k=k+1;
      listIndices(j,i)=k;  % lookup
      angleList(k,:)=[alpha beta];
  end;
end;
angleList=angleList(1:k,:);  % truncate the list.
