function aliImgs=TransformImages(imgs, xytfr)
% function aliImgs=TransformImages(imgs, xytfr);
% Transform the image stack imgs according to the shifts, rotates and flips
% in the array xytfr.
% xytfr contains the coordinate system operations needed to return the
% image to align to the reference.

[n ny nim]=size(imgs);
nim2=size(xytfr,1);

if nim2<nim
    error('Not enough elements in xytfr');
end;

aliImgs=single(zeros(n,n,nim));
parfor j=1:nim
    pars=xytfr(j,:);  % necessary for parfor operation.  Has to be float
    theta=pars(3);
    iflip=pars(4);
    im=shiftf(double(imgs(:,:,j)),-pars(1:2));
    if iflip
        im=circshift(flipud(im),[mod(n+1,2) 0]);  % flip before rotate.
    end;
    if theta~=0
        im=grotate(im,-theta);
    end;
    aliImgs(:,:,j)=im;
end;
