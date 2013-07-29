function rotImages=rsRotateImage(m,alphas,kernelsize)
%  Create a stack of rotated images given the input image or image stackm.
%  The size of m must be either n x n or n x n x numel(alphas).  n must be
%  divisible by 4.
%  Each image is rotated cw by an angle alpha, given in degrees.
%  This is basically grotate() modified to process a stack.  For a large
%  stack it runs about 3x faster than grotate.
% m is the input 2d image, n x n with n divisible by 4.
% m is assumed to be strictly band-limited, i.e. containing no spatial
% frequencies above 0.5-1/n in any direction.
% kernelsize=3 gives about 1% power variation under rotation.
% kernelsize=5 gives ~0.1% but takes 50% longer.
% n=256, takes 75 ms per image with kernelsize=5 on i7 quad.
% n=40, takes 2 ms per image.

if nargin<3
    kernelsize=5;
end;

[n n1 nim]=size(m);
nAlphas=numel(alphas);
if n~=n1
    error('Images must be square');
end;
if nim~=1 && nim~=nAlphas
    error('Image stack must contain either 1 or numel(alphas) images.');
end;

rotImages=single(zeros(n,n,nAlphas));

% % old code
% for i=1:nAlphas
%     rotImages(:,:,i)=grotate(m,-alphas(i)*pi/180,kernelsize);
% end;
% return

nw=kernelsize;

[w1 ov]=gridMakeKaiserTable(kernelsize,'grid');
w1=w1';

wcomp=gridMakePreComp(n,nw);
P2=gridMakeNullFT(n,2);

np=P2.np;
np1=P2.np1;
np1ctr=P2.np1ctr;
sp1=P2.sp1;
ovctr=ov/2+1;

% Vectorized code.
nw2=floor(kernelsize/2);  % half-width of kernel (=(nw-1)/2).

% Source coordinates
[is,js]=ndgrid(-np/2:np/2-1);
is=reshape(is,np*np,1);
js=reshape(js,np*np,1);

if nim==1
    P=gridMakePaddedFT(m,'grid',wcomp);
end;

for ind=1:nAlphas
    if nim>1
%         This could be sped up if coded here...
        P=gridMakePaddedFT(m(:,:,ind),'grid',wcomp);
    end;
    s=sind(alphas(ind));
    c=cosd(alphas(ind));
    plane=zeros(np*np,1);  % vector to receive the data
    
    % Object coordinates
    ip=c*is+s*js+np1ctr;
    ip=min(sp1+np,max(ip,sp1+1));  % prevent coords from going out of bounds
    ipint=round(ip);
    ipfrac=floor((ip-ipint)*ov)+ovctr;
    
    jp=-s*is+c*js+np1ctr;
    jp=min(sp1+np,max(jp,sp1+1));
    jpint=round(jp);
    jpfrac=floor((jp-jpint)*ov)+ovctr;
    
    addrs=ipint+np1*(jpint-1);  % 1-dim addresses in the PadFT array.
    for i1=-nw2:nw2
        for j1=-nw2:nw2
            plane=plane+w1(ipfrac,i1+nw2+1).*w1(jpfrac,j1+nw2+1)...
                .*P.PadFT(addrs+i1+(np1*j1));
        end;
    end;
    
    % % Convert the vector to a 2d array
    P2.PadFT(sp1+1:sp1+np,sp1+1:sp1+np)=reshape(plane,np,np);
    rotImages(:,:,ind)=gridRecoverRealImage(P2);
end;
