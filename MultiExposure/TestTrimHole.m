% TestTrimHole.m
% Experiment with auto-trimming of thick carbon

mainPath='/Volumes/TetraData/EMWork/Hideki/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Merge/';
filename='DE_20120122_220714_719m.mrc';
mainPath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/Merged/';

mainPath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/Merged/';
infoPath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/Info/';
filename='004_sq02_1_04';
postscr='m.mrc';

% d=dir(mainPath);
% dirIndex=6;
% filename=d(dirIndex).name

figure(1);
SetGrayscale;

[m pixA]=ReadEMFile([mainPath filename postscr]);
load([infoPath filename 'mi.mat']);

% m=rot90(m,3);  % old file is not rotated.
%%
n0=256;
n=size(m);
fm=fftn(m);

%%

f=RadiusNorm(n)/pixA;  % Frequency in inverse A

% Butterworth bandpass
hexp=8;
fh=.005;
fl=pixA/30;
% fh=.01;
fl=.2;
H=1./(1+(fh./f).^hexp+(f./fl).^hexp);
imacs(H);
normWidth=.3;
binFactor=round(n(1)/n0);
frac=.5;
nw=n/binFactor;
width=20*normWidth;
erw=width;

% Compute the local variance, binned
A=BinImage(abs(ifftn(fm.*ifftshift(H)).^2),binFactor);
subplot(231);
imacs(A.^.2);
subplot(233);
imacs(BinImage(m,binFactor/4));
q=A(:);
nq=numel(q);
% h=hist(q.^.2,1000);
% semilogy(h);
subplot(234);
qsrt=sort(q);
thrMin=qsrt(round(nq*.25));  % get the 1% value
thrMax=qsrt(round(nq*.92));  % get the 99% value
% thr=qsrt(round(nq*frac));
thr=thrMin+frac*(thrMax-thrMin);

msk=A>thr;
imacs(msk);
subplot(235);
expmask=GaussFiltDCT(msk,.133/width)>.1;
imacs(expmask);
subplot(236);
em2=GaussFiltDCT(expmask,.133/erw)>.9;
imacs(em2);
subplot(232);
imacs(GaussFilt(Downsample(1-em2,n/4),.05).*BinImage(m,4));

%% Extend to edge

emsk=1-em2;   % =1 in interior
nmin=40;
mr=emsk;
for irot=1:4
    for j=1:nw(2)  % scan each row
        if mr(1,j)==1  % a border pixel is enabled
            q=find(mr(:,j)==0,1);
            if numel(q)>0 && q<nmin
                mr(1:q,j)=0;
            end;
        end;
    end;
    mr=rot90(mr);
end;
imacs(mr);
emsk=mr;
%% Save the mask in the mi structure

mask.merge='AND';
mask.encoding='RLE';  % RLE
mask.data=RunLengthEncode(emsk);

maskIndex=3;  % This is index for automasking
if ~isfield(mi,'mask')
    mi.mask=struct('merge',[],'encoding',[],'data',[]);
end;
    mi.mask(maskIndex)=mask;

save([infoPath filename 'mi.mat'],'mi');


%%
emsk=RunLengthDecode(mi.mask(3).data);
subplot(1,1,1);
emskx=ExpandImage(single(emsk),size(m,1)/(4*size(emsk,1)));
cMsk=repmat(rot90(emskx),[1 1 3]);
cMsk(:,:,3)=max(cMsk(:,:,3),.7);
cMsk(:,:,1:2)=max(cMsk(:,:,1:2),.3);
mscl=imscale(Downsample(m,n/4),256,1e-3)/256;
cImg=repmat(rot90(mscl),[1 1 3]);  % grayscale image
image(cMsk.*cImg);
