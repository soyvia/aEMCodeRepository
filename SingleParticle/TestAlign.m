% TestAlign

figure(1);
SetGrayscale;
load /Users/fred/matlabWork/EMProcessing/testimg  % loads img
n=size(img,1);
nim=10;
refs=single(zeros(n,n,2));
msk=fuzzymask(n,2,n*0.4,n*0.05);
refs(:,:,1)=img.*msk;
refs(:,:,2)=circshift(flipud(img),[1 0]).*msk;  % flipped reference

msk2=fuzzymask(n,2,n*.45,n*.05);
ref=img;
imgs=single(zeros(n,n,nim));
sigma=0;
dt=5*pi/180;  % 10 degrees per step.
sx=0;
sy=0;

for i=1:nim
    theta=(i-1)*dt;
%     sx=i/4;
%     sy=-i/3;
    imgs(:,:,i)=shiftf(grotate(img,theta)+sigma*randn(n,n),[-sx -sy]).*msk2;  % cw rotation
end;
sp=spCreateStackParams(imgs,2.9);

% imovie(imgs,.1);

mode=3;

[sp aliImgs]=spAlignRot(imgs,sp,refs(:,:,1));
vals=[sp.rot sp.flip]

% [sp aliImgs]=spTransformStack(imgs, sp);
figure(1);
imovie(aliImgs,.1);