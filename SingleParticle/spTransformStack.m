function [sp, newStack]=spTransformStack(stack,sp)
n=size(stack,1);  % get the stack dimension
ds=sp.boxSize/n;  % stack downsampling ratio
npar=size(sp.trans,1);
xytfr=[sp.trans/ds sp.rot*pi/180 single(sp.flip)...
    single(sp.class) zeros(npar,1)];
newStack=TransformImages(stack, xytfr);

% mark that we've done the translation and rotation
sp.trans=0*sp.trans;
sp.rot=0*sp.rot;
sp.flip=0*sp.flip;
